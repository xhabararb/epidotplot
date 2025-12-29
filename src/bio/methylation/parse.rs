use crate::bio::MethylationPath;
use crate::bio::dna::domain::SequenceRegion;
use crate::bio::error::{IoError, ParseError};
use crate::bio::methylation::bed::{
    BedSource, RecordSource, parse_methylation_table, parse_single_global_slice, parse_single_id,
    parse_single_id_slice,
};
use crate::bio::methylation::bedgraph::BedGraphSource;
use crate::bio::methylation::domain::MethylationValue;
use crate::bio::methylation::parse::RecordError::PositionRangeNotPlusOne;
use crate::bio::methylation::{MethylationSites, SingleMethylation};
use crate::config::MethylationThreshold;
use crate::error::{CoreError, MalformedInputError};
use std::borrow::Cow;
use std::fmt::Formatter;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::{Range, RangeInclusive};
use std::path::Path;
use std::rc::Rc;
use std::{error, fmt};
use thiserror::Error;

/// Minimal methylation record.
///
/// This is the reduced data actually needed downstream:
/// - `chrom`: chromosome identifier
/// - `start` / `end`: genomic interval (usually single-base)
/// - `methylation`: fraction or percentage of methylated reads
///
/// All other BED/bedGraph fields are ignored here.
pub struct MinimalRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub percent_modified: MethylationValue,
}

/// Trait implemented by pluggable methylation parsers (BED, bedGraph).
///
/// Each parser turns input files into a `MethylationTable`.
pub trait MethylationParser: fmt::Debug {
    fn parse(
        &self,
        path: &MethylationPath,
        region: Option<SequenceRegion>,
        ordered_sequences: &[(&str, usize)],
        println: &dyn Fn(Cow<'_, str>),
    ) -> Result<Option<SingleMethylation>, CoreError> {
        match path {
            MethylationPath::Bed(p) => self.parse_bed(p, region, ordered_sequences, println),
            MethylationPath::BedGraph(p) => {
                self.parse_bedgraph(p, region, ordered_sequences, println)
            }
        }
    }
    fn parse_bed(
        &self,
        path: &Path,
        region: Option<SequenceRegion>,
        ordered_sequences: &[(&str, usize)],
        println: &dyn Fn(Cow<'_, str>),
    ) -> Result<Option<SingleMethylation>, CoreError>;
    fn parse_bedgraph(
        &self,
        path: &Path,
        region: Option<SequenceRegion>,
        ordered_sequences: &[(&str, usize)],
        println: &dyn Fn(Cow<'_, str>),
    ) -> Result<Option<SingleMethylation>, CoreError>;
}

#[derive(Debug)]
pub struct NativeMethylationParser {
    pub min_threshold: MethylationThreshold,
}

impl MethylationParser for NativeMethylationParser {
    fn parse_bed(
        &self,
        path: &Path,
        region: Option<SequenceRegion>,
        ordered_sequences: &[(&str, usize)],
        println: &dyn Fn(Cow<'_, str>),
    ) -> Result<Option<SingleMethylation>, CoreError> {
        parse_methylation::<BedSource>(path, region, ordered_sequences, println)
    }

    fn parse_bedgraph(
        &self,
        path: &Path,
        region: Option<SequenceRegion>,
        ordered_sequences: &[(&str, usize)],
        println: &dyn Fn(Cow<'_, str>),
    ) -> Result<Option<SingleMethylation>, CoreError> {
        parse_methylation::<BedGraphSource>(path, region, ordered_sequences, println)
    }
}

pub(in crate::bio::methylation) struct MethylationContext {
    min_value: Option<MethylationValue>,
    max_value: Option<MethylationValue>,
    /// If true, ranges other than [n, n+1] are leniently ignored, as opposed to failing.
    forgive_non_plus_one: bool,
    range: Option<RangeInclusive<MethylationValue>>,
}

impl MethylationContext {
    #[inline]
    pub const fn new_with_estimate(forgive_non_plus_one: bool) -> Self {
        Self {
            min_value: None,
            max_value: None,
            forgive_non_plus_one,
            range: None,
        }
    }
}

#[derive(Debug, Error)]
pub enum MethylationError {
    Record(RecordError),
    Other(Box<dyn error::Error>),
}

impl fmt::Display for MethylationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Methylation error: {}",
            match self {
                Self::Record(err) => err,
                Self::Other(err) => err.as_ref(),
            }
        )
    }
}

#[derive(Debug, Error)]
pub enum RecordError {
    ValueOutOfRange(MethylationValue, RangeInclusive<MethylationValue>),
    PositionRangeNotPlusOne(RangeInclusive<usize>),
}

impl fmt::Display for RecordError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::ValueOutOfRange(value, range) => {
                write!(
                    f,
                    "value {value} is out of range (expected [{}, {}])",
                    range.start(),
                    range.end()
                )
            }
            PositionRangeNotPlusOne(range) => {
                write!(
                    f,
                    "invalid interval [{}, {}]: end must equal start + 1",
                    range.start(),
                    range.end()
                )
            }
        }
    }
}

impl MethylationContext {
    pub fn check_record(&mut self, record: &MinimalRecord) -> Result<bool, RecordError> {
        if record.end != record.start + 1 {
            if self.forgive_non_plus_one {
                return Ok(false);
            }
            return Err(PositionRangeNotPlusOne(record.start..=record.end));
        }

        let methylation_value = record.percent_modified;

        /*if methylation_value < *self.min_threshold {
            return Ok(false);
        }*/

        if let Some(range) = &self.range
            && !range.contains(&methylation_value)
        {
            return Err(RecordError::ValueOutOfRange(
                methylation_value,
                range.clone(),
            ));
        }
        // range is not specified, so we store min and max for bound estimation
        self.min_value = Some(
            self.min_value
                .map_or(methylation_value, |t| t.min(methylation_value)),
        );
        self.max_value = Some(
            self.max_value
                .map_or(methylation_value, |t| t.max(methylation_value)),
        );

        Ok(true)
    }
}

fn parse_methylation<S: RecordSource>(
    path: &Path,
    region: Option<SequenceRegion>,
    sequences_ranges: &[(&str, usize)],
    println: &dyn Fn(Cow<'_, str>),
) -> Result<Option<SingleMethylation>, CoreError> {
    let mut ctx = MethylationContext::new_with_estimate(false);
    let res = match region {
        None => {
            let table = parse_methylation_table::<S>(&mut ctx, path, sequences_ranges)?;
            Ok::<Option<MethylationSites>, CoreError>(Some(MethylationSites::Multiple(table)))
        }

        Some(SequenceRegion::Id { id }) => {
            Ok(parse_single_id::<S>(&mut ctx, path, Rc::clone(&id))?.map(MethylationSites::Single))
        }

        Some(SequenceRegion::IdSlice { id, range }) => {
            Ok(
                parse_single_id_slice::<S>(&mut ctx, path, id.as_str(), range)?
                    .map(MethylationSites::Single),
            )
        }

        Some(SequenceRegion::Slice { range }) => {
            Ok(
                parse_single_global_slice::<S>(&mut ctx, path, range, sequences_ranges)?
                    .map(MethylationSites::Single),
            )
        }
    }?;

    // methylation values are loaded as %
    match res {
        Some(res) => {
            if let Some(min) = ctx.min_value
                && let Some(max) = ctx.max_value
            {
                let min = min.as_percent();
                let max = max.as_percent();

                // all values fall in [0,1], so they are quite probably fractions rather than % -> multiply by 100
                // e.g., 0.1 is 10%, not 0.1%
                let f_range = 0.0..=1.0;
                if f_range.contains(&min) && f_range.contains(&max) {
                    println(Cow::Borrowed(
                        "All methylation values fall into [0,1], so they are treated as fractions, not %",
                    ));
                    return Ok(res.collapse(sequences_ranges).map(|mut t| {
                        t.correct_percents_into_fractions();
                        t
                    }));
                }

                // the values are valid %, so no transformation is needed
                let p_range = 0.0..=100.0;
                if p_range.contains(&min) && p_range.contains(&max) {
                    return Ok(res.collapse(sequences_ranges));
                }

                // values fit neither percent range [0,100] nor fraction range [0,1]
                return Err(MethylationError::Other(
                    "cannot infer methylation scale (not in [0,1] or [0,100])".into(),
                )
                .into());
            }
            Ok(res.collapse(sequences_ranges))
        }
        None => Ok(res.and_then(|t| t.collapse(sequences_ranges))),
    }
}

pub(in crate::bio::methylation) fn try_advance_field<'a, I, T, F, E>(
    iter: &mut I,
    parse: F,
    missing_column_err: &'static str,
) -> Result<T, ParseError>
where
    I: Iterator<Item = &'a str>,
    F: Fn(&'a str) -> Result<T, E>,
    E: Into<ParseError>,
{
    let s = iter
        .next()
        .ok_or_else(|| MalformedInputError(missing_column_err.into()))?;
    parse(s).map_err(Into::into)
}

pub(in crate::bio::methylation) fn parse_chrom<'a, I: Iterator<Item = &'a str>>(
    mut columns: &mut I,
) -> Result<String, ParseError> {
    try_advance_field(
        &mut columns,
        |t| Ok::<_, ParseError>(t.to_owned()),
        "missing Chrom column",
    )
}

pub(in crate::bio::methylation) fn parse_start<'a, I: Iterator<Item = &'a str>>(
    mut columns: &mut I,
) -> Result<usize, ParseError> {
    try_advance_field(
        &mut columns,
        |t| {
            t.parse().map_err(|err| {
                MalformedInputError(format!("invalid ChromStart value ({err})").into())
            })
        },
        "missing ChromStart column",
    )
}

pub(in crate::bio::methylation) fn parse_end<'a, I: Iterator<Item = &'a str>>(
    mut columns: &mut I,
) -> Result<usize, ParseError> {
    try_advance_field(
        &mut columns,
        |t| {
            t.parse().map_err(|err| {
                MalformedInputError(format!("invalid ChromEnd value ({err})").into())
            })
        },
        "missing ChromEnd column",
    )
}

#[allow(clippy::type_complexity)]
pub(in crate::bio::methylation) fn get_records<P: AsRef<Path>, F>(
    path: P,
    region: Option<SequenceRegion>,
    line_parser: F,
) -> Result<Box<dyn Iterator<Item = Result<Option<MinimalRecord>, CoreError>>>, CoreError>
where
    F: Fn(
            &str,
            Option<&Rc<String>>,
            Option<&Range<usize>>,
        ) -> Result<Option<MinimalRecord>, ParseError>
        + 'static,
{
    let path_buf = path.as_ref().to_path_buf();
    let file = File::open(path).map_err(|err| IoError {
        err,
        message: Some(format!("at {}", path_buf.display())),
    })?;

    let (region_id, region_range) = region.map_or((None, None), SequenceRegion::into_parts);

    let reader = BufReader::new(file);
    Ok(Box::new(reader.lines().map(move |line_res| {
        match line_res {
            Ok(line) => line_parser(&line, region_id.as_ref(), region_range.as_ref())
                .map_err(|e| MalformedInputError(e.into()).into()),
            Err(err) => Err(IoError {
                err,
                message: Some(format!(
                    "Failed to read a line from an opened file at {}",
                    path_buf.display()
                )),
            }
            .into()),
        }
    })))
}
