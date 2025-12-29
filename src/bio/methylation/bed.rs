use crate::bio::dna::SequenceChunkId;
use crate::bio::dna::domain::SequenceRegion;
use crate::bio::error::ParseError;
use crate::bio::methylation::domain::MethylationValue;
use crate::bio::methylation::parse::{
    MethylationContext, MethylationError, MinimalRecord, get_records, parse_chrom, parse_end,
    parse_start, try_advance_field,
};
use crate::bio::methylation::{SequenceMethylation, SingleMethylation};
use crate::error::{CoreError, InvalidInputError, MalformedInputError};
use std::collections::{HashMap, HashSet};
use std::ops::Range;
use std::path::Path;
use std::rc::Rc;

pub(in crate::bio::methylation) fn parse_methylation_table<S: RecordSource>(
    ctx: &mut MethylationContext,
    path: &Path,
    sequences_ranges: &[(&str, usize)],
) -> Result<HashMap<Box<str>, SequenceMethylation>, CoreError> {
    let iter = S::records(path, None)?;
    let mut chroms: HashMap<String, HashMap<usize, MethylationValue>> = HashMap::new();

    let allowed_ids_set = sequences_ranges.iter().map(|t| t.0).collect::<HashSet<_>>();

    for rec in iter {
        let Some(rec) = rec? else { continue };
        if !ctx.check_record(&rec).map_err(MethylationError::Record)? {
            continue;
        }

        if !allowed_ids_set.contains(rec.chrom.as_str()) {
            continue;
        }

        chroms
            .entry(rec.chrom)
            .or_default()
            .insert(rec.start, rec.percent_modified);
    }

    let chroms = chroms
        .into_iter()
        .map(|(id, sites)| (id.into_boxed_str(), SequenceMethylation { sites }))
        .collect();

    Ok(chroms)
}

pub(in crate::bio::methylation) fn parse_single_id<S: RecordSource>(
    ctx: &mut MethylationContext,
    path: &Path,
    id: Rc<String>,
) -> Result<Option<SingleMethylation>, CoreError> {
    if id.is_empty() {
        return Err(InvalidInputError("empty id".into()).into());
    }

    let iter = S::records(path, Some(SequenceRegion::Id { id: Rc::clone(&id) }))?;
    let mut sites = HashMap::new();

    let mut max_index = 0;
    for rec in iter {
        let rec = rec?;
        let Some(rec) = rec else { continue };
        if rec.chrom != id.as_str() || !ctx.check_record(&rec).map_err(MethylationError::Record)? {
            continue;
        }

        sites.insert(rec.start, rec.percent_modified);
        max_index = max_index.max(rec.start);
    }

    Ok(sites_to_single_methylation(sites, SequenceChunkId::Id(id)))
}

pub(in crate::bio::methylation) fn parse_single_id_slice<S: RecordSource>(
    ctx: &mut MethylationContext,
    path: &Path,
    id: impl Into<String>,
    local_range: Range<usize>,
) -> Result<Option<SingleMethylation>, CoreError> {
    let id = id.into();
    if id.is_empty() || local_range.is_empty() {
        let mut problems = Vec::with_capacity(2);
        if id.is_empty() {
            problems.push("empty id".to_owned());
        }
        if local_range.is_empty() {
            problems.push(format!(
                "invalid range: [{}, {}]",
                local_range.start, local_range.end
            ));
        }
        return Err(InvalidInputError(problems.join(", ").into()).into());
    }

    let iter = S::records(path, None)?;
    let mut sites = HashMap::new();

    for rec in iter {
        let rec = rec?;
        let Some(rec) = rec else { continue };
        if !local_range.contains(&rec.start)
            || rec.chrom != id.as_str()
            || !ctx.check_record(&rec).map_err(MethylationError::Record)?
        {
            continue;
        }
        sites.insert(rec.start - local_range.start, rec.percent_modified);
    }

    Ok(sites_to_single_methylation(
        sites,
        SequenceChunkId::Id(Rc::new(id)),
    ))
}

pub(in crate::bio::methylation) fn parse_single_global_slice<S: RecordSource>(
    ctx: &mut MethylationContext,
    path: &Path,
    global_range: Range<usize>,
    ordered_sequences: &[(&str, usize)],
) -> Result<Option<SingleMethylation>, CoreError> {
    if global_range.is_empty() {
        return Err(InvalidInputError(
            format!(
                "invalid range: [{}, {}]",
                global_range.start, global_range.end
            )
            .into(),
        )
        .into());
    }

    let mut offsets = HashMap::with_capacity(ordered_sequences.len());
    let mut cum = 0;
    for &(id, len) in ordered_sequences {
        offsets.insert(id, cum);
        cum += len;
    }
    let iter = S::records(path, None)?;
    let mut sites = HashMap::new();
    let mut max_index = 0;

    for rec in iter {
        let rec = rec?;
        let Some(rec) = rec else { continue };
        let Some(off) = offsets.get(rec.chrom.as_str()) else {
            continue;
        };
        let global_pos = off + rec.start;

        if !global_range.contains(&global_pos)
            || !ctx.check_record(&rec).map_err(MethylationError::Record)?
        {
            continue;
        }

        sites.insert(global_pos - global_range.start, rec.percent_modified);
        max_index = max_index.max(rec.start);
    }

    Ok(sites_to_single_methylation(
        sites,
        SequenceChunkId::Global(global_range),
    ))
}

fn sites_to_single_methylation(
    sites: HashMap<usize, MethylationValue>,
    id: SequenceChunkId,
) -> Option<SingleMethylation> {
    if sites.is_empty() {
        None
    } else {
        Some(SingleMethylation {
            id,
            methylation: SequenceMethylation { sites },
        })
    }
}

pub(in crate::bio::methylation) trait RecordSource {
    type Iter<'a>: Iterator<Item = Result<Option<MinimalRecord>, CoreError>>
    where
        Self: 'a;
    fn records(path: &Path, region: Option<SequenceRegion>) -> Result<Self::Iter<'_>, CoreError>;
}

pub(in crate::bio::methylation) struct BedSource;

impl RecordSource for BedSource {
    type Iter<'a> = Box<dyn Iterator<Item = Result<Option<MinimalRecord>, CoreError>>>;

    fn records(path: &Path, region: Option<SequenceRegion>) -> Result<Self::Iter<'_>, CoreError> {
        get_records(path, region, parse_line_bed_minimal)
    }
}

/// Parse one BED line into a full `BedFull` record.
///
/// Currently unused by the `MinimalRecord` pipeline.
#[allow(unused)]
fn parse_line_bed_minimal(
    line: &str,
    region_id: Option<&Rc<String>>,
    region_range: Option<&Range<usize>>,
) -> Result<Option<MinimalRecord>, ParseError> {
    let mut columns = line.split_whitespace();

    let chrom = parse_chrom(&mut columns)?;

    if !region_id.is_none_or(|id| id.as_str() == chrom.as_str()) {
        return Ok(None);
    }

    let start = parse_start(&mut columns)?;

    let end = parse_end(&mut columns)?;

    if !region_range.is_none_or(|range| range.contains(&start)) {
        return Ok(None);
    }

    try_advance_field(
        &mut columns,
        |t| Ok::<_, ParseError>(()),
        "missing Modified Base Code (Name) column",
    )?;

    try_advance_field(
        &mut columns,
        |t| Ok::<_, ParseError>(()),
        "missing Score column",
    )?;

    try_advance_field(
        &mut columns,
        |t| Ok::<_, ParseError>(()),
        "missing Strand column",
    )?;

    try_advance_field(
        &mut columns,
        |t| Ok::<_, ParseError>(()),
        "missing ThickStart column",
    )?;

    try_advance_field(
        &mut columns,
        |t| Ok::<_, ParseError>(()),
        "missing ThickEnd column",
    )?;

    try_advance_field(
        &mut columns,
        |t| Ok::<_, ParseError>(()),
        "missing ItemRgb column",
    )?;

    try_advance_field(
        &mut columns,
        |t| Ok::<_, ParseError>(()),
        "missing Valid Coverage column",
    )?;

    let percent_modified = try_advance_field(
        &mut columns,
        |t| {
            t.parse().map_err(|err| {
                MalformedInputError(format!("invalid Percent Modified value ({err})").into())
            })
        },
        "missing Percent Modified column",
    )
    .map(MethylationValue::from_percent)?;

    Ok(Some(MinimalRecord {
        chrom,
        start,
        end,
        percent_modified,
    }))
}

#[cfg(test)]
mod bed_tests {
    use super::*;

    #[test]
    fn bed_empty_line() {
        let r = parse_line_bed_minimal("", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Chrom column");
    }

    #[test]
    fn bed_whitespace_line() {
        let r = parse_line_bed_minimal("   ", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Chrom column");
    }

    #[test]
    fn bed_missing_start() {
        let r = parse_line_bed_minimal("chr1", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing ChromStart column");
    }

    #[test]
    fn bed_missing_end() {
        let r = parse_line_bed_minimal("chr1 100", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing ChromEnd column");
    }

    #[test]
    fn bed_missing_modified_code() {
        let r = parse_line_bed_minimal("chr1 100 200", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Modified Base Code (Name) column");
    }

    #[test]
    fn bed_missing_score() {
        let r = parse_line_bed_minimal("chr1 100 200 X", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Score column");
    }

    #[test]
    fn bed_missing_strand() {
        let r = parse_line_bed_minimal("chr1 100 200 X 0", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Strand column");
    }

    #[test]
    fn bed_missing_thick_start() {
        let r = parse_line_bed_minimal("chr1 100 200 X 0 +", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing ThickStart column");
    }

    #[test]
    fn bed_missing_thick_end() {
        let r = parse_line_bed_minimal("chr1 100 200 X 0 + 100", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing ThickEnd column");
    }

    #[test]
    fn bed_missing_item_rgb() {
        let r = parse_line_bed_minimal("chr1 100 200 X 0 + 100 200", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing ItemRgb column");
    }

    #[test]
    fn bed_missing_valid_coverage() {
        let r = parse_line_bed_minimal("chr1 100 200 X 0 + 100 200 255,0,0", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Valid Coverage column");
    }

    #[test]
    fn bed_missing_percent_modified() {
        let r = parse_line_bed_minimal("chr1 100 200 X 0 + 100 200 255,0,0 5", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Percent Modified column");
    }

    #[test]
    fn bed_start_not_number() {
        let r = parse_line_bed_minimal("chr1 xyz 200 X 0 + 100 200 255,0,0 5 10.0", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("ChromStart"));
    }

    #[test]
    fn bed_end_not_number() {
        let r = parse_line_bed_minimal("chr1 100 xyz X 0 + 100 200 255,0,0 5 10.0", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("ChromEnd"));
    }

    #[test]
    fn bed_percent_modified_not_float() {
        let r = parse_line_bed_minimal("chr1 100 200 X 0 + 100 200 255,0,0 5 nope", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("Percent Modified"));
    }

    #[test]
    fn bed_line_ok() {
        let r = parse_line_bed_minimal("chr1 100 200 m 0 + 100 200 255,0,0 5 7.5", None, None)
            .unwrap()
            .unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert_eq!(r.percent_modified, MethylationValue::from_percent(7.5));
    }

    #[test]
    fn bed_negative_start() {
        let r = parse_line_bed_minimal("chr1 -100 200 X 0 + 100 200 255,0,0 5 10.0", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("ChromStart"));
    }

    #[test]
    fn bed_negative_end() {
        let r = parse_line_bed_minimal("chr1 100 -200 X 0 + 100 200 255,0,0 5 10.0", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("ChromEnd"));
    }
}
