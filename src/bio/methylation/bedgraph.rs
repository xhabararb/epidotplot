use crate::bio::dna::domain::SequenceRegion;
use crate::bio::error::ParseError;
use crate::bio::methylation::bed::RecordSource;
use crate::bio::methylation::domain::MethylationValue;
use crate::bio::methylation::parse::{
    MinimalRecord, get_records, parse_chrom, parse_end, parse_start, try_advance_field,
};
use crate::error::{CoreError, MalformedInputError};
use std::ops::Range;
use std::path::Path;
use std::rc::Rc;

pub(in crate::bio::methylation) struct BedGraphSource;

impl RecordSource for BedGraphSource {
    type Iter<'a> = Box<dyn Iterator<Item = Result<Option<MinimalRecord>, CoreError>>>;

    fn records(path: &Path, region: Option<SequenceRegion>) -> Result<Self::Iter<'_>, CoreError> {
        get_records(path, region, parse_line_bedgraph_minimal)
    }
}

/// Parse one bedGraph line into a `MinimalRecord`.
fn parse_line_bedgraph_minimal(
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

    let percent_modified = try_advance_field(
        &mut columns,
        |t| {
            t.parse().map_err(|err| {
                MalformedInputError(format!("invalid methylation DataValue ({err})").into())
            })
        },
        "missing methylation DataValue column",
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
mod bedgraph_tests {
    use super::*;

    #[test]
    fn empty_line() {
        let r = parse_line_bedgraph_minimal("", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Chrom column");
    }

    #[test]
    fn whitespace_line() {
        let r = parse_line_bedgraph_minimal("  ", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert_eq!(format!("{e}"), "missing Chrom column");
    }

    #[test]
    fn missing_start() {
        let r = parse_line_bedgraph_minimal("chr1", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("missing ChromStart column"));
    }

    #[test]
    fn missing_end() {
        let r = parse_line_bedgraph_minimal("chr1 100", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("missing ChromEnd column"));
    }

    #[test]
    fn missing_percent_modified() {
        let r = parse_line_bedgraph_minimal("chr1 100 200", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("missing methylation DataValue column"));
    }

    #[test]
    fn start_not_number() {
        let r = parse_line_bedgraph_minimal("chr1 xyz 200 0.5", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("ChromStart"));
    }

    #[test]
    fn end_not_number() {
        let r = parse_line_bedgraph_minimal("chr1 100 xyz 0.5", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("ChromEnd"));
    }

    #[test]
    fn value_not_float() {
        let r = parse_line_bedgraph_minimal("chr1 100 200 nope", None, None);
        assert!(matches!(r, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(e))) = r else {
            unreachable!()
        };
        assert!(format!("{e}").contains("methylation DataValue"));
    }

    #[test]
    fn line_ok() {
        let line = "chr1 100 200 0.5";
        let record = parse_line_bedgraph_minimal(line, None, None)
            .unwrap()
            .unwrap();
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.start, 100);
        assert_eq!(record.end, 200);
        assert_eq!(record.percent_modified, MethylationValue::from_percent(0.5));
    }

    #[test]
    fn negative_start() {
        let line = "chr1 -100 200 0.5";
        let record = parse_line_bedgraph_minimal(line, None, None);
        assert!(record.is_err());
        assert!(matches!(record, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(err))) = record else {
            unreachable!()
        };
        assert!(format!("{err}").contains("ChromStart"));
    }

    #[test]
    fn negative_end() {
        let line = "chr1 100 -200 0.5";
        let record = parse_line_bedgraph_minimal(line, None, None);
        assert!(record.is_err());
        assert!(matches!(record, Err(ParseError::MalformedInput(_))));
        let Err(ParseError::MalformedInput(MalformedInputError(err))) = record else {
            unreachable!()
        };
        assert!(format!("{err}").contains("ChromEnd"));
    }
}

#[cfg(test)]
mod methyl_tests {
    use super::*;
    use crate::bio::MethylationPath;
    use crate::bio::dna::SequenceChunkId;
    use crate::bio::dna::domain::SequenceRegion;
    use crate::bio::methylation::parse::{MethylationParser, NativeMethylationParser};
    use crate::error::CoreError;
    use std::rc::Rc;
    use std::{fs::File, io::Write, path::PathBuf};

    fn write_tmp(name: &str, content: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!("{}_test.bed", name));
        let mut f = File::create(&p).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        p
    }

    fn mv(p: f32) -> MethylationValue {
        MethylationValue::from_percent(p)
    }

    #[test]
    fn bed_multi_chromosome_collapses_in_ordered_sequences_order() {
        let bed = "\
chrX\t0\t1\tm\t0\t+\t0\t1\t255,0,0\t8\t10.0\t7\t1\t0\t0\t4\t0\t2\n\
chr1\t0\t1\tm\t0\t+\t0\t1\t255,0,0\t8\t20.0\t7\t1\t0\t0\t4\t0\t2\n\
chrX\t1\t2\tm\t0\t+\t1\t2\t255,0,0\t8\t30.0\t7\t1\t0\t0\t4\t0\t2\n\
chr3\t1\t2\tm\t0\t+\t1\t2\t255,0,0\t8\t40.0\t7\t1\t0\t0\t4\t0\t2\n\
chr1\t1\t2\tm\t0\t+\t1\t2\t255,0,0\t8\t40.0\t7\t1\t0\t0\t4\t0\t2\n\
";
        let path = write_tmp("bed_multi", bed);
        let ordered = vec![("chr1", 2), ("chrX", 2)];
        let p = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let r = p
            .parse(&MethylationPath::Bed(path), None, &ordered, &|_| ())
            .unwrap()
            .unwrap();

        assert_eq!(r.id, SequenceChunkId::Global(0..4));
        let s = r.methylation.sites;

        assert_eq!(s[&0], mv(20.0));
        assert_eq!(s[&1], mv(40.0));
        assert_eq!(s[&2], mv(10.0));
        assert_eq!(s[&3], mv(30.0));
    }

    #[test]
    fn bed_non_plus_one_interval_errors() {
        let bed = "\
chr1\t0\t2\tm\t0\t+\t0\t2\t255,0,0\t8\t10.0\t7\t1\t0\t0\t4\t0\t2\n\
";
        let path = write_tmp("np1_err", bed);
        let ordered = vec![("chr1", 3)];
        let p = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let err = p.parse(&MethylationPath::Bed(path), None, &ordered, &|_| ());
        assert!(matches!(err, Err(CoreError::Methylation(_))));
    }

    #[test]
    fn bed_global_slice_shifts_positions() {
        let bed = "\
a\t0\t1\tm\t0\t+\t0\t1\t255,0,0\t8\t10.0\t7\t1\t0\t0\t4\t0\t2\n\
a\t1\t2\tm\t0\t+\t1\t2\t255,0,0\t8\t20.0\t7\t1\t0\t0\t4\t0\t2\n\
b\t0\t1\tm\t0\t+\t0\t1\t255,0,0\t8\t30.0\t7\t1\t0\t0\t4\t0\t2\n\
b\t1\t2\tm\t0\t+\t1\t2\t255,0,0\t8\t40.0\t7\t1\t0\t0\t4\t0\t2\n\
";
        let path = write_tmp("g_slice", bed);
        let ordered = vec![("a", 2), ("b", 2)];
        let p = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let f = SequenceRegion::Slice { range: 1..3 };
        let r = p
            .parse(&MethylationPath::Bed(path), Some(f), &ordered, &|_| ())
            .unwrap()
            .unwrap();

        let s = r.methylation.sites;
        assert_eq!(s.len(), 2);
        assert_eq!(s[&0], mv(20.0));
        assert_eq!(s[&1], mv(30.0));
    }

    #[test]
    fn bed_id_region_keeps_only_requested_chrom() {
        let bed = "\
b\t2\t3\tm\t0\t+\t0\t1\t255,0,0\t8\t20.0\t7\t1\t0\t0\t4\t0\t2\n\
a\t0\t1\tm\t0\t+\t0\t1\t255,0,0\t8\t10.0\t7\t1\t0\t0\t4\t0\t2\n\
b\t0\t1\tm\t0\t+\t0\t1\t255,0,0\t8\t2.0\t7\t1\t0\t0\t4\t0\t2\n\
";
        let path = write_tmp("id", bed);
        let ordered = vec![("a", 1), ("b", 2)];
        let p = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let f = SequenceRegion::Id {
            id: Rc::new("b".into()),
        };
        let r = p
            .parse(&MethylationPath::Bed(path), Some(f), &ordered, &|_| ())
            .unwrap()
            .unwrap();

        let s = r.methylation.sites;
        assert_eq!(s.len(), 2);
        assert_eq!(s[&0], mv(2.0));
        assert_eq!(s[&2], mv(20.0));
    }

    #[test]
    fn bed_threshold_filters_out_low_values() {
        let bed = "\
chr1\t0\t1\tm\t0\t+\t0\t1\t255,0,0\t8\t5.0\t7\t1\t0\t0\t4\t0\t2\n\
chr1\t1\t2\tm\t0\t+\t1\t2\t255,0,0\t8\t50.0\t7\t1\t0\t0\t4\t0\t2\n\
";
        let path = write_tmp("thr", bed);
        let ordered = vec![("chr1", 2)];
        let p = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(10.0).into(),
        };

        let mut r = p
            .parse(&MethylationPath::Bed(path), None, &ordered, &|_| ())
            .unwrap()
            .unwrap();
        r.filter_by_methylation_threshold(p.min_threshold);
        let s = r.methylation.sites;
        assert_eq!(s.len(), 1);
        assert_eq!(s[&1], mv(50.0));
    }
}
