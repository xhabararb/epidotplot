use crate::bio::methylation::domain::MethylationValue;
use std::io::Write;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

#[allow(dead_code)]
fn write_tmp(name: &str, ext: &str, content: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let p = std::env::temp_dir().join(format!("{name}_{nanos}.{ext}"));
    let mut f = std::fs::File::create(&p).unwrap();
    f.write_all(content.as_bytes()).unwrap();
    p
}

#[allow(dead_code)]
fn mv(p: f32) -> MethylationValue {
    MethylationValue::from_percent(p)
}

#[allow(dead_code)]
fn bed_line(id: &str, start: usize, end: usize, percent: f32) -> String {
    format!(
        "{id}\t{start}\t{end}\tm\t0\t+\t{start}\t{end}\t255,0,0\t8\t{percent:.1}\t7\t1\t0\t0\t4\t0\t2\n"
    )
}

#[cfg(test)]
mod tests {
    use super::{bed_line, mv, write_tmp};
    use crate::bio::dna::domain::SequenceRegion;
    use crate::bio::dna::parse::{NativeSequenceParser, SequenceParser};
    use crate::bio::methylation::domain::MethylationValue;
    use crate::bio::methylation::parse::{MethylationParser, NativeMethylationParser};
    use crate::bio::{MethylationPath, SequencePath};
    use std::rc::Rc;

    #[test]
    fn integration_multifasta_bed_global_collapse() {
        let fasta = "\
>a
ACGT
AC
>b

>c
TT
TTGG
>d
AAA
";

        let bed = [
            bed_line("a", 1, 2, 10.0),
            bed_line("c", 3, 4, 20.0),
            bed_line("d", 0, 1, 30.0),
        ]
        .join("");

        let fasta_p = write_tmp("int_fa", "fa", fasta);
        let bed_p = write_tmp("int_bed", "bed", &bed);

        let seq_parser = NativeSequenceParser;
        let mut ordered = Vec::new();
        let seq = seq_parser
            .parse(&SequencePath::Fasta(fasta_p.clone()), None, &mut ordered)
            .unwrap()
            .unwrap();

        let ordered_pairs: Vec<(&str, usize)> = ordered
            .iter()
            .map(|(id, len)| (id.as_str(), *len))
            .collect();

        assert_eq!(ordered_pairs, vec![("a", 6), ("b", 0), ("c", 6), ("d", 3),]);

        assert_eq!(
            seq.data,
            b"ACGTAC"
                .iter()
                .chain(b"".iter())
                .chain(b"TTTTGG".iter())
                .chain(b"AAA".iter())
                .copied()
                .collect::<Vec<u8>>()
        );

        let mp = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let m = mp
            .parse(
                &MethylationPath::Bed(bed_p.clone()),
                None,
                &ordered_pairs,
                &|_| (),
            )
            .unwrap()
            .unwrap();

        let offs_a = 0;
        let offs_b = offs_a + 6;
        let offs_c = offs_b + 0;
        let offs_d = offs_c + 6;

        let s = m.methylation.sites;

        println!("{s:?}");
        assert_eq!(s.len(), 3);
        assert_eq!(s[&(offs_a + 1)], mv(10.0));
        assert_eq!(s[&(offs_c + 3)], mv(20.0));
        assert_eq!(s[&(offs_d + 0)], mv(30.0));
    }

    #[test]
    fn integration_multifasta_bed_global_slice() {
        let fasta = "\
>a
AAAAAA
>b
CCGG
CCGG
>c
TTCCGGTT
>d
GGCCGGCCGG
>e
AAAACAAA
";

        let bed = [
            bed_line("b", 1, 2, 10.0),
            bed_line("b", 6, 7, 12.0),
            bed_line("c", 2, 3, 20.0),
            bed_line("c", 4, 5, 25.0),
            bed_line("d", 3, 4, 30.0),
            bed_line("e", 3, 4, 99.0),
        ]
        .join("");

        let fasta_p = write_tmp("slice_fa", "fa", fasta);
        let bed_p = write_tmp("slice_bed", "bed", &bed);

        let seq_parser = NativeSequenceParser;
        let mut ordered = Vec::new();

        let region = SequenceRegion::Slice { range: 9..25 };

        let seq = seq_parser
            .parse(
                &SequencePath::Fasta(fasta_p.clone()),
                Some(region.clone()),
                &mut ordered,
            )
            .unwrap()
            .unwrap();

        let ordered_pairs: Vec<(&str, usize)> = ordered
            .iter()
            .map(|(id, len)| (id.as_str(), *len))
            .collect();

        assert_eq!(ordered_pairs, vec![("a", 6), ("b", 8), ("c", 8), ("d", 10)]);

        assert_eq!(&seq.data[..], b"GCCGGTTCCGGTTGGC");

        let mp = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let m = mp
            .parse(
                &MethylationPath::Bed(bed_p.clone()),
                Some(region),
                &ordered_pairs,
                &|_| (),
            )
            .unwrap()
            .unwrap();

        let s = m.methylation.sites;

        println!("{s:?}");
        assert_eq!(s.len(), 3);
        assert_eq!(s[&3], mv(12.0));
        assert_eq!(s[&7], mv(20.0));
        assert_eq!(s[&9], mv(25.0));
    }

    #[test]
    fn integration_multifasta_bed_local_slice_misaligned() {
        let fasta = "\
>a
AAAAAA
>b
CCGG
CCGG
>c
TTCCGGTT
>d
GGCCGGCCGG
>e
AAAACAAA
";

        let bed = [
            bed_line("b", 1, 2, 10.0),
            bed_line("b", 6, 7, 12.0),
            bed_line("c", 2, 3, 20.0),
            bed_line("c", 4, 5, 25.0),
            bed_line("c", 7, 8, 12.0), // edge of an interval, expected to be ignored
            bed_line("d", 3, 4, 30.0),
            bed_line("e", 3, 4, 99.0),
        ]
        .join("");

        let fasta_p = write_tmp("local_fa", "fa", fasta);
        let bed_p = write_tmp("local_bed", "bed", &bed);

        let region = SequenceRegion::IdSlice {
            id: Rc::new("c".into()),
            range: 2..7,
        };

        let seq_parser = NativeSequenceParser;
        let mut ordered = Vec::new();
        let seq = seq_parser
            .parse(
                &SequencePath::Fasta(fasta_p.clone()),
                Some(region.clone()),
                &mut ordered,
            )
            .unwrap()
            .unwrap();

        let ordered_pairs: Vec<(&str, usize)> = ordered
            .iter()
            .map(|(id, len)| (id.as_str(), *len))
            .collect();

        assert_eq!(ordered_pairs, vec![("a", 6), ("b", 8), ("c", 8),]);

        // last 2 of b + c + first 2 of d
        assert_eq!(&seq.data[..], b"CCGGT");

        let mp = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let m = mp
            .parse(
                &MethylationPath::Bed(bed_p.clone()),
                Some(region),
                &ordered_pairs,
                &|_| (),
            )
            .unwrap()
            .unwrap();

        let s = &m.methylation.sites;

        println!("{s:?}");
        assert_eq!(s.len(), 2);
        assert_eq!(s[&0], mv(20.0));
        assert_eq!(s[&2], mv(25.0));

        let res = m.check_against_sequence(&seq);
        println!("{res:?}");
        assert!(matches!(res, Err(2)));
    }

    #[test]
    fn integration_multifasta_bed_local_slice_aligned() {
        let fasta = "\
>a
AAAAAA
>b
CCGG
CCGG
>c
TTCCCGTT
>d
GGCCGGCCGG
>e
AAAACAAA
";

        let bed = [
            bed_line("b", 1, 2, 10.0),
            bed_line("b", 6, 7, 12.0),
            bed_line("c", 2, 3, 20.0),
            bed_line("c", 4, 5, 25.0),
            bed_line("c", 7, 8, 12.0), // edge of an interval, expected to be ignored
            bed_line("d", 3, 4, 30.0),
            bed_line("e", 3, 4, 99.0),
        ]
        .join("");

        let fasta_p = write_tmp("local_fa", "fa", fasta);
        let bed_p = write_tmp("local_bed", "bed", &bed);

        let region = SequenceRegion::IdSlice {
            id: Rc::new("c".into()),
            range: 2..7,
        };

        let seq_parser = NativeSequenceParser;
        let mut ordered = Vec::new();
        let seq = seq_parser
            .parse(
                &SequencePath::Fasta(fasta_p.clone()),
                Some(region.clone()),
                &mut ordered,
            )
            .unwrap()
            .unwrap();

        let ordered_pairs: Vec<(&str, usize)> = ordered
            .iter()
            .map(|(id, len)| (id.as_str(), *len))
            .collect();

        assert_eq!(ordered_pairs, vec![("a", 6), ("b", 8), ("c", 8),]);

        // last 2 of b + c + first 2 of d
        assert_eq!(&seq.data[..], b"CCCGT");

        let mp = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let m = mp
            .parse(
                &MethylationPath::Bed(bed_p.clone()),
                Some(region),
                &ordered_pairs,
                &|_| (),
            )
            .unwrap()
            .unwrap();

        let s = &m.methylation.sites;

        println!("{s:?}");
        assert_eq!(s.len(), 2);
        assert_eq!(s[&0], mv(20.0));
        assert_eq!(s[&2], mv(25.0));

        let res = m.check_against_sequence(&seq);
        println!("{res:?}");
        assert!(res.is_ok());
    }

    #[test]
    fn integration_multifasta_bed_global_slice_threshold() {
        let fasta = "\
>a
AAAAAA
>b
CCGG
CCGG
>c
TTCCGGTT
>d
GGCCGGCCGG
>e
AAAACAAA
";

        let bed = [
            bed_line("b", 1, 2, 10.0),
            bed_line("b", 6, 7, 12.0),
            bed_line("c", 2, 3, 20.0),
            bed_line("c", 4, 5, 25.0),
            bed_line("d", 3, 4, 30.0),
            bed_line("e", 3, 4, 99.0),
        ]
        .join("");

        let fasta_p = write_tmp("slice_fa", "fa", fasta);
        let bed_p = write_tmp("slice_bed", "bed", &bed);

        let region = SequenceRegion::Slice { range: 9..25 };

        let seq_parser = NativeSequenceParser;
        let mut ordered = Vec::new();
        let seq = seq_parser
            .parse(
                &SequencePath::Fasta(fasta_p.clone()),
                Some(region.clone()),
                &mut ordered,
            )
            .unwrap()
            .unwrap();

        let ordered_pairs: Vec<(&str, usize)> = ordered
            .iter()
            .map(|(id, len)| (id.as_str(), *len))
            .collect();

        assert_eq!(
            ordered_pairs,
            vec![("a", 6), ("b", 8), ("c", 8), ("d", 10),]
        );

        assert_eq!(&seq.data[..], b"GCCGGTTCCGGTTGGC");

        let mp = NativeMethylationParser {
            min_threshold: MethylationValue::from_percent(0.0).into(),
        };

        let m = mp
            .parse(
                &MethylationPath::Bed(bed_p.clone()),
                Some(region),
                &ordered_pairs,
                &|_| (),
            )
            .unwrap()
            .unwrap();

        let s = &m.methylation.sites;

        println!("{s:?}");
        assert_eq!(s.len(), 3);
        assert_eq!(s[&3], mv(12.0));
        assert_eq!(s[&7], mv(20.0));
        assert_eq!(s[&9], mv(25.0));

        println!("{:?}", m.check_against_sequence(&seq));
        assert!(matches!(m.check_against_sequence(&seq), Err(3)));
    }
}
