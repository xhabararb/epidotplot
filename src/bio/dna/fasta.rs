use crate::bio::dna::SequenceChunkId;
use crate::bio::dna::domain::SequenceBytes;
use crate::bio::error::IoError;
use crate::error::{CoreError, InvalidInputError, MalformedInputError};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;
use std::rc::Rc;

fn open_lines(path: &Path) -> Result<impl Iterator<Item = Result<String, CoreError>>, CoreError> {
    let file = File::open(path).map_err(|err| IoError {
        err,
        message: Some(format!("at {}", path.display())),
    })?;
    Ok(BufReader::new(file).lines().map(|l| {
        l.map_err(|err| {
            IoError {
                err,
                message: Some("failed to read FASTA line".into()),
            }
            .into()
        })
    }))
}

#[inline]
fn validate_dna_symbol(b: u8) -> Result<u8, MalformedInputError> {
    let b = b.to_ascii_uppercase();
    match b {
        b'A' | b'C' | b'G' | b'T' | b'N' => Ok(b),
        b'-' => Err(MalformedInputError("unexpected gap symbol '-'".into())),
        _ => Err(MalformedInputError(
            format!(
                "invalid DNA sequence code '{}' (0x{b:02X}); allowed: A, C, G, T, N (uppercase/lowercase)",
                b as char
            )
            .into(),
        )),
    }
}

#[inline]
pub fn contains_n(kmer: &[u8]) -> bool {
    kmer.contains(&b'N')
}

#[inline]
fn parse_id_line(line: &str) -> Result<String, CoreError> {
    if !line.starts_with('>') {
        return Err(MalformedInputError("FASTA header must start with '>'".into()).into());
    }
    line[1..]
        .split_whitespace()
        .next()
        .filter(|s| !s.is_empty())
        .map(str::to_owned)
        .ok_or(MalformedInputError("missing identifier".into()).into())
}

pub(in crate::bio::dna) fn read_all<P: AsRef<Path>>(
    path: P,
    id_order_out: &mut Vec<(Rc<String>, usize)>,
) -> Result<Option<SequenceBytes>, CoreError> {
    let mut data = Vec::with_capacity(1 << 20);
    let mut id: Option<String> = None;
    let mut record_start = 0usize;

    for line in open_lines(path.as_ref())? {
        let line = line?;
        // skip comments
        if line.starts_with(';') {
            continue;
        }

        if line.starts_with('>') {
            if let Some(i) = id.take() {
                id_order_out.push((Rc::new(i), data.len() - record_start));
                record_start = data.len();
            }
            id = Some(parse_id_line(&line)?);
            continue;
        }

        for b in line.bytes() {
            data.push(validate_dna_symbol(b)?);
        }
    }

    let Some(id) = id else {
        return Err(MalformedInputError("no id starting with '>' found in FASTA".into()).into());
    };
    let id = Rc::new(id);

    id_order_out.push((Rc::clone(&id), data.len() - record_start));

    if data.is_empty() {
        Ok(None)
    } else {
        let range = 0..data.len();
        Ok(Some(SequenceBytes {
            id: if id_order_out.len() > 1 {
                SequenceChunkId::Global(0..data.len())
            } else {
                SequenceChunkId::Id(id)
            },
            data,
            range,
            multi_fasta: Some(id_order_out.len() > 1),
        }))
    }
}
pub(in crate::bio::dna) fn read_id<P: AsRef<Path>>(
    path: P,
    target: &Rc<String>,
    local: Option<Range<usize>>,
    ordered_out: &mut Vec<(Rc<String>, usize)>,
) -> Result<Option<SequenceBytes>, CoreError> {
    let mut seq = Vec::with_capacity(1 << 14);
    let mut in_target = false;
    let mut count = 0usize;
    let mut current_id: Option<String> = None;

    for line in open_lines(path.as_ref())? {
        let line = line?;
        // skip comments
        if line.starts_with(';') {
            continue;
        }
        if line.starts_with('>') {
            if let Some(id) = current_id.take() {
                ordered_out.push((Rc::new(id), count));
                count = 0;
            }

            let name = parse_id_line(&line)?;
            if in_target {
                break;
            }
            in_target = name.as_str() == target.as_str();
            current_id = Some(name);
            continue;
        }

        let line = line.trim_end();

        for b in line.bytes() {
            validate_dna_symbol(b)?;
        }
        if in_target {
            seq.extend(line.bytes());
        }
        count += line.len();
    }

    if seq.is_empty() {
        return Ok(None);
    }

    let (range, id_variant) = if let Some(r) = local {
        if r.end > seq.len() {
            return Ok(None);
        }
        let len = r.end - r.start;
        seq.copy_within(r.clone(), 0);
        seq.truncate(len);
        (0..len, SequenceChunkId::Local(Rc::clone(target), r))
    } else {
        (0..seq.len(), SequenceChunkId::Id(Rc::clone(target)))
    };

    if ordered_out
        .last()
        .as_ref()
        .is_none_or(|t| t.0.as_str() != target.as_str())
    {
        ordered_out.push((Rc::clone(target), seq.len()));
    }

    Ok(Some(SequenceBytes {
        id: id_variant,
        data: seq,
        range,
        multi_fasta: Some(ordered_out.len() > 1),
    }))
}

pub(in crate::bio::dna) fn read_global<P: AsRef<Path>>(
    path: P,
    global: Range<usize>,
    ordered_out: &mut Vec<(Rc<String>, usize)>,
) -> Result<Option<SequenceBytes>, CoreError> {
    if global.is_empty() {
        return Err(InvalidInputError(
            format!("invalid range: [{}, {}]", global.start, global.end).into(),
        )
        .into());
    }

    let mut data = Vec::with_capacity(global.len());
    let mut id: Option<String> = None;
    let mut offset = 0usize;
    let mut record_len = 0usize;

    for line in open_lines(path.as_ref())? {
        let line = line?;

        // skip comments
        if line.starts_with(';') {
            continue;
        }

        if line.starts_with('>') {
            if let Some(i) = id.take() {
                //if offset < global.end && offset + record_len > global.start {
                ordered_out.push((Rc::new(i), record_len));
                //}
                offset += record_len;
                record_len = 0;
                if offset >= global.end {
                    break;
                }
            }
            id = Some(parse_id_line(&line)?);
            continue;
        }

        for b in line.bytes() {
            let base = validate_dna_symbol(b)?;
            let abs = offset + record_len;
            if abs >= global.start && abs < global.end {
                data.push(base);
            }
            record_len += 1;
        }
    }

    if let Some(i) = id
    /*&& offset < global.end
    && offset + record_len > global.start*/
    {
        ordered_out.push((Rc::new(i), record_len));
    }

    if data.is_empty() {
        Ok(None)
    } else {
        Ok(Some(SequenceBytes {
            id: SequenceChunkId::Global(global.clone()),
            data,
            range: global,
            multi_fasta: Some(ordered_out.len() > 1),
        }))
    }
}

#[cfg(test)]
mod fasta_tests {
    use super::*;
    use std::{fs::File, io::Write, path::PathBuf};

    fn write_tmp_fasta(name: &str, content: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!("{}_test.fasta", name));
        let mut f = File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn global_across_irregular_multi_contigs() {
        let a = "ACGTACGTACCGT"; // 0..=12
        let b = "TTTTGGTTTTTTT"; // 13..=25
        let c = "GGG"; // 26..=28
        let content = format!(">a\n{a}\n>b\n{b}\n>c\n{c}\n");
        let path = write_tmp_fasta("multi_irregular", &content);
        let mut ordered = Vec::new();

        let out = read_global(&path, 10..25, &mut ordered).unwrap().unwrap();

        println!("{}", String::from_utf8_lossy(&out.data));
        println!("{:?}", ordered);

        assert_eq!(out.data, b"CGTTTTTGGTTTTTT");
        assert_eq!(ordered.len(), 2);
        assert_eq!(ordered[0].0.as_str(), "a");
        assert_eq!(ordered[0].1, 13);
        assert_eq!(ordered[1].0.as_str(), "b");
        assert_eq!(ordered[1].1, 13);
    }
    #[test]
    fn global_header_without_identifier_errors() {
        let path = write_tmp_fasta(
            "global_empty_id",
            ">a\nAAA\n>\nCCC\n", // the second header is empty
        );
        let mut ordered = Vec::new();

        let err = read_global(&path, 0..10, &mut ordered);
        let Err(err) = err else {
            panic!();
        };
        assert!(
            matches!(err, CoreError::MalformedInput(err) if format!("{err}").contains("missing identifier"))
        );
    }

    #[test]
    fn global_handles_missing_final_newline() {
        let path = write_tmp_fasta(
            "global_no_eof_newline",
            ">a\nAAA\n>b\nCCC", // no \n at the end
        );
        let mut ordered = Vec::new();

        let out = read_global(&path, 1..5, &mut ordered).unwrap().unwrap();

        assert_eq!(out.data, b"AACC");
        assert_eq!(ordered.len(), 2);
        assert_eq!(ordered[0].0.as_str(), "a");
        assert_eq!(ordered[0].1, 3);
        assert_eq!(ordered[1].0.as_str(), "b");
        assert_eq!(ordered[1].1, 3);
    }

    #[test]
    fn global_tolerates_blank_lines_everywhere() {
        let path = write_tmp_fasta(
            "global_blank_lines",
            "\n\n>a\nAAA\n\n>b\nCCC\n\n>c\nGGG\n\n",
        );
        let mut ordered = Vec::new();

        let out = read_global(&path, 1..8, &mut ordered).unwrap().unwrap();

        assert_eq!(out.data, b"AACCCGG");
        assert_eq!(ordered.len(), 3);
        assert_eq!(ordered[0].0.as_str(), "a");
        assert_eq!(ordered[1].0.as_str(), "b");
        assert_eq!(ordered[2].0.as_str(), "c");
    }

    #[test]
    fn global_partial_span_across_records() {
        let path = write_tmp_fasta("global_cut", ">a\nAAA\n>b\nCCC\n");
        let mut ordered = Vec::new();
        let out = read_global(&path, 2..5, &mut ordered).unwrap().unwrap();
        assert_eq!(out.data, b"AACCC"[1..4]);
        assert_eq!(ordered.len(), 2);
        assert_eq!(ordered[0].0.as_str(), "a");
        assert_eq!(ordered[0].1, 3);
        assert_eq!(ordered[1].0.as_str(), "b");
        assert_eq!(ordered[1].1, 3);
    }

    #[test]
    fn global_inside_one_contig() {
        let path = write_tmp_fasta("inside", ">a\nACGTACGT\n>b\nTTTT\n");
        let mut ordered = Vec::new();
        let out = read_global(&path, 2..5, &mut ordered).unwrap().unwrap();
        assert_eq!(out.data, b"GTA");
        assert_eq!(ordered.len(), 1);
        assert_eq!(ordered[0].0.as_str(), "a");
        assert_eq!(ordered[0].1, 8);
    }

    #[test]
    fn global_exact_boundaries_touching_excluded() {
        let path = write_tmp_fasta("exact_bound", ">a\nAAA\n>b\nTTT\n>c\nGGG\n");
        let mut ordered = Vec::new();
        let out = read_global(&path, 3..6, &mut ordered).unwrap().unwrap();
        assert_eq!(out.data, b"TTT");
        assert_eq!(ordered.len(), 2);
        assert_eq!(ordered[0].0.as_str(), "a");
        assert_eq!(ordered[0].1, 3);

        assert_eq!(ordered[1].0.as_str(), "b");
        assert_eq!(ordered[1].1, 3);
    }

    #[test]
    fn global_span_multiple_contigs() {
        let path = write_tmp_fasta("span_multi", ">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n");
        let mut ordered = Vec::new();
        let out = read_global(&path, 2..9, &mut ordered).unwrap().unwrap();
        assert_eq!(out.data, b"AAAACCCCGGGG"[2..9]);
        assert_eq!(ordered.len(), 3);
        assert_eq!(ordered[0].1, 4);
        assert_eq!(ordered[1].1, 4);
        assert_eq!(ordered[2].1, 4);
    }

    #[test]
    fn bruh() {
        let path = write_tmp_fasta("span_multi", ">a\nAAA\nA\n>b\nCC\nCC\n>c\nG\nGGG\n");
        let mut ordered = Vec::new();
        let out = read_global(&path, 2..9, &mut ordered).unwrap().unwrap();
        assert_eq!(out.data, b"AAAACCCCGGGG"[2..9]);
    }

    #[test]
    fn global_no_overlap_before_all() {
        let path = write_tmp_fasta("before_all", ">a\nAAA\n>b\nTTT\n");
        let mut ordered = Vec::new();
        let out = read_global(&path, 0..0, &mut ordered);
        assert!(matches!(out, Err(CoreError::InvalidInput(_))));
    }

    #[test]
    fn global_no_overlap_after_all() {
        let path = write_tmp_fasta("after_all", ">a\nAA\n>b\nTT\n");
        let mut ordered = Vec::new();
        let out = read_global(&path, 10..12, &mut ordered).unwrap();
        assert!(out.is_none());
        println!("{:?}", ordered);
        assert_eq!(ordered.len(), 2);
    }

    #[test]
    fn blank_sequence_still_counted_if_intersects() {
        let path = write_tmp_fasta("blank", ">a\nAAA\n>b\n>c\nCCC\n");
        let mut ordered = Vec::new();
        let out = read_global(&path, 2..7, &mut ordered).unwrap().unwrap();
        assert_eq!(out.data, b"AAACCC"[2..6]);
        assert_eq!(ordered.len(), 3);
        assert_eq!(ordered[0].0.as_str(), "a");
        assert_eq!(ordered[1].0.as_str(), "b");
        assert_eq!(ordered[2].0.as_str(), "c");
    }

    #[test]
    fn empty_file_returns_none() {
        let path = write_tmp_fasta("empty", "");
        let out = read_all(&path, &mut Vec::new());
        assert!(matches!(out, Err(CoreError::MalformedInput(_))));
    }

    #[test]
    fn id_stores_counts_for_all_contigs_before_target() {
        let path = write_tmp_fasta("id_counts_all", ">a\nAAA\n>b\nTTATTTGG\n>c\nGGG\n");
        let mut ordered = Vec::new();

        let out = read_id(&path, &Rc::new("b".to_string()), None, &mut ordered)
            .unwrap()
            .unwrap();

        assert_eq!(out.data, b"TTATTTGG");
        assert_eq!(ordered.len(), 2);
        assert_eq!(ordered[0], (Rc::new("a".into()), 3));
        assert_eq!(ordered[1], (Rc::new("b".into()), 8));
    }

    #[test]
    fn local_reads_subrange_and_preserves_all_counts() {
        let path = write_tmp_fasta("local_counts_all", ">x\nAAA\n>y\nTTATTTGG\n>z\nGG\n");
        let mut ordered = Vec::new();

        let out = read_id(&path, &Rc::new("y".to_string()), Some(2..5), &mut ordered)
            .unwrap()
            .unwrap();

        assert_eq!(out.data, b"ATT");
        assert_eq!(ordered.len(), 2);
        assert_eq!(ordered[0], (Rc::new("x".into()), 3));
        assert_eq!(ordered[1], (Rc::new("y".into()), 8));
    }

    #[test]
    fn global_counts_and_order_across_multiple_contigs() {
        let path = write_tmp_fasta("global_counts_all", ">a\nAAA\n>b\nTTTT\n>c\nGGTAA\n");
        let mut ordered = Vec::new();

        let out = read_global(&path, 2..9, &mut ordered).unwrap().unwrap();

        assert_eq!(out.data, b"AAATTTTGGTAA"[2..9]);
        assert_eq!(ordered.len(), 3);
        assert_eq!(ordered[0], (Rc::new("a".into()), 3));
        assert_eq!(ordered[1], (Rc::new("b".into()), 4));
        assert_eq!(ordered[2], (Rc::new("c".into()), 5));
    }

    #[test]
    fn all_reads_and_counts_every_contig_in_order() {
        let path = write_tmp_fasta("all_counts_all", ">a\nAAA\n>b\nTTTT\n>c\nGGG\n");
        let mut ordered = Vec::new();

        let out = read_all(&path, &mut ordered).unwrap().unwrap();

        assert_eq!(out.data, b"AAATTTTGGG");
        assert_eq!(ordered.len(), 3);
        assert_eq!(ordered[0], (Rc::new("a".into()), 3));
        assert_eq!(ordered[1], (Rc::new("b".into()), 4));
        assert_eq!(ordered[2], (Rc::new("c".into()), 3));
    }
}

/*#[cfg(test)]
mod tests {
    use super::*;
    use std::{io::Write, path::PathBuf};

    fn write_tmp_fasta(name: &str, content: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!("{}_test.fasta", name));
        let mut file = File::create(&path).expect("failed to create tmp fasta");
        file.write_all(content.as_bytes()).expect("failed to write tmp fasta");
        path
    }

    #[test]
    fn empty_fasta() {
        let path = write_tmp_fasta("empty", "");
        let records = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect("failed to parse empty fasta");
        assert!(records.is_empty());
    }

    #[test]
    fn one_sequence() {
        let path = write_tmp_fasta("one", ">seq1\nACTG\n");
        let records = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect("failed to parse one-sequence fasta");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id.as_str(), "seq1");
        assert_eq!(records[0].seq, b"ACTG");
    }

    #[test]
    fn two_sequences() {
        let path = write_tmp_fasta("two", ">seq1\nACTG\n>seq2\nTTAA\n");
        let records = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect("failed to parse two-sequence fasta");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id.as_str(), "seq1");
        assert_eq!(records[0].seq, b"ACTG");
        assert_eq!(records[1].id.as_str(), "seq2");
        assert_eq!(records[1].seq, b"TTAA");
    }

    #[test]
    fn header_with_multiple_words() {
        let path = write_tmp_fasta("multi", ">seq1 some description here\nACTG\n");
        let records = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect("failed to parse fasta with multi-word header");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id.as_str(), "seq1");
        assert_eq!(records[0].seq, b"ACTG");
    }

    #[test]
    fn header_without_id() {
        let path = write_tmp_fasta("noid", ">\nACTG\n");
        let err = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect_err("expected malformed input error");
        matches!(err, CoreError::MalformedInput(_));
    }

    #[test]
    fn header_without_sequence() {
        let path = write_tmp_fasta("emptyseq", ">seq1\n>seq2\nACTG\n");
        let records = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect("failed to parse fasta with empty sequence");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id.as_str(), "seq1");
        assert!(records[0].seq.is_empty());
        assert_eq!(records[1].id.as_str(), "seq2");
        assert_eq!(records[1].seq, b"ACTG");
    }

    #[test]
    fn lowercase_is_uppercased() {
        let path = write_tmp_fasta("lowercase", ">seq1\nactg\n");
        let records = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect("failed to parse fasta with lowercase");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq, b"ACTG");
    }

    #[test]
    fn sequence_with_n() {
        let path = write_tmp_fasta("with_n", ">seq1\nACNGTN\n");
        let records = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect("failed to parse fasta with N bases");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq, b"ACNGTN");
    }

    #[test]
    fn sequence_with_gap_rejected() {
        let path = write_tmp_fasta("with_gap", ">seq1\nAC-TG\n");
        let err = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect_err("expected malformed input error (gap not allowed)");
        matches!(err, CoreError::MalformedInput(_));
    }

    #[test]
    fn sequence_with_gap_accepted() {
        let mut ctx = SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated);
        ctx.accept_gap = true;
        let path = write_tmp_fasta("with_gap_ok", ">seq1\nAC-TG\n");
        let records =
            read_fasta(&ctx, &path, None).expect("failed to parse fasta with accepted gap");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq, b"AC-TG");
    }

    #[test]
    fn sequence_with_invalid_char() {
        let path = write_tmp_fasta("invalid", ">seq1\nACXG\n");
        let err = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect_err("expected malformed input error");
        matches!(err, CoreError::MalformedInput(_));
    }

    #[test]
    fn sequence_with_whitespace_rejected() {
        let path = write_tmp_fasta("whitespace", ">seq1\nAC TG\n");
        let err = read_fasta(
            &SequenceParserContext::default_with_coord_mode(CoordinateMode::Isolated),
            &path,
            None,
        )
        .expect_err("expected malformed input error (whitespace not allowed)");
        matches!(err, CoreError::MalformedInput(_));
    }
}
*/
