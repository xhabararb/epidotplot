use crate::bio::dna::domain::{SequenceBytes, SequenceRegion};
use crate::bio::methylation::SingleMethylation;
use crate::bio::{MethylationPath, SequencePath};
use crate::config::{Config, PerAxis, SequenceMode};
use crate::core::run::{Reporter, process_sequences};
use crate::error::{CoreError, InvalidInputError};
use crate::util::to_num_pretty;
use std::rc::Rc;

pub fn handle_single(
    cfg: &mut Config,
    sequence_path: &SequencePath,
    methylation_path: Option<&MethylationPath>,
) -> Result<(), CoreError> {
    let region = cfg.input.region.fst.clone();

    let reporter = Reporter::default();

    let (mut sequence, out) = parse_sequence(cfg, sequence_path, &reporter)?;

    {
        let word_len = *cfg.word_len;
        let len = sequence.data.len();
        if len < word_len {
            return Err(InvalidInputError(
                format!("sequence is shorter ({len}) than word length ({word_len}), aborting...")
                    .into(),
            )
            .into());
        }
    }

    let mut methylation = if let Some(mut methylation) = methylation_path
        .as_ref()
        .map(|path| parse_methylation(cfg, path, region, &out, &reporter))
        .transpose()?
        .flatten()
    {
        if cfg.input.strict_methylation_base_match {
            methylation.check_against_sequence(&sequence).map_err(|pos| {
                InvalidInputError(format!(
                    "methylation at 0-based position {pos} is inconsistent with the sequence. To permit this, run without --pedantic.").into()
                )
            })?;
        } else {
            let (old_len, new_len) = methylation.prune_mismatches(&sequence);
            if old_len > new_len {
                println!(
                    "removed {} inconsistent methylation entries (mC on non-C positions)",
                    to_num_pretty(&(old_len - new_len))
                );
            }
        }
        Some(methylation)
    } else {
        None
    };

    match (cfg.input.fst_sequence_mode, cfg.input.snd_sequence_mode) {
        (SequenceMode::ReverseComplement, SequenceMode::ReverseComplement) => {
            println!("transforming sequence data into reverse complement...");
            sequence.rev_complement();
            if let Some(meth) = methylation.as_mut() {
                meth.rev(0..sequence.len());
            }
        }
        (SequenceMode::Original, SequenceMode::Original) => {}
        _ => return Err(CoreError::Other {
            err: "a single-sequence input dotplot cannot reverse only one axis; reverse-complementing must be applied to both axes or neither".into(),
            msg: None })
    }

    let seq_rc = Rc::new(sequence);
    let meth_rc = methylation.map(|m| {
        let m_rc = Rc::new(m);
        PerAxis {
            fst: Rc::clone(&m_rc),
            snd: m_rc,
        }
    });
    process_sequences(
        &PerAxis {
            fst: Rc::clone(&seq_rc),
            snd: seq_rc,
        },
        meth_rc,
        cfg,
        &reporter,
    )
}

#[allow(clippy::type_complexity)]
fn parse_sequence(
    cfg: &Config,
    path: &SequencePath,
    reporter: &Reporter,
) -> Result<(SequenceBytes, Vec<(Rc<String>, usize)>), CoreError> {
    let region = cfg.input.region.fst.clone();
    let mut out = Vec::new();

    let seq_spinner = reporter.create_spinner("parsing sequence data...");

    let seq = match cfg.sequence_parser.parse(path, region, &mut out)? {
        Some(seq) if !seq.is_empty() => seq,
        _ => {
            return Err(InvalidInputError("sequence: no data".into()).into());
        }
    };

    seq_spinner.finish_with_message("finished parsing sequence data");

    Ok((seq, out))
}

fn parse_methylation(
    cfg: &Config,
    path: &MethylationPath,
    region: Option<SequenceRegion>,
    sequences_ids: &[(Rc<String>, usize)],
    reporter: &Reporter,
) -> Result<Option<SingleMethylation>, CoreError> {
    let meth_spinner = reporter.create_spinner("parsing methylation data...");

    let seq_ids = sequences_ids
        .iter()
        .map(|(id, size)| (id.as_str(), *size))
        .collect::<Vec<_>>();

    let methylation = cfg
        .methylation_parser
        .parse(path, region, seq_ids.as_slice(), &|str| {
            meth_spinner.println(str);
        })?
        .map(|mut t| {
            t.filter_by_methylation_threshold(cfg.methylation_threshold);
            t
        });

    meth_spinner.finish_with_message("finished parsing methylation data");

    if methylation.is_none() {
        println!("there is no methylation data for the sequence, skipping...");
    }

    Ok(methylation)
}
