use crate::config::{Config, PerAxis, SequenceMode};

use crate::plot::scale::TileScale;

use crate::bio::dna::domain::{SequenceBytes, SequenceRegion};
use crate::bio::methylation::SingleMethylation;
use crate::bio::{MethylationPath, SequencePath};
use crate::core::run::{Reporter, process_sequences};
use crate::error::{CoreError, InvalidInputError};
use crate::util::to_num_pretty;
use image::{ImageBuffer, Rgb, RgbImage};
use std::rc::Rc;

/// Handle a distinct pair of sequences (cross-dotplot).
pub fn handle_pairwise(
    cfg: &mut Config,
    sequences_paths: PerAxis<SequencePath>,
    methylations_paths: Option<PerAxis<MethylationPath>>,
) -> Result<(), CoreError> {
    let regions = cfg.input.region.map_ref(Clone::clone);

    let reporter = Reporter::default();

    let PerAxis {
        fst: (mut fst_sequence, fst_out),
        snd: (mut snd_sequence, snd_out),
    } = parse_sequences(cfg, sequences_paths, &reporter)?;

    {
        let word_len = *cfg.word_len;
        let fst_len = fst_sequence.data.len();
        let snd_len = snd_sequence.data.len();
        if fst_len < word_len {
            return Err(InvalidInputError(format!(
                "first sequence is shorter ({fst_len}) than word length ({word_len}), aborting..."
            ).into()).into());
        }

        if snd_len < word_len {
            return Err(InvalidInputError(format!(
                "second sequence is shorter ({snd_len}) than word length ({word_len}), aborting..."
            ).into()).into());
        }
    }

    let mut methylation: Option<PerAxis<SingleMethylation>> = {
        if let Some(PerAxis { mut fst, mut snd }) = methylations_paths
            .map(|paths| {
                parse_methylation(
                    cfg,
                    paths,
                    regions,
                    PerAxis {
                        fst: &fst_out,
                        snd: &snd_out,
                    },
                    &reporter,
                )
            })
            .transpose()?
            .flatten()
            .map(From::from)
        {
            if cfg.input.strict_methylation_base_match {
                let msg = |pos: usize| {
                    format!(
                        "methylation at 0-based position {pos} is inconsistent with the sequence. To permit this, run without --pedantic."
                    )
                };
                fst.check_against_sequence(&fst_sequence).map_err(|pos| {
                    InvalidInputError(format!("first track: {}", msg(pos)).into())
                })?;
                snd.check_against_sequence(&snd_sequence).map_err(|pos| {
                    InvalidInputError(format!("second track: {}", msg(pos)).into())
                })?;
            } else {
                let (old_len, new_len) = fst.prune_mismatches(&fst_sequence);
                if old_len > new_len {
                    println!(
                        "first methylation track: removed {} inconsistent methylation entries (mC on non-C positions)",
                        to_num_pretty(&(old_len - new_len))
                    );
                }
                let (old_len, new_len) = snd.prune_mismatches(&snd_sequence);
                if old_len > new_len {
                    println!(
                        "second methylation track: removed {} inconsistent methylation entries (mC on non-C positions)",
                        to_num_pretty(&(old_len - new_len))
                    );
                }
            }

            Some(PerAxis { fst, snd })
        } else {
            None
        }
    };

    if cfg.input.fst_sequence_mode == SequenceMode::ReverseComplement {
        println!("transforming first sequence data into reverse complement...");
        fst_sequence.rev_complement();
        if let Some(m) = methylation.as_mut() {
            m.fst.rev(0..fst_sequence.len());
        }
    }
    if cfg.input.snd_sequence_mode == SequenceMode::ReverseComplement {
        println!("transforming second sequence data into reverse complement...");
        snd_sequence.rev_complement();
        if let Some(m) = methylation.as_mut() {
            m.snd.rev(0..snd_sequence.len());
        }
    }

    process_sequences(
        &PerAxis {
            fst: Rc::new(fst_sequence),
            snd: Rc::new(snd_sequence),
        },
        methylation.map(|t| t.map(Rc::new)),
        cfg,
        &reporter,
    )
}

fn parse_methylation(
    cfg: &Config,
    paths: PerAxis<MethylationPath>,
    regions: PerAxis<Option<SequenceRegion>>,
    sequences_ids: PerAxis<&[(Rc<String>, usize)]>,
    reporter: &Reporter,
) -> Result<Option<(SingleMethylation, SingleMethylation)>, CoreError> {
    let meth_spinner = reporter.create_spinner("parsing methylation data...");

    let parse = |path, region, ordered_out| {
        cfg.methylation_parser
            .parse(path, region, ordered_out, &|str| meth_spinner.println(str))
    };

    let PerAxis {
        fst: fst_path,
        snd: snd_path,
    } = paths;
    let PerAxis {
        fst: fst_region,
        snd: snd_region,
    } = regions;

    let seq_ids = sequences_ids.map(|t| {
        t.iter()
            .map(|(id, size)| (id.as_str(), *size))
            .collect::<Vec<_>>()
    });

    let fst_methylation = parse(&fst_path, fst_region, seq_ids.fst.as_slice())?.map(|mut t| {
        t.filter_by_methylation_threshold(cfg.methylation_threshold);
        t
    });
    let snd_methylation = parse(&snd_path, snd_region, seq_ids.snd.as_slice())?.map(|mut t| {
        t.filter_by_methylation_threshold(cfg.methylation_threshold);
        t
    });

    meth_spinner.finish_with_message("finished parsing methylation data");

    if fst_methylation.is_none() || snd_methylation.is_none() {
        if fst_methylation.is_none() {
            println!("there is no methylation data for the first sequence.");
        }

        if snd_methylation.is_none() {
            println!("there is no methylation data for the second sequence.");
        }

        println!("skipping methylation data...");
    }
    Ok(fst_methylation.zip(snd_methylation))
}

#[allow(clippy::type_complexity)]
fn parse_sequences(
    cfg: &mut Config,
    sequences_paths: PerAxis<SequencePath>,
    reporter: &Reporter,
) -> Result<PerAxis<(SequenceBytes, Vec<(Rc<String>, usize)>)>, CoreError> {
    #[allow(clippy::type_complexity)]
    let parse = |path, region| -> Result<(SequenceBytes, Vec<(Rc<String>, usize)>), CoreError> {
        let mut ordered_out = Vec::new();
        Ok((
            match cfg.sequence_parser.parse(path, region, &mut ordered_out)? {
                Some(seq) if !seq.is_empty() => seq,
                _ => {
                    return Err(InvalidInputError("first sequence: no data".into()).into());
                }
            },
            ordered_out,
        ))
    };

    let (fst_region, snd_region) = (cfg.input.region.fst.take(), cfg.input.region.snd.take());
    let PerAxis {
        fst: fst_seq_path,
        snd: snd_seq_path,
    } = sequences_paths;

    let seq_spinner = reporter.create_spinner("parsing sequence data...");
    let ((fst_sequence, fst_out), (snd_sequence, snd_out)) = (
        parse(&fst_seq_path, fst_region)?,
        parse(&snd_seq_path, snd_region)?,
    );
    seq_spinner.finish_with_message("finished parsing sequence data");

    Ok(PerAxis {
        fst: (fst_sequence, fst_out),
        snd: (snd_sequence, snd_out),
    })
}

pub struct Dotplot {
    pub buf: Vec<u8>,
    pub scale: TileScale,
}

impl Dotplot {
    pub const WHITE: u8 = 255;
    pub const BLACK: u8 = 0;
    pub fn new(scale: TileScale) -> Self {
        Self {
            buf: vec![Self::WHITE; scale.out_w * scale.out_h * 3],
            scale,
        }
    }

    pub fn into_rgb_image(self) -> Option<ImageBuffer<Rgb<u8>, Vec<u8>>> {
        RgbImage::from_raw(self.scale.out_w as u32, self.scale.out_h as u32, self.buf)
    }

    pub fn enlarge_black_dots(&mut self, dot_size: u32) {
        let r = dot_size as usize;

        let w = self.scale.out_w;
        let h = self.scale.out_h;

        for y in 0..h {
            for x in 0..w {
                let i = (y * w + x) * 3;
                let is_black = self.buf[i] == 0 && self.buf[i + 1] == 0 && self.buf[i + 2] == 0;
                if !is_black {
                    continue;
                }

                for dy in y.saturating_sub(r)..=y {
                    for dx in x.saturating_sub(r)..=x {
                        let j = (dy * w + dx) * 3;
                        self.buf[j] = 0;
                        self.buf[j + 1] = 0;
                        self.buf[j + 2] = 0;
                    }
                }
            }
        }
    }

    pub fn enlarge_non_bw_dots(&mut self, dot_size: u32) {
        let r = dot_size as usize;

        let w = self.scale.out_w;
        let h = self.scale.out_h;

        for y in 0..h {
            for x in 0..w {
                let i = (y * w + x) * 3;
                let r0 = self.buf[i];
                let g0 = self.buf[i + 1];
                let b0 = self.buf[i + 2];

                // dont extend black and white
                if (r0, g0, b0) == (0, 0, 0) || (r0, g0, b0) == (255, 255, 255) {
                    continue;
                }

                for dy in y.saturating_sub(r)..=y {
                    for dx in x.saturating_sub(r)..=x {
                        let j = (dy * w + dx) * 3;

                        let a = self.buf[j];
                        let b = self.buf[j + 1];
                        let c = self.buf[j + 2];

                        // methylation colors are a layer below the black dotplot
                        if (a, b, c) == (0, 0, 0) {
                            continue;
                        }

                        if (a, b, c) == (255, 255, 255) {
                            // copied into white (as opposed to dilution), otherwise averaged
                            self.buf[j] = r0;
                            self.buf[j + 1] = g0;
                            self.buf[j + 2] = b0;
                        } else {
                            let a = self.buf[j];
                            let b = self.buf[j + 1];
                            let c = self.buf[j + 2];

                            self.buf[j] = ((u16::from(a) + u16::from(r0)) >> 1) as u8;
                            self.buf[j + 1] = ((u16::from(b) + u16::from(g0)) >> 1) as u8;
                            self.buf[j + 2] = ((u16::from(c) + u16::from(b0)) >> 1) as u8;
                        }
                    }
                }
            }
        }
    }
}
