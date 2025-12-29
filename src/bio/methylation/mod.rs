use crate::bio::dna::SequenceChunkId;
use crate::bio::dna::domain::SequenceBytes;
use crate::bio::methylation::domain::MethylationValue;
use crate::config::MethylationThreshold;
use crate::util::to_num_pretty;
use itertools::Itertools;
use std::collections::HashMap;
use std::fmt;
use std::ops::Range;

pub mod bed;
pub mod bedgraph;
pub mod domain;
pub mod parse;
pub mod plot;

pub struct SingleMethylation {
    pub id: SequenceChunkId,
    pub methylation: SequenceMethylation,
}

impl SingleMethylation {
    /// Reverses methylation site coordinates within the given interval.
    /// The caller must supply an interval that fully covers all sites.
    ///
    /// # Panics
    /// Panics if any methylation site lies outside `coords`. This indicates a
    /// violated post-processing invariant and is considered unreachable in
    /// correct usage.
    pub fn rev(&mut self, coords: Range<usize>) {
        let Range { start, end } = coords;
        let mut rev = HashMap::with_capacity(self.methylation.sites.len());

        for (&pos, val) in &self.methylation.sites {
            assert!(
                coords.contains(&pos),
                "methylation sites are out of specified bounds"
            );
            rev.insert(start + (end - 1 - pos), *val);
        }

        self.methylation.sites = rev;
    }

    pub fn check_against_sequence(&self, sequence: &SequenceBytes) -> Result<(), usize> {
        self.methylation
            .sites
            .keys()
            .sorted()
            .find(|&&pos| !sequence.is_cytosine_at(pos))
            .map_or(Ok(()), |&err_pos| Err(err_pos))
    }

    // old, new
    pub fn prune_mismatches(&mut self, sequence: &SequenceBytes) -> (usize, usize) {
        let old_len = self.methylation.sites.len();
        self.methylation
            .sites
            .retain(|&pos, _| sequence.is_cytosine_at(pos));
        let new_len = self.methylation.sites.len();
        (old_len, new_len)
    }

    pub fn correct_percents_into_fractions(&mut self) {
        self.methylation.correct_percents_into_fractions();
    }
    pub fn filter_by_methylation_threshold(&mut self, threshold: MethylationThreshold) {
        self.methylation
            .sites
            .retain(|_, meth| meth.as_percent() >= threshold.as_percent());
    }
}

pub enum MethylationSites {
    Single(SingleMethylation),
    Multiple(HashMap<Box<str>, SequenceMethylation>),
}

impl MethylationSites {
    pub fn collapse(self, ordered_sequences: &[(&str, usize)]) -> Option<SingleMethylation> {
        if ordered_sequences.is_empty() {
            return None;
        }

        Some(match self {
            Self::Single(single) => single,
            Self::Multiple(table) => {
                let (sites, _max_index) = Self::collapse_global(table, ordered_sequences)?;
                let end = ordered_sequences.iter().map(|t| t.1).sum::<usize>();
                SingleMethylation {
                    id: SequenceChunkId::Global(0..end),
                    methylation: SequenceMethylation { sites },
                }
            }
        })
    }

    pub fn collapse_global(
        mut chroms: HashMap<Box<str>, SequenceMethylation>,
        ordered_sequences: &[(&str, usize)],
    ) -> Option<(HashMap<usize, MethylationValue>, usize)> {
        if ordered_sequences.is_empty() {
            return None;
        }

        let mut global = HashMap::new();
        let mut offset = 0usize;
        let mut max_index = 0usize;

        for (id, len) in ordered_sequences {
            if let Some(seq_meth) = chroms.remove(*id) {
                for (&pos, val) in &seq_meth.sites {
                    let shifted = offset + pos;
                    global.insert(shifted, *val);
                    max_index = max_index.max(shifted);
                }
            }
            offset += len;
        }

        Some((global, max_index))
    }
}

/// Methylation data for a single chromosome.
///
/// Stores site-level information, precomputed site positions,
/// and the maximum genomic index present.
#[derive(Default)]
//fixme better name like segment methylation
pub struct SequenceMethylation {
    pub sites: HashMap<usize, MethylationValue>,
}

impl SequenceMethylation {
    pub fn sites(&self) -> Vec<(usize, MethylationValue)> {
        self.sites.iter().map(|s| (*s.0, *s.1)).collect()
    }
    pub fn correct_percents_into_fractions(&mut self) {
        // e.g., fraction 0.1 read as 0.1% instead of 10%, which is 100x less
        for val in &mut self.sites.values_mut() {
            *val = MethylationValue::from_percent(val.as_percent() * 100.0);
        }
    }
}

impl fmt::Debug for SequenceMethylation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ChromosomeMethylation")
            .field("number_of_sites", &to_num_pretty(&self.sites.len()))
            .finish()
    }
}
