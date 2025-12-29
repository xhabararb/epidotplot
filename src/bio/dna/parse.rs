use crate::bio::SequencePath;
use crate::bio::dna::domain::{SequenceBytes, SequenceRegion};
use crate::bio::dna::fasta;
use crate::error::CoreError;
use std::fmt;
use std::path::Path;
use std::rc::Rc;

pub struct NativeSequenceParser;

impl fmt::Debug for NativeSequenceParser {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "native sequence parser for FASTA files")
    }
}

impl SequenceParser for NativeSequenceParser {
    fn parse_fasta(
        &self,
        path: &Path,
        region: Option<SequenceRegion>,
        ordered_out: &mut Vec<(Rc<String>, usize)>,
    ) -> Result<Option<SequenceBytes>, CoreError> {
        match region {
            None => fasta::read_all(path, ordered_out),
            Some(SequenceRegion::Id { id }) => fasta::read_id(path, &id, None, ordered_out),
            Some(SequenceRegion::IdSlice { id, range }) => {
                fasta::read_id(path, &id, Some(range), ordered_out)
            }
            Some(SequenceRegion::Slice { range }) => fasta::read_global(path, range, ordered_out),
        }
    }
}

pub trait SequenceParser: fmt::Debug {
    fn parse_fasta(
        &self,
        path: &Path,
        region: Option<SequenceRegion>,
        ordered_out: &mut Vec<(Rc<String>, usize)>,
    ) -> Result<Option<SequenceBytes>, CoreError>;
    fn parse(
        &self,
        path: &SequencePath,
        region: Option<SequenceRegion>,
        ordered_out: &mut Vec<(Rc<String>, usize)>,
    ) -> Result<Option<SequenceBytes>, CoreError> {
        match path {
            SequencePath::Fasta(path) => self.parse_fasta(path, region, ordered_out),
        }
    }
}
