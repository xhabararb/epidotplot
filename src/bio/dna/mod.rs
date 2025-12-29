use crate::util::to_num_pretty;
use bincode::{Decode, Encode};
use std::hash::{Hash, Hasher};
use std::mem;
use std::ops::Range;
use std::rc::Rc;

pub mod domain;
pub mod fasta;
pub mod parse;

#[derive(Debug, Eq, Clone, Encode, Decode)]
pub enum SequenceChunkId {
    Global(Range<usize>),
    Local(Rc<String>, Range<usize>),
    Id(Rc<String>),
}

impl SequenceChunkId {
    pub fn to_region_specifier_str_pretty(&self) -> String {
        match self {
            Self::Global(r) => {
                format!("{}-{}", to_num_pretty(&r.start), to_num_pretty(&r.end))
            }
            Self::Local(id, r) => {
                format!("{id}:{}-{}", to_num_pretty(&r.start), to_num_pretty(&r.end))
            }
            Self::Id(id) => id.as_str().to_owned(),
        }
    }

    pub fn to_filename_str(&self) -> String {
        match self {
            Self::Global(r) => {
                format!("{}-{}", r.start, r.end)
            }
            Self::Local(id, r) => {
                format!("{}_{}-{}", id.as_str(), r.start, r.end)
            }
            Self::Id(id) => id.to_string(),
        }
    }
}

impl PartialEq for SequenceChunkId {
    fn eq(&self, other: &Self) -> bool {
        use SequenceChunkId::{Global, Id, Local};

        if mem::discriminant(self) != mem::discriminant(other) {
            return false;
        }

        match (self, other) {
            (Global(r1), Global(r2)) => r1 == r2,
            (Local(n1, r1), Local(n2, r2)) => n1 == n2 && r1 == r2,
            (Id(n1), Id(n2)) => n1 == n2,
            _ => false,
        }
    }
}

impl Hash for SequenceChunkId {
    fn hash<H: Hasher>(&self, state: &mut H) {
        mem::discriminant(self).hash(state);
        match self {
            Self::Global(r) => r.hash(state),
            Self::Local(name, r) => {
                name.hash(state);
                r.hash(state);
            }
            Self::Id(name) => name.hash(state),
        }
    }
}
