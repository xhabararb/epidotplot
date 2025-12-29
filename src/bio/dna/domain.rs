use crate::bio::dna::SequenceChunkId;
use bincode::{Decode, Encode};
use std::fmt;
use std::ops::Range;
use std::rc::Rc;

pub trait Complement {
    fn complement(self) -> Self;
    fn complement_ref(&self) -> Self;
}

impl Complement for u8 {
    fn complement(self) -> Self {
        match self {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            _ => self,
        }
    }

    fn complement_ref(&self) -> Self {
        self.complement()
    }
}

#[derive(Encode, Decode)]
pub struct SequenceBytes {
    pub id: SequenceChunkId,
    pub data: Vec<u8>,
    pub range: Range<usize>,
    pub multi_fasta: Option<bool>,
}

impl SequenceBytes {
    pub const fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
    pub fn rev_complement(&mut self) {
        self.data.reverse();
        for b in &mut self.data {
            *b = Complement::complement(*b);
        }
    }
    #[inline]
    pub const fn len(&self) -> usize {
        self.data.len()
    }
    #[inline]
    pub fn is_cytosine_at(&self, pos: usize) -> bool {
        self.data.get(pos).is_some_and(|&t| Self::is_cytosine(t))
    }

    #[inline]
    const fn is_cytosine(base: u8) -> bool {
        base == b'C'
    }
}

#[derive(Clone)]
pub enum SequenceRegion {
    /// A plain sequence identifier (e.g. `"chr1"` or `"contig-abc"`).
    Id { id: Rc<String> },

    /// A genomic interval in absolute (global) coordinates, without an explicit reference.
    ///
    /// Example: `"100-500"`.
    /// Used when positions refer to a shared or concatenated coordinate system.
    Slice { range: Range<usize> },

    /// A genomic interval associated with a specific reference sequence.
    ///
    /// Example: `"chr1:100-500"`.
    /// The coordinates are **relative to that reference sequence** rather than global.
    IdSlice { id: Rc<String>, range: Range<usize> },
}

impl PartialEq for SequenceRegion {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Id { id }, Self::Id { id: other_id }) => id.as_str() == other_id.as_str(),
            (Self::Slice { range }, Self::Slice { range: other_range }) => *range == *other_range,
            (
                Self::IdSlice { id, range },
                Self::IdSlice {
                    id: other_id,
                    range: other_range,
                },
            ) => id.as_str() == other_id.as_str() && *range == *other_range,
            _ => false,
        }
    }
}

impl SequenceRegion {
    /// Parses a textual region specification into a [`SequenceRegion`].
    ///
    /// # Rules
    /// - `<id>:<start>-<end>` → [`IdSlice`] if the interval is valid (`start ≤ end`)
    /// - `<start>-<end>` → [`Slice`] if valid (same as above)
    /// - otherwise → [`Id`]
    ///
    /// Only the **rightmost `:`** or **`-`** is considered a delimiter,
    /// and they take effect **only when the interval parses successfully**;
    /// otherwise, the string is treated as a plain identifier.
    ///
    /// # Returns
    /// `None` if the identifier part is empty (e.g. `":100-200"`, `":"`, or an empty string).
    /// In all other cases, a valid [`SequenceFilter`] is returned.
    pub fn parse(s: &str) -> Option<Self> {
        fn try_parse_interval(s: &str) -> Option<Range<usize>> {
            let (l, h) = s.rsplit_once('-')?;
            let (l, h) = (l.parse::<usize>().ok()?, h.parse::<usize>().ok()?);
            (l < h).then_some(l..h)
        }

        let get_id = || SequenceRegion::Id {
            id: Rc::new(s.into()),
        };

        if let Some((id, tail)) = s.rsplit_once(':') {
            return (!id.is_empty()).then(|| {
                try_parse_interval(tail).map_or_else(get_id, |range| Self::IdSlice {
                    id: Rc::new(id.into()),
                    range,
                })
            });
        }

        (!s.is_empty())
            .then(|| try_parse_interval(s).map_or_else(get_id, |range| Self::Slice { range }))
    }

    #[inline]
    pub fn into_parts(self) -> (Option<Rc<String>>, Option<Range<usize>>) {
        match self {
            Self::Id { id } => (Some(id), None),
            Self::Slice { range } => (None, Some(range)),
            Self::IdSlice { id, range } => (Some(id), Some(range)),
        }
    }
}

impl From<SequenceChunkId> for SequenceRegion {
    fn from(value: SequenceChunkId) -> Self {
        match value {
            SequenceChunkId::Global(range) => Self::Slice { range },
            SequenceChunkId::Local(id, range) => Self::IdSlice { id, range },
            SequenceChunkId::Id(id) => Self::Id { id },
        }
    }
}

impl fmt::Debug for SequenceRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Id { id: name } => write!(f, "{name}"),
            Self::Slice {
                range: Range { start, end },
            } => write!(f, "{start}-{end}"),
            Self::IdSlice {
                id: name,
                range: Range { start, end, .. },
            } => {
                write!(f, "{name}:{start}-{end}")
            }
        }
    }
}

#[cfg(test)]
mod sequence_region_tests {
    use super::*;

    #[test]
    fn id_slice_ok() {
        assert_eq!(
            SequenceRegion::parse("chr1:100-500"),
            Some(SequenceRegion::IdSlice {
                id: Rc::new("chr1".into()),
                range: 100..500
            })
        );
    }

    #[test]
    fn id_slice_decreasing() {
        assert_eq!(
            SequenceRegion::parse("chr1:500-100"),
            Some(SequenceRegion::Id {
                id: Rc::new("chr1:500-100".into())
            })
        );
    }

    #[test]
    fn id_slice_empty() {
        assert_eq!(Rc::new("ab"), Rc::new("ab".into()));
        assert_eq!(
            SequenceRegion::parse("chr1:12-12"),
            Some(SequenceRegion::Id {
                id: Rc::new("chr1:12-12".into())
            })
        );
    }

    #[test]
    fn slice_ok() {
        assert_eq!(
            SequenceRegion::parse("5-10"),
            Some(SequenceRegion::Slice { range: 5..10 })
        );
    }

    #[test]
    fn id_ok() {
        assert_eq!(
            SequenceRegion::parse("chrX"),
            Some(SequenceRegion::Id {
                id: Rc::new("chrX".into())
            })
        );
    }

    #[test]
    fn empty_id_rejected() {
        for s in ["", ":10-20", ":", ":-"] {
            assert_eq!(SequenceRegion::parse(s), None);
        }
    }

    #[test]
    fn bad_intervals() {
        for s in ["-", "-10", "10-", "-10-20", "10-20-30", "-10--10", "10--10"] {
            assert_eq!(
                SequenceRegion::parse(s),
                Some(SequenceRegion::Id {
                    id: Rc::new(s.into())
                })
            );
        }

        for s in [
            "id:-",
            "id:-10",
            "id:10-",
            "id:-10-20",
            "id:10-20-30",
            "id:-10--10",
            "id:10--10",
        ] {
            assert_eq!(
                SequenceRegion::parse(s),
                Some(SequenceRegion::Id {
                    id: Rc::new(s.into())
                })
            );
        }
    }

    #[test]
    fn id_slice_invalid_interval_falls_back() {
        assert_eq!(
            SequenceRegion::parse("chr1:abc-def"),
            Some(SequenceRegion::Id {
                id: Rc::new("chr1:abc-def".into())
            })
        );
        assert_eq!(
            SequenceRegion::parse("chr1:10-20-30"),
            Some(SequenceRegion::Id {
                id: Rc::new("chr1:10-20-30".into())
            })
        );
    }

    #[test]
    fn rightmost_colon_and_hyphen_used() {
        assert_eq!(
            SequenceRegion::parse("ref:part:10-20"),
            Some(SequenceRegion::IdSlice {
                id: Rc::new("ref:part".into()),
                range: 10..20
            })
        );
        assert_eq!(
            SequenceRegion::parse("ref:10-30:10-20"),
            Some(SequenceRegion::IdSlice {
                id: Rc::new("ref:10-30".into()),
                range: 10..20
            })
        );
    }
}
