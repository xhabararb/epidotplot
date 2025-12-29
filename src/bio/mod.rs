use std::fmt;
use std::path::{Path, PathBuf};

pub mod dna;
pub mod error;
pub mod methylation;
mod tests;

#[derive(Debug, Clone)]
pub enum SequencePath {
    Fasta(PathBuf),
}

#[derive(Debug, Clone)]
pub enum MethylationPath {
    Bed(PathBuf),
    BedGraph(PathBuf),
}

impl AsRef<Path> for SequencePath {
    fn as_ref(&self) -> &Path {
        match self {
            Self::Fasta(path) => path.as_ref(),
        }
    }
}

impl PathLike for SequencePath {}

impl AsRef<Path> for MethylationPath {
    fn as_ref(&self) -> &Path {
        match self {
            Self::BedGraph(path) | Self::Bed(path) => path.as_ref(),
        }
    }
}

impl PathLike for MethylationPath {}

#[derive(Debug, Clone, Copy)]
pub enum PathError {
    NoPath,
    NotUnicode,
}

impl std::error::Error for PathError {}

impl fmt::Display for PathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoPath => f.write_str("no path"),
            Self::NotUnicode => f.write_str("path is not a unicode string"),
        }
    }
}

/// Abstraction over path-like values.
pub trait PathLike: AsRef<Path> + fmt::Debug {
    /// Returns the underlying path reference.
    fn path(&self) -> &Path {
        self.as_ref()
    }
    /// Returns the filename as UTF-8.
    ///
    /// # Errors
    /// Returns an error if the path has no filename, or it is not valid Unicode.
    fn filename_str(&self) -> Result<&str, PathError> {
        self.path()
            .file_name()
            .ok_or(PathError::NoPath)?
            .to_str()
            .ok_or(PathError::NotUnicode)
    }
}
