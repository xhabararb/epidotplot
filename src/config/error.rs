use crate::bio::error::IoError;
use std::{error, fmt, fmt::Formatter, path::PathBuf, str};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ConfigError {
    #[error(transparent)]
    Other(#[from] OtherError),
    #[error(transparent)]
    Cli(#[from] clap::Error),
    #[error(transparent)]
    Io(#[from] IoError),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
    #[error("unsupported file kind: {kind}")]
    UnsupportedFileKind { kind: String },
    #[error("unknown format of file {0}")]
    UnknownFormat(PathBuf),
    #[error("unsupported combination: {detail}")]
    UnsupportedCombination { detail: String },
    #[error("invalid value for *{field}* ({reason})")]
    InvalidValue { field: &'static str, reason: String },
}

#[derive(Debug, Error)]
pub struct OtherError(pub Box<dyn error::Error>);

impl OtherError {
    pub fn from_str_repr(s: &str) -> Self {
        Self(s.into())
    }

    pub fn from_general(err: impl Into<Box<dyn error::Error>>) -> Self {
        Self(err.into())
    }
}

impl fmt::Display for OtherError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
