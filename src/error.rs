use crate::bio::error::{IoError, ParseError};
use crate::bio::methylation::parse::MethylationError;
use crate::plot::render::RenderError;
use std::{error, error::Error, fmt};
use thiserror::Error;

#[derive(Debug)]
pub struct OffsetPosition {
    pub line: usize,
    pub column: usize,
}

#[derive(Debug, Error)]
pub enum CoreError {
    Io(#[from] IoError),
    InvalidInput(#[from] InvalidInputError),
    MalformedInput(#[from] MalformedInputError),
    Other {
        #[source]
        err: Box<dyn Error>,
        msg: Option<String>,
    },
    Parse(#[from] ParseError),
    Render(#[from] RenderError),
    Methylation(#[from] MethylationError),
}

#[derive(Debug, Error)]
pub struct MalformedInputError(pub Box<dyn error::Error>);
impl fmt::Display for MalformedInputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Malformed input: {}", self.0)
    }
}

#[derive(Debug, Error)]
pub struct InvalidInputError(pub Box<dyn error::Error>);

impl fmt::Display for InvalidInputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Invalid input: {}", self.0)
    }
}

impl fmt::Display for CoreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(s) => write!(f, "{s}"),
            Self::MalformedInput(s) => write!(f, "{s}"),
            Self::Other { err, msg: Some(m) } => write!(f, "Error: {err} ({m})"),
            Self::Other { err, msg: _ } => write!(f, "Error: {err}"),
            Self::Render(err) => write!(f, "{err}"),
            Self::Methylation(err) => write!(f, "{err}"),
            Self::Parse(err) => write!(f, "{err}"),
            Self::Io(err) => write!(f, "{err}"),
        }
    }
}
