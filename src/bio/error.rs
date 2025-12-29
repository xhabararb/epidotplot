use std::{fmt, io};
use thiserror::Error;

use crate::error::{InvalidInputError, MalformedInputError};

#[derive(Debug, Error)]
pub enum ParseError {
    #[error(transparent)]
    Io(#[from] IoError),
    #[error(transparent)]
    MalformedInput(#[from] MalformedInputError),
    #[error(transparent)]
    InvalidInput(#[from] InvalidInputError),
}

#[derive(Debug, Error)]
pub struct IoError {
    #[source]
    pub err: io::Error,
    pub message: Option<String>,
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.message {
            Some(message) => write!(f, "IO Error: {} ({message})", self.err),
            None => write!(f, "IO Error: {}", self.err),
        }
    }
}
