pub mod cli;
pub mod error;
mod file;
pub mod parse;

use crate::bio::dna::parse::SequenceParser;

use crate::bio::dna::domain::SequenceRegion;
use crate::bio::methylation::domain::MethylationValue;
use crate::bio::methylation::parse::MethylationParser;
use crate::bio::{MethylationPath, PathError, PathLike, SequencePath};
use crate::config::cli::Cli;
use crate::config::file::load_json;
use crate::config::{error::ConfigError, parse::PartialConfig};
use crate::error::InvalidInputError;
use crate::plot::PlotStyle;
use crate::util::to_num_pretty;
use clap::Parser;
use std::ops::Deref;
use std::path::PathBuf;

/// A pair of values bound to two plot axes.
///
/// This type deliberately avoids naming the axes `x`/`y` to prevent
/// accidental assumptions about orientation. Use when you need
/// to carry around “axis-paired” data without privileging one side.
#[derive(Debug)]
pub struct PerAxis<T> {
    pub fst: T,
    pub snd: T,
}

/// Sequence input specification.
///
/// - `Single` — one sequence used for both axes (self-dotplot).
/// - `Pair` — two sequences, one per axis.
#[derive(Debug, Clone)]
pub enum SequenceInput {
    Single(SequencePath),
    Pair(PerAxis<SequencePath>),
}

impl<T> From<(T, T)> for PerAxis<T> {
    fn from((fst, snd): (T, T)) -> Self {
        Self { fst, snd }
    }
}

impl<T> PerAxis<T> {
    #[inline]
    pub fn map<U>(self, f: impl Fn(T) -> U) -> PerAxis<U> {
        let Self { fst, snd } = self;
        PerAxis {
            fst: f(fst),
            snd: f(snd),
        }
    }
    #[inline]
    pub fn map_ref<U>(&self, f: impl Fn(&T) -> U) -> PerAxis<U> {
        PerAxis {
            fst: f(&self.fst),
            snd: f(&self.snd),
        }
    }
    pub fn into_inner(self) -> (T, T) {
        let Self { fst, snd } = self;
        (fst, snd)
    }
}

impl<T: Copy> PerAxis<T> {
    pub const fn copied(&self) -> Self {
        let Self { fst, snd } = *self;
        Self { fst, snd }
    }
}

impl<T: Clone> Clone for PerAxis<T> {
    fn clone(&self) -> Self {
        Self {
            fst: self.fst.clone(),
            snd: self.snd.clone(),
        }
    }
}

impl SequenceInput {
    pub fn two_filenames(&self) -> Result<(String, String), PathError> {
        Ok(match self {
            Self::Single(single) => {
                let r = single.filename_str()?.to_owned();
                (r.clone(), r)
            }
            Self::Pair(PerAxis { fst, snd }) => (
                fst.filename_str()?.to_owned(),
                snd.filename_str()?.to_owned(),
            ),
        })
    }
}

/// Methylation input specification.
///
/// - `Single` — one methylation track used for both axes.
/// - `Pair` — two methylation tracks, one per axis.
#[derive(Debug, Clone)]
pub enum MethylationInput {
    Single(MethylationPath),
    Pair(PerAxis<MethylationPath>),
}

/// Threshold for calling a site methylated.
#[derive(Debug, Clone, Copy)]
pub struct MethylationThreshold(MethylationValue);

impl Default for MethylationThreshold {
    fn default() -> Self {
        Self(MethylationValue::from_percent(0.0))
    }
}

impl From<MethylationValue> for MethylationThreshold {
    fn from(value: MethylationValue) -> Self {
        Self(value)
    }
}

impl Deref for MethylationThreshold {
    type Target = MethylationValue;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Word length (`k`) for k-mer matching.
#[derive(Debug, Clone, Copy)]
pub struct WordLength(usize);

impl Default for WordLength {
    fn default() -> Self {
        Self(21)
    }
}

impl Deref for WordLength {
    type Target = usize;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceMode {
    Original,
    ReverseComplement,
}

impl Default for SequenceMode {
    fn default() -> Self {
        Self::Original
    }
}

#[derive(Debug, Clone)]
pub struct Input {
    /// Sequence input specification (FASTA, single or pair).
    pub sequence_path: SequenceInput,
    /// Methylation input specification (BED/bedGraph, single or pair).
    pub methylation_path: Option<MethylationInput>,

    pub fst_sequence_mode: SequenceMode,
    pub snd_sequence_mode: SequenceMode,

    pub strict_methylation_base_match: bool,
    pub region: PerAxis<Option<SequenceRegion>>,
}

/// Global configuration resolved from CLI, ENV and file (now unsupported).
///
/// Holds input sources, thresholds, style settings,
/// parser contexts and output options.
#[derive(Debug)]
pub struct Config {
    pub input: Input,

    pub parallel: bool,

    /// Output directory for generated plots and data.
    pub output_dir: Option<PathBuf>,

    /// Threshold for calling a site methylated.
    pub methylation_threshold: MethylationThreshold,
    /// Word length (`k`) used in k-mer matching.
    pub word_len: WordLength,

    /// Plot style configuration (e.g., colors, dot sizes).
    pub style: PlotStyle,

    /// Parser for sequence input.
    pub sequence_parser: Box<dyn SequenceParser>,
    /// Parser for methylation input.
    pub methylation_parser: Box<dyn MethylationParser>,
}

impl Config {
    pub fn validate(&self) -> Result<(), InvalidInputError> {
        let methyl_threshold = self.methylation_threshold.as_percent();
        if !(0.0..=100.0).contains(&methyl_threshold) {
            return Err(InvalidInputError(
                format!(
                    "methylation threshold % must be within [0,100], but is {methyl_threshold:.1}%"
                )
                .into(),
            ));
        }

        let allowed_plot_side = 350..=12_500;
        if !allowed_plot_side.contains(&self.style.plot_side) {
            return Err(InvalidInputError(
                format!(
                    "the plot's longer dimension must lie within [{}, {}] range",
                    to_num_pretty(allowed_plot_side.start()),
                    to_num_pretty(allowed_plot_side.end())
                )
                .into(),
            ));
        }

        if *self.word_len < 1 {
            return Err(InvalidInputError(
                "word (k-mer) length must be positive".into(),
            ));
        }

        Ok(())
    }
}

pub fn load_config() -> Result<Config, ConfigError> {
    let cli = Cli::try_parse().map_err(ConfigError::Cli)?;

    let cli::Commands::Plot(ref args) = cli.command;
    let precedence = if let Some(file) = &args.config_file {
        println!("loading config from file at {}", file.display());
        let file = load_json(file.as_path())?;
        vec![file, PartialConfig::from(cli)]
    } else {
        vec![PartialConfig::from(cli)]
    };

    PartialConfig::from(precedence).try_into()
}
