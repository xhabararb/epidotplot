use crate::config::{PerAxis, SequenceMode, WordLength, parse::PartialConfig};
use clap::{Args, Parser, Subcommand};
use std::borrow::Cow;
use std::fmt::{Debug, Display, Formatter};
use std::path::PathBuf;
use std::{error, num};

#[derive(Debug, Parser)]
#[command(subcommand_precedence_over_arg = false)]
pub(in crate::config) struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub(in crate::config) enum Commands {
    /// Generate a dotplot from the given sequence and optional methylation input.
    #[command(
        override_usage = "epidotplot plot <--sequence <PATH> | --fst-sequence <PATH>  --snd-sequence <PATH>> [OPTIONS]"
    )]
    Plot(CommonArgs),
}

#[derive(Debug, Parser)]
pub(in crate::config) struct CommonArgs {
    #[command(flatten)]
    pub input: Input,

    #[command(flatten)]
    pub regions: Regions,

    /// Enable parallel computation. Disabled by default.
    #[arg(short, long, action = clap::ArgAction::SetTrue)]
    pub parallel: bool,

    /// Output directory for rendered plots and other data. Default is the current directory.
    #[arg(short, long, value_name = "DIR")]
    pub output_dir: Option<PathBuf>,

    /// Percent threshold (0–100) for considering a site methylated, lower values are ignored. Default is 0.
    #[arg(short, long, value_name = "FLOAT")]
    pub methylation_threshold: Option<f32>,

    /// Word length (k-mer size) for search and exact matching within sequences. Must be positive. Default is 21.
    #[arg(short, long, value_name = "UINT")]
    pub word_length: Option<usize>,

    /// Load configuration from JSON (CLI flags take precedence).
    #[arg(short, long = "config", value_name = "FILE")]
    pub config_file: Option<PathBuf>,

    /// Enlarge the plot when its computed size falls below the target. Usually useful for small plots, though larger plots may be slightly smaller due to rounding too.
    #[arg(long, action = clap::ArgAction::SetTrue)]
    pub enlarge_small: bool,
    /// Pixel size of a single dotplot dot. Default is 1.
    #[arg(short, long, value_name = "UINT")]
    pub dot_size: Option<u32>,
    /// Maximum plot dimension; the longer axis is scaled to this size. Default is 4000.
    #[arg(long, value_name = "UINT")]
    pub plot_side: Option<u32>,
}

#[derive(Debug, Args)]
pub(in crate::config) struct Input {
    // *------------------------
    // Sequence input
    // *------------------------
    /// Path to the x-axis input DNA sequence in FASTA/Multi-FASTA (pairwise mode).
    #[arg(
        long,
        value_name = "PATH",
        requires = "snd_sequence",
        conflicts_with = "sequence"
    )]
    pub fst_sequence: Option<PathBuf>,
    /// Path to the y-axis input DNA sequence in FASTA/Multi-FASTA (pairwise mode).
    #[arg(
        long,
        value_name = "PATH",
        requires = "fst_sequence",
        conflicts_with = "sequence"
    )]
    pub snd_sequence: Option<PathBuf>,
    /// Path to the input DNA sequence in FASTA/Multi-FASTA (single-sequence mode: both x- and y-axis).
    #[arg(long, value_name = "PATH", conflicts_with_all = ["fst_sequence", "fst_sequence"])]
    pub sequence: Option<PathBuf>,

    /// Reverse-complement x-axis sequence (usable for both pairwise and single-sequence modes).
    #[arg(long, action = clap::ArgAction::SetTrue)]
    pub fst_rev: bool,
    /// Reverse-complement y-axis sequence (usable for both pairwise and single-sequence modes).
    #[arg(long, action = clap::ArgAction::SetTrue)]
    pub snd_rev: bool,

    // *------------------------
    // Methylation input
    // *------------------------
    /// Path to the x-axis input DNA methylation site data in bedMethyl/BedGraph format (pairwise mode).
    #[arg(long, value_name = "PATH")]
    pub fst_methylation: Option<PathBuf>,
    /// Path to the y-axis input DNA methylation site data in bedMethyl/BedGraph format (pairwise mode).
    #[arg(long, value_name = "PATH")]
    pub snd_methylation: Option<PathBuf>,
    /// Path to the input DNA methylation site data in bedMethyl/BedGraph format (single-sequence mode: both x- and y-axis).
    #[arg(long, value_name = "PATH")]
    pub methylation: Option<PathBuf>,

    /// Ignore methylation sites that don't correspond to cytosines in the sequence, printing a warning instead of aborting.
    #[arg(long = "forgiving", action = clap::ArgAction::SetTrue)]
    pub forgiving_methylation_base_match: bool,
}

#[derive(Debug, Args)]
pub(in crate::config) struct Regions {
    /// Selects only a specific region of the input sequence(s) to be used for plotting. REGION may be:
    ///
    /// <ID>
    ///
    /// – Entire sequence with this ID (same as in the FASTA/methylation input)
    ///
    /// <START-END>
    ///
    /// – Global interval on the concatenated sequences (concatenation follows MultiFASTA file order).
    ///
    /// - (0-based, half-open: start inclusive, end exclusive)
    ///
    /// <ID:START-END> – interval local to the sequence ID
    ///
    /// (0-based, half-open: start inclusive, end exclusive; coordinates local to that sequence)
    #[arg(long, value_name = "REGION")]
    pub region: Option<String>,

    /// Region for the x-axis. Same REGION syntax as --region.
    #[arg(long, value_name = "REGION")]
    pub fst_region: Option<String>,

    /// Region for the y-axis. Same REGION syntax as --region.
    #[arg(long, value_name = "REGION")]
    pub snd_region: Option<String>,
}

impl From<Cli> for PartialConfig {
    fn from(cli: Cli) -> Self {
        let Commands::Plot(args) = cli.command;
        let Regions {
            region,
            fst_region,
            snd_region,
        } = args.regions;

        Self {
            sequence: args.input.sequence,
            sequence_fst: args.input.fst_sequence,
            sequence_snd: args.input.snd_sequence,

            fst_sequence_mode: Some(if args.input.fst_rev {
                SequenceMode::ReverseComplement
            } else {
                SequenceMode::default()
            }),
            snd_sequence_mode: Some(if args.input.snd_rev {
                SequenceMode::ReverseComplement
            } else {
                SequenceMode::default()
            }),

            region: match (region, fst_region, snd_region) {
                (Some(f), None, None) => Some(PerAxis {
                    fst: Some(f.clone()),
                    snd: Some(f),
                }),
                (None, fst, snd) => Some(PerAxis { fst, snd }),
                _ => unreachable!(),
            },

            methylation: args.input.methylation,
            methylation_fst: args.input.fst_methylation,
            methylation_snd: args.input.snd_methylation,

            output_dir: args.output_dir,
            methylation_threshold: args.methylation_threshold,
            word_length: args.word_length.map(WordLength),
            forgiving_methylation_base_match: args.input.forgiving_methylation_base_match,
            parallel: args.parallel,

            enlarge_small: Some(args.enlarge_small),
            dot_size: args.dot_size,
            plot_side: args.plot_side,
        }
    }
}

#[derive(Debug)]
pub enum IntervalError {
    InvalidFormat(Option<Cow<'static, str>>),
    Parse(num::ParseIntError),
}

impl Display for IntervalError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidFormat(Some(msg)) => write!(f, "invalid interval: {msg}"),
            Self::InvalidFormat(None) => write!(f, "invalid interval"),
            Self::Parse(err) => write!(f, "failed to parse interval: {err}"),
        }
    }
}

impl error::Error for IntervalError {}
