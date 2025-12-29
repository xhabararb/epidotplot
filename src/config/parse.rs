use crate::config::cli::IntervalError;
use crate::config::error::OtherError;
use crate::config::{
    Input, MethylationInput, MethylationThreshold, PerAxis, SequenceInput, SequenceMode, WordLength,
};

use crate::bio::dna::domain::SequenceRegion;
use crate::bio::dna::parse::NativeSequenceParser;

use crate::bio::methylation::domain::MethylationValue;
use crate::bio::methylation::parse::NativeMethylationParser;
use crate::bio::{MethylationPath, SequencePath};
use crate::config::{Config, PlotStyle, error::ConfigError};
use std::ops::Range;
use std::path::Path;
use std::path::PathBuf;

#[derive(Debug, Default, Clone)]
pub(in crate::config) struct PartialConfig {
    pub sequence: Option<PathBuf>,
    pub sequence_fst: Option<PathBuf>,
    pub sequence_snd: Option<PathBuf>,

    pub region: Option<PerAxis<Option<String>>>,

    pub parallel: bool,

    pub fst_sequence_mode: Option<SequenceMode>,
    pub snd_sequence_mode: Option<SequenceMode>,

    pub forgiving_methylation_base_match: bool,

    pub methylation: Option<PathBuf>,
    pub methylation_fst: Option<PathBuf>,
    pub methylation_snd: Option<PathBuf>,

    pub output_dir: Option<PathBuf>,

    pub methylation_threshold: Option<f32>,
    pub word_length: Option<WordLength>,

    pub enlarge_small: Option<bool>,
    pub dot_size: Option<u32>,
    pub plot_side: Option<u32>,
}

pub(in crate::config) fn parse_sequence_os(s: &Path) -> Result<SequencePath, ConfigError> {
    parse_path_os(s, &[(&["fasta", "fa", "fna"], |p| SequencePath::Fasta(p))])
}

pub(in crate::config) fn parse_methylation_os(s: &Path) -> Result<MethylationPath, ConfigError> {
    parse_path_os(
        s,
        &[
            (&["bed"], |p| MethylationPath::Bed(p)),
            (&["bedgraph"], |p| MethylationPath::BedGraph(p)),
        ],
    )
}

#[allow(clippy::type_complexity)]
fn parse_path_os<T>(s: &Path, formats: &[(&[&str], fn(PathBuf) -> T)]) -> Result<T, ConfigError>
where
    T: Sized,
{
    let ext = s
        .extension()
        .ok_or_else(|| OtherError::from_str_repr("file extension missing"))?
        .to_str()
        .ok_or_else(|| OtherError::from_str_repr("file extension is not valid Unicode"))?;

    for (extensions, constructor) in formats {
        if extensions.iter().any(|f| f.eq_ignore_ascii_case(ext)) {
            return Ok(constructor(PathBuf::from(s)));
        }
    }

    let allowed_formats = formats
        .iter()
        .flat_map(|t| t.0.iter().copied())
        .collect::<Vec<_>>()
        .join(", ");

    Err(OtherError::from_str_repr(&format!(
        "unsupported file format; please use one of: {allowed_formats}"
    ))
    .into())
}

macro_rules! prefer_right_opt {
    ($self:ident, $rhs:ident, $field:ident) => {
        $self.$field = match ($self.$field.take(), $rhs.$field) {
            (None, Some(r)) => Some(r),
            (Some(_), Some(r)) => Some(r),
            (Some(l), None) => Some(l),
            (None, None) => None,
        };
    };
}

macro_rules! prefer_right_bool {
    ($self:ident, $rhs:ident, $field:ident) => {
        $self.$field = if $rhs.$field {
            $rhs.$field
        } else {
            $self.$field
        };
    };
}

impl PartialConfig {
    fn merge_right(mut self, rhs: Self) -> Self {
        prefer_right_opt!(self, rhs, sequence);
        prefer_right_opt!(self, rhs, sequence_fst);
        prefer_right_opt!(self, rhs, sequence_snd);

        prefer_right_opt!(self, rhs, fst_sequence_mode);
        prefer_right_opt!(self, rhs, snd_sequence_mode);

        prefer_right_opt!(self, rhs, region);

        prefer_right_bool!(self, rhs, forgiving_methylation_base_match);

        prefer_right_opt!(self, rhs, methylation);
        prefer_right_opt!(self, rhs, methylation_fst);
        prefer_right_opt!(self, rhs, methylation_snd);

        prefer_right_bool!(self, rhs, parallel);

        prefer_right_opt!(self, rhs, output_dir);

        prefer_right_opt!(self, rhs, word_length);

        prefer_right_opt!(self, rhs, methylation_threshold);

        prefer_right_opt!(self, rhs, enlarge_small);
        prefer_right_opt!(self, rhs, dot_size);
        prefer_right_opt!(self, rhs, plot_side);
        self
    }
}

impl<I> From<I> for PartialConfig
where
    I: IntoIterator<Item = Self>,
{
    fn from(it: I) -> Self {
        it.into_iter().fold(Self::default(), Self::merge_right)
    }
}

impl TryFrom<PartialConfig> for Config {
    type Error = ConfigError;

    fn try_from(mut raw: PartialConfig) -> Result<Self, Self::Error> {
        let sequence_path = match (raw.sequence, raw.sequence_fst, raw.sequence_snd) {
            (Some(single), None, None) => {
                SequenceInput::Single(parse_sequence_os(single.as_path())?)
            }
            (None, Some(fst), Some(snd)) => SequenceInput::Pair(PerAxis {
                fst: parse_sequence_os(fst.as_path())?,
                snd: parse_sequence_os(snd.as_path())?,
            }),
            (None, None, None) => {
                return Err(OtherError::from_str_repr("missing sequence input").into());
            }
            _ => {
                return Err(ConfigError::UnsupportedCombination {
                    detail: "cannot specify both single and pairwise sequence inputs".to_string(),
                });
            }
        };

        let methylation_path = match (raw.methylation, raw.methylation_fst, raw.methylation_snd) {
            (Some(single), None, None) => Some(MethylationInput::Single(parse_methylation_os(
                single.as_path(),
            )?)),
            (None, Some(fst), Some(snd)) => Some(MethylationInput::Pair(PerAxis {
                fst: parse_methylation_os(fst.as_path())?,
                snd: parse_methylation_os(snd.as_path())?,
            })),
            (None, None, None) => None,
            _ => {
                return Err(ConfigError::UnsupportedCombination {
                    detail: "cannot specify both single and pairwise sequence paths".to_string(),
                });
            }
        };

        let empty_id_err = || ConfigError::InvalidValue {
            field: "region",
            reason: "empty identifier in region specification".to_owned(),
        };

        let methylation_threshold = raw
            .methylation_threshold
            .map(|f| MethylationThreshold(MethylationValue::from_percent(f)))
            .unwrap_or_default();

        Ok(Self {
            input: Input {
                sequence_path,
                methylation_path,
                fst_sequence_mode: raw.fst_sequence_mode.unwrap_or_default(),
                snd_sequence_mode: raw.snd_sequence_mode.unwrap_or_default(),
                strict_methylation_base_match: raw.forgiving_methylation_base_match,
                region: {
                    match raw.region.take() {
                        None => Ok::<_, ConfigError>(PerAxis {
                            fst: None,
                            snd: None,
                        }),
                        Some(PerAxis { fst, snd }) => Ok(PerAxis {
                            fst: fst
                                .map(|s| SequenceRegion::parse(&s).ok_or_else(empty_id_err))
                                .transpose()?,
                            snd: snd
                                .map(|s| SequenceRegion::parse(&s).ok_or_else(empty_id_err))
                                .transpose()?,
                        }),
                    }?
                },
            },
            output_dir: raw.output_dir,
            methylation_threshold,
            word_len: raw.word_length.unwrap_or_default(),
            style: PlotStyle::new_from_options(raw.dot_size, raw.plot_side, raw.enlarge_small),
            sequence_parser: Box::new(NativeSequenceParser),
            methylation_parser: Box::new(NativeMethylationParser {
                min_threshold: methylation_threshold,
            }),
            parallel: raw.parallel,
        })
    }
}

fn _parse_interval(str: &str) -> Result<Range<usize>, IntervalError> {
    let str = str.chars().filter(|c| *c != ',').collect::<String>();
    let (l, h) = str.split_once('-').ok_or_else(|| {
        IntervalError::InvalidFormat(Some(
            "missing dash delimiter (-), expected: low-high".into(),
        ))
    })?;
    let (l, h) = (
        l.parse::<usize>().map_err(IntervalError::Parse)?,
        h.parse::<usize>().map_err(IntervalError::Parse)?,
    );
    Ok(l..h)
}
