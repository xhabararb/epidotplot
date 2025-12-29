use std::path::{Path, PathBuf};

use crate::bio::error::IoError;
use crate::config::{PerAxis, SequenceMode, WordLength, error::ConfigError, parse::PartialConfig};
use serde::Deserialize;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct JsonConfig {
    pub sequence: Option<String>,
    pub sequence_fst: Option<String>,
    pub sequence_snd: Option<String>,

    pub fst_rev: Option<bool>,
    pub snd_rev: Option<bool>,

    pub region: Option<String>,
    pub fst_region: Option<String>,
    pub snd_region: Option<String>,

    pub methylation: Option<String>,
    pub methylation_fst: Option<String>,
    pub methylation_snd: Option<String>,

    pub forgiving: Option<bool>,
    pub parallel: Option<bool>,

    pub output_dir: Option<String>,

    pub methylation_threshold: Option<f32>,
    pub word_length: Option<usize>,

    pub enlarge_small: Option<bool>,
    pub dot_size: Option<u32>,
    pub plot_side: Option<u32>,
}

impl TryFrom<JsonConfig> for PartialConfig {
    type Error = ConfigError;

    fn try_from(json: JsonConfig) -> Result<Self, Self::Error> {
        let region = match (json.region, json.fst_region, json.snd_region) {
            (Some(f), None, None) => Some(PerAxis {
                fst: Some(f.clone()),
                snd: Some(f),
            }),
            (None, fst, snd) => Some(PerAxis { fst, snd }),
            _ => {
                return Err(ConfigError::UnsupportedCombination {
                    detail: "cannot specify both pairwise and single input regions at once"
                        .to_owned(),
                });
            }
        };

        Ok(Self {
            sequence: json.sequence.map(PathBuf::from),
            sequence_fst: json.sequence_fst.map(PathBuf::from),
            sequence_snd: json.sequence_snd.map(PathBuf::from),
            region,
            parallel: json.parallel.unwrap_or(false),
            fst_sequence_mode: Some(if json.fst_rev.unwrap_or(false) {
                SequenceMode::ReverseComplement
            } else {
                SequenceMode::default()
            }),
            snd_sequence_mode: Some(if json.snd_rev.unwrap_or(false) {
                SequenceMode::ReverseComplement
            } else {
                SequenceMode::default()
            }),
            forgiving_methylation_base_match: json.forgiving.unwrap_or(false),
            methylation: json.methylation.map(PathBuf::from),
            methylation_fst: json.methylation_fst.map(PathBuf::from),
            methylation_snd: json.methylation_snd.map(PathBuf::from),
            output_dir: json.output_dir.map(PathBuf::from),
            methylation_threshold: json.methylation_threshold,
            word_length: json.word_length.map(WordLength),
            enlarge_small: json.enlarge_small,
            dot_size: json.dot_size,
            plot_side: json.plot_side,
        })
    }
}

pub fn load_json(path: &Path) -> Result<PartialConfig, ConfigError> {
    let data = std::fs::read_to_string(path)
        .map_err(|err| ConfigError::Io(IoError { err, message: None }))?;
    let cfg: JsonConfig = serde_json::from_str(&data).map_err(ConfigError::Json)?;
    PartialConfig::try_from(cfg)
}
