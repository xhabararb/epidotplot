use epidotplot::config::error::ConfigError;
use epidotplot::config::{Config, MethylationInput, PerAxis, SequenceInput, load_config};
use epidotplot::error::{CoreError, InvalidInputError};
use epidotplot::plot::render::{global_font, init_font};
use epidotplot::run_pairwise::handle_pairwise;
use epidotplot::run_single::handle_single;
use std::process::ExitCode;

fn main() -> ExitCode {
    let usage_err = "See `./epidotplot --help` for usage.";
    let mut cfg = match load_config() {
        Ok(cfg) => cfg,
        Err(ConfigError::Cli(err)) => err.exit(),
        Err(ConfigError::Json(err)) => {
            eprintln!("failed to load JSON configuration: {err}");
            eprintln!("{usage_err}");
            return ExitCode::FAILURE;
        }
        Err(err) => {
            eprintln!("Configuration error: {err}");
            eprintln!("{usage_err}");
            return ExitCode::FAILURE;
        }
    };
    init_font();
    _ = global_font(); // panics when uninitialized

    if let Err(err) = cfg.validate() {
        eprintln!("Invalid configuration: {err}");
        return ExitCode::FAILURE;
    }

    if let Err(err) = run(&mut cfg) {
        eprintln!("{err}");
        return ExitCode::FAILURE;
    }

    ExitCode::SUCCESS
}

/// Top-level entry point for generating dotplots.
///
/// Dispatches to either identity mode (1 seq + 1 meth) or
/// distinct mode (2 seqs + 2 meths) based on config.
pub fn run(cfg: &mut Config) -> Result<(), CoreError> {
    let seq_in = cfg.input.sequence_path.clone();
    let meth_in = cfg.input.methylation_path.clone();
    let modes_eq = cfg.input.fst_sequence_mode == cfg.input.snd_sequence_mode;

    match (seq_in, meth_in) {
        (SequenceInput::Single(s), Some(MethylationInput::Single(m))) if modes_eq => {
            handle_single(cfg, &s, Some(&m))
        }
        (SequenceInput::Single(s), Some(MethylationInput::Single(m))) if !modes_eq => {
            handle_pairwise(
                cfg,
                PerAxis {
                    fst: s.clone(),
                    snd: s,
                },
                Some(PerAxis {
                    fst: m.clone(),
                    snd: m,
                }),
            )
        }

        (SequenceInput::Single(s), None) if modes_eq => handle_single(cfg, &s, None),
        (SequenceInput::Single(s), None) if !modes_eq => handle_pairwise(
            cfg,
            PerAxis {
                fst: s.clone(),
                snd: s,
            },
            None,
        ),

        (SequenceInput::Pair(s), Some(MethylationInput::Pair(m))) => {
            handle_pairwise(cfg, s, Some(m))
        }

        (SequenceInput::Pair(s), None) => handle_pairwise(cfg, s, None),

        _ => Err(InvalidInputError(
            "expected (1 sequence, 1 methylation) or (2 sequences, 2 methylations)".into(),
        )
        .into()),
    }?;

    Ok(())
}
