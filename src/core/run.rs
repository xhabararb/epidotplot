use crate::bio::dna::domain::SequenceBytes;
use crate::bio::dna::fasta::contains_n;
use crate::bio::error::IoError;
use crate::bio::methylation::SingleMethylation;
use crate::bio::methylation::plot::{MethylationCoordinate, add_methylation};
use crate::config::{Config, PerAxis, SequenceMode};
use crate::core::sais;
use crate::error::CoreError;
use crate::plot::io::save_png;
use crate::plot::render::annotate_image;
use crate::plot::scale::TileScale;
use crate::plot::{Annotations, AxisAnnotation};
use crate::run_pairwise::Dotplot;
use rayon::ThreadPoolBuilder;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSliceMut;
use std::borrow::Cow;
use std::cell::RefCell;
use std::ops::DerefMut;
use std::rc::Rc;
use std::sync::{Arc, Condvar, Mutex};
use std::time::Duration;

pub struct Reporter {
    bar_style: indicatif::ProgressStyle,
    spinner_style: indicatif::ProgressStyle,
}

impl Default for Reporter {
    fn default() -> Self {
        Self {
            bar_style: Self::create_bar_style(),
            spinner_style: Self::create_spinner_style(),
        }
    }
}

impl Reporter {
    fn create_bar_style() -> indicatif::ProgressStyle {
        indicatif::ProgressStyle::with_template(
            "[{elapsed_precise}] [{bar:40.green/white}] {pos}/{len} {msg}",
        )
        .expect("failed to create progress bar style")
        .progress_chars("## ")
    }

    fn create_spinner_style() -> indicatif::ProgressStyle {
        indicatif::ProgressStyle::with_template("[{elapsed_precise}] {spinner:.green} {msg}")
            .expect("failed to create spinner style")
    }

    pub fn create_bar(
        &self,
        total: usize,
        label: impl Into<Cow<'static, str>>,
    ) -> indicatif::ProgressBar {
        indicatif::ProgressBar::new(total as u64)
            .with_style(self.bar_style.clone())
            .with_message(label)
    }

    pub fn create_spinner(&self, label: impl Into<Cow<'static, str>>) -> indicatif::ProgressBar {
        let spinner = indicatif::ProgressBar::new_spinner()
            .with_style(self.spinner_style.clone())
            .with_message(label);
        spinner.enable_steady_tick(Duration::from_millis(150));
        spinner
    }
}

#[allow(clippy::too_many_arguments)]
pub fn compute_methylated_dotplot(
    config: &Config,
    sequences: &PerAxis<&[u8]>,
    methylations: Option<PerAxis<Rc<SingleMethylation>>>,
    fst_sa: &[u32],
    fst_lcp: &[u32],
    scale: TileScale,
    ctx: &Reporter,
) -> Dotplot {
    let out_w = scale.out_w;
    let mut dotplot = Dotplot::new(scale);

    let buf = dotplot.buf.as_mut_slice();
    let mut write_meth_dots = |c: MethylationCoordinate| {
        let MethylationCoordinate { x, y, w } = c;
        let rgb = (config.style.colormap_raw)(w);
        let i = y * out_w * 3 + x * 3;
        buf[i..i + 3].copy_from_slice(&rgb);
    };

    if let Some(methylations) = methylations {
        let PerAxis { fst, snd } = methylations;
        let fst = fst
            .methylation
            .sites
            .iter()
            .map(|t| (*t.0, *t.1))
            .collect::<Vec<_>>();
        let snd = snd
            .methylation
            .sites
            .iter()
            .map(|t| (*t.0, *t.1))
            .collect::<Vec<_>>();
        add_methylation(
            &PerAxis {
                fst: &fst,
                snd: &snd,
            },
            &dotplot.scale,
            &mut write_meth_dots,
            ctx,
        );
    }

    let dot_size = config.style.dot_size;
    if dot_size > 1 {
        dotplot.enlarge_non_bw_dots(config.style.dot_size);
    }

    run_stripes(
        sequences,
        fst_sa,
        fst_lcp,
        &dotplot.scale,
        dotplot.buf.as_mut_slice(),
        *config.word_len,
        ctx,
        config.parallel,
    );
    if dot_size > 1 {
        dotplot.enlarge_black_dots(config.style.dot_size);
    }

    dotplot
}

#[allow(clippy::too_many_arguments)]
fn run_stripes(
    sequences: &PerAxis<&[u8]>,
    fst_sa: &[u32],
    fst_lcp: &[u32],
    scale: &TileScale,
    buf: &mut [u8],
    k: usize,
    ctx: &Reporter,
    parallel: bool,
) {
    let stripes = scale.out_h;

    let stride = 3 * scale.out_w;
    let total_kmers_y = sequences.snd.len().saturating_sub(k) + 1;

    if parallel {
        let thread_count = 64;
        let (acquire, release) =
            create_acquire_release_for_stripe_buf_reuse_parallel(thread_count, scale.out_w);
        let pool = ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build()
            .expect("failed to create a local thread pool");

        pool.install(|| {
            let bar = Arc::new(Mutex::new(
                ctx.create_bar(stripes, "calculating similarity between k-mers (parallel)"),
            ));
            bar.lock().expect("progress bar mutex poisoned").inc(0);
            buf.par_chunks_mut(stride)
                .enumerate()
                .for_each(|(stripe_pos_y, stripe_out)| {
                    let b = Arc::clone(&bar);
                    process_stripe(
                        stripe_pos_y,
                        stripe_out,
                        sequences,
                        fst_sa,
                        fst_lcp,
                        k,
                        scale,
                        total_kmers_y,
                        &acquire,
                        &release,
                        &|| (*b).lock().expect("progress bar mutex poisoned").inc(1),
                    );
                });
            (*bar).lock().expect("progress bar mutex poisoned").finish();
        });
    } else {
        let bar = ctx.create_bar(
            stripes,
            "calculating similarity between k-mers (sequential)",
        );
        let inc = || bar.inc(1);
        buf.chunks_mut(stride)
            .enumerate()
            .for_each(|(stripe_pos_y, stripe_out)| {
                let b = RefCell::new(Some(vec![0usize; scale.out_w]));
                let acquire = || b.replace(None).unwrap_or_else(|| vec![0usize; scale.out_w]);
                let release = |mut buf: Vec<usize>| {
                    buf.fill(0);
                    _ = b.replace(Some(buf));
                };
                process_stripe(
                    stripe_pos_y,
                    stripe_out,
                    sequences,
                    fst_sa,
                    fst_lcp,
                    k,
                    scale,
                    total_kmers_y,
                    acquire,
                    release,
                    &inc,
                );
            });
    }
}

#[allow(clippy::too_many_arguments)]
fn process_stripe<T: DerefMut<Target = [usize]>>(
    stripe_pos_y: usize,
    stripe_out: &mut [u8],
    sequences: &PerAxis<&[u8]>,
    fst_sa: &[u32],
    fst_lcp: &[u32],
    k: usize,
    scale: &TileScale,
    total_kmers_y: usize,
    acquire_stripe_buf: impl FnOnce() -> T,
    release_stripe_buf: impl FnOnce(T),
    report_progress: &dyn Fn(),
) {
    let stripe_start = stripe_pos_y * scale.tile_h;
    let stripe_end = ((stripe_pos_y + 1) * scale.tile_h).min(total_kmers_y);
    if stripe_start >= total_kmers_y {
        return;
    }

    let mut tile_hits = acquire_stripe_buf();

    (stripe_start..stripe_end).for_each(|y_pos| {
        let kmer = &sequences.snd[y_pos..y_pos + k];
        if let Some(range) = sais::find_kmer_range(sequences.fst, fst_sa, fst_lcp, kmer, contains_n)
            && let Some(slice) = fst_sa.get(range)
        {
            for &pos_x in slice {
                let tile_x = (pos_x as usize / scale.tile_w).min(scale.out_w - 1);
                tile_hits[tile_x] += 1;
            }
        }
    });

    for (tile_x, &hits) in tile_hits.iter().enumerate() {
        if hits > 0 {
            let i = tile_x * 3;
            stripe_out[i..i + 3].fill(Dotplot::BLACK);
        }
    }

    release_stripe_buf(tile_hits);

    report_progress();
}

#[allow(clippy::type_complexity)]
fn create_acquire_release_for_stripe_buf_reuse_parallel(
    thread_count: usize,
    buf_size: usize,
) -> (
    Box<dyn Fn() -> Vec<usize> + Send + Sync>,
    Box<dyn Fn(Vec<usize>) + Send + Sync>,
) {
    let pool = Arc::new((
        Mutex::new(
            (0..thread_count)
                .map(|_| vec![0usize; buf_size])
                .collect::<Vec<_>>(),
        ),
        Condvar::new(),
    ));
    let acquire = {
        let pool = Arc::clone(&pool);
        Box::new(move || {
            let (lock, cv) = &*pool;
            let mut guard = lock
                .lock()
                .expect("poisoned stripe buffer mutex while acquiring");
            loop {
                if let Some(buf) = guard.pop() {
                    break buf;
                }
                guard = cv
                    .wait(guard)
                    .expect("wait on poisoned stripe buffer mutex");
            }
        })
    };

    let release = {
        Box::new(move |mut buf: Vec<usize>| {
            // reinitialize counts in the buffer to zeros
            buf.fill(0);
            let (lock, cv) = &*pool;
            lock.lock()
                .expect("poisoned stripe buffer mutex while releasing")
                .push(buf);
            cv.notify_one();
        })
    };

    (acquire, release)
}

pub fn process_sequences(
    input: &PerAxis<Rc<SequenceBytes>>,
    methylations: Option<PerAxis<Rc<SingleMethylation>>>,
    config: &mut Config,
    reporter: &Reporter,
) -> Result<(), CoreError> {
    let (fst_seq, snd_seq) = (input.fst.data.as_slice(), input.snd.data.as_slice());
    let (sa, lcp) = {
        let spinner = reporter.create_spinner(format!(
            "building suffix array and LCP array{}",
            if config.parallel { " (parallel)" } else { "" }
        ));
        let res = sais::build_sa_lcp(
            fst_seq,
            config.parallel,
            |s| spinner.println(s),
            |s| spinner.println(s),
        )
        .inspect_err(|_| {
            spinner.finish_and_clear();
        })?;
        spinner.finish_with_message("finished building suffix array and LCP array");
        res
    };

    let scale = TileScale::new_with_max_side(
        config.style.plot_side as usize,
        fst_seq.len(),
        snd_seq.len(),
    );

    let dotplot = compute_methylated_dotplot(
        config,
        &PerAxis {
            fst: fst_seq,
            snd: snd_seq,
        },
        methylations,
        &sa,
        &lcp,
        scale,
        reporter,
    );

    annotate_and_save(config, dotplot, input)
}

fn annotate_and_save(
    config: &Config,
    dotplot: Dotplot,
    input: &PerAxis<Rc<SequenceBytes>>,
) -> Result<(), CoreError> {
    let SequenceBytes {
        id: fst_id,
        multi_fasta: fst_multi_fasta,
        ..
    } = input.fst.as_ref();
    let SequenceBytes {
        id: snd_id,
        multi_fasta: snd_multi_fasta,
        ..
    } = input.snd.as_ref();

    let fst_range = 0..input.fst.data.len();
    let snd_range = 0..input.snd.data.len();

    let (fst_filename, snd_filename) =
        config
            .input
            .sequence_path
            .two_filenames()
            .map_err(|err| CoreError::Other {
                err: Box::new(err),
                msg: None,
            })?;

    let fst_axis = AxisAnnotation::new_fasta(
        fst_filename,
        fst_id.to_region_specifier_str_pretty(),
        fst_range,
        config.input.fst_sequence_mode == SequenceMode::ReverseComplement,
        fst_multi_fasta.unwrap_or(false),
    );

    let snd_axis = AxisAnnotation::new_fasta(
        snd_filename,
        snd_id.to_region_specifier_str_pretty(),
        snd_range,
        config.input.snd_sequence_mode == SequenceMode::ReverseComplement,
        snd_multi_fasta.unwrap_or(false),
    );

    let img = dotplot.into_rgb_image().ok_or(CoreError::Other {
        err: "failed to create image from raw buffer".into(),
        msg: None,
    })?;

    let annotations = Annotations {
        axes: Annotations::distinct_axes(fst_axis, snd_axis),
        word_len: config.word_len,
        methylation_threshold: Some(config.methylation_threshold),
    };

    println!("rendering image...");
    let img = annotate_image(
        img,
        &annotations,
        config.style.enlarge_small,
        config.style.plot_side,
    );

    let filename = format!("{}__{}", fst_id.to_filename_str(), snd_id.to_filename_str());
    save_png(&img, config.output_dir.as_deref(), Some(&filename))
        .map(|res| println!("plot saved as PNG to {}", res.display()))
        .map_err(|err| {
            IoError {
                err,
                message: Some("failed to write PNG".into()),
            }
            .into()
        })
}
