use image::RgbImage;
use png::{Compression, Encoder, FilterType};
use std::ffi::OsStr;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::time::Duration;
use std::{fs, io, thread};

/// Ensures that the given path is unique. If the file already exists, this
/// function attempts to generate a non-conflicting filename by appending a
/// fallback suffix (e.g., timestamp). Up to three attempts are made.
///
/// # Errors
/// Returns an error if the target directory cannot be created, if I/O
/// operations fail, or if a unique filename cannot be generated within
/// the allowed number of attempts.
fn ensure_unique_path(
    path: PathBuf,
    default_filename: &str,
    default_ext: &str,
    fallback_suffix: impl Fn() -> String,
) -> io::Result<PathBuf> {
    let dir = if let Some(parent) = path.parent() {
        if !parent.exists() {
            fs::create_dir_all(parent)?;
        }
        parent
    } else {
        Path::new(".")
    };

    if !path.exists() {
        return Ok(path);
    }

    println!(
        "file {} already exists, appending timestamp...",
        path.to_str().unwrap_or("<invalid UTF-8>")
    );

    let stem = path
        .file_stem()
        .and_then(OsStr::to_str)
        .unwrap_or(default_filename);
    let ext = path
        .extension()
        .and_then(OsStr::to_str)
        .unwrap_or(default_ext);

    for _ in 0..3 {
        thread::sleep(Duration::from_millis(100));
        let candidate = dir.join(format!("{}__{}.{}", stem, fallback_suffix(), ext));
        if !candidate.exists() {
            return Ok(candidate);
        }
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        "file already exists even after 3 attempts",
    ))
}

/// Prepares an output path by combining an optional directory, an optional
/// filename, and the required extension. If the resulting path already exists,
/// a unique variant is produced via `ensure_unique_path`.
///
/// # Errors
/// Propagates any I/O errors encountered during directory creation or from
/// `ensure_unique_path`.
fn prepare_path(dir: Option<&Path>, filename: Option<&str>, ext: &str) -> io::Result<PathBuf> {
    let filename = format!("{}.{ext}", filename.unwrap_or("output"));
    ensure_unique_path(
        dir.unwrap_or_else(|| Path::new(".")).join(filename),
        "output",
        "out",
        timestamp,
    )
}

pub fn save_png(img: &RgbImage, dir: Option<&Path>, filename: Option<&str>) -> io::Result<PathBuf> {
    use png;

    let path = prepare_path(dir, filename, "png")?;
    let file = File::create(&path)?;

    let mut encoder = Encoder::new(file, img.width(), img.height());
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);

    encoder.set_compression(Compression::Best);
    encoder.set_filter(FilterType::NoFilter);

    let mut writer = encoder.write_header()?;
    writer.write_image_data(img.as_raw().as_slice())?;

    Ok(path)
}

#[inline]
fn timestamp() -> String {
    chrono::Local::now().format("%Y-%m-%d_%H-%M-%S").to_string()
}
