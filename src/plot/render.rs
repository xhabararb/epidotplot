use crate::bio::methylation::domain::MethylationValue;
use crate::config::MethylationThreshold;
use crate::plot::Annotations;
use crate::util::color::colormap_weighted;
use crate::util::to_num_pretty;
use ab_glyph::{Font, FontArc, PxScale, ScaleFont};
use image::imageops::FilterType;
use image::{ImageBuffer, Rgb, RgbImage, imageops};
use imageproc::drawing::{draw_filled_rect_mut, draw_text_mut};
use imageproc::rect::Rect;
use std::num::NonZeroUsize;
use std::ops::Range;
use std::sync::OnceLock;
use std::{fmt, io};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum RenderError {
    Io(#[from] io::Error),
    DimensionsOutOfRange {
        width: usize,
        height: usize,
        max: Option<(NonZeroUsize, NonZeroUsize)>,
    },
    ExternalImageCrate(#[from] image::ImageError),
}

impl fmt::Display for RenderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let msg = match self {
            Self::DimensionsOutOfRange { width, height, max } => {
                if let Some((mw, mh)) = max {
                    format!("image dimensions {width}×{height} exceed maximum {mw}×{mh}")
                } else {
                    format!("image dimensions {width}×{height} out of range")
                }
            }
            Self::ExternalImageCrate(e) => format!("image processing ({e})"),
            Self::Io(io) => format!("IO ({io})"),
        };
        write!(f, "Render error: {msg}")
    }
}

static FONT: OnceLock<FontArc> = OnceLock::new();

/// Initializes the global font. The use of panicky `expect` is justified, as valid
/// embedded font data is a mandatory program invariant.
///
/// # Panics
/// If the embedded font cannot be parsed.
#[deny(dead_code)]
pub fn init_font() {
    FONT.get_or_init(|| {
        let font_data = include_bytes!("../../public/DejaVuSans.ttf");
        FontArc::try_from_slice(font_data as &[u8]).expect("failed to parse font")
    });
}

/// # Panics
/// If the font has not been initialized.
pub fn global_font() -> &'static FontArc {
    FONT.get().expect("font not initialized")
}

pub fn annotate_image(
    img: ImageBuffer<Rgb<u8>, Vec<u8>>,
    annotations: &Annotations,
    enlarge_small: bool,
    plot_side: u32,
) -> RgbImage {
    let mut lines = vec![
        format!(
            "{} vs. {}",
            &annotations.axes.fst.filename_and_type(),
            &annotations.axes.snd.filename_and_type(),
        ),
        String::new(),
        format!("Word length: {}", *annotations.word_len),
    ];

    if let Some(threshold) = annotations.methylation_threshold {
        lines.push(format!(
            "Methylation threshold: {:.1}%",
            threshold.as_percent()
        ));
    }

    if enlarge_small && img.width() < plot_side && img.height() < plot_side {
        lines.push(String::new());
        lines.push(format!(
            "(plot was enlarged to fit {}px resolution)",
            to_num_pretty(&plot_side)
        ));
    }

    let sq = plot_side;
    let pad = sq / 10;

    let sidebar_w = (pad as f32 * 1.2).round() as u32;
    let sidebar_h = (sq as f32 * 0.75).round() as u32;

    let line_h = pad / 5;
    // largest expected tick label
    let tick_label_space =
        (text_width(global_font(), PxScale::from(line_h as f32), "3,000,000,000") as f32 * 1.2)
            .round() as u32;
    let pad_h = pad + tick_label_space;
    let line_space = line_h / 4;
    let pad_without_lines = pad;
    let pad_v = pad_without_lines + (lines.len() as u32 + 2) * (line_space + line_h);
    let border_thickness = (12.0 * (f64::from(plot_side) / 4_000.0)).round() as u32; // 12 px ~ 4000 px side

    let img = {
        let sq = f64::from(sq);
        let (w, h) = (f64::from(img.width()), f64::from(img.height()));
        if enlarge_small && w < sq && h < sq {
            if w < h {
                let f = f64::from(plot_side) / h;
                let large_w = f * w;
                imageops::resize(&img, large_w as u32, sq as u32, FilterType::Nearest)
            } else {
                let f = f64::from(plot_side) / w;
                let large_h = f * h;
                imageops::resize(&img, sq as u32, large_h as u32, FilterType::Nearest)
            }
        } else {
            img
        }
    };

    let width = img.width();
    let height = img.height();

    let canvas_w = sq + 2 * pad_h + sidebar_w;
    let canvas_h = sq + 2 * pad_v + tick_label_space;

    let scale = PxScale::from(line_h as f32);

    // the plot is at most 4k by 4k
    // 4000×4000 white background
    let mut canvas = RgbImage::from_pixel(canvas_w, canvas_h, Rgb([255, 255, 255]));

    let (plot_start_x, plot_start_y) = {
        let inner_pad_h = (sq - width) / 2;
        let inner_pad_v = (sq - height) / 2;
        let x = pad_h + inner_pad_h;
        let y = pad_v + inner_pad_v;
        imageops::overlay(&mut canvas, &img, i64::from(x), i64::from(y));
        draw_borders(
            &mut canvas,
            x,
            y,
            width,
            height,
            border_thickness,
            Rgb([0, 0, 0]),
        );
        (x, y)
    };

    draw_axes_labels(
        &mut canvas,
        plot_start_x,
        plot_start_y,
        width,
        height,
        line_h,
        &annotations.axes.fst.name,
        &annotations.axes.snd.name,
    );

    draw_lines(
        &mut canvas,
        line_h,
        Rgb([0, 0, 0]),
        pad_without_lines,
        pad_h,
        &lines,
        line_space,
    );

    /*draw_borders(
        &mut canvas,
        pad_h,
        pad_v,
        sq,
        sq,
        border_thickness,
        Rgb([0, 0, 0]),
    );*/

    draw_sidebar(
        &mut canvas,
        pad_h,
        pad_v,
        sq,
        sidebar_w,
        sidebar_h,
        colormap_weighted,
        &[(0.0, "0%"), (0.3, "30%"), (0.5, "50%"), (1.0, "100%")],
        scale,
        border_thickness,
        annotations.methylation_threshold,
        plot_side as f32,
    );

    // misaligned ticks are noticeable on a small dataset
    let adjust_ticks = move |range: &Range<usize>| -> u32 {
        if range.len() >= 100 {
            return 6;
        }
        for segments in [5u32, 4, 3, 2, 1] {
            if range.len().is_multiple_of(segments as usize) {
                return segments + 1;
            }
        }
        2
    };

    let tick_count = adjust_ticks(&annotations.axes.snd.range);

    let tick_h = (line_h as f32 * 0.7).round() as u32;
    // if the plot is too small, ticks would overlap
    if tick_count * line_h < height {
        draw_ticks_vertical(
            &mut canvas,
            plot_start_x,
            plot_start_y,
            width,
            height,
            scale,
            tick_h,
            border_thickness,
            tick_count,
            Rgb([0, 0, 0]),
            annotations.axes.snd.range.start,
            annotations.axes.snd.range.end,
        );
    }

    let tick_count = adjust_ticks(&annotations.axes.fst.range);

    if tick_count * line_h < width {
        draw_ticks_horizontal(
            &mut canvas,
            plot_start_x,
            plot_start_y,
            width,
            height,
            tick_h,
            scale,
            border_thickness,
            tick_count,
            Rgb([0, 0, 0]),
            annotations.axes.fst.range.start,
            annotations.axes.fst.range.end,
        );
    }

    canvas
}

#[allow(clippy::too_many_arguments)]
fn draw_axes_labels(
    canvas: &mut ImageBuffer<Rgb<u8>, Vec<u8>>,
    plot_start_x: u32,
    plot_start_y: u32,
    plot_w: u32,
    plot_h: u32,
    line_h: u32,
    x_label: &str,
    y_label: &str,
) {
    let scale = PxScale::from(line_h as f32);
    let font = global_font();

    {
        let x_width = text_width(font, scale, x_label);
        if x_width <= plot_w {
            let x_start = plot_start_x + (plot_w - x_width) / 2;
            let y_start = plot_start_y - line_h - (scale.y / 2.0) as u32;

            draw_text_mut(
                canvas,
                Rgb([0, 0, 0]),
                x_start as i32,
                y_start as i32,
                scale,
                font,
                x_label,
            );
        }
    }

    let y_width = text_width(font, scale, y_label);
    if y_width <= plot_h {
        let x_start = plot_start_x + plot_w + (scale.y / 2.0) as u32;
        let y_start = plot_start_y + (plot_h - y_width) / 2;

        let mut buf = RgbImage::from_pixel(y_width, scale.y as u32, Rgb([255, 255, 255]));
        draw_text_mut(&mut buf, Rgb([0, 0, 0]), 0, 0, scale, font, y_label);
        let rotated = imageops::rotate90(&buf);

        imageops::overlay(canvas, &rotated, x_start as i64, y_start as i64);
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_ticks_vertical(
    canvas: &mut ImageBuffer<Rgb<u8>, Vec<u8>>,
    plot_start_x: u32,
    plot_start_y: u32,
    _plot_w: u32,
    plot_h: u32,
    scale: PxScale,
    tick_h: u32,
    tick_w: u32,
    tick_count: u32,
    tick_color: Rgb<u8>,
    start: usize,
    end: usize,
) {
    let mut draw = |x, y, text: String| {
        draw_text_mut(
            canvas,
            tick_color,
            x - text_width(global_font(), scale, &text) as i32 - (tick_h / 2) as i32,
            y - (scale.y / 2.0) as i32,
            scale,
            global_font(),
            &text,
        );
        draw_filled_rect_mut(canvas, Rect::at(x, y).of_size(tick_h, tick_w), tick_color);
    };

    let length = end - start;
    let value_step = length / (tick_count - 1) as usize;
    let mut value = start;

    let x_step = plot_h / (tick_count - 1);

    let x = (plot_start_x - tick_h) as i32;
    let mut y = plot_start_y;
    for i in 0..tick_count {
        let t = if i == tick_count - 1 { end } else { value };

        draw(x, y as i32, to_num_pretty(&t));
        y += x_step;
        value += value_step;
    }
    //draw(x, (plot_start_y + plot_h) as i32, to_num_pretty(&end));
}

#[allow(clippy::too_many_arguments)]
fn draw_ticks_horizontal(
    canvas: &mut ImageBuffer<Rgb<u8>, Vec<u8>>,
    plot_start_x: u32,
    plot_start_y: u32,
    plot_w: u32,
    plot_h: u32,
    tick_h: u32,
    scale: PxScale,
    tick_w: u32,
    tick_count: u32,
    tick_color: Rgb<u8>,
    start: usize,
    end: usize,
) {
    debug_assert!(tick_count > 1);

    let mut draw = |x: i32, y, text: String| {
        let text_w = text_width(global_font(), scale, &text);
        let mut buf = RgbImage::from_pixel(text_w, scale.y as u32, Rgb([255, 255, 255]));
        draw_text_mut(&mut buf, tick_color, 0, 0, scale, global_font(), &text);
        let rotated = imageops::rotate270(&buf);

        imageops::overlay(
            canvas,
            &rotated,
            i64::from(x) - (scale.y / 2.0) as i64,
            i64::from(y) + (tick_h as f32 * 1.5) as i64,
        );

        draw_filled_rect_mut(canvas, Rect::at(x, y).of_size(tick_w, tick_h), tick_color);
    };

    let length = end - start;
    let value_step = length / (tick_count - 1) as usize;
    let mut value = start;

    let x_step = plot_w / (tick_count - 1);
    let y = (plot_start_y + plot_h) as i32;
    let mut x = plot_start_x;

    for _ in 0..tick_count {
        draw(x as i32, y, to_num_pretty(&value));
        x += x_step;
        value += value_step;
    }
    // the last tick is fixedly set to prevent rounding errors
    draw((plot_start_x + plot_w) as i32, y, to_num_pretty(&end));
}

fn text_width(font: &FontArc, scale: PxScale, s: &str) -> u32 {
    let sf = font.as_scaled(scale);
    s.chars()
        .map(|ch| sf.h_advance(font.glyph_id(ch)))
        .sum::<f32>() as u32
}

fn draw_lines(
    canvas: &mut ImageBuffer<Rgb<u8>, Vec<u8>>,
    line_h: u32,
    text_color: Rgb<u8>,
    mut pad_v: u32,
    pad_h: u32,
    lines: &[String],
    line_space: u32,
) {
    let scale = PxScale::from(line_h as f32);
    let font = global_font();

    for line in lines {
        draw_text_mut(
            canvas,
            text_color,
            pad_h as i32,
            pad_v as i32,
            scale,
            font,
            line,
        );

        pad_v += line_space + line_h;
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_sidebar(
    canvas: &mut ImageBuffer<Rgb<u8>, Vec<u8>>,
    pad_h: u32,
    pad_v: u32,
    sq: u32,
    sidebar_w: u32,
    sidebar_h: u32,
    color_fn: fn(MethylationValue) -> [u8; 3],
    stops: &[(f32, &str)], // in fraction format: [0.0,1.0]
    scale: PxScale,
    border_thickness: u32,
    methylation_threshold: Option<MethylationThreshold>,
    plot_side: f32,
) {
    let space_h = (canvas.width() - pad_h - sq - sidebar_w) / 2;

    let start_h = pad_h + sq + space_h;
    let pad_v = (sq - sidebar_h) / 2 + pad_v;

    for y in 0..sidebar_h {
        let f = 1.0 - (y as f32 / sidebar_h as f32);
        let color = color_fn(MethylationValue::from_fraction(f));
        for x in 0..sidebar_w {
            canvas.put_pixel(start_h + x, pad_v + y, Rgb(color));
        }
    }

    for (f, label) in stops {
        let f = 1.0 - f;
        draw_text_mut(
            canvas,
            Rgb([0, 0, 0]),
            (start_h + sidebar_w + (35.0 * (plot_side / 4000.0)).round() as u32) as i32,
            pad_v as i32 + (sidebar_h as f32).mul_add(f, -(scale.y / 2.0)) as i32,
            scale,
            global_font(),
            label,
        );
    }

    if let Some(threshold) = methylation_threshold {
        let f = threshold.as_fraction();
        if (0.0..=1.0).contains(&f) {
            let y = pad_v as i32 + (sidebar_h as f32 * (1.0 - f)) as i32;
            draw_filled_rect_mut(
                canvas,
                Rect::at(start_h as i32, y).of_size(sidebar_w, border_thickness),
                Rgb([0, 0, 0]),
            );
        }
    }
}

fn draw_borders(
    canvas: &mut ImageBuffer<Rgb<u8>, Vec<u8>>,
    pad_h: u32,
    pad_v: u32,
    w: u32,
    h: u32,
    border: u32,
    color: Rgb<u8>,
) {
    let h_start = pad_h;
    let v_start = pad_v;
    let h_end = pad_h + w;
    let v_end = pad_v + h;

    // top horizontal line
    draw_filled_rect_mut(
        canvas,
        Rect::at(h_start as i32, v_start as i32).of_size(w, border),
        color,
    );
    // bottom horizontal line
    draw_filled_rect_mut(
        canvas,
        Rect::at(h_start as i32, v_end as i32).of_size(w, border),
        color,
    );
    // left vertical line
    draw_filled_rect_mut(
        canvas,
        Rect::at(h_start as i32, v_start as i32).of_size(border, h),
        color,
    );
    // right vertical line
    draw_filled_rect_mut(
        canvas,
        Rect::at(h_end as i32, v_start as i32).of_size(border, h),
        color,
    );
}
