use crate::bio::methylation::domain::MethylationValue;

#[inline]
fn srgb_to_linear(c: u8) -> f32 {
    let s = f32::from(c) / 255.0;
    if s <= 0.04045 {
        s / 12.92
    } else {
        ((s + 0.055) / 1.055).powf(2.4)
    }
}

#[inline]
fn linear_to_srgb(c: f32) -> u8 {
    let s = if c <= 0.003_130_8 {
        c * 12.92
    } else {
        1.055f32.mul_add(c.powf(1.0 / 2.4), -0.055)
    };
    (s * 255.0).round().clamp(0.0, 255.0) as u8
}

#[inline]
pub fn lerp_linear_channel(a: u8, b: u8, t: f32) -> u8 {
    let la = srgb_to_linear(a);
    let lb = srgb_to_linear(b);
    let lc = (lb - la).mul_add(t, la);
    linear_to_srgb(lc)
}

#[inline]
pub fn lerp_rgb_linear(a: [u8; 3], b: [u8; 3], t: f32) -> [u8; 3] {
    [
        lerp_linear_channel(a[0], b[0], t),
        lerp_linear_channel(a[1], b[1], t),
        lerp_linear_channel(a[2], b[2], t),
    ]
}

const fn lighten_rgb(rgb: [u8; 3], factor: f32) -> [u8; 3] {
    let [r, g, b] = rgb;
    let f = factor.clamp(0.0, 1.0);
    [
        (r as f32 + (255.0 - r as f32) * f) as u8,
        (g as f32 + (255.0 - g as f32) * f) as u8,
        (b as f32 + (255.0 - b as f32) * f) as u8,
    ]
}

pub const BLUE: [u8; 3] = [4, 55, 242];
pub const PURPLE: [u8; 3] = [127, 0, 255];
pub const RED: [u8; 3] = [255, 0, 0];
pub const ORAN: [u8; 3] = [255, 68, 51];
const LIGHTER: f32 = 0.5;
const STOPS: [[u8; 3]; 3] = [
    lighten_rgb(BLUE, LIGHTER),
    lighten_rgb(PURPLE, LIGHTER),
    lighten_rgb(RED, LIGHTER),
];
const W: [f32; STOPS.len()] = [0.0, 0.3, 0.5];

#[allow(clippy::const_is_empty)]
/// # Panics
/// Only if `W`/`STOPS` are empty or have mismatched lengths (internal invariant violation).
pub fn colormap_weighted(v: MethylationValue) -> [u8; 3] {
    debug_assert!(!W.is_empty());
    debug_assert!(!STOPS.is_empty());

    let v = v.as_fraction();
    let v = v.clamp(0.0, 1.0);

    for i in 0..W.len() - 1 {
        let start = W[i];
        let end = W[i + 1];
        if v >= start && v <= end {
            let t = (v - start) / (end - start);
            return lerp_rgb_linear(STOPS[i], STOPS[i + 1], t);
        }
    }

    *STOPS.last().unwrap()
}
