use crate::bio::methylation::domain::MethylationValue;
use crate::config::{MethylationThreshold, PerAxis, WordLength};
use crate::util::color::colormap_weighted;
use crate::util::str::truncate_ellipsis;
use std::borrow::Cow;
use std::ops::Range;

pub mod io;
pub mod render;
pub mod scale;

#[derive(Debug, Clone)]
pub struct AxisAnnotation {
    pub filename: String,
    pub filetype: Option<&'static str>,
    pub name: String,
    pub range: Range<usize>,
    pub is_rev_comp: bool,
}

impl AxisAnnotation {
    pub fn new_fasta(
        filename: String,
        name: String,
        range: Range<usize>,
        is_rev_comp: bool,
        multi: bool,
    ) -> Self {
        Self {
            filename,
            filetype: Some(if multi { "Multi FASTA" } else { "FASTA" }),
            name,
            range,
            is_rev_comp,
        }
    }
    #[inline]
    pub fn axis_label(&self) -> Cow<'_, str> {
        if self.is_rev_comp {
            Cow::Owned(format!("{} (rev. comp.)", self.name))
        } else {
            Cow::Borrowed(self.name.as_str())
        }
    }

    pub fn filename_and_type(&self) -> String {
        let filename_trunc = truncate_ellipsis(self.filename.as_str(), 40, true).into_owned();
        if let Some(filetype) = self.filetype {
            format!("{filetype} ({filename_trunc})")
        } else {
            filename_trunc
        }
    }
}

#[derive(Clone)]
pub struct Annotations<'a> {
    pub axes: PerAxis<Cow<'a, AxisAnnotation>>,
    pub word_len: WordLength,
    pub methylation_threshold: Option<MethylationThreshold>,
}

impl<'a> Annotations<'a> {
    #[inline]
    pub fn distinct_axes<T: Into<AxisAnnotation>>(
        fst: T,
        snd: T,
    ) -> PerAxis<Cow<'a, AxisAnnotation>> {
        PerAxis {
            fst: Cow::Owned(fst.into()),
            snd: Cow::Owned(snd.into()),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PlotStyle {
    pub dot_size: u32,
    pub colormap_raw: fn(MethylationValue) -> [u8; 3],
    pub plot_side: u32,
    pub enlarge_small: bool,
}

impl PlotStyle {
    pub fn new_from_options(
        dot_size: Option<u32>,
        plot_side: Option<u32>,
        enlarge_small: Option<bool>,
    ) -> Self {
        Self {
            dot_size: dot_size.unwrap_or(1),
            plot_side: plot_side.unwrap_or(4_000),
            enlarge_small: enlarge_small.unwrap_or(false),
            colormap_raw: colormap_weighted,
        }
    }
}
