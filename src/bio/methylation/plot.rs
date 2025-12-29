use crate::bio::methylation::domain::MethylationValue;
use crate::config::PerAxis;
use crate::core::run::Reporter;
use crate::plot::scale::TileScale;
use std::borrow::Cow;
use std::collections::HashMap;

pub struct MethylationCoordinate {
    pub x: usize,
    pub y: usize,
    pub w: MethylationValue,
}

pub fn add_methylation(
    methylations: &PerAxis<&[(usize, MethylationValue)]>,
    scale: &TileScale,
    write_dot: &mut dyn FnMut(MethylationCoordinate),
    ctx: &Reporter,
) {
    let tile_stats = |fst: bool,
                      meth: &[(usize, MethylationValue)],
                      tile_dim: usize,
                      border: usize|
     -> HashMap<usize, (usize, f64)> {
        let ord = if fst { "first" } else { "second" };
        let bar = ctx.create_bar(
            meth.len(),
            Cow::Owned(format!("{ord} methylation: binning dots")),
        );

        let mut tiles = HashMap::new();
        for &(pos, fraction) in meth {
            bar.inc(1);
            let tile = pos / tile_dim;
            if tile >= border {
                println!(
                    "warning: tile coord (axis: {ord}, pos: {pos}) is out of bounds ({tile} >= {border})"
                );
                continue;
            }

            let entry = tiles.entry(tile).or_insert((0, 0.0));
            entry.0 += 1;
            entry.1 += f64::from(fraction.as_percent());
        }
        bar.finish();

        tiles
    };

    let PerAxis {
        fst: fst_methylation,
        snd: snd_methylation,
    } = methylations;

    let fst_tiles = tile_stats(true, fst_methylation, scale.tile_w, scale.out_w);
    let snd_tiles = tile_stats(false, snd_methylation, scale.tile_h, scale.out_h);

    let bar = ctx.create_bar(fst_tiles.len() * snd_tiles.len(), "averaging dots");

    for (&tx, &(c1, s1)) in &fst_tiles {
        for (&ty, &(c2, s2)) in &snd_tiles {
            bar.inc(1);

            if c1 == 0 || c2 == 0 {
                continue;
            }

            let avg = (s1 + s2) / (c1 + c2) as f64;
            let avg = MethylationValue::from_percent(avg as f32);

            write_dot(MethylationCoordinate {
                x: tx,
                y: ty,
                w: avg,
            });
        }
    }
    bar.finish();
}
