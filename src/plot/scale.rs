#[derive(Debug, Clone)]
pub struct TileScale {
    pub tile_w: usize,
    pub tile_h: usize,
    pub tiles_x: usize,
    pub tiles_y: usize,
    pub out_w: usize,
    pub out_h: usize,
}

impl TileScale {
    fn new(original_w: usize, original_h: usize, out_w: usize, out_h: usize) -> Self {
        let tile_w = (original_w as f64 / out_w as f64).ceil() as usize;
        let tile_h = (original_h as f64 / out_h as f64).ceil() as usize;

        let tiles_x = original_w.div_ceil(tile_w).max(1);
        let tiles_y = original_h.div_ceil(tile_h).max(1);

        Self {
            tile_w,
            tile_h,
            tiles_x,
            tiles_y,
            out_w: tiles_x,
            out_h: tiles_y,
        }
    }

    pub fn new_with_max_side(max_side: usize, original_w: usize, original_h: usize) -> Self {
        if original_w < max_side && original_h < max_side {
            return Self::new(original_w, original_h, original_w, original_h);
        }

        let (out_w, out_h) = {
            let max_side = max_side as f64;
            let (w, h) = (original_w as f64, original_h as f64);
            if w >= h {
                (max_side as usize, (h * max_side / w).ceil() as usize)
            } else {
                ((w * max_side / h).ceil() as usize, max_side as usize)
            }
        };

        Self::new(original_w, original_h, out_w, out_h)
    }
}
