use core::fmt;
use std::ops::{Add, AddAssign, Div};

#[derive(Clone, Copy, PartialOrd, PartialEq)]
pub struct MethylationValue(f32);

impl fmt::Debug for MethylationValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MethylationValue({}%)", self.0)
    }
}

impl fmt::Display for MethylationValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}%", self.0)
    }
}

impl MethylationValue {
    pub fn from_fraction(f: f32) -> Self {
        MethylationValue(f * 100.0)
    }

    pub const fn from_percent(p: f32) -> Self {
        MethylationValue(p)
    }

    pub const fn as_fraction(&self) -> f32 {
        self.0 / 100.0
    }

    pub const fn as_percent(&self) -> f32 {
        self.0
    }

    pub const fn min(self, other: Self) -> Self {
        Self(self.0.min(other.0))
    }
    pub const fn max(self, other: Self) -> Self {
        Self(self.0.max(other.0))
    }
}

impl AddAssign for MethylationValue {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl Add for MethylationValue {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Div for MethylationValue {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(self.0 / rhs.0)
    }
}
