use std::ops::{Index, IndexMut};
use super::{Wigner3jm, Wigner6j};
use super::internal::{phase, sort3, sort4};

/// Regge square for Wigner 3-jm symbols, arranged in row-major order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Regge3jm(pub [i32; 9]);

impl Index<(usize, usize)> for Regge3jm {
    type Output = i32;
    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.0[i * 3 + j]
    }
}

impl IndexMut<(usize, usize)> for Regge3jm {
    #[inline]
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        &mut self.0[i * 3 + j]
    }
}

impl From<Wigner3jm> for Regge3jm {
    #[inline]
    fn from(Wigner3jm { tj1, tm1, tj2, tm2, tj3, tm3 }: Wigner3jm) -> Self {
        Regge3jm([
            (-tj1 + tj2 + tj3) / 2,
            (tj1 - tj2 + tj3) / 2,
            (tj1 + tj2 - tj3) / 2,
            (tj1 - tm1) / 2,
            (tj2 - tm2) / 2,
            (tj3 - tm3) / 2,
            (tj1 + tm1) / 2,
            (tj2 + tm2) / 2,
            (tj3 + tm3) / 2,
        ])
    }
}

impl Regge3jm {
    #[inline]
    pub fn magic_sum(&self) -> i32 {
        self[(0, 0)] + self[(0, 1)] + self[(0, 2)]
    }

    #[inline]
    pub fn swap(&mut self, (i1, j1): (usize, usize), (i2, j2): (usize, usize)) {
        let x = self[(i1, j1)];
        self[(i1, j1)] = self[(i2, j2)];
        self[(i2, j2)] = x;
    }

    #[inline]
    pub fn swap_rows(&mut self, i1: usize, i2: usize, parity: &mut bool) {
        for j in 0 .. 3 {
            self.swap((i1, j), (i2, j));
        }
        *parity ^= i1 != i2;
    }

    #[inline]
    pub fn swap_cols(&mut self, j1: usize, j2: usize, parity: &mut bool) {
        for i in 0 .. 3 {
            self.swap((i, j1), (i, j2));
        }
        *parity ^= j1 != j2;
    }

    #[inline]
    pub fn transpose(&mut self) {
        self.swap((0, 1), (1, 0));
        self.swap((0, 2), (2, 0));
        self.swap((1, 2), (2, 1));
    }

    /// Canonicalize the Regge square using the ordering specified in [Tuzun
    /// et al (1998)](https://doi.org/10.1016/S0010-4655(98)00065-4).
    #[inline]
    pub fn canonicalize(&mut self) -> (CanonicalRegge3jm, i32) {
        let magic_sum = self.magic_sum();
        let parity = &mut false;

        let mut min_i = 0;
        let mut min_j = 0;
        let mut s = self[(0, 0)];
        let mut l = self[(0, 0)];
        for i in 0 .. 3 {
            for j in 0 .. 3 {
                if self[(i, j)] < s {
                    s = self[(i, j)];
                    min_i = i;
                    min_j = j;
                }
                if self[(i, j)] > l {
                    l = self[(i, j)];
                }
            }
        }

        // move S to (0, 0)
        self.swap_rows(0, min_i, parity);
        self.swap_cols(0, min_j, parity);

        // move L to (0, 1)
        if self[(0, 2)] == l {
            self.swap_cols(1, 2, parity);
        } else if self[(1, 0)] == l {
            self.transpose();
        } else if self[(2, 0)] == l {
            self.transpose();
            self.swap_cols(1, 2, parity);
        }

        // canonicalize last two rows
        if (self[(1, 1)], self[(1, 2)]) > (self[(2, 1)], self[(2, 2)]) {
            self.swap_rows(1, 2, parity);
        }

        (CanonicalRegge3jm {
            l: self[(0, 1)] as _,
            x: self[(1, 0)] as _,
            t: self[(2, 2)] as _,
            b: self[(1, 1)] as _,
            s: self[(0, 0)] as _,
        }, if *parity {
            phase(magic_sum)
        } else {
            1
        })
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct CanonicalRegge3jm {
    pub l: u8,
    pub x: u8,
    pub t: u8,
    pub b: u8,
    pub s: u8,
}

impl CanonicalRegge3jm {
    /// Index into a table ordered according to [Rasch and Yu
    /// (2003)](https://doi.org/10.1137/S1064827503422932).
    #[inline]
    pub fn index(self) -> usize {
        let l = self.l as usize;
        let x = self.x as usize;
        let t = self.t as usize;
        let b = self.b as usize;
        let s = self.s as usize;
        debug_assert!(l >= x);
        debug_assert!(x >= t);
        debug_assert!(t >= b);
        debug_assert!(b >= s);
        l * (24 + l * (50 + l * (35 + l * (10 + l)))) / 120
            + x * (6 + x * (11 + x * (6 + x))) / 24
            + t * (2 + t * (3 + t)) / 6
            + b * (b + 1) / 2
            + s
    }

    #[inline]
    pub fn len(tj_max: i32) -> usize {
        Self {
            l: tj_max as u8 + 1,        // j_max = j_max + j_max - 0
            .. Default::default()
        }.index()
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct CanonicalRegge6j {
    pub e: u8,
    pub l: u8,
    pub x: u8,
    pub t: u8,
    pub b: u8,
    pub s: u8,
}

impl From<Wigner6j> for CanonicalRegge6j {
    #[inline]
    fn from(Wigner6j { tj1, tj2, tj3, tj4, tj5, tj6 }: Wigner6j) -> Self {
        // there's a typo in Rasch and Yu (2003):
        // (3.12) should say α3 ≥ α2 ≥ α1
        let (alpha1, alpha2, alpha3) = sort3(
            (tj1 + tj2 + tj4 + tj5) / 2,
            (tj1 + tj3 + tj4 + tj6) / 2,
            (tj2 + tj3 + tj5 + tj6) / 2,
        );
        let (beta4, beta3, beta2, beta1) = sort4(
            (tj1 + tj2 + tj3) / 2,
            (tj1 + tj5 + tj6) / 2,
            (tj2 + tj4 + tj6) / 2,
            (tj3 + tj4 + tj5) / 2,
        );
        Self {
            s: (alpha1 - beta1) as _,
            b: (alpha1 - beta2) as _,
            t: (alpha1 - beta3) as _,
            x: (alpha1 - beta4) as _,
            l: (alpha2 - beta4) as _,
            e: (alpha3 - beta4) as _,
        }
    }
}

impl CanonicalRegge6j {
    /// Index into a table ordered according to [Rasch and Yu
    /// (2003)](https://doi.org/10.1137/S1064827503422932).
    #[inline]
    pub fn index(self) -> usize {
        debug_assert!(self.e >= self.l);
        let e = self.e as usize;
        e * (120 + e * (274 + e * (225 + e * (85 + e * (15 + e))))) / 720
            + CanonicalRegge3jm {
                l: self.l,
                x: self.x,
                b: self.b,
                t: self.t,
                s: self.s,
            }.index()
    }

    #[inline]
    pub fn len(tj_max: i32) -> usize {
        Self {
            // α − β will always cancel at least two of the j's, with two j's
            // remaining
            e: tj_max as u8 + 1,
            .. Default::default()
        }.index()
    }
}
