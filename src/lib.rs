#[cfg(test)]
extern crate md5;
extern crate rug;

pub mod internal;

use std::cmp::Ordering;
use rug::{Integer, Rational};
use rug::ops::Pow;
use internal::*;

/// Represents a mathematical expression of the form `s √(n / d)` where
///
/// - `s` is a sign (`-1`, `0`, or `+1`),
/// - `n` is a nonnegative numerator, and
/// - `d` is a positive denominator.
///
/// This can be converted to a floating-point number via `f64::from(…)`.
///
/// Defaults to zero.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SignedSqrtRational(Rational);

impl SignedSqrtRational {
    /// Construct a `SignedSqrtRational` equal to `c √r`.
    #[inline]
    pub fn new(c: Integer, r: Rational) -> Self {
        let sign = Rational::from(ordering_to_i32(c.sign()));
        let radical = Rational::from(c.pow(2)) * r;
        SignedSqrtRational(sign * radical)
    }

    /// Equivalent to `self.cmp(&Self::zero())`.
    #[inline]
    pub fn sign(&self) -> Ordering {
        self.0.sign()
    }

    /// Returns the squared value of the expression.
    #[inline]
    pub fn sq(self) -> Rational {
        self.signed_sq().abs()
    }

    /// Returns the squared value of the expression, with the sign adjusted to
    /// match the sign of the original expression.
    #[inline]
    pub fn signed_sq(self) -> Rational {
        self.0
    }
}

impl From<SignedSqrtRational> for f32 {
    #[inline]
    fn from(s: SignedSqrtRational) -> Self {
        let sign = ordering_to_i32(s.sign()) as f32;
        let radical = s.sq().to_f32().sqrt();
        sign * radical
    }
}

impl From<SignedSqrtRational> for f64 {
    #[inline]
    fn from(s: SignedSqrtRational) -> Self {
        let sign = ordering_to_i32(s.sign()) as f64;
        let radical = s.sq().to_f64().sqrt();
        sign * radical
    }
}

/// Calculate Clebsch-Gordan coefficients.
///
/// ```text
/// ⟨j1 j2 m1 m2|j1 j2 j12 m12⟩
/// ```
pub fn clebsch_gordan(
    tj1: i32,
    tm1: i32,
    tj2: i32,
    tm2: i32,
    tj12: i32,
    tm12: i32,
) -> SignedSqrtRational
{
    let SignedSqrtRational(z) = wigner_3jm_raw_c(tj1, tm1, tj2, tm2, tj12, -tm12);
    SignedSqrtRational(z * Rational::from(tj12 + 1))
}

/// Calculates Wigner 3-jm symbols.
///
/// ```text
/// ⎛j1 j2 j3⎞
/// ⎝m1 m2 m3⎠
/// ```
pub fn wigner_3jm(
    tj1: i32,
    tm1: i32,
    tj2: i32,
    tm2: i32,
    tj3: i32,
    tm3: i32,
) -> SignedSqrtRational
{
    let sign = Rational::from(phase((tj1 - tj2 - tm3) / 2));
    let SignedSqrtRational(radical) =
        wigner_3jm_raw_c (tj1, tm1, tj2, tm2, tj3, tm3);
    SignedSqrtRational(sign * radical)
}

/// Calculates Wigner 6-j symbols.
///
/// ```text
/// ⎧j1 j2 j3⎫
/// ⎩j4 j5 j6⎭
/// ```
pub fn wigner_6j(
    tj1: i32,
    tj2: i32,
    tj3: i32,
    tj4: i32,
    tj5: i32,
    tj6: i32,
) -> SignedSqrtRational
{
    if triangle_condition(tj1, tj2, tj3) &&
        triangle_condition(tj1, tj5, tj6) &&
        triangle_condition(tj4, tj2, tj6) &&
        triangle_condition(tj4, tj5, tj3)
    {
        wigner_6j_raw(tj1, tj2, tj3, tj4, tj5, tj6)
    } else {
        Default::default()
    }
}

/// Calculates Wigner 9-j symbols.
///
/// ```text
/// ⎧ja jb jc⎫
/// ⎨jd je jf⎬
/// ⎩jg jh ji⎭
/// ```
pub fn wigner_9j(
    tja: i32,
    tjb: i32,
    tjc: i32,
    tjd: i32,
    tje: i32,
    tjf: i32,
    tjg: i32,
    tjh: i32,
    tji: i32,
) -> SignedSqrtRational
{
    if triangle_condition(tja, tjb, tjc) &&
        triangle_condition(tjd, tje, tjf) &&
        triangle_condition(tjg, tjh, tji) &&
        triangle_condition(tja, tjd, tjg) &&
        triangle_condition(tjb, tje, tjh) &&
        triangle_condition(tjc, tjf, tji)
    {
        wigner_9j_raw(tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji)
    } else {
        Default::default()
    }
}
