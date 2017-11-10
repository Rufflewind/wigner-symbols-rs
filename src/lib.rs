extern crate rug;

use std::cmp::Ordering;
use rug::{Integer, Rational};
use rug::ops::Pow;

/// Reinterpret ordering as a sign:
///
/// ```text
/// Less => -1
/// Equal => 0
/// Greater => +1
/// ```
#[inline]
pub fn ordering_to_i32(ordering: Ordering) -> i32 {
    match ordering {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

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
        self.0.abs()
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

/// Calculate the binomial coefficient `C(n, k)`.
#[inline]
fn binomial(n: i32, k: i32) -> Integer {
    Integer::from(n).binomial(k as u32)
}

/// Calculate the falling factorial, i.e. the product of the integers `[n, n - k)`.
#[inline]
fn falling_factorial(n: i32, k: i32) -> Integer {
    let mut r = Integer::from(1);
    for i in n - k .. n {
        r *= Integer::from(i + 1);
    }
    r
}

/// Calculate the factorial `n!`.
#[inline]
fn factorial(n: i32) -> Integer {
    let mut i = Integer::default();
    i.assign_factorial(n as u32);
    i
}

#[inline]
fn phase(phi: i32) -> i32 {
    if phi % 2 == 0 {
        1
    } else {
        -1
    }
}

/// Check `|j1 − j2| ≤ j3 ≤ j1 + j2` and `j1 + j2 + j3 ∈ ℤ`.
#[inline]
pub fn triangle_condition(tj1: i32, tj2: i32, tj3: i32) -> bool {
    let d = tj1 + tj2 - tj3;
    d >= 0 && d % 2 == 0 && tj3 - (tj1 - tj2).abs() >= 0
}


/// Calculate a Clebsch-Gordan coefficient:
///
/// ```text
/// ⟨j1 j2 m1 m2|j1 j2 j12 m12⟩
/// ```
#[inline]
pub fn clebsch_gordan(
    tj1: i32,
    tm1: i32,
    tj2: i32,
    tm2: i32,
    tj12: i32,
    tm12: i32,
) -> SignedSqrtRational
{
    let SignedSqrtRational(z) = wigner3j_raw_c(tj1, tm1, tj2, tm2, tj12, -tm12);
    SignedSqrtRational(z * Rational::from(tj12 + 1))
}

/// Calculate a Wigner 3-j symbol:
///
/// ```text
/// ⎛j1 j2 j3⎞
/// ⎝m1 m2 m3⎠
/// ```
#[inline]
pub fn wigner3j(
    tj1: i32,
    tm1: i32,
    tj2: i32,
    tm2: i32,
    tj3: i32,
    tm3: i32,
) -> SignedSqrtRational
{
    let sign = Rational::from(phase((tj1 - tj2 - tm3) / 2));
    let SignedSqrtRational(radical) = wigner3j_raw_c (tj1, tm1, tj2, tm2, tj3, tm3);
    SignedSqrtRational(sign * radical)
}

/// Calculate the Wigner 3-j symbol times `(−1) ^ (j1 − j2 − m3)`.
#[inline]
pub fn wigner3j_raw_c(
    tj1: i32,
    tm1: i32,
    tj2: i32,
    tm2: i32,
    tj3: i32,
    tm3: i32,
) -> SignedSqrtRational
{
    let jm1 = (tj1 + tm1) / 2;
    let jmr1 = (tj1 + tm1) % 2;
    let jm2 = (tj2 + tm2) / 2;
    let jmr2 = (tj2 + tm2) % 2;
    if
        tm1 + tm2 + tm3 == 0 &&
        tm1.abs() <= tj1 &&
        tm2.abs() <= tj2 &&
        tm3.abs() <= tj3 &&
        jmr1 == 0 &&
        jmr2 == 0 &&
        triangle_condition(tj1, tj2, tj3)
    {
        wigner3j_raw(jm1, jm2, tj1, tm1, tj2, tm2, tj3, tm3)
    } else {
        Default::default()
    }
}

/// Calculate the Wigner 3-j symbol times `(−1) ^ (j1 − j2 − m3)`.
/// The selection rules are not checked.
///
/// ```text
/// j1 + m1, j2 + m2,
/// tj1, tm1, tj2, tm2, tj3, tm3
/// ```
pub fn wigner3j_raw(
    jm1: i32,
    jm2: i32,
    tj1: i32,
    tm1: i32,
    tj2: i32,
    tm2: i32,
    tj3: i32,
    tm3: i32,
) -> SignedSqrtRational
{
    let jjj1 = (tj1 - tj2 + tj3) / 2;
    let jjj2 = (tj2 - tj3 + tj1) / 2;
    let jjj3 = (tj3 - tj1 + tj2) / 2;
    let jjj  = (tj1 + tj2 + tj3) / 2 + 1;
    let jsm1 = (tj1 - tm1) / 2;
    let jm3  = (tj3 + tm3) / 2;
    let kmin = sort3(0, tj1 - tj3 + tm2, tj2 - tj3 - tm1).2 / 2;
    let kmax = sort3(jjj2, jsm1, jm2).0;
    let z1 = Rational::from((
        binomial(tj1, jjj1) * binomial(tj2, jjj2) * binomial(tj3, jjj3),
        binomial(tj1, jm1) * binomial(tj2, jm2) * binomial(tj3, jm3),
    )) * triangular_factor_raw(jjj, jjj1, jjj2, jjj3);
    let z2 = if kmin > kmax {
        Integer::default()
    } else {
        let c0 = Integer::from(phase(kmin))
            * binomial(jjj2, kmin)
            * binomial(jjj1, jsm1 - kmin)
            * binomial(jjj3, jm2 - kmin);
        let mut s = c0.clone();
        let mut c = c0;
        for k in kmin + 1 .. kmax + 1 {
            c = c
                * Integer::from(jjj2 - k + 1) / Integer::from(k)
                * Integer::from(jsm1 - k + 1) / Integer::from(jjj1 - (jsm1 - k))
                * Integer::from(jm2  - k + 1) / Integer::from(jjj3 - (jm2  - k));
            s -= &c;
        }
        s
    };
    SignedSqrtRational::new(z2, z1)
}

/// Calculate a Wigner 6-j symbol:
///
/// ```text
/// ⎧j1 j2 j3⎫
/// ⎩j4 j5 j6⎭
/// ```
pub fn wigner6j(
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
        wigner6j_raw(tj1, tj2, tj3, tj4, tj5, tj6)
    } else {
        Default::default()
    }
}

/// Calculate the Wigner 6-j symbol.  The selection rules are not checked.
pub fn wigner6j_raw(
    tja: i32,
    tjb: i32,
    tjc: i32,
    tjd: i32,
    tje: i32,
    tjf: i32,
) -> SignedSqrtRational
{
    let z1 =
        triangular_factor(tja, tje, tjf)
        * triangular_factor(tjd, tjb, tjf)
        * triangular_factor(tjd, tje, tjc)
        / triangular_factor(tja, tjb, tjc);
    let z2 = tetrahedral_sum(tja, tje, tjf, tjd, tjb, tjc);
    SignedSqrtRational::new(z2, z1)
}

/// Calculate a Wigner 9-j symbol:
///
/// ```text
/// ⎧ja jb jc⎫
/// ⎨jd je jf⎬
/// ⎩jg jh ji⎭
/// ```
pub fn wigner9j(
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
        wigner9j_raw(tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji)
    } else {
        Default::default()
    }
}

/// Calculate the Wigner 9-j symbol.  The selection rules are not checked.
pub fn wigner9j_raw(
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
    let tkmin = *[
        (tjh - tjd).abs(),
        (tjb - tjf).abs(),
        (tja - tji).abs(),
    ].iter().max().unwrap();
    let tkmax = *[
        tjh + tjd,
        tjb + tjf,
        tja + tji,
    ].iter().min().unwrap();
    let z2 = (0 .. (tkmax - tkmin) / 2 + 1).map(|i| {
        let tk = tkmin + i * 2;
        Integer::from(minus_one_pow(tk * (tk + 1)))
            * tetrahedral_sum(tja, tjb, tjc, tjf, tji, tk)
            * tetrahedral_sum(tjf, tjd, tje, tjh, tjb, tk)
            * tetrahedral_sum(tjh, tji, tjg, tja, tjd, tk)
    }).sum();
    let z1 =
        triangular_factor(tja, tjb, tjc) *
        triangular_factor(tjd, tje, tjf) *
        triangular_factor(tjg, tjh, tji) *
        triangular_factor(tja, tjd, tjg) *
        triangular_factor(tjb, tje, tjh) *
        triangular_factor(tjc, tjf, tji);
    SignedSqrtRational::new(z2, z1)
}

#[inline]
pub fn sort2<T: Ord>(a: T, b: T) -> (T, T) {
    if b < a {
        (b, a)
    } else {
        (a, b)
    }
}

#[inline]
pub fn sort3<T: Ord>(a: T, b: T, c: T) -> (T, T, T) {
    let (a, b) = sort2(a, b);
    if c < a {
        (c, a, b)
    } else if c < b {
        (a, c, b)
    } else {
        (a, b, c)
    }
}

/// Calculate the triangular factor:
///
/// ```text
/// Δ(j1, j2, j3) = (−j1 + j2 _ j3)! (j1 − j2 + j3)! (j1 + j2 − j3)!
///               / (j1 + j2 + j3 + 1)!
/// ```
///
#[inline]
pub fn triangular_factor(tj1: i32, tj2: i32, tj3: i32) -> Rational {
    let jjja = (tj3 - tj1 + tj2) / 2;
    let jjjb = (tj1 - tj2 + tj3) / 2;
    let jjjc = (tj2 - tj3 + tj1) / 2;
    let jjj = (tj1 + tj2 + tj3) / 2 + 1;
    triangular_factor_raw(jjj, jjja, jjjb, jjjc)
}

/// Calculate `ja! jb! jc! / jd!`.
#[inline]
pub fn triangular_factor_raw(jd: i32, ja: i32, jb: i32, jc: i32) -> Rational {
    let (ju, jv, jw) = sort3(ja, jb, jc);
    Rational::from((
        factorial(ju) * factorial(jv),
        falling_factorial(jd, jd - jw),
    ))
}

/// Calculate the symbol in the paper by L. Wei that is enclosed in square
/// brackets:
///
/// ```text
/// ⎡j11 j12 j13⎤
/// ⎣j21 j22 j23⎦
/// ```
///
/// This is essentially a Wigner 6-j symbol without the triangular factors,
/// although the ordering of the arguments is a bit funky here.
pub fn tetrahedral_sum(
    tja: i32,
    tje: i32,
    tjf: i32,
    tjd: i32,
    tjb: i32,
    tjc: i32,
) -> Integer
{
    let jjja = (tjc - tja + tjb) / 2;
    let jjjb = (tja - tjb + tjc) / 2;
    let jjjc = (tjb - tjc + tja) / 2;
    let jabc = (tja + tjb + tjc) / 2;
    let jaef = (tja + tje + tjf) / 2;
    let jdbf = (tjd + tjb + tjf) / 2;
    let jdec = (tjd + tje + tjc) / 2;
    let kmin = *[jabc, jdec, jdbf, jaef].iter().max().unwrap();
    let kmax = *[
        tja + tjd + tjb + tje,
        tjb + tje + tjc + tjf,
        tja + tjd + tjc + tjf,
    ].iter().max().unwrap() / 2;
    (kmin .. kmax + 1).map(|k| {
        Integer::from(minus_one_pow(k))
            * binomial(k + 1, k - jabc)
            * binomial(jjja, k - jaef)
            * binomial(jjjb, k - jdbf)
            * binomial(jjjc, k - jdec)
    }).sum()
}

#[inline]
pub fn minus_one_pow(e: i32) -> i32 {
    if e % 2 == 0 {
        1
    } else {
        -1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        assert_eq!(f64::from(SignedSqrtRational::default()), 0.0);
        assert_eq!(f64::from(SignedSqrtRational::default()), 0.0);

        assert_eq!(wigner3j(0, 0, 0, 0, 0, 0),
                   SignedSqrtRational::new(1.into(), 1.into()));
        assert_eq!(clebsch_gordan(5, -3, 1, 1, 4, -2),
                   SignedSqrtRational::new((-1).into(), (2, 3).into()));
        // TODO: Add more extensive testing
    }
}
