//! Contents of this module are subject to change.

use std::cmp::Ordering;
use std::ops::Range;
use rug::{Integer, Rational};
use super::{SignedSqrt, Wigner3jm, Wigner6j, Wigner9j};

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

#[inline]
pub fn sort4<T: Ord>(a: T, b: T, c: T, d: T) -> (T, T, T, T) {
    let (a, b) = sort2(a, b);
    let (c, d) = sort2(c, d);
    if c < a {
        if d < a {
            (c, d, a, b)
        } else if d < b {
            (c, a, d, b)
        } else {
            (c, a, b, d)
        }
    } else if c < b {
        if d < b {
            (a, c, d, b)
        } else {
            (a, c, b, d)
        }
    } else {
        (a, b, c, d)
    }
}

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

/// Calculate the binomial coefficient `C(n, k)`.
#[inline]
pub fn binomial(n: i32, k: i32) -> Integer {
    Integer::from(n).binomial(k as u32)
}

/// Calculate the falling factorial, i.e. the product of the integers `[n, n - k)`.
#[inline]
pub fn falling_factorial(n: i32, k: i32) -> Integer {
    let mut r = Integer::from(1);
    for i in n - k .. n {
        r *= Integer::from(i + 1);
    }
    r
}

/// Calculate the factorial `n!`.
#[inline]
pub fn factorial(n: i32) -> Integer {
    Integer::factorial(n as u32).into()
}

#[inline]
pub fn phase(phi: i32) -> i32 {
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

/// Calculate the Wigner 3-jm symbol times `(−1) ^ (j1 − j2 − m3)`.
pub fn wigner_3jm_raw_c(this: Wigner3jm) -> SignedSqrt {
    let Wigner3jm { tj1, tm1, tj2, tm2, tj3, tm3 } = this;
    let jmr1 = (tj1 + tm1) % 2;
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
        wigner_3jm_raw(this)
    } else {
        Default::default()
    }
}

/// Calculate the Wigner 3-jm symbol times `(−1) ^ (j1 − j2 − m3)`.
/// The selection rules are not checked.
pub fn wigner_3jm_raw(this: Wigner3jm) -> SignedSqrt {
    let Wigner3jm { tj1, tm1, tj2, tm2, tj3, tm3 } = this;
    let jjj1 = (tj1 - tj2 + tj3) / 2;
    let jjj2 = (tj2 - tj3 + tj1) / 2;
    let jjj3 = (tj3 - tj1 + tj2) / 2;
    let jjj  = (tj1 + tj2 + tj3) / 2 + 1;
    let jm1 = (tj1 + tm1) / 2;
    let jm2 = (tj2 + tm2) / 2;
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
            c = -c
                * Integer::from(jjj2 - k + 1) / Integer::from(k)
                * Integer::from(jsm1 - k + 1) / Integer::from(jjj1 - (jsm1 - k))
                * Integer::from(jm2 - k + 1) / Integer::from(jjj3 - (jm2  - k));
            s += &c;
        }
        s
    };
    SignedSqrt::new(z2, z1)
}

/// Calculate the Wigner 6-j symbol.  The selection rules are not checked.
pub fn wigner_6j_raw(this: Wigner6j) -> SignedSqrt {
    let Wigner6j { tj1, tj2, tj3, tj4, tj5, tj6 } = this;
    let z1 =
        triangular_factor(tj1, tj5, tj6)
        * triangular_factor(tj4, tj2, tj6)
        * triangular_factor(tj4, tj5, tj3)
        / triangular_factor(tj1, tj2, tj3);
    let z2 = tetrahedral_sum(tj1, tj5, tj6, tj4, tj2, tj3);
    SignedSqrt::new(z2, z1)
}

/// Calculate the Wigner 9-j symbol.  The selection rules are not checked.
pub fn wigner_9j_raw(this: Wigner9j) -> SignedSqrt {
    let Wigner9j { tj1, tj2, tj3, tj4, tj5, tj6, tj7, tj8, tj9 } = this;
    let tkmin = sort3(
        (tj8 - tj4).abs(),
        (tj2 - tj6).abs(),
        (tj1 - tj9).abs(),
    ).2;
    let tkmax = sort3(
        tj8 + tj4,
        tj2 + tj6,
        tj1 + tj9,
    ).0;
    let z2 = (0 .. (tkmax - tkmin) / 2 + 1).map(|i| {
        let tk = tkmin + i * 2;
        Integer::from(phase(tk) * (tk + 1))
            * tetrahedral_sum(tj1, tj2, tj3, tj6, tj9, tk)
            * tetrahedral_sum(tj6, tj4, tj5, tj8, tj2, tk)
            * tetrahedral_sum(tj8, tj9, tj7, tj1, tj4, tk)
    }).sum();
    let z1 =
        triangular_factor(tj1, tj2, tj3) *
        triangular_factor(tj4, tj5, tj6) *
        triangular_factor(tj7, tj8, tj9) *
        triangular_factor(tj1, tj4, tj7) *
        triangular_factor(tj2, tj5, tj8) *
        triangular_factor(tj3, tj6, tj9);
    SignedSqrt::new(z2, z1)
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
        Integer::from(phase(k))
            * binomial(k + 1, k - jabc)
            * binomial(jjja, k - jaef)
            * binomial(jjjb, k - jdbf)
            * binomial(jjjc, k - jdec)
    }).sum()
}

#[inline]
pub fn intersect_ranges(a: Range<i32>, b: Range<i32>) -> Range<i32> {
    a.start.max(b.start) .. a.end.min(b.end)
}

pub struct Step<I> {
    pub iter: I,
    pub step: usize,
}

impl<I: Iterator> Iterator for Step<I> {
    type Item = I::Item;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let item = self.iter.next();
        if item.is_some() && self.step >= 2 {
            self.iter.nth(self.step - 2);
        }
        item
    }
}

#[inline]
pub fn get_triangular_tjs(tj_max: i32, tj1: i32, tj2: i32) -> Step<Range<i32>> {
    Step {
        iter: (tj1 - tj2).abs() .. tj_max.min(tj1 + tj2) + 1,
        step: 2,
    }
}

#[inline]
pub fn get_bitriangular_tjs(
    tj_max: i32,
    tj1: i32,
    tj2: i32,
    tj3: i32,
    tj4: i32,
) -> Step<Range<i32>>
{
    Step {
        iter: if (tj1 + tj2 + tj3 + tj4) % 2 != 0 {
            0 .. 0
        } else {
            intersect_ranges(
                get_triangular_tjs(tj_max, tj1, tj2).iter,
                get_triangular_tjs(tj_max, tj3, tj4).iter,
            )
        },
        step: 2,
    }
}

/// Get all projection quantum numbers that lie within the multiplet of `j`.
#[inline]
pub fn get_tms(tj: i32) -> Step<Range<i32>> {
    Step { iter: -tj .. tj + 1, step: 2 }
}

/// Get all possible arguments of the Wigner 3-j symbol that satisfy the
/// selection rules up to a maximum of `j_max`.
pub fn get_3tjms(
    tj_max: i32,
    callback: &mut FnMut(Wigner3jm),
) {
    for tj1 in 0 .. tj_max + 1 {
    for tj2 in 0 .. tj_max + 1 {
    for tj3 in get_triangular_tjs(tj_max, tj1, tj2) {
    for tm1 in get_tms(tj1) {
    for tm2 in get_tms(tj2) {
        let tm3 = -(tm1 + tm2);
        if tm3.abs() > tj3 {
            continue;
        }
        callback(Wigner3jm { tj1, tm1, tj2, tm2, tj3, tm3 });
    }
    }
    }
    }
    }
}

/// Get all possible arguments of the Wigner 6-j symbol that satisfy the
/// selection rules up to a maximum of `j_max`.
pub fn get_6tjs(
    tj_max: i32,
    callback: &mut FnMut(Wigner6j),
) {
    for tj1 in 0 .. tj_max + 1 {
    for tj2 in 0 .. tj_max + 1 {
    for tj3 in get_triangular_tjs(tj_max, tj1, tj2) {
    for tj4 in 0 .. tj_max + 1 {
    for tj5 in get_triangular_tjs(tj_max, tj4, tj3) {
    for tj6 in get_bitriangular_tjs(tj_max, tj1, tj5, tj4, tj2) {
        callback(Wigner6j { tj1, tj2, tj3, tj4, tj5, tj6 });
    }
    }
    }
    }
    }
    }
}

/// Get all possible arguments of the Wigner 9-j symbol that satisfy the
/// selection rules up to a maximum of `j_max`.
pub fn get_9tjs(
    tj_max: i32,
    callback: &mut FnMut(Wigner9j),
) {
    for tj1 in 0 .. tj_max + 1 {
    for tj2 in 0 .. tj_max + 1 {
    for tj3 in get_triangular_tjs(tj_max, tj1, tj2) {
    for tj4 in 0 .. tj_max + 1 {
    for tj5 in 0 .. tj_max + 1 {
    for tj6 in get_triangular_tjs(tj_max, tj4, tj5) {
    for tj7 in get_triangular_tjs(tj_max, tj1, tj4) {
    for tj8 in get_triangular_tjs(tj_max, tj2, tj5) {
    for tj9 in get_bitriangular_tjs(tj_max, tj7, tj8, tj3, tj6) {
        callback(Wigner9j { tj1, tj2, tj3, tj4, tj5, tj6, tj7, tj8, tj9 });
    }
    }
    }
    }
    }
    }
    }
    }
    }
}
