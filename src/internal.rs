//! Contents of this module are subject to change.

use std::cmp::Ordering;
use std::ops::Range;
use rug::{Integer, Rational};
use super::SignedSqrtRational;

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
    let mut i = Integer::default();
    i.assign_factorial(n as u32);
    i
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
#[inline]
pub fn wigner_3jm_raw_c(
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
        wigner_3jm_raw(jm1, jm2, tj1, tm1, tj2, tm2, tj3, tm3)
    } else {
        Default::default()
    }
}

/// Calculate the Wigner 3-jm symbol times `(−1) ^ (j1 − j2 − m3)`.
/// The selection rules are not checked.
///
/// ```text
/// j1 + m1, j2 + m2,
/// tj1, tm1, tj2, tm2, tj3, tm3
/// ```
pub fn wigner_3jm_raw(
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
            c = -c
                * Integer::from(jjj2 - k + 1) / Integer::from(k)
                * Integer::from(jsm1 - k + 1) / Integer::from(jjj1 - (jsm1 - k))
                * Integer::from(jm2 - k + 1) / Integer::from(jjj3 - (jm2  - k));
            s += &c;
        }
        s
    };
    SignedSqrtRational::new(z2, z1)
}

/// Calculate the Wigner 6-j symbol.  The selection rules are not checked.
pub fn wigner_6j_raw(
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

/// Calculate the Wigner 9-j symbol.  The selection rules are not checked.
pub fn wigner_9j_raw(
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
    let tkmin = sort3(
        (tjh - tjd).abs(),
        (tjb - tjf).abs(),
        (tja - tji).abs(),
    ).2;
    let tkmax = sort3(
        tjh + tjd,
        tjb + tjf,
        tja + tji,
    ).0;
    let z2 = (0 .. (tkmax - tkmin) / 2 + 1).map(|i| {
        let tk = tkmin + i * 2;
        Integer::from(phase(tk) * (tk + 1))
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
        Integer::from(phase(k))
            * binomial(k + 1, k - jabc)
            * binomial(jjja, k - jaef)
            * binomial(jjjb, k - jdbf)
            * binomial(jjjc, k - jdec)
    }).sum()
}

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
    callback: &mut FnMut(i32, i32, i32, i32, i32, i32),
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
        callback(tj1, tm1, tj2, tm2, tj3, tm3);
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
    callback: &mut FnMut(i32, i32, i32, i32, i32, i32),
) {
    for tja in 0 .. tj_max + 1 {
    for tjb in 0 .. tj_max + 1 {
    for tjc in get_triangular_tjs(tj_max, tja, tjb) {
    for tjd in 0 .. tj_max + 1 {
    for tje in get_triangular_tjs(tj_max, tjd, tjc) {
    for tjf in get_bitriangular_tjs(tj_max, tja, tje, tjd, tjb) {
        callback(tja, tjb, tjc, tjd, tje, tjf);
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
    callback: &mut FnMut(i32, i32, i32, i32, i32, i32, i32, i32, i32),
) {
    for tja in 0 .. tj_max + 1 {
    for tjb in 0 .. tj_max + 1 {
    for tjc in get_triangular_tjs(tj_max, tja, tjb) {
    for tjd in 0 .. tj_max + 1 {
    for tje in 0 .. tj_max + 1 {
    for tjf in get_triangular_tjs(tj_max, tjd, tje) {
    for tjg in get_triangular_tjs(tj_max, tja, tjd) {
    for tjh in get_triangular_tjs(tj_max, tjb, tje) {
    for tji in get_bitriangular_tjs(tj_max, tjg, tjh, tjc, tjf) {
        callback(tja, tjb, tjc, tjd, tje, tjf, tjg, tjh, tji);
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

#[cfg(test)]
mod tests {
    use std::fmt;
    use std::io::Write;
    use md5;
    use super::*;
    use super::super::*;

    const CG_HASHES: &[(i32, &str)] = &[
        (5, "e74c501299b456a6cb29e4f5714e9061"), // 681
        (10, "b6d0770101f4ebdaa9a55d94f07b001f"), // 11487
        (15, "9192023f26dae0eebcce11afa7372eb6"), // 69272
        (20, "75ef56391b61e1bb2336e36ac7834216"), // 259523
        (25, "5901128892a264b73b5479b70b331fd0"), // 737113
        (30, "75ef56391b61e1bb2336e36ac7834216"), // 1747984
        (40, "2f9b936ea977249c1fea8a22d190a4cf"), // 6931995
    ];

    const W6J_HASHES: &[(i32, &str)] = &[
        (5, "26c24e568fc96f1732ebb3130a46f22a"), // 1479
        (10, "f892f4b466e0558179ca870941d0a456"), // 42393
        (15, "f50b0163194cef1699727b7064760ec0"), // 363196
        (20, "e1b5dad0f1469cc54b6139533f982815"), // 1766270
        (25, "f326bf6e12a94120d2f46582e95e92f8"), // 6171698
    ];

    const W9J_HASHES: &[(i32, &str)] = &[
        (3, "4005ef20e2ed8c789917dce99d027bc4"), // 1616
        (4, "92cfc13320e7fd6a34b3970ebef58e06"), // 9060
        (5, "d596fa3960aafae148754b6f3274507d"), // 38031
        (7, "7b338708ef3aa4ba0a4f5bd8c8b4e6aa"), // 401899
        (10, "479c0a020eaceff5539e2dda2200c1ab"), // 5898846
    ];

    fn lookup<'a, K: Eq, V>(table: &'a [(K, V)], key: &K) -> Option<&'a V> {
        table.iter().find(|&&(ref k, _)| k == key).map(|x| &x.1)
    }

    struct RenderValue<'a>(&'a SignedSqrtRational);

    impl<'a> fmt::Display for RenderValue<'a> {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let r = self.0.clone().signed_sq();
            write!(f, "{}/{}", r.numer(), r.denom())
        }
    }

    #[test]
    fn test_clebsch_gordan_and_wigner_3jm() {
        let tj_max = 25;
        let mut f = md5::Context::new();
        get_3tjms(tj_max, &mut |tj1, tm1, tj2, tm2, tj3, tm3| {
            let c = clebsch_gordan(tj1, tm1, tj2, tm2, tj3, -tm3);
            let w = wigner_3jm(tj1, tm1, tj2, tm2, tj3, tm3);
            assert_eq!(
                c.clone().signed_sq(),
                Rational::from((tj3 + 1) * phase((tj1 - tj2 - tm3) / 2))
                    * w.signed_sq()
            );
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                tj1, tm1, tj2, tm2, tj3, -tm3,
                RenderValue(&c),
            ).unwrap();
        });
        assert_eq!(&format!("{:x}", f.compute()),
                   *lookup(CG_HASHES, &tj_max).expect("hash not available"));
    }

    #[test]
    fn test_wigner_6j() {
        let tj_max = 15;
        let mut f = md5::Context::new();
        get_6tjs(tj_max, &mut |tj1, tj2, tj3, tj4, tj5, tj6| {
            let w = wigner_6j(tj1, tj2, tj3, tj4, tj5, tj6);
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                tj1, tj2, tj3, tj4, tj5, tj6,
                RenderValue(&w),
            ).unwrap();
        });
        assert_eq!(&format!("{:x}", f.compute()),
                   *lookup(W6J_HASHES, &tj_max).expect("hash not available"));
    }

    #[test]
    fn test_wigner_9j() {
        let tj_max = 7;
        let mut f = md5::Context::new();
        get_9tjs(tj_max, &mut |tj1, tj2, tj3, tj4, tj5, tj6, tj7, tj8, tj9| {
            let w = wigner_9j(tj1, tj2, tj3, tj4, tj5, tj6, tj7, tj8, tj9);
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                tj1, tj2, tj3, tj4, tj5, tj6, tj7, tj8, tj9,
                RenderValue(&w),
            ).unwrap();
        });
        assert_eq!(&format!("{:x}", f.compute()),
                   *lookup(W9J_HASHES, &tj_max).expect("hash not available"));
    }

    #[test]
    fn test_signed_sqrt_rational() {
        assert_eq!(f64::from(SignedSqrtRational::default()), 0.0);
        assert_eq!(f64::from(SignedSqrtRational::new(
            10.into(),
            (1, 4).into(),
        )), 5.0);
    }

    #[test]
    fn test_sort3() {
        assert_eq!(sort3(1, 2, 3), (1, 2, 3));
        assert_eq!(sort3(2, 1, 3), (1, 2, 3));
        assert_eq!(sort3(2, 3, 1), (1, 2, 3));
        assert_eq!(sort3(3, 2, 1), (1, 2, 3));
        assert_eq!(sort3(3, 1, 2), (1, 2, 3));
        assert_eq!(sort3(1, 3, 2), (1, 2, 3));
    }
}
