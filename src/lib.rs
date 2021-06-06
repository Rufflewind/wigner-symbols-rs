extern crate rug;

pub mod internal;
pub mod regge;

use std::cmp::Ordering;
use std::ops::Mul;
use rug::{Integer, Rational};
use rug::ops::Pow;

/// Signed square root of a rational number
///
/// This represents a mathematical expression of the form `s √(n / d)` where
///
/// - `s` is a sign (`-1`, `0`, or `+1`),
/// - `n` is a nonnegative numerator, and
/// - `d` is a positive denominator.
///
/// Internally, it is represented by the rational number `s × n / d`.
///
/// This can be converted to a floating-point number via `f64::from(…)`.
///
/// Defaults to zero.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SignedSqrt(pub Rational);

impl SignedSqrt {
    /// Construct a `SignedSqrt` equal to `c √r`.
    #[inline]
    pub fn new(c: Integer, r: Rational) -> Self {
        let sign = Rational::from(internal::ordering_to_i32(c.cmp0()));
        let radical = Rational::from(c.pow(2)) * r;
        SignedSqrt(sign * radical)
    }

    /// Equivalent to `self.cmp(&Self::from(0))`.
    #[inline]
    pub fn sign(&self) -> Ordering {
        self.0.cmp0()
    }

    /// Returns the square of the expression.
    #[inline]
    pub fn sq(self) -> Rational {
        self.signed_sq().abs()
    }

    /// Returns the square of the expression, but with the sign adjusted to
    /// match the sign of the original expression.
    #[inline]
    pub fn signed_sq(self) -> Rational {
        self.0
    }
}

impl Mul<SignedSqrt> for SignedSqrt {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        SignedSqrt(self.signed_sq() * other.signed_sq())
    }
}

impl Mul<i32> for SignedSqrt {
    type Output = Self;
    fn mul(self, other: i32) -> Self::Output {
        self * Self::from(other)
    }
}

impl Mul<SignedSqrt> for i32 {
    type Output = SignedSqrt;
    fn mul(self, other: SignedSqrt) -> Self::Output {
        SignedSqrt::from(self) * other
    }
}

impl From<i32> for SignedSqrt {
    #[inline]
    fn from(s: i32) -> Self {
        i64::from(s).into()
    }
}

impl From<i64> for SignedSqrt {
    #[inline]
    fn from(s: i64) -> Self {
        Self::new(s.into(), 1.into())
    }
}

impl From<SignedSqrt> for f32 {
    #[inline]
    fn from(s: SignedSqrt) -> Self {
        let sign = internal::ordering_to_i32(s.sign()) as f32;
        let radical = s.sq().to_f32().sqrt();
        sign * radical
    }
}

impl From<SignedSqrt> for f64 {
    #[inline]
    fn from(s: SignedSqrt) -> Self {
        let sign = f64::from(internal::ordering_to_i32(s.sign()));
        let radical = s.sq().to_f64().sqrt();
        sign * radical
    }
}

/// Clebsch-Gordan coefficient
///
/// ```text
/// ⟨j1 j2 m1 m2|j1 j2 j12 m12⟩
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ClebschGordan {
    pub tj1: i32,
    pub tm1: i32,
    pub tj2: i32,
    pub tm2: i32,
    pub tj12: i32,
    pub tm12: i32,
}

impl From<Wigner3jm> for ClebschGordan {
    fn from(this: Wigner3jm) -> Self {
        let Wigner3jm { tj1, tm1, tj2, tm2, tj3, tm3 } = this;
        Self { tj1, tm1, tj2, tm2, tj12: tj3, tm12: -tm3 }
    }
}

impl ClebschGordan {
    pub fn value(self) -> SignedSqrt {
        SignedSqrt((self.tj12 + 1).into())
            * internal::wigner_3jm_raw_c(self.into())
    }
}

/// Wigner 3-jm symbol
///
/// ```text
/// ⎛j1 j2 j3⎞
/// ⎝m1 m2 m3⎠
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Wigner3jm {
    pub tj1: i32,
    pub tm1: i32,
    pub tj2: i32,
    pub tm2: i32,
    pub tj3: i32,
    pub tm3: i32,
}

impl From<ClebschGordan> for Wigner3jm {
    fn from(this: ClebschGordan) -> Self {
        let ClebschGordan { tj1, tm1, tj2, tm2, tj12, tm12 } = this;
        Self { tj1, tm1, tj2, tm2, tj3: tj12, tm3: -tm12 }
    }
}

impl Wigner3jm {
    pub fn value(self) -> SignedSqrt {
        internal::phase((self.tj1 - self.tj2 - self.tm3) / 2)
            * internal::wigner_3jm_raw_c(self)
    }
}

/// Wigner 6-j symbol
///
/// ```text
/// ⎧j1 j2 j3⎫
/// ⎩j4 j5 j6⎭
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Wigner6j {
    pub tj1: i32,
    pub tj2: i32,
    pub tj3: i32,
    pub tj4: i32,
    pub tj5: i32,
    pub tj6: i32,
}

impl Wigner6j {
    pub fn value(self) -> SignedSqrt {
        if internal::triangle_condition(self.tj1, self.tj2, self.tj3) &&
            internal::triangle_condition(self.tj1, self.tj5, self.tj6) &&
            internal::triangle_condition(self.tj4, self.tj2, self.tj6) &&
            internal::triangle_condition(self.tj4, self.tj5, self.tj3)
        {
            internal::wigner_6j_raw(self)
        } else {
            Default::default()
        }
    }
}

/// Wigner 9-j symbol
///
/// ```text
/// ⎧ja jb jc⎫
/// ⎨jd je jf⎬
/// ⎩jg jh ji⎭
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Wigner9j {
    pub tj1: i32,
    pub tj2: i32,
    pub tj3: i32,
    pub tj4: i32,
    pub tj5: i32,
    pub tj6: i32,
    pub tj7: i32,
    pub tj8: i32,
    pub tj9: i32,
}

impl Wigner9j {
    pub fn value(self) -> SignedSqrt {
        if internal::triangle_condition(self.tj1, self.tj2, self.tj3) &&
            internal::triangle_condition(self.tj4, self.tj5, self.tj6) &&
            internal::triangle_condition(self.tj7, self.tj8, self.tj9) &&
            internal::triangle_condition(self.tj1, self.tj4, self.tj7) &&
            internal::triangle_condition(self.tj2, self.tj5, self.tj8) &&
            internal::triangle_condition(self.tj3, self.tj6, self.tj9)
        {
            internal::wigner_9j_raw(self)
        } else {
            Default::default()
        }
    }
}

/// Symmetrized Wigner 12-j symbol of the second kind
///
/// ```text
/// ⎧j1  j2  j3  j4⎫
/// |j5  j6  j7  j8|
/// ⎩j9 j10 j11 j12⎭
/// ```
///
/// See Yutsis et al (1962), page 62, equation (19.3).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Wigner12jSecond {
    pub tj1: i32,
    pub tj2: i32,
    pub tj3: i32,
    pub tj4: i32,
    pub tj5: i32,
    pub tj6: i32,
    pub tj7: i32,
    pub tj8: i32,
    pub tj9: i32,
    pub tj10: i32,
    pub tj11: i32,
    pub tj12: i32,
}

impl Wigner12jSecond {
    pub fn value(self) -> SignedSqrt {
        if internal::triangle_condition(self.tj1, self.tj5, self.tj9) &&
            internal::triangle_condition(self.tj1, self.tj6, self.tj11) &&
            internal::triangle_condition(self.tj2, self.tj7, self.tj9) &&
            internal::triangle_condition(self.tj2, self.tj8, self.tj11) &&
            internal::triangle_condition(self.tj3, self.tj5, self.tj10) &&
            internal::triangle_condition(self.tj3, self.tj6, self.tj12) &&
            internal::triangle_condition(self.tj4, self.tj7, self.tj10) &&
            internal::triangle_condition(self.tj4, self.tj8, self.tj12)
        {
            internal::wigner_12j_second_raw(self)
        } else {
            Default::default()
        }
    }
}
