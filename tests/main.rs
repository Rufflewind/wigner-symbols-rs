extern crate md5;
extern crate permutohedron;
extern crate rug;
extern crate wigner_symbols;

use std::{cmp, fmt, hash};
use std::collections::HashMap;
use std::io::Write;
use rug::Rational;
use wigner_symbols::*;
use wigner_symbols::internal::*;
use wigner_symbols::regge::*;

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

struct RenderValue<'a>(&'a SignedSqrt);

impl<'a> fmt::Display for RenderValue<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let r = self.0.clone().signed_sq();
        write!(f, "{}/{}", r.numer(), r.denom())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Entry<K, V>(pub K, pub V);

impl<K: PartialEq, V> PartialEq for Entry<K, V> {
    fn eq(&self, other: &Entry<K, V>) -> bool {
        self.0.eq(&other.0)
    }
}

impl<K: Eq, V> Eq for Entry<K, V> {}

impl<K: hash::Hash, V> hash::Hash for Entry<K, V> {
    fn hash<H: hash::Hasher>(&self, hasher: &mut H) {
        self.0.hash(hasher)
    }
}

impl<K: PartialOrd, V> PartialOrd for Entry<K, V> {
    fn partial_cmp(&self, other: &Entry<K, V>) -> Option<cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<K: Ord, V> Ord for Entry<K, V> {
    fn cmp(&self, other: &Entry<K, V>) -> cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

#[derive(Debug)]
pub struct TestPermutations<T> {
    pub slice: Vec<T>,
    pub first: bool,
}

impl<T: Clone> TestPermutations<T> {
    pub fn new(slice: &[T]) -> Self {
        Self { slice: slice.to_owned(), first: true }
    }
}

impl<T: Ord + Clone> Iterator for TestPermutations<T> {
    type Item = Vec<Entry<T, usize>>;
    fn next(&mut self) -> Option<Self::Item> {
        use permutohedron::LexicalPermutation;
        if self.first {
            self.first = false;
        } else if !self.slice.next_permutation() {
            return None;
        }
        Some(self.slice.iter().cloned().enumerate()
             .map(|(i, x)| Entry(x, i)).collect())
    }
}

#[test]
fn test_sort2() {
    for start in &[
        [1, 2],
        [1, 1],
    ] {
        for mut xs in TestPermutations::new(start) {
            let (x0, x1) = sort2(xs[0], xs[1]);
            xs.sort();
            assert_eq!(&xs, &[x0, x1]);
        }
    }
}

#[test]
fn test_sort3() {
    for start in &[
        [1, 2, 3],
        [1, 1, 2],
        [1, 1, 1],
    ] {
        for mut xs in TestPermutations::new(start) {
            let (x0, x1, x2) = sort3(xs[0], xs[1], xs[2]);
            xs.sort();
            assert_eq!(&xs, &[x0, x1, x2]);
        }
    }
}

#[test]
fn test_sort4() {
    for start in &[
        [1, 2, 3, 4],
        [1, 1, 2, 3],
        [1, 1, 2, 2],
        [1, 1, 1, 2],
        [1, 1, 1, 1],
    ] {
        for mut xs in TestPermutations::new(start) {
            let (x0, x1, x2, x3) = sort4(xs[0], xs[1], xs[2], xs[3]);
            xs.sort();
            assert_eq!(&xs, &[x0, x1, x2, x3]);
        }
    }
}

#[test]
fn test_clebsch_gordan_and_wigner_3jm() {
    let tj_max = 25;
    let mut f = md5::Context::new();
    get_3tjms(tj_max, &mut |w3jm| {
        let cg = ClebschGordan::from(w3jm);
        let c = cg.value();
        assert_eq!(
            c.clone().signed_sq(),
            Rational::from(
                (w3jm.tj3 + 1)
                    * phase((w3jm.tj1 - w3jm.tj2 - w3jm.tm3) / 2)
            ) * w3jm.value().signed_sq()
        );
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            cg.tj1,
            cg.tm1,
            cg.tj2,
            cg.tm2,
            cg.tj12,
            cg.tm12,
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
    get_6tjs(tj_max, &mut |w6j| {
        let w = w6j.value();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            w6j.tj1, w6j.tj2, w6j.tj3,
            w6j.tj4, w6j.tj5, w6j.tj6,
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
    get_9tjs(tj_max, &mut |w9j| {
        let w = w9j.value();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            w9j.tj1, w9j.tj2, w9j.tj3,
            w9j.tj4, w9j.tj5, w9j.tj6,
            w9j.tj7, w9j.tj8, w9j.tj9,
            RenderValue(&w),
        ).unwrap();
    });
    assert_eq!(&format!("{:x}", f.compute()),
               *lookup(W9J_HASHES, &tj_max).expect("hash not available"));
}

#[test]
fn test_signed_sqrt_rational() {
    assert_eq!(f64::from(SignedSqrt::default()), 0.0);
    assert_eq!(f64::from(SignedSqrt::new(10.into(), (1, 4).into())), 5.0);
}

#[test]
fn test_regge3jm() {
    let tj_max = 25;
    let mut map = HashMap::new();
    let n = CanonicalRegge3jm::len(tj_max);
    let mut vec = vec![Default::default(); n];
    internal::get_3tjms(tj_max, &mut |w3jm| {
        let value = w3jm.value();
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        let canon_value = phase * value;
        assert_eq!(
            *map.entry(regge).or_insert(canon_value.clone()),
            canon_value.clone()
        );
        vec[regge.index()] = (regge, canon_value);
    });
    for (regge, value) in map {
        assert_eq!(vec[regge.index()], (regge, value));
    }
}

#[test]
fn test_regge6j() {
    let tj_max = 15;
    let mut map = HashMap::new();
    let n = CanonicalRegge6j::len(tj_max);
    let mut vec = vec![Default::default(); n];
    internal::get_6tjs(tj_max, &mut |w6j| {
        let value = w6j.value();
        let regge = CanonicalRegge6j::from(w6j);
        assert_eq!(
            *map.entry(regge).or_insert(value.clone()),
            value.clone()
        );
        vec[regge.index()] = (regge, value);
    });
    for (regge, value) in &map {
        assert_eq!(vec[regge.index()], (*regge, value.clone()));
    }
}
