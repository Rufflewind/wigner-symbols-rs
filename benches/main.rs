#![feature(test)]

extern crate fnv;
extern crate rand;
extern crate test;
extern crate wigner_symbols;

use fnv::FnvHashMap;
use rand::distributions::IndependentSample;
use test::Bencher;
use wigner_symbols::{internal, ClebschGordan, Wigner3jm};
use wigner_symbols::regge::{CanonicalRegge3jm, Regge3jm};

fn wigner_3jm_arg_table(tj_min: i32, tj_max: i32) -> Vec<Wigner3jm> {
    let mut table = Vec::default();
    internal::get_3tjms(tj_max, &mut |w3jm| {
        if internal::sort3(w3jm.tj1, w3jm.tj2, w3jm.tj3).2 >= tj_min {
            table.push(w3jm);
        }
    });
    table
}

fn clebsch_gordan_bench(b: &mut Bencher, tj_min: i32, tj_max: i32) {
    let table = wigner_3jm_arg_table(tj_min, tj_max);
    let range = rand::distributions::Range::new(0, table.len());
    let mut rng = rand::weak_rng();
    b.iter(|| {
        let w3jm = table[range.ind_sample(&mut rng)];
        test::black_box(ClebschGordan::from(w3jm).value());
    });
}

#[bench]
fn bench_clebsch_gordan_00_05(b: &mut Bencher) {
    clebsch_gordan_bench(b, 0, 5);
}

#[bench]
fn bench_clebsch_gordan_05_10(b: &mut Bencher) {
    clebsch_gordan_bench(b, 5, 10);
}

#[bench]
fn bench_clebsch_gordan_10_15(b: &mut Bencher) {
    clebsch_gordan_bench(b, 10, 15);
}

#[bench]
fn bench_clebsch_gordan_15_20(b: &mut Bencher) {
    clebsch_gordan_bench(b, 15, 20);
}

#[bench]
fn bench_clebsch_gordan_20_25(b: &mut Bencher) {
    clebsch_gordan_bench(b, 20, 25);
}

// this one just uses a plain FnvHashMap, which requires storing both keys and
// values
#[bench]
fn bench_regge3jm_map_25(b: &mut Bencher) {
    let tj_max = 25;
    let table = wigner_3jm_arg_table(0, tj_max);

    let mut map = FnvHashMap::default();
    for &w3jm in &table {
        let value = w3jm.value();
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        map.entry(regge).or_insert(f64::from(phase * value));
    }

    let range = rand::distributions::Range::new(0, table.len());
    let mut rng = rand::weak_rng();
    b.iter(|| {
        let w3jm = table[range.ind_sample(&mut rng)];
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        test::black_box(phase as f64 * map.get(&regge).unwrap());
    });
}

// this one just uses the ordering scheme of Rasch and Yu, which only requires
// storing values but has a density of about 50%
#[bench]
fn bench_regge3jm_vec_25(b: &mut Bencher) {
    let tj_max = 25;
    let table = wigner_3jm_arg_table(0, tj_max);

    let n = CanonicalRegge3jm::len(tj_max);
    let mut vec = vec![Default::default(); n];
    for &w3jm in &table {
        let value = w3jm.value();
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        vec[regge.index()] = f64::from(phase * value);
    }

    let range = rand::distributions::Range::new(0, table.len());
    let mut rng = rand::weak_rng();
    b.iter(|| {
        let w3jm = table[range.ind_sample(&mut rng)];
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        test::black_box(phase as f64 * vec[regge.index()]);
    });
}
