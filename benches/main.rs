#![feature(test)]

extern crate fnv;
extern crate rand;
extern crate test;
extern crate wigner_symbols;

use fnv::FnvHashMap;
use rand::Rng;
use test::Bencher;
use wigner_symbols::{internal, ClebschGordan, Wigner3jm, Wigner6j, Wigner9j};
use wigner_symbols::regge::{CanonicalRegge3jm, CanonicalRegge6j, Regge3jm};

fn wigner_3jm_arg_table(tj_min: i32, tj_max: i32) -> Vec<Wigner3jm> {
    let mut table = Vec::default();
    internal::get_3tjms(tj_max, &mut |w3jm| {
        if [
            w3jm.tj1, w3jm.tj2, w3jm.tj3,
        ].iter().max().unwrap() >= &tj_min {
            table.push(w3jm);
        }
    });
    table
}

fn wigner_6j_arg_table(tj_min: i32, tj_max: i32) -> Vec<Wigner6j> {
    let mut table = Vec::default();
    internal::get_6tjs(tj_max, &mut |w6j| {
        if [
            w6j.tj1, w6j.tj2, w6j.tj3,
            w6j.tj4, w6j.tj5, w6j.tj6,
        ].iter().max().unwrap() >= &tj_min {
            table.push(w6j);
        }
    });
    table
}

fn wigner_9j_arg_table(tj_min: i32, tj_max: i32) -> Vec<Wigner9j> {
    let mut table = Vec::default();
    internal::get_9tjs(tj_max, &mut |w9j| {
        if [
            w9j.tj1, w9j.tj2, w9j.tj3,
            w9j.tj4, w9j.tj5, w9j.tj6,
            w9j.tj7, w9j.tj8, w9j.tj9,
        ].iter().max().unwrap() >= &tj_min {
            table.push(w9j);
        }
    });
    table
}

fn clebsch_gordan_bench(b: &mut Bencher, tj_min: i32, tj_max: i32) {
    let mut table = wigner_3jm_arg_table(tj_min, tj_max);
    rand::thread_rng().shuffle(&mut table);
    let mut i = 0;
    b.iter(|| {
        test::black_box(ClebschGordan::from(table[i % table.len()]).value());
        i += 1;
    });
}

fn wigner_6j_bench(b: &mut Bencher, tj_min: i32, tj_max: i32) {
    let mut table = wigner_6j_arg_table(tj_min, tj_max);
    rand::thread_rng().shuffle(&mut table);
    let mut i = 0;
    b.iter(|| {
        test::black_box(table[i % table.len()].value());
        i += 1;
    });
}

// don't go all the way to 15 if you don't want to lose all your RAM
fn wigner_9j_bench(b: &mut Bencher, tj_min: i32, tj_max: i32) {
    let mut table = wigner_9j_arg_table(tj_min, tj_max);
    rand::thread_rng().shuffle(&mut table);
    let mut i = 0;
    b.iter(|| {
        test::black_box(table[i % table.len()].value());
        i += 1;
    });
}

#[bench]
fn bench_clebsch_gordan_00_04(b: &mut Bencher) {
    clebsch_gordan_bench(b, 0, 4);
}

#[bench]
fn bench_clebsch_gordan_05_09(b: &mut Bencher) {
    clebsch_gordan_bench(b, 5, 9);
}

#[bench]
fn bench_clebsch_gordan_10_14(b: &mut Bencher) {
    clebsch_gordan_bench(b, 10, 14);
}

#[bench]
fn bench_clebsch_gordan_15_19(b: &mut Bencher) {
    clebsch_gordan_bench(b, 15, 19);
}

#[bench]
fn bench_clebsch_gordan_20_24(b: &mut Bencher) {
    clebsch_gordan_bench(b, 20, 24);
}

#[bench]
fn bench_wigner_6j_00_04(b: &mut Bencher) {
    wigner_6j_bench(b, 0, 4);
}

#[bench]
fn bench_wigner_6j_05_09(b: &mut Bencher) {
    wigner_6j_bench(b, 5, 9);
}

#[bench]
fn bench_wigner_6j_10_14(b: &mut Bencher) {
    wigner_6j_bench(b, 10, 14);
}

#[bench]
fn bench_wigner_6j_15_19(b: &mut Bencher) {
    wigner_6j_bench(b, 15, 19);
}

#[bench]
fn bench_wigner_6j_20_24(b: &mut Bencher) {
    wigner_6j_bench(b, 20, 24);
}

#[bench]
fn bench_wigner_9j_00_02(b: &mut Bencher) {
    wigner_9j_bench(b, 0, 2);
}

#[bench]
fn bench_wigner_9j_03_05(b: &mut Bencher) {
    wigner_9j_bench(b, 3, 5);
}

#[bench]
fn bench_wigner_9j_06_08(b: &mut Bencher) {
    wigner_9j_bench(b, 6, 8);
}

// this one just uses a plain FnvHashMap, which requires storing both keys and
// values
#[bench]
fn bench_regge3jm_map_25(b: &mut Bencher) {
    let tj_max = 25;
    let mut table = wigner_3jm_arg_table(0, tj_max);
    rand::thread_rng().shuffle(&mut table);

    let mut map = FnvHashMap::default();
    for &w3jm in &table {
        let value = w3jm.value();
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        map.entry(regge).or_insert_with(|| f64::from(phase * value));
    }

    let mut i = 0;
    b.iter(|| {
        let w3jm = table[i % table.len()];
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        test::black_box(f64::from(phase) * map[&regge]);
        i += 1;
    });
}

// this one just uses the ordering scheme of Rasch and Yu, which only requires
// storing values but has a density of about 50%
#[bench]
fn bench_regge3jm_vec_25(b: &mut Bencher) {
    let tj_max = 25;
    let mut table = wigner_3jm_arg_table(0, tj_max);
    rand::thread_rng().shuffle(&mut table);

    let n = CanonicalRegge3jm::len(tj_max);
    let mut vec = vec![Default::default(); n];
    for &w3jm in &table {
        let value = w3jm.value();
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        vec[regge.index()] = f64::from(phase * value);
    }

    let mut i = 0;
    b.iter(|| {
        let w3jm = table[i % table.len()];
        let (regge, phase) = Regge3jm::from(w3jm).canonicalize();
        test::black_box(f64::from(phase) * vec[regge.index()]);
        i += 1;
    });
}

#[bench]
fn bench_regge6j_map_25(b: &mut Bencher) {
    let tj_max = 25;
    let mut table = wigner_6j_arg_table(0, tj_max);
    rand::thread_rng().shuffle(&mut table);

    let mut map = FnvHashMap::default();
    for &w6j in &table {
        let value = w6j.value();
        let regge = CanonicalRegge6j::from(w6j);
        map.entry(regge).or_insert_with(|| f64::from(value));
    }

    let mut i = 0;
    b.iter(|| {
        let w6j = table[i % table.len()];
        let regge = CanonicalRegge6j::from(w6j);
        test::black_box(map[&regge]);
        i += 1;
    });
}

#[bench]
fn bench_regge6j_vec_25(b: &mut Bencher) {
    let tj_max = 25;
    let mut table = wigner_6j_arg_table(0, tj_max);
    rand::thread_rng().shuffle(&mut table);

    let n = CanonicalRegge6j::len(tj_max);
    let mut vec = vec![Default::default(); n];
    for &w6j in &table {
        let value = w6j.value();
        let regge = CanonicalRegge6j::from(w6j);
        vec[regge.index()] = f64::from(value);
    }

    let mut i = 0;
    b.iter(|| {
        let w6j = table[i % table.len()];
        let regge = CanonicalRegge6j::from(w6j);
        test::black_box(vec[regge.index()]);
        i += 1;
    });
}
