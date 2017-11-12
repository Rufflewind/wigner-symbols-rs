#![feature(test)]

extern crate rand;
extern crate test;
extern crate wigner_symbols;

use rand::distributions::IndependentSample;
use test::Bencher;
use wigner_symbols::{internal, clebsch_gordan};

fn clebsch_gordan_bench(b: &mut Bencher, tj_min: i32, tj_max: i32) {
    let mut args = Vec::default();
    internal::get_3tjms(tj_max, &mut |tj1, tm1, tj2, tm2, tj3, tm3| {
        if internal::sort3(tj1, tj2, tj3).2 >= tj_min {
            args.push((tj1, tm1, tj2, tm2, tj3, tm3));
        }
    });
    let range = rand::distributions::Range::new(0, args.len());
    let mut rng = rand::thread_rng();
    b.iter(|| {
        let (tj1, tm1, tj2, tm2, tj3, tm3) = args[range.ind_sample(&mut rng)];
        test::black_box(clebsch_gordan(tj1, tm1, tj2, tm2, tj3, -tm3));
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
