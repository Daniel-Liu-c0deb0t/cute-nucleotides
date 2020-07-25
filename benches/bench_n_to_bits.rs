use std::ptr;
use criterion::*;
use rust_fast::n_to_bits::*;
use rust_fast::n_to_bits2::*;

// Note: memory allocation takes a nontrivial amount of time!
// For fair comparison, all functions must allocate memory for its output data.

fn bench_n_to_bits(c: &mut Criterion) {
    let n = black_box(get_nucleotides(10000));

    let mut group = c.benchmark_group("n_to_bits");

    group.bench_function("n_to_bits_lut", |b| b.iter(|| n_to_bits_lut(&n)));
    group.bench_function("n_to_bits_pext", |b| b.iter(|| n_to_bits_pext(&n)));
    group.bench_function("n_to_bits_shift", |b| b.iter(|| n_to_bits_shift(&n)));
    group.bench_function("n_to_bits_movemask", |b| b.iter(|| n_to_bits_movemask(&n)));
    group.bench_function("n_to_bits_mul", |b| b.iter(|| n_to_bits_mul(&n)));
    group.bench_function("memcpy", |b| b.iter(|| unsafe {let mut dest = vec![0u8; n.len()]; ptr::copy_nonoverlapping(n.as_ptr(), dest.as_mut_ptr(), n.len()); dest}));

    group.finish();
}

fn bench_n_to_bits2(c: &mut Criterion) {
    let n = black_box(get_nucleotides_undetermined(8000));

    let mut group = c.benchmark_group("n_to_bits2");

    group.bench_function("n_to_bits2_lut", |b| b.iter(|| n_to_bits2_lut(&n)));
    group.bench_function("n_to_bits2_pext", |b| b.iter(|| n_to_bits2_pext(&n)));

    group.finish();
}

fn bench_bits_to_n(c: &mut Criterion) {
    let bits = black_box(get_bits(10000));
    let len = black_box(4 * 10000);

    let mut group = c.benchmark_group("bits_to_n");

    group.bench_function("bits_to_n_lut", |b| b.iter(|| bits_to_n_lut(&bits, len)));
    group.bench_function("bits_to_n_shuffle", |b| b.iter(|| bits_to_n_shuffle(&bits, len)));
    group.bench_function("bits_to_n_pdep", |b| b.iter(|| bits_to_n_pdep(&bits, len)));
    group.bench_function("bits_to_n_clmul", |b| b.iter(|| bits_to_n_clmul(&bits, len)));

    group.finish();
}

fn bench_bits_to_n2(c: &mut Criterion) {
    let bits = black_box(get_bits_undetermined(8000));
    let len = black_box(5 * 8000);

    let mut group = c.benchmark_group("bits_to_n2");

    group.bench_function("bits_to_n2_lut", |b| b.iter(|| bits_to_n2_lut(&bits, len)));
    group.bench_function("bits_to_n2_pdep", |b| b.iter(|| bits_to_n2_pdep(&bits, len)));

    group.finish();
}

criterion_group!(benches, bench_n_to_bits, bench_bits_to_n, bench_n_to_bits2, bench_bits_to_n2);
criterion_main!(benches);

fn get_nucleotides(repeat: usize) -> Vec<u8> {
    b"ATCG".repeat(repeat)
}

fn get_nucleotides_undetermined(repeat: usize) -> Vec<u8> {
    b"ATCGN".repeat(repeat)
}

fn get_bits(repeat: usize) -> Vec<u64> {
    n_to_bits_lut(&(b"ATCG".repeat(repeat)))
}

fn get_bits_undetermined(repeat: usize) -> Vec<u64> {
    n_to_bits2_lut(&(b"ATCGN".repeat(repeat)))
}
