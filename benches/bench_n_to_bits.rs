use std::ptr;
use criterion::*;
use rust_fast::n_to_bits::*;

// Note: memory allocation takes a nontrivial amount of time!
// For fair comparison, all functions must allocate memory for its output data.

fn bench_n_to_bits(c: &mut Criterion) {
    let n = black_box(get_nucleotides(1000));

    let mut group = c.benchmark_group("n_to_bits");

    group.bench_function("n_to_bits_lut", |b| b.iter(|| n_to_bits_lut(&n)));
    group.bench_function("n_to_bits_pext", |b| b.iter(|| n_to_bits_pext(&n)));
    group.bench_function("n_to_bits_mul", |b| b.iter(|| n_to_bits_mul(&n)));
    group.bench_function("memcpy", |b| b.iter(|| unsafe {let mut dest = vec![0u8; n.len()]; ptr::copy_nonoverlapping(n.as_ptr(), dest.as_mut_ptr(), n.len()); dest}));

    group.finish();
}

fn bench_bits_to_n(c: &mut Criterion) {
    let bits = black_box(get_bits(1000));
    let len = black_box(4 * 1000);

    let mut group = c.benchmark_group("bits_to_n");

    group.bench_function("bits_to_n_lut", |b| b.iter(|| bits_to_n_lut(&bits, len)));
    group.bench_function("bits_to_n_shuffle", |b| b.iter(|| bits_to_n_shuffle(&bits, len)));
    group.bench_function("bits_to_n_pdep", |b| b.iter(|| bits_to_n_pdep(&bits, len)));
    group.bench_function("memcpy", |b| b.iter(|| unsafe {let mut dest = vec![0u64; bits.len()]; ptr::copy_nonoverlapping(bits.as_ptr(), dest.as_mut_ptr(), bits.len()); dest}));

    group.finish();
}

criterion_group!(benches, bench_n_to_bits, bench_bits_to_n);
criterion_main!(benches);

fn get_nucleotides(repeat: usize) -> Vec<u8> {
    b"ATCG".repeat(repeat)
}

fn get_bits(repeat: usize) -> Vec<u64> {
    let bits = 0b11011000u64;

    let mut curr = 0u64;

    for i in 0..8 {
        curr |= bits << (i << 3);
    }

    vec![curr].repeat((repeat >> 3) + if repeat & 7 == 0 {0} else {1})
}
