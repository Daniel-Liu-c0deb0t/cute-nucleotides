# rust-fast
A collection of cute tests and benchmarks in Rust for various tasks.

These programs should be run on x86 CPUs that supports AVX2 and BMI2 instructions.

**Warning: there is a lot of unsafe code! Your eyes may make you think that the code is written
in C. No, it is (disappointingly) 100% organic Rust. Read it at your own risk.**

## Cute speedups for converting nucleotides to bits
### Motivation


### Nucleotides to bits conversion algorithms
* `pext` algorithm (requires AVX2 and BMI2)
* magic multiplication algorithm (requires AVX2)

### Bits to nucleotides conversion algorithms


### Cuter speedups for converting undetermined nucleotides to bits
#### More motivation


#### Nucleotides to bits conversion algorithms


#### Bits to nucleotides conversion algorithms


### Cute speedups in cute micro-benchmarks

| Algorithm                        | Throughput   |
|----------------------------------|--------------|
| n_to_bits/n_to_bits_lut          | 1.9777 GiB/s |
| n_to_bits/n_to_bits_pext         | 19.869 GiB/s |
| n_to_bits/n_to_bits_shift        | 23.052 GiB/s |
| n_to_bits/n_to_bits_movemask     | 28.962 GiB/s |
| n_to_bits/n_to_bits_mul          | 25.286 GiB/s |
| n_to_bits/memcpy                 | 23.599 GiB/s |
| bits_to_n/bits_to_n_lut          | 1.8257 GiB/s |
| bits_to_n/bits_to_n_shuffle      | 30.224 GiB/s |
| bits_to_n/bits_to_n_pdep         | 14.216 GiB/s |
| bits_to_n/bits_to_n_clmul        | 11.819 GiB/s |
| n_to_bits2/n_to_bits2_lut        | 1.7887 GiB/s |
| n_to_bits2/n_to_bits2_pext       | 11.787 GiB/s |
| bits_to_n2/bits_to_n2_lut        | 1.3751 GiB/s |
| bits_to_n2/bits_to_n2_pdep       | 10.175 GiB/s |
