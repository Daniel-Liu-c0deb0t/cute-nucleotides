# cute-nucleotides
Cute tricks for efficiently encoding and decoding nucleotides.

This should be run on x86 CPUs that support AVX2 and BMI2 instructions (so modern Intel and AMD CPUs).
Many functions are not written in a cross-platform way, since I'm too lazy to do so for mere test programs.

**Warning: there is a lot of unsafe code! Your eyes may trick you into thinking that the code is written
in C. No, it is (unfortunately) 100% organic Rust. Read it at your own risk.**

## Motivation
It turns out that when dealing with millions and billions of sequenced DNA/RNA data, a large amount
of memory space is necessary (gasp!). Therefore, it is important to compress the DNA nucleotides
`{A, T, C, G}` or RNA nucleotides `{A, U, C, G}` to save precious space. There are general-purpose
compression techniques like Gzip, but they are not suitable for specific cases, like fast, in-memory
compression/decompression, or when comparisons/calculations are done over the compressed genomic data.

The alternative is a simple compression scheme translating a string of nucleotides (`{A, T/U, C, G}`)
to a string of bits (`{00, 10, 01, 11}`). Usual text formats take up one byte per nucleotide, while this
only takes two bits, for a four times reduction in space. At the character level, this is obviously optimal,
and on random strings of nucleotide characters, it is difficult to do much better. To save space, these bits
are packed into 64-bit integer words.

We get the following benefits from this encoding technique:
* Less space compared to using byte strings.
* Works well with short fragments of strings, and does not require decoding large blocks like in most
compression schemes. Random access is ok.
* Likely fast enough to not be a bottleneck in larger software tools.
* Many operations (like Hamming distance) can be done directly on the bit strings without decoding.

These benefits make this encoding technique very commonly used. The natural question is: how fast can
encoding and decoding be done if we squeeze every last cycle out of a CPU? Spoiler: we can go faster
than `memcpy` (`std::ptr::copy_nonoverlapping` in Rust).

## Nucleotides to bits conversion
The goal here is fast, case-insensitive conversion of a string of nucleotides to pairs of bits:
`{A, T/U, C, G} -> {00, 10, 01, 11}`. We do not have to validate the data for whether other characters are
present, and we can avoid all branches in the hot loops.

Here, I give a brief overview of every method that I implemented:

* **n_to_bits_lut**. This is just a straightforward loop through each input nucleotide
and using a lookup table to map each byte to two bits. Some bit manipulation magic is used to pack them
into 64-bit integers. Nothing special here. After this, we are stepping into the domain of vectorized,
platform-specific instructions. Most notably, the AVX instruction set allows us to manipulate vector registers
of 32 bytes (256 bits) efficiently.

* **n_to_bits_pext (AVX2 and BMI2)**. It turns out that we don't need a lookup table. In fact, if you stare
at the ASCII table for a bit, you will realize that the second and third bits of each nucleotide character
is exactly the two bits that we want. Additionally, the second and third bits for the `U` and `T` characters
are the same! Now, if only there is a way to extract those bits and pack them together! Well, there is a pretty
cool BMI2 instruction called `_pext_u64`, which extracts bits that are selected based on a mask, and packs
them into a new 64-bit integer. Exactly what we need. Four `_pext_u64` instructions are used to build up one
64-bit integer. On modern CPUs, this instruction is pretty fast, with a throughput of one per cycle.

* **n_to_bits_shift (AVX2)**. If we use the good ol' vector `_mm256_and_si256` instruction to only keep the
second and third bits in each byte character, then it becomes pretty obvious that a bunch of vector shifts and
`_mm256_or_si256` instructions can be used to pack four nucleotides (2 bits per nucleotide) into a single byte
in parallel. Doing so cleverly, by first combining pairs of 8-bit chunks, and then pairs of 16-bit chunks,
results in a fast vectorized algorithm without the `_pext_u64` instruction. Notice that the difficult part of
the problem of converting nucleotides to binary is actually packing those pairs of bits together.

* **n_to_bits_movemask (AVX2 and BMI2)**. There aren't many bit-level shuffle/scatter/gather instructions in
the AVX2 repertoire, so the `_mm256_movemask_epi8` instruction is one of the few. This instruction is very
powerful: it extracts the most significant bit out of each byte in a vector, and packs them into a 32-bit
integer. We can do this twice, along with shifting the bits, to extract two bits out of every single byte.
The only problem is that we have to interleave the bits from the two 32-bit integers we get. This
interleaving can be done using the BMI2 `_pdep_u64` instruction, which is used to deposit bits at every
other position in the final 64-bit integer. Both `_mm256_movemask_epi8` and `_pdep_u64` have throughputs of
one per cycle on modern CPUs, so they are quite fast.

* **n_to_bits_mul (AVX2)**. Surprisingly, multiplication can be used as a fast bit-level shuffle operation.
Crazy, right? Here's an example, where we want to create the binary number `ab00ab` from just `0000ab` (`a` and
`b` are any two bits) by simply multiplying the magic binary number `010001`:
```
ab00ab = 0000ab * 010001

        0000ab
*       010001
--------------
        0000ab
+   0000ab
-------------
        ab00ab
```
Multiplication by a constant is just shifts to where the bits are turned on in the constant and vertically adding.
If there aren't any potential overlaps between one bits after shifting, then there won't be any carries to mess up
the result. Therefore, multiplication by a specific constant number with selected bits turned on is a quick bit-level
shuffle, and it can be used to pack four pairs of bits into a single byte with just a single 32-bit multiply. AVX2
provides a way for eight multiplications to happen simultaneously! Compared to the other methods, this method can be
easily implemented without vector instructions, on just 32-bit words, and it does not require BMI2 support.

## Bits to nucleotides conversion


## Cuter algorithms for converting undetermined nucleotides to bits
### More motivation


### Nucleotides to bits conversion


### Bits to nucleotides conversion


## Cute speedups in cute micro-benchmarks
All benchmarks were ran on a Intel Core i9-9880H (Coffee Lake-H) CPU with a clock rate of 2.3 Ghz. Run times were measured on byte strings with
40,000 nucleotides. For shorter strings, the speedup of vectorized methods should be less noticeable due to overhead.
Recall that the lookup table (`lut`) methods are the naive, scalar algorithms.

First, we compare the throughputs of the algorithms for converting nucleotides to binary strings:

| Nucleotides to bits    | Throughput   |
|------------------------|--------------|
| n_to_bits_lut          | 1.9777 GiB/s |
| n_to_bits_pext         | 19.869 GiB/s |
| n_to_bits_shift        | 23.052 GiB/s |
| **n_to_bits_movemask** | 28.962 GiB/s |
| n_to_bits_mul          | 25.286 GiB/s |
| memcpy                 | 23.599 GiB/s |

The winner is `n_to_bits_movemask`, beating even `memcpy` (`std::ptr::copy_nonoverlapping`). Note that all methods
(including `memcpy`) allocate memory for the output first, which actually takes a nontrivial amount of time.

The decoding operation of converting bits to nucleotides shows a more dramatic difference between the winner and the other
methods:

| Bits to nucleotides    | Throughput   |
|------------------------|--------------|
| bits_to_n_lut          | 1.8257 GiB/s |
| **bits_to_n_shuffle**  | 30.224 GiB/s |
| bits_to_n_pdep         | 14.216 GiB/s |
| bits_to_n_clmul        | 11.819 GiB/s |

Giving these results, it seems like `bits_to_n_clmul`, which uses SSE and PCLMULQDQ instructions to handle 128-bit vectors,
won't be able to beat `bits_to_n_shuffle` even if `bits_to_n_shuffle` translated to handle 128-bit vectors instead of 256-bit
vectors.

| Undetermined nucleotides | Throughput   |
|--------------------------|--------------|
| n_to_bits2_lut           | 1.7887 GiB/s |
| **n_to_bits2_pext**      | 11.787 GiB/s |
| bits_to_n2_lut           | 1.3751 GiB/s |
| **bits_to_n2_pdep**      | 10.175 GiB/s |

For encoding and decoding undetermined nucleotides, the vectorized algorithms provide clear speedups over the lookup-table-based
naive methods.
