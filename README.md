# cute-nucleotides
Cute tricks for efficiently encoding and decoding nucleotides.

## Running
To run the tests, use cargo with a special flag telling it to target your CPU, for maximum efficiency.
```
RUSTFLAGS="-C target-cpu=native" cargo test
```
You can also run the benchmarks:
```
RUSTFLAGS="-C target-cpu=native" cargo bench
```

These should all run on x86 CPUs that support AVX2 and BMI2 instructions (so modern Intel and AMD CPUs).
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
encoding and decoding be done if we squeeze every last cycle out of my CPU? Spoiler: we can go faster
than `memcpy` (`std::ptr::copy_nonoverlapping` in Rust).

In general, we can avoid costly branches, reduce data dependencies, use undiscernible magic numbers,
and exploit special vectorized instructions for maximum speed.

## Nucleotides to bits conversion
The goal here is fast, case-insensitive conversion of a string of nucleotides to pairs of bits:
`{A, T/U, C, G} -> {00, 10, 01, 11}`. Here, I give a brief overview of every method that I implemented:

* **n_to_bits_lut**. This is just a straightforward loop through each input nucleotide
and using a lookup table to map each byte to two bits. Some bit manipulation magic is used to pack them
into 64-bit integers. Nothing special here. After this, we are stepping into the domain of vectorized,
platform-specific instructions. Most notably, the AVX instruction set allows us to manipulate vector registers
of 32 bytes (256 bits) efficiently with SIMD (Single Instruction Multiple Data) operations.

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
Crazy, right? Here's an example, where we want to create the binary number `dcba00` from `dc00ba` (`a`, `b`, `c`,
and `b` are any bits) by simply multiplying by a magic binary number and masking out extraneous bits (in practice,
truncation and shuffles are used instead of a mask):

    ```
    dcba00 = dc00ba * 000101

            dc00ba
    *       000101
    --------------
            dc00ba
    +     dc00ba    shift and add
    --------------
          dcdcbaba
    &     00111100
    --------------
            dcba00
    ```

    Multiplication by a constant is just shifts to where the bits are turned on in the constant and vertically adding.
    If there aren't any potential overlaps between one bits after shifting, then there won't be any carries to mess up
    the result. Therefore, multiplication by a specific constant number with selected bits turned on is a quick bit-level
    shuffle, and it can be used to pack four pairs of bits into a single byte with just a single 32-bit multiply. Additionally,
    AVX2 provides a way for eight multiplications to happen simultaneously! Compared to the other methods, this method can be
    easily implemented without vector instructions, on just 32-bit words, and it does not require BMI2 support.

## Bits to nucleotides conversion
When converting bits to nucleotides, the goal is to unpack each byte with four nucleotides of two bits each into 32 bits.
We need an extra length parameter, since we do not know when the sequence of bits terminates. Here is an overview of
the implemented methods:

* **bits_to_n_lut**. The most basic implementation. Two bits in, one byte character out. A small lookup table is used to
convert bits to characters.

Now, let's look at vectorized methods with scary magic numbers.

* **bits_to_n_shuffle (AVX2)**. The main goal of each bits-to-nucleotide conversion algorithm is to distribute the two bits
of each nucleotide to the low half of each byte. Then, a `_mm256_shuffle_epi8` instruction can use those two bits to lookup
the full byte character. The simplest way to distribute the two bits to each byte is to first do a byte-granularity shuffle
to duplicate each byte of four nucleotides four times, and then use shifts to move pairs of bits to the start of each byte.
Naively, this requires four shifts, because there is no single AVX2 instruction for separately shifting the bits within
every byte by different amounts. However, we can reduce the number of shifts by noticing that we do not need to get the bits
to the starting two bit positions of each byte. We just need to get them to be within the first four bits, since the
`_mm256_shuffle_epi8` lookup table can use the first four bits to index into a 16-byte lookup table. We only need to shift
every other 16-bit chunk to ensure that all four groups of nucleotide bits land in the first 4 bits of every byte. For any
two bits `a` and `b`, the lookup table will map `00ba` and `ba00` to the same nucleotide byte character.

* **bits_to_n_pdep (AVX2 and BMI2)**. The `_pdep_u64` instruction, which distributes bits to specific locations based on a
mask is the perfect instruction for the task of depositing two bits that represent one nucleotide to the start of each byte.
After depositing the bits in chunks of 64-bits, the last step is to do a shuffle lookup to convert those bits to byte
characters.

* **bits_to_n_clmul (SSSE3 and PCLMULQDQ)**. It is also possible to use multiplication to distribute bits to the start of
each byte. However, in this case, it is difficult to distribute each byte (that contains four nucleotides) to its correct
locations without portions of each byte overlapping. When these portions overlap, during the addition step of multiplication,
there will be carries that cause more significant bits (bits to the left) to be wrong. The solution to this is carry-less
multiplication (with a throughput of one per cycle on recent CPUs), but it only works for multiplying 64-bit integers from
128-bit vectors. This means that it is slower than the other methods, and we need an extra step to combine 64-bit products
using the `_mm_movelh_ps` instruction. Regardless, this method is pretty cool. Here is an example of carry-less multiplication
in action to shuffle `0000dcba` into `00dc00ba` (`a`, `b`, `c`, and `d` are any bits):

    ```
    00dc00ba = clmul(0000dcba, 00000101)

              0000dcba
    *         00000101
    ------------------
              0000dcba
    ^       0000dcba    XOR instead of addition!
    ------------------
          000000dc??ba  we don't care about '?' bits
    &     000000110011  they won't affect other bits
    ------------------
              00dc00ba
    ```

## Cuter algorithms for converting undetermined nucleotides to bits
### More motivation


### Nucleotides to bits conversion


### Bits to nucleotides conversion


## Cute speedups in cute micro-benchmarks
All benchmarks were ran on a Intel Core i9-9880H (Coffee Lake-H) CPU with a clock rate of 2.3 Ghz. Run times were measured
on byte strings with 40,000 nucleotides. For shorter strings, the speedup of vectorized methods should be less noticeable
due to overhead. Recall that the lookup table (`lut`) methods are the naive, scalar algorithms.

First, we compare the throughputs of the algorithms for converting nucleotides to binary strings:

| Nucleotides to bits    | Throughput       |
|------------------------|------------------|
| n_to_bits_lut          | 1.9777 GiB/s     |
| n_to_bits_pext         | 19.869 GiB/s     |
| n_to_bits_shift        | 23.052 GiB/s     |
| **n_to_bits_movemask** | **28.962 GiB/s** |
| n_to_bits_mul          | 25.286 GiB/s     |
| memcpy                 | 23.599 GiB/s     |

The winner is `n_to_bits_movemask`, beating even `memcpy` (`std::ptr::copy_nonoverlapping`). Note that all methods
(including `memcpy`) allocate memory for the output first, which actually takes a nontrivial amount of time. In practice,
`n_to_bits_mul` is probably the most flexible method for different platforms that is fast enough.

The decoding operation of converting bits to nucleotides shows a more dramatic difference between the winner and the other
methods:

| Bits to nucleotides    | Throughput       |
|------------------------|------------------|
| bits_to_n_lut          | 1.8257 GiB/s     |
| **bits_to_n_shuffle**  | **30.224 GiB/s** |
| bits_to_n_pdep         | 14.216 GiB/s     |
| bits_to_n_clmul        | 11.819 GiB/s     |

Giving these results, it seems like `bits_to_n_clmul`, which uses SSE and PCLMULQDQ instructions to handle 128-bit vectors,
won't be able to beat `bits_to_n_shuffle` even if `bits_to_n_shuffle` is translated to use 128-bit vectors instead of 256-bit
vectors.

| Undetermined nucleotides | Throughput       |
|--------------------------|------------------|
| n_to_bits2_lut           | 1.7887 GiB/s     |
| **n_to_bits2_pext**      | **11.787 GiB/s** |
| bits_to_n2_lut           | 1.3751 GiB/s     |
| **bits_to_n2_pdep**      | 10.175 GiB/s     |

For encoding and decoding undetermined nucleotides, the vectorized algorithms provide clear speedups over the lookup-table-based
naive methods.

This was an insightful experiment that unvealed many bit twiddling tricks. If these tricks were implemented into software tools,
then encoding/decoding DNA/RNA sequences should not be a bottleneck. In fact, it should lead to nontrivial memory savings, if
encoded binary strings are converted to nucleotides on the fly when they are needed. Unfortunately, its probably not possible to
reach the lofty 30 GiB/s throughput in practice, due to memory access bottlenecks.
