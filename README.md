# cute nucleotides 🧬 💻
Cute tricks for vectorized binary encoding and decoding of nucleotides.

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

These benefits make this encoding technique very commonly used in software tools. The natural question is: how fast
can encoding and decoding be done if we squeeze every last cycle out of my CPU? Spoiler: we can go faster
than `memcpy` (`std::ptr::copy_nonoverlapping` in Rust). In general, we avoid costly branches, reduce data dependencies,
use indiscernible magic numbers, and exploit special vectorized instructions for maximum efficiency. In the experiments, I compare
many different techniques for encoding and decoding.

In addition to an empirical investigation on the best methods for encoding/decoding DNA, I also hope that this could
be a good resource for some bit-twiddling micro-optimization tricks. Helpful resources that inspired me
are Daniel Lemire's [blog](https://lemire.me/blog/), Geoff Langdale's [blog](https://branchfree.org/), literally all
of Peter Cordes's [StackOverflow](https://stackoverflow.com/users/224132/peter-cordes) answers, the
[Bit Twiddling Hacks](https://graphics.stanford.edu/~seander/bithacks.html) page by Sean Eron Anderson, and the *very*
detailed instruction tables and optimization guide on Agner Fog's [website](https://www.agner.org/optimize/).

This experiment was inspired by Pierre Marijon's Optimization Friday [work](https://github.com/natir/fbio) and our discussion
on Twitter.

## Encoding 🧬 ➡️ 🔟
The goal here is fast, case-insensitive conversion of a string of nucleotides to pairs of bits:
`{A, T/U, C, G} -> {00, 10, 01, 11}`. For example:
```
Nucleotides =         A        T        G
ASCII       =  01000001 01010100 01000111
Encoding    =        00       10       11
            =                      001011
```
Here, I give a brief overview of every method that I implemented:

* **n_to_bits_lut**. This is just a straightforward loop through each input nucleotide
and using a lookup table to map each byte to two bits. Some bit manipulation magic is used to pack them
into 64-bit integers. Nothing special here. After this, we are stepping into the domain of vectorized,
platform-specific instructions. Most notably, the AVX instruction set allows us to manipulate vector registers
of 32 bytes (256 bits) efficiently with SIMD (Single Instruction Multiple Data) operations.

* **n_to_bits_pext (AVX2 and BMI2)**. It turns out that we don't need a lookup table. In fact, if you stare
at the ASCII table for a bit, you will realize that the second and third bits of each nucleotide character
is exactly the two bits that we want:

    ```
        ASCII      mask       bits
    A = 01000001 & 00000110 = 00000 00 0
    T = 01010100 & 00000110 = 00000 10 0
    U = 01010101 & 00000110 = 00000 10 0
    C = 01000011 & 00000110 = 00000 01 0
    G = 01000111 & 00000110 = 00000 11 0
    ```

    Additionally, the second and third bits for the `U` and `T` characters are the same! Now, if only there is
    a way to extract those bits and pack them together! Well, there is a pretty cool BMI2 instruction called
    `_pext_u64`, which extracts bits that are selected based on a mask, and packs them into a new 64-bit integer.
    Exactly what we need. Four `_pext_u64` instructions are used to build up one 64-bit integer. On modern CPUs,
    this instruction is pretty fast, with a throughput of one per cycle.

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
interleaving can be done using the BMI2 `_pdep_u64` instruction (described in a blog
[post](https://lemire.me/blog/2018/01/08/how-fast-can-you-bit-interleave-32-bit-integers/) by Daniel Lemire), which
is used to deposit bits at every other position in the final 64-bit integer. Both `_mm256_movemask_epi8` and
`_pdep_u64` have throughputs of one per cycle on modern CPUs, so they are quite fast. If the `_pdep_u64` instruction
is not available, and a lookup-table-based bit interleave algorithm is used, then this method would probably be slower.

* **n_to_bits_mul (AVX2)**. Surprisingly, multiplication can be used as a fast bit-level shuffle operation.
Crazy, right? Here's an example, where we want to create the binary number `dcba0000` from `00dc00ba` (`a`, `b`, `c`,
and `b` are any bits) by simply multiplying by a magic binary number and masking out extraneous bits (in practice,
truncation and shuffles are used instead of a mask):

    ```
    dcba0000 = 00dc00ba * 00010100

           00dc00ba
    *      00010100
    ---------------
         00dc00ba
    +  00dc00ba      multiply = shift and add
    ---------------
       00dcdcbaba
    &  000011110000  mask out extraneous bits
    ---------------
           dcba0000
    ```

    Multiplication by a constant is just shifts to where the bits are turned on in the constant and vertically adding.
    If there aren't any potential overlaps between one bits after shifting, then there won't be any carries to mess up
    the result. Therefore, multiplication by a specific constant number with selected bits turned on is a quick bit-level
    shuffle, and it can be used to pack four pairs of bits into a single byte with just a single 32-bit multiply. On modern CPUs,
    multiplications should be very cheap. Additionally, AVX2 provides a way for eight multiplications to happen simultaneously!
    Compared to the other methods, this method can be easily implemented without vector instructions, on just 32-bit words,
    and it does not require BMI2 support.

## Decoding 🔟 ➡️ 🧬
When converting bits to nucleotides, the goal is to unpack each byte with four nucleotides of two bits each into 32 bits.
We need an extra length parameter, since we do not know when the sequence of bits terminates. Example:
```
Encoding    =                      001011
            =        00       10       11
ASCII       =  01000001 01010100 01000111
Nucleotides =         A        T        G
```
The decoding methods use similar techniques to the encoding methods. Here is an overview of the implemented methods:

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
characters. This should be pretty fast, but this involves an extra step to convert four 64-bit integers into a 256-bit vector.

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
          000000dc??ba  we don't care about '?' bits, as
    &     000000110011  they won't affect other bits
    ------------------
              00dc00ba
    ```

    Carry-less multiplication was brought to my attention by Geoff Langdale's [post](https://branchfree.org/2019/03/06/code-fragment-finding-quote-pairs-with-carry-less-multiply-pclmulqdq/)
    on using it for fast parallel prefix XOR, and I have been itching to apply carry-less multiplication to *something*.
    I have to admit, multiplication in general is one of my favorite bit-twiddling algorithm. Normal multiplication can
    also be used for fast prefix sum.

## Cuter algorithms for converting undetermined nucleotides to bits
### More motivation
In some cases, it is important to keep track of an extra "undetermined nucleotide" represented with `N`, in addition to the
typical `{A, T, C, G}` nucleotides. Usually, this arises from ambiguities in sequencing, and the exact nucleotide cannot be
determined. In some sequence alignment applications, it may be beneficial to keep track of a gap character represented with
`-`. Regardless, it is important to keep track of an extra character in addition to the four common nucleotides. This encoding
will give much of the same benefits as the two-bit encoding, since it is pretty simple.

In the implementations, we will consider the problem of encoding and decoding `{A, T/U, C, G, N}` (case-insensitive) to/from bits.

### Encoding 🧬 ➡️ 🔟
Unlike the easy encoding step with four nucleotides, the five nucleotide encoding is slightly more involved. I worked on this
way back in 10th grade, but cute vector instructions made me revisit this. Anyways, the simple naive method is to encode each
nucleotide to three bits instead of two bits. However, this is a bit of a waste given that an extra bit is necessary to
represent a single extra nucleotide. A better way is to encode three nucleotides into seven bits. The seven bit encoding is
just the number after the three nucleotides are mapped to numbers from 0 to 4 and used as digits in a base-5 number. For example,
given the three nucleotides `a`, `b`, and `c`, we get `encoding = c * 5^2 + b * 5^1 + a * 5^0`. This gives a number from 0 to 124,
which fits into 7 bits. Its possible to use the rest of the numbers 125-127 to indicate the end of an encoded sequence, but we take
the easy route and just separately save the length of the encoded sequence. In summary, we convert triplets of nucleotides to 7
bits each, and then pack 9 of those 7-bit chunks into a single 64-bit integer. For every three nucleotides, this takes only
7 bits, instead of 9 bits with the naive method, where each nucleotide is encoded as three bits. Here is an example:
```
Nucleotides =        A         N         G
Digits      =        0         4         3
Base-5      =                          043
Base-10     =  0 * 5^2 + 4 * 5^1 + 3 * 5^0
            =                           23
Base-2      =                      0010111
```
We will now discuss the implemented algorithms:

* **n_to_bits2_lut**. This is a straightforward scalar implementation: read in three nucleotides, use the equation
`encoding = c * 5^2 + b * 5^1 + a * 5^0`, then write out 7 bits to the current 64-bit integer.

* **n_to_bits2_pext (AVX2 and BMI2)**. Alright, this gets a bit complicated. First of all, we map 27 nucleotide ASCII characters
to integers between 0 and 4. After this, the easiest method is to do some sort of transpose, separating each nucleotide from each
triplet into different vectors. We would have three vectors of 9 nucleotides each, and we can multiply each vector by `5^2`, `5^1`,
or `5^0`, and then vertically add them all to calculate the encoding based on the formula. However, there is a triplet of nucleotides
that cross the 128-bit lanes, so we cannot directly shuffle bytes to line up the three nucleotides in each vector for vertical
addition, as the `_mm256_shuffle_epi8` instruction cannot cross lanes. This means that we need a lane-crossing permute with
`_mm256_permutevar8x32_epi32` to move bytes across the lane boundary. Finally, we need to ensure that each nucleotide is padded to
16 bits, because there is no straightforward 8-bit vector multiply instruction. Once we get our encoded 7 bits for every three nucleotides,
we can use `_pext_u64` to pack 9 of those 7-bit chunks into a single 64-bit integer.

    Sounds pretty cool! This seems pretty much optimal, right? Wrong. I realized that the `_mm256_addubs_epi16` instruction existed, which
    does vertical multiplies and then adds pairs of adjacent products horizontally. This means that instead of splitting each 27 byte vector
    into three vectors of 9 bytes each, we can split it into two and use `_mm256_addubs_epi16` to multiply one vector by both `5^2` and `5^1`
    and add horizontally (`c * 5^2 + b * 5^1`) in a single step, and then vertically adding the other vector (`a * 5^0`). For a simplified example,
    consider the following, which denotes each vector as a list of bytes:

    ```
    First,

        shuffle [f, e, d, c, b, a]
     -> [0, 0, f, e, c, b] and [0, 0, 0, d, 0, a]

    Next,

        _mm256_addubs_epi16([0, 0, f, e, c, b], [0, 0, 5^2, 5^1, 5^2, 5^1])
      = [0, 0, f * 5^2, e * 5^1, c * 5^2, b * 5^1]
     -> [0, 0, 0, f * 5^2 + e * 5^1, 0, c * 5^2 + b * 5^1]

    Finally,

        [0, 0, 0, f * 5^2 + e * 5^1, 0, c * 5^2 + b * 5^1] + [0, 0, 0, d, 0, a]
      = [0, 0, 0, f * 5^2 + e * 5^1 + d, 0, c * 5^2 + b * 5^1 + a]
    ```

    Each character represents an arbitrary byte. This technique leads to a nice 11% speedup over transposing into three vectors of 9
    nucleotide bytes each. Unfortunately, since we are dealing with with chunks of 27 bytes with AVX2 vectors of at most 32 bytes,
    we have to use awkward unaligned loads. On modern CPUs, unaligned vector loads/stores to memory locations not a multiple of 32
    shouldn't result in massive speed penalties. It might be possible to use multiplication or bit shifts instead of `_pext_u64`, but
    they may only yield speedups when `_pext_u64` is very slow on a certain CPU.

### Decoding 🔟 ➡️ 🧬
Here is an example of decoding a triplet of nucleotides:
```
Base-2      =                      0010111
Base-10     =                           23
            =  0 * 5^2 + 4 * 5^1 + 3 * 5^0
Base-5      =                          043
Digits      =        0         4         3
Nucleotides =        A         N         G
```
When decoding a triplet of nucleotides, we need integer divisions and modulos:
```
encoding = c * 5^2 + b * 5^1 + a * 5^0
a = encoding % 5
b = (encoding / 5^1) % 5
c = encoding / 5^2
```
Absolutely not fun at all. Divisions and modulos are very, very slow compared to addition or multiplication. The AVX2 instruction set
does not even have integer division operations. To make use of floating point division, the each integer has to be converted into 32-bit
floats, which also isn't fun. There is a solution, which we will see very soon. Our implemented algorithms are the following:

* **bits_to_n2_lut**. Very basic implementation, where 7 bits are read in and we do division/modulos to get three nucleotides out. A
simple lookup table is used to convert the integers (from 0 to 4) to characters for each nucleotide.

* **bits_to_n2_pdep (AVX2 and BMI2)**. First, we use pdep to distribute each 7-bit chunk to 8-bit chunks. This will allow byte-level shuffles to work.
We place five 8-bit chunks in the high lane and four 8-bit chunks in the low lane of a 256-bit vector. This will help us avoid
lane-crossing operations. Additionally, we zero-extend 8-bit chunks to 16-bit chunks. Now, it is time to somehow do vector divisions and
modulos. For division, the idea is to use multiplication by a reciprocal (represented as a fixed-point number) instead of division:

    ```
    n / d = high_16_bits(n * d'), where
    d' = 2^16 / d + 1
    ```

    Note that `d'` is a constant value, since in this case `d = 5^1` or `d = 5^2`, and multiplying two 16-bit integers results in a 32-bit
    product. For division, we can multiply the quotient by the divisor and subtract that product from the original `n`:

    ```
    n % d = n - high_16_bits(n * d') * d
    ```

    However, Daniel Lemire [proposed](https://lemire.me/blog/2019/02/08/faster-remainders-when-the-divisor-is-a-constant-beating-compilers-and-libdivide/)
    a faster and more direct method, which uses the low bits of the product after multiplying by the reciprocal of the divisor in order to divide:

    ```
    n % d = high_16_bits(low_16_bits(n * d') * d)
    ```

    AVX2 conveniently provides the `_mm256_mullo_epi16` and `_mm256_mulhi_epu16` instructions, which allow for easy truncation of the high/low
    bits of the product. Using Daniel Lemire's method provides a ~5% speedup. The last step is to interleave the three vectors that were created
    through these division/modulo operations, use `_mm256_permutevar8x32_epi32` to pack the bytes into a contiguous 27-byte chunk, and map the
    nucleotide integers from 0 to 4 to actual ASCII characters.

## Cute speedups in cute micro-benchmarks 📈
All benchmarks were ran on a Intel Core i9-9880H (Coffee Lake-H) CPU with a clock rate of 2.3 GHz (4.8 GHz turbo boost). Run times were measured
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

The winner is `n_to_bits_movemask`, beating even `memcpy` (using `std::ptr::copy_nonoverlapping` to just make a copy of the entire byte string).
Note that all methods (including `memcpy`) allocate memory for the output first, which actually takes a nontrivial amount of time. This is slightly unfair
to `memcpy`, since it must allocate more memory to fit the entire byte string. In practice, `n_to_bits_mul` is probably the most flexible method for different
platforms that is fast enough.

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
vectors. Oh well.

| Undetermined nucleotides | Throughput       |
|--------------------------|------------------|
| n_to_bits2_lut           | 1.7887 GiB/s     |
| **n_to_bits2_pext**      | **11.787 GiB/s** |
| bits_to_n2_lut           | 1.3751 GiB/s     |
| **bits_to_n2_pdep**      | **10.175 GiB/s** |

For encoding and decoding undetermined nucleotides, the vectorized algorithms provide clear speedups over the lookup-table-based
naive methods.

## Discussion
Let's do a quick back-of-the-envelop calculation: processing data at ~30 GiB/s means that we only spend `1s / (30GiB / 32B) = ~1 ns` for
every 32-byte chunk of data. One nanosecond is 2-5 cycles for a 2.3 GHz clock rate (4.8 GHz turbo boost) CPU. This
should represent the time taken for one iteration of the `bits_to_n_shuffle` [loop](src/n_to_bits.rs#L280-L298). If we make the overly simplistic assumption
that one instruction takes one cycle to execute in sequence, and everything is in the CPU cache, then we expect each iteration to take maybe 9-10+ cycles, due to all
the operations we are doing! How could this be?

A high-level answer is that since each iteration in our loop operates on independent data, and there are a large number of iterations, we only
need to care about throughput (how many cycles it takes before we are able to launch a new instruction), instead of latency (how many cycles it
takes to get the output data of some instruction). The CPU can execute other instructions while waiting for some instructions to finish and output
their result. CPUs can also magically start multiple instructions per cycle. This allows it to be super fast.

More specifically, modern CPUs have many execution ports that handle different micro-operations (each instruction can be one or more micro-operations).
Some common instructions, like bitwise AND, OR, or addition, are a single micro-operation that can be executed by one of many ports, while vector shuffles are bottlenecked
on a single port on recent Intel CPUs. Each execution port can start a micro-operation every cycle, but multiple micro-operations can be launched in the same clock
cycle if they use different ports. This means that is may be possible to do, for example, four vector additions in the same cycle on recent Intel CPUs,
if there are no data dependencies. However, there is also a limit on the number of micro-operations that can be executed each cycle. On recent
Intel CPUs, around four micro-operations total can be executed each cycle. Also, it is possible for combinations of instructions, like load-and-then-vector-add, to be fused into
one micro-operation. You get the idea: the reason why a piece of code is fast or slow in practice is very complicated. However, we can reasonably
assume that `bits_to_n_shuffle`'s iterations are fast due to the CPU's out-of-order execution capability across multiple iterations of the loop.
In my code, I did not focus too heavily on optimizing CPU out-of-order-execution, but instead I focused more on general algorithms.

Usually, you do not need to worry about these details (or other more subtle CPU magic). Cache misses, unpredictable branches, long dependency chains,
lots of (not inlined) function calls, wasteful memory allocations, subpar algorithms, etc. are much more significant bottlenecks in most programs.

Note: In general, looking at the Rust SIMD code won't do us any good, since LLVM might output different assembly instructions and apply different
optimizations.

## Conclusion
This was an insightful experiment that unvealed many bit twiddling tricks. These tricks should lead to nontrivial memory savings
without a significant speed penalty. Unfortunately, its probably not possible to reach the lofty 30 GiB/s throughput in practice,
due to memory access bottlenecks.

Note that there is another interesting encoding technique that I have explored before: mapping `{A, T, C, G} -> {110, 011, 000, 101}`. These
encoding bits are *pairwise equidistant* under bitwise Hamming distance (every pair of encodings differ by two bits). This allows
nucleotide Hamming distance (number of different nucleotides) calculations to be done with just three instructions: an XOR, a popcount,
and a right shift to divide by two. For short strings, this is highly efficient when many Hamming distance calculations are needed.
