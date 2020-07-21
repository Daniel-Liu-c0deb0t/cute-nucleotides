#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use std::alloc;

static BYTE_LUT: [u8; 128] = {
    let mut lut = [0u8; 128];
    lut[b'a' as usize] = 0b00;
    lut[b't' as usize] = 0b10;
    lut[b'c' as usize] = 0b01;
    lut[b'g' as usize] = 0b11;
    lut[b'A' as usize] = 0b00;
    lut[b'T' as usize] = 0b10;
    lut[b'C' as usize] = 0b01;
    lut[b'G' as usize] = 0b11;
    lut
};

static BITS_LUT: [u8; 4] = {
    let mut lut = [0u8; 4];
    lut[0b00] = b'A';
    lut[0b10] = b'T';
    lut[0b01] = b'C';
    lut[0b11] = b'G';
    lut
};

pub fn n_to_bits_lut(n: &[u8]) -> Vec<u64> {
    let mut res = vec![0u64; (n.len() >> 5) + if n.len() & 31 == 0 {0} else {1}];

    unsafe {
        for i in 0..n.len() {
            let offset = i >> 5;
            let shift = (i & 31) << 1;
            *res.get_unchecked_mut(offset) = *res.get_unchecked(offset)
                | ((*BYTE_LUT.get_unchecked(*n.get_unchecked(i) as usize) as u64) << shift);
        }
    }

    res
}

pub fn bits_to_n_lut(bits: &[u64], len: usize) -> Vec<u8> {
    if len > (bits.len() << 5) {
        panic!("The length is greater than the number of nucleotides!");
    }

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(len, 1);
        let res_ptr = alloc::alloc(layout);

        for i in 0..len {
            let offset = i >> 5;
            let shift = (i & 31) << 1;
            let curr = *bits.get_unchecked(offset);
            *res_ptr.offset(i as isize) = *BITS_LUT.get_unchecked(((curr >> shift) & 0b11) as usize);
        }

        Vec::from_raw_parts(res_ptr, len, len)
    }
}

union AlignedArray {
    v: __m256i,
    a: [u64; 4]
}

pub fn n_to_bits_pext(n: &[u8]) -> Vec<u64> {
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5;
    let len = end_idx + if n.len() & 31 == 0 {0} else {1};

    let ascii_mask = 0x0606060606060606; // 0b...00000110

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(len << 3, 8);
        let res_ptr = alloc::alloc(layout) as *mut u64;

        let mut arr = [AlignedArray{v: _mm256_undefined_si256()}, AlignedArray{v: _mm256_undefined_si256()}];

        for i in 0..end_idx as isize {
            let arr_idx = (i as usize) & 1;
            // fast conversion of unaligned data to aligned
            (*arr.get_unchecked_mut(arr_idx)).v = _mm256_loadu_si256(ptr.offset(i));

            // ascii_mask uses a special property of ATCG ASCII characters in binary
            // hide latency
            let a = _pext_u64((*arr.get_unchecked(arr_idx)).a[0], ascii_mask);
            let b = _pext_u64((*arr.get_unchecked(arr_idx)).a[1], ascii_mask);
            let c = _pext_u64((*arr.get_unchecked(arr_idx)).a[2], ascii_mask);
            let d = _pext_u64((*arr.get_unchecked(arr_idx)).a[3], ascii_mask);

            // combine low 16 bits in each 64 bit chunk
            *res_ptr.offset(i) = a | (b << 16) | (c << 32) | (d << 48);
        }

        if n.len() & 31 > 0 {
            *res_ptr.offset(end_idx as isize) = *n_to_bits_lut(&n[(end_idx << 5)..]).get_unchecked(0);
        }

        Vec::from_raw_parts(res_ptr, len, len)
    }
}

pub fn n_to_bits_mul(n: &[u8]) -> Vec<u64> {
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5;
    let len = end_idx + if n.len() & 31 == 0 {0} else {1};

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(len << 3, 8);
        let res_ptr = alloc::alloc(layout) as *mut u64;

        let ascii_mask = _mm256_set1_epi8(0b00000110);
        let mul_mask = {
            let mut m = 0u32;
            // m |= 1 << (length - input byte offset + output bit offset - 1 LSB to ignore);
            m |= 1 << (32 -  8 + 0 - 1);
            m |= 1 << (32 - 16 + 2 - 1);
            m |= 1 << (32 - 24 + 4 - 1);
            m |= 1 << (32 - 32 + 6 - 1);
            _mm256_set1_epi32(m as i32)
        };
        let shuffle_mask = _mm256_set_epi32(-1, -1, -1, 0x0F0B0703, -1, -1, -1, 0x0F0B0703);

        let mut arr = [AlignedArray{v: _mm256_undefined_si256()}, AlignedArray{v: _mm256_undefined_si256()}];

        for i in 0..end_idx as isize {
            let v = _mm256_loadu_si256(ptr.offset(i));
            // mask out unimportant bits
            let v = _mm256_and_si256(v, ascii_mask);
            // multiply to pack left exactly 4 nucleotides (8 bits)
            let v = _mm256_mullo_epi32(v, mul_mask);
            let arr_idx = (i as usize) & 1;
            // extract last 8 bits of every 32 bit integer
            (*arr.get_unchecked_mut(arr_idx)).v = _mm256_shuffle_epi8(v, shuffle_mask);
            // combine first 32 bits from both lanes
            *res_ptr.offset(i) = (*arr.get_unchecked(arr_idx)).a[0] | ((*arr.get_unchecked(arr_idx)).a[2] << 32);
        }

        if n.len() & 31 > 0 {
            *res_ptr.offset(end_idx as isize) = *n_to_bits_lut(&n[(end_idx << 5)..]).get_unchecked(0);
        }

        Vec::from_raw_parts(res_ptr, len, len)
    }
}

pub fn bits_to_n_shuffle(bits: &[u64], len: usize) -> Vec<u8> {
    if len > (bits.len() << 5) {
        panic!("The length is greater than the number of nucleotides!");
    }

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(bits.len() << 5, 32);
        let ptr = alloc::alloc(layout) as *mut __m256i;

        let shuffle_mask = _mm256_set_epi32(0x07070707, 0x06060606, 0x05050505, 0x04040404, 0x03030303, 0x02020202, 0x01010101, 0x00000000);
        let lo_mask = _mm256_set1_epi16(0b0000110000000011);
        let lut_i32 = (b'A' as i32) | ((b'C' as i32) << 8) | ((b'T' as i32) << 16) | ((b'G' as i32) << 24);
        let lut = _mm256_set_epi32(b'G' as i32, b'T' as i32, b'C' as i32, lut_i32, b'G' as i32, b'T' as i32, b'C' as i32, lut_i32);

        for i in 0..bits.len() {
            let curr = *bits.get_unchecked(i) as i64;
            let v = _mm256_set1_epi64x(curr);
            // duplicate each byte four times
            let v1 = _mm256_shuffle_epi8(v, shuffle_mask);
            // separately right shift each 16-bit chunk by 0 or 4 bits
            let v2 = _mm256_srli_epi16(v1, 4);
            // merge together shifted chunks
            let v = _mm256_blend_epi16(v1, v2, 0b10101010i32);
            // only keep two bits in each byte
            // either 0b0011 or 0b1100
            let v = _mm256_and_si256(v, lo_mask);
            // use lookup table to convert nucleotide bits to bytes
            let v = _mm256_shuffle_epi8(lut, v);
            _mm256_store_si256(ptr.offset(i as isize), v);
        }

        Vec::from_raw_parts(ptr as *mut u8, len, bits.len() << 5)
    }
}

pub fn bits_to_n_pdep(bits: &[u64], len: usize) -> Vec<u8> {
    if len > (bits.len() << 5) {
        panic!("The length is greater than the number of nucleotides!");
    }

    let scatter_mask = 0x0303030303030303u64;

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(bits.len() << 5, 32);
        let ptr = alloc::alloc(layout) as *mut __m256i;

        let lut_i32 = (b'A' as i32) | ((b'C' as i32) << 8) | ((b'T' as i32) << 16) | ((b'G' as i32) << 24);
        let lut = _mm256_set_epi32(0, 0, 0, lut_i32, 0, 0, 0, lut_i32);

        for i in 0..bits.len() {
            let curr = *bits.get_unchecked(i);
            // spread out nucleotide bits to first 2 bits of each byte
            let a = _pdep_u64(curr, scatter_mask) as i64;
            let b = _pdep_u64(curr >> 16, scatter_mask) as i64;
            let c = _pdep_u64(curr >> 32, scatter_mask) as i64;
            let d = _pdep_u64(curr >> 48, scatter_mask) as i64;
            let v = _mm256_set_epi64x(d, c, b, a);
            // lookup table from nucleotide bits to bytes
            let v = _mm256_shuffle_epi8(lut, v);
            _mm256_store_si256(ptr.offset(i as isize), v);
        }

        Vec::from_raw_parts(ptr as *mut u8, len, bits.len() << 5)
    }
}

pub fn bits_to_n_clmul(bits: &[u64], len: usize) -> Vec<u8> {
    if len > (bits.len() << 5) {
        panic!("The length is greater than the number of nucleotides!");
    }

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(bits.len() << 5, 16);
        let ptr = alloc::alloc(layout) as *mut __m128i;

        let lo_shuffle_mask = _mm_set_epi32(0xFFFFFF03u32 as i32, 0xFFFFFF02u32 as i32, 0xFFFFFF01u32 as i32, 0xFFFFFF00u32 as i32);
        let hi_shuffle_mask = _mm_set_epi32(0xFFFFFF07u32 as i32, 0xFFFFFF06u32 as i32, 0xFFFFFF05u32 as i32, 0xFFFFFF04u32 as i32);
        let mul_mask = {
            let mut m = 0u64;
            // m |= 1 << (byte offset - bit offset);
            m |= 1 << ( 0 - 0);
            m |= 1 << ( 8 - 2);
            m |= 1 << (16 - 4);
            m |= 1 << (24 - 6);
            _mm_set_epi64x(0, m as i64)
        };
        let lo_mask = _mm_set1_epi8(0b00000011);
        let lut_i32 = (b'A' as i32) | ((b'C' as i32) << 8) | ((b'T' as i32) << 16) | ((b'G' as i32) << 24);
        let lut = _mm_set1_epi32(lut_i32);

        for i in 0..bits.len() {
            let curr = *bits.get_unchecked(i) as i64;
            let v = _mm_set1_epi64x(curr);
            // spread out bytes to the low 8 bits of each 32 bit chunk
            let lo_v = _mm_shuffle_epi8(v, lo_shuffle_mask);
            let hi_v = _mm_shuffle_epi8(v, hi_shuffle_mask);
            // multiply by mask to shift to correct positions
            // carry-less multiply will ensure that separate bytes do not interfere with each other
            // handle 64 bit chunks separately
            let lo_v1 = _mm_clmulepi64_si128(lo_v, mul_mask, 0x00);
            let lo_v2 = _mm_clmulepi64_si128(lo_v, mul_mask, 0x0F);
            let hi_v1 = _mm_clmulepi64_si128(hi_v, mul_mask, 0x00);
            let hi_v2 = _mm_clmulepi64_si128(hi_v, mul_mask, 0x0F);
            // combine the two low chunks of 64 bits into 128 bit vectors
            // casts are free
            let lo_v = _mm_castps_si128(_mm_movelh_ps(_mm_castsi128_ps(lo_v1), _mm_castsi128_ps(lo_v2)));
            let hi_v = _mm_castps_si128(_mm_movelh_ps(_mm_castsi128_ps(hi_v1), _mm_castsi128_ps(hi_v2)));
            // only keep low bits
            let lo_v = _mm_and_si128(lo_v, lo_mask);
            let hi_v = _mm_and_si128(hi_v, lo_mask);
            // use lookup table to convert nucleotide bits to bytes
            let lo_v = _mm_shuffle_epi8(lut, lo_v);
            let hi_v = _mm_shuffle_epi8(lut, hi_v);
            _mm_store_si128(ptr.offset((i << 1) as isize), lo_v);
            _mm_store_si128(ptr.offset(((i << 1) + 1) as isize), hi_v);
        }

        Vec::from_raw_parts(ptr as *mut u8, len, bits.len() << 5)
    }
}

// A = 00, T = 10, C = 01, G = 11

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n_to_bits_lut() {
        assert_eq!(n_to_bits_lut(b"ATCGATCGATCGATCGATCGATCGATCGATCG"),
                vec![0b1101100011011000110110001101100011011000110110001101100011011000]);
        assert_eq!(n_to_bits_lut(b"ATCG"), vec![0b11011000]);
    }

    #[test]
    fn test_bits_to_n_lut() {
        assert_eq!(bits_to_n_lut(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                "ATCGATCGATCGATCGATCGATCGATCGATCG".as_bytes());
    }

    #[test]
    fn test_n_to_bits_pext() {
        assert_eq!(n_to_bits_pext(b"ATCGATCGATCGATCGATCGATCGATCGATCG"),
                vec![0b1101100011011000110110001101100011011000110110001101100011011000]);
        assert_eq!(n_to_bits_pext(b"ATCG"), vec![0b11011000]);
    }

    #[test]
    fn test_n_to_bits_mul() {
        assert_eq!(n_to_bits_mul(b"ATCGATCGATCGATCGATCGATCGATCGATCG"),
                vec![0b1101100011011000110110001101100011011000110110001101100011011000]);
        assert_eq!(n_to_bits_mul(b"ATCG"), vec![0b11011000]);
    }

    #[test]
    fn test_bits_to_n_shuffle() {
        assert_eq!(bits_to_n_shuffle(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                "ATCGATCGATCGATCGATCGATCGATCGATCG".as_bytes());
    }

    #[test]
    fn test_bits_to_n_pdep() {
        assert_eq!(bits_to_n_pdep(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                "ATCGATCGATCGATCGATCGATCGATCGATCG".as_bytes());
    }

    #[test]
    fn test_bits_to_n_clmul() {
        assert_eq!(bits_to_n_clmul(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                "ATCGATCGATCGATCGATCGATCGATCGATCG".as_bytes());
    }
}
