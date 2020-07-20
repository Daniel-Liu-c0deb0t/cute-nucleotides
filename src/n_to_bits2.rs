#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use std::alloc;

static BYTE_LUT: [u8; 128] = {
    let mut lut = [0u8; 128];
    lut[b'a' as usize] = 0b000;
    lut[b'c' as usize] = 0b001;
    lut[b'g' as usize] = 0b010;
    lut[b'n' as usize] = 0b011;
    lut[b't' as usize] = 0b100;
    lut[b'A' as usize] = 0b000;
    lut[b'C' as usize] = 0b001;
    lut[b'G' as usize] = 0b010;
    lut[b'N' as usize] = 0b011;
    lut[b'T' as usize] = 0b100;
    lut
};

/*static BITS_LUT: [u8; 4] = {
    let mut lut = [0u8; 4];
    lut[0b00] = b'A';
    lut[0b10] = b'T';
    lut[0b01] = b'C';
    lut[0b11] = b'G';
    lut
};*/

pub fn n_to_bits2_lut(n: &[u8]) -> Vec<u64> {
    let mut res = vec![0u64; (n.len() / 27) + if n.len() % 27 == 0 {0} else {1}];
    let len = n.len() / 3;

    unsafe {
        for i in 0..len {
            let idx = i * 3;
            let res_offset = i / 9;
            let res_shift = (i % 9) * 7;
            let a = *BYTE_LUT.get_unchecked(*n.get_unchecked(idx) as usize);
            let b = (*BYTE_LUT.get_unchecked(*n.get_unchecked(idx + 1) as usize)) * 5;
            let c = (*BYTE_LUT.get_unchecked(*n.get_unchecked(idx + 2) as usize)) * 25;
            let encoding = (a + b + c) as u64;
            *res.get_unchecked_mut(res_offset) = *res.get_unchecked(res_offset) | (encoding << res_shift);
        }

        let leftover = n.len() % 3;

        if leftover > 0 {
            let idx = len * 3;
            let res_offset = len / 9;
            let res_shift = (len % 9) * 7;
            let a = *BYTE_LUT.get_unchecked(*n.get_unchecked(idx) as usize);
            let b = if leftover >= 2 {(*BYTE_LUT.get_unchecked(*n.get_unchecked(idx + 1) as usize)) * 5} else {0};
            let encoding = (a + b) as u64;
            *res.get_unchecked_mut(res_offset) = *res.get_unchecked(res_offset) | (encoding << res_shift);
        }
    }

    res
}

/*pub fn bits_to_n_lut(bits: &[u64], len: usize) -> Vec<u8> {
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
}*/

union AlignedArray {
    v: __m256i,
    a: [u64; 4]
}

pub fn n_to_bits2_pext(n: &[u8]) -> Vec<u64> {
    let mut ptr = n.as_ptr();
    let end_idx = if n.len() < 5 {0} else {(n.len() - 5) / 27};
    let len = (n.len() / 27) + if n.len() % 27 == 0 {0} else {1};

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(len << 3, 8);
        let res_ptr = alloc::alloc(layout) as *mut u64;

        let lut = {
            let mut lut = 0;
            lut |= 0b000 << (((b'A' as i64) & 0b111) << 3);
            lut |= 0b001 << (((b'C' as i64) & 0b111) << 3);
            lut |= 0b010 << (((b'G' as i64) & 0b111) << 3);
            lut |= 0b011 << (((b'N' as i64) & 0b111) << 3);
            lut |= 0b100 << (((b'T' as i64) & 0b111) << 3);
            _mm256_set1_epi64x(lut)
        };
        let permute_mask = _mm256_set_epi32(6, 5, 4, 3, 3, 2, 1, 0);
        let shuffle_mask1 = _mm256_set_epi16(-1, -1, -1, -1, 28, 25, 22, 19, -1, -1, -1, 12,  9,  6,  3,  0);
        let shuffle_mask2 = _mm256_set_epi16(-1, -1, -1, -1, 29, 26, 23, 20, -1, -1, -1, 13, 10,  7,  4,  1);
        let shuffle_mask3 = _mm256_set_epi16(-1, -1, -1, -1, 30, 27, 24, 21, -1, -1, -1, 14, 11,  8,  5,  2);
        let mul5 = _mm256_set1_epi16(5);
        let mul25 = _mm256_set1_epi16(25);
        let pack_right_mask = 0x007F007F007F007Fu64; // 0b...0000000001111111

        let mut arr = [AlignedArray{v: _mm256_undefined_si256()}, AlignedArray{v: _mm256_undefined_si256()}];

        for i in 0..end_idx as isize {
            let v = _mm256_loadu_si256(ptr as *const __m256i);

            let v = _mm256_shuffle_epi8(lut, v);
            let v = _mm256_permutevar8x32_epi32(v, permute_mask);

            let a = _mm256_shuffle_epi8(v, shuffle_mask1);
            let b = _mm256_shuffle_epi8(v, shuffle_mask2);
            let c = _mm256_shuffle_epi8(v, shuffle_mask3);

            let b = _mm256_mullo_epi16(b, mul5);
            let c = _mm256_mullo_epi16(c, mul25);

            let ab = _mm256_add_epi16(a, b);
            let arr_idx = (i as usize) & 1;
            (*arr.get_unchecked_mut(arr_idx)).v = _mm256_add_epi16(ab, c);

            let a = _pext_u64((*arr.get_unchecked(arr_idx)).a[0], pack_right_mask);
            let b = (*arr.get_unchecked(arr_idx)).a[1];
            let c = _pext_u64((*arr.get_unchecked(arr_idx)).a[2], pack_right_mask);

            // combine a, b, and c into a 63 bit chunk
            *res_ptr.offset(i) = a | (b << 28) | (c << 35);

            ptr = ptr.offset(27);
        }

        if end_idx < len {
            let end = n_to_bits2_lut(&n[(end_idx * 27)..]);

            for i in 0..end.len() {
                *res_ptr.offset((end_idx + i) as isize) = *end.get_unchecked(i);
            }
        }

        Vec::from_raw_parts(res_ptr, len, len)
    }
}

/*pub fn bits_to_n_shuffle(bits: &[u64], len: usize) -> Vec<u8> {
    if len > (bits.len() << 5) {
        panic!("The length is greater than the number of nucleotides!");
    }

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(bits.len() << 5, 32);
        let ptr = alloc::alloc(layout) as *mut __m256i;

        let shuffle_mask = _mm256_set_epi32(0x07070707, 0x06060606, 0x05050505, 0x04040404, 0x03030303, 0x02020202, 0x01010101, 0x00000000);
        let shift1 = _mm256_set1_epi32(0);
        let shift2 = _mm256_set1_epi32(2);
        let shift3 = _mm256_set1_epi32(4);
        let shift4 = _mm256_set1_epi32(6);
        let blend_mask8 = _mm256_set1_epi16(0xFF00u16 as i16);
        let lo_mask = _mm256_set1_epi8(0b00000011);
        let lut_i32 = (b'A' as i32) | ((b'C' as i32) << 8) | ((b'T' as i32) << 16) | ((b'G' as i32) << 24);
        let lut = _mm256_set_epi32(0, 0, 0, lut_i32, 0, 0, 0, lut_i32);

        for i in 0..bits.len() {
            let curr = *bits.get_unchecked(i) as i64;
            let v = _mm256_set1_epi64x(curr);
            // duplicate each byte four times
            let v = _mm256_shuffle_epi8(v, shuffle_mask);
            // separately right shift each 32-bit chunk by 0, 2, 4, and 6 bits
            let a = _mm256_srlv_epi32(v, shift1);
            let b = _mm256_srlv_epi32(v, shift2);
            let c = _mm256_srlv_epi32(v, shift3);
            let d = _mm256_srlv_epi32(v, shift4);
            // merge together shifted chunks
            let ab = _mm256_blendv_epi8(a, b, blend_mask8);
            let cd = _mm256_blendv_epi8(c, d, blend_mask8);
            let abcd = _mm256_blend_epi16(ab, cd, 0b10101010i32);
            // only keep the two low bits per byte
            let v = _mm256_and_si256(abcd, lo_mask);
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
}*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n_to_bits2_lut() {
        assert_eq!(n_to_bits2_lut(b"ATCGNATCGNATCGNATCGNATCGNATCGNATCGN"),
                vec![0b110011101110110010001010110110101101100111011101100100010101101, 0b1000101011011010110]);
        assert_eq!(n_to_bits2_lut(b"ATCGN"), vec![0b100010101101]);
    }

    /*#[test]
    fn test_bits_to_n_lut() {
        assert_eq!(bits_to_n_lut(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                b"ATCGATCGATCGATCGATCGATCGATCGATCG");
    }*/

    #[test]
    fn test_n_to_bits2_pext() {
        assert_eq!(n_to_bits2_pext(b"ATCGNATCGNATCGNATCGNATCGNATCGNATCGN"),
                vec![0b110011101110110010001010110110101101100111011101100100010101101, 0b1000101011011010110]);
        assert_eq!(n_to_bits2_pext(b"ATCGN"), vec![0b100010101101]);
    }

    /*#[test]
    fn test_n_to_bits_mul() {
        assert_eq!(n_to_bits_mul(b"ATCGATCGATCGATCGATCGATCGATCGATCG"),
                vec![0b1101100011011000110110001101100011011000110110001101100011011000]);
        assert_eq!(n_to_bits_mul(b"ATCG"), vec![0b11011000]);
    }

    #[test]
    fn test_bits_to_n_shuffle() {
        assert_eq!(bits_to_n_shuffle(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                b"ATCGATCGATCGATCGATCGATCGATCGATCG");
    }

    #[test]
    fn test_bits_to_n_pdep() {
        assert_eq!(bits_to_n_pdep(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                b"ATCGATCGATCGATCGATCGATCGATCGATCG");
    }

    #[test]
    fn test_bits_to_n_clmul() {
        assert_eq!(bits_to_n_clmul(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                b"ATCGATCGATCGATCGATCGATCGATCGATCG");
    }*/
}
