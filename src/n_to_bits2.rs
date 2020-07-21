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
        let shuffle_mask1 = _mm256_set_epi16(-1, -1, -1, -1, 0xFF1Cu16 as i16, 0xFF19u16 as i16, 0xFF16u16 as i16, 0xFF13u16 as i16,
                -1, -1, -1, 0xFF0Cu16 as i16, 0xFF09u16 as i16, 0xFF06u16 as i16, 0xFF03u16 as i16, 0xFF00u16 as i16);
        let shuffle_mask2 = _mm256_set_epi16(-1, -1, -1, -1, 0xFF1Du16 as i16, 0xFF1Au16 as i16, 0xFF17u16 as i16, 0xFF14u16 as i16,
                -1, -1, -1, 0xFF0Du16 as i16, 0xFF0Au16 as i16, 0xFF07u16 as i16, 0xFF04u16 as i16, 0xFF01u16 as i16);
        let shuffle_mask3 = _mm256_set_epi16(-1, -1, -1, -1, 0xFF1Eu16 as i16, 0xFF1Bu16 as i16, 0xFF18u16 as i16, 0xFF15u16 as i16,
                -1, -1, -1, 0xFF0Eu16 as i16, 0xFF0Bu16 as i16, 0xFF08u16 as i16, 0xFF05u16 as i16, 0xFF02u16 as i16);
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

pub fn bits_to_n2_pdep(bits: &[u64], len: usize) -> Vec<u8> {
    if len > (bits.len() * 27) {
        panic!("The length is greater than the number of nucleotides!");
    }

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(bits.len() * 27 + 5, 32);
        let res_ptr = alloc::alloc(layout);
        let mut ptr = res_ptr as *mut __m256i;

        let deposit_mask = 0x7F7F7F7F7F7F7F7Fu64; // 0b...01111111
        let shuffle_mask = _mm256_set_epi16(-1, -1, -1, 0xFF04u16 as i16, 0xFF03u16 as i16, 0xFF02u16 as i16, 0xFF01u16 as i16, 0xFF00u16 as i16,
                -1, -1, -1, -1, 0xFF03u16 as i16, 0xFF02u16 as i16, 0xFF01u16 as i16, 0xFF00u16 as i16);
        let mul5 = _mm256_set1_epi16(5);
        let div5 = _mm256_set1_epi16((((i16::MAX as u32) + 1) / 5u32) as i16);
        let div25 = _mm256_set1_epi16((((i16::MAX as u32) + 1) / 25u32) as i16);
        let ab_shuffle_mask = _mm256_set_epi64x(0xFFFF0908FF0706FFu64 as i64, 0x0504FF0302FF0100u64 as i64,
                0xFFFF0908FF0706FFu64 as i64, 0x0504FF0302FF0100u64 as i64);
        let c_shuffle_mask = _mm256_set_epi64x(0xFF08FFFF06FFFF04u64 as i64, 0xFFFF02FFFF00FFFFu64 as i64,
                0xFF08FFFF06FFFF04u64 as i64, 0xFFFF02FFFF00FFFFu64 as i64);
        let permute_mask = _mm256_set_epi32(7, 7, 6, 5, 4, 2, 1, 0);
        let lut = {
            let mut lut = 0;
            lut |= (b'A' as i64) <<  0;
            lut |= (b'C' as i64) <<  8;
            lut |= (b'G' as i64) << 16;
            lut |= (b'N' as i64) << 24;
            lut |= (b'T' as i64) << 32;
            _mm256_set1_epi64x(lut)
        };

        for i in 0..bits.len() {
            let curr = *bits.get_unchecked(i) as i64;
            let a = _pdep_u64(curr as u64, deposit_mask) as i64;
            let b = ((curr >> 56) << 32) | (a >> 32);

            let v = _mm256_set_epi64x(0, b, 0, a);

            let v = _mm256_shuffle_epi8(v, shuffle_mask);

            let v_div1 = v;
            let v_div5 = _mm256_mulhrs_epi16(v_div1, div5);
            let v_div25 = _mm256_mulhrs_epi16(v_div1, div25);
            // v[i] = c[i] * 5^2 + b[i] * 5^1 + a[i] * 5^0
            // a[i] = (c[i] * 5^2 + b[i] * 5^1 + a[i] * 5^0) - (c[i] * 5^2 + b[i] * 5^1)
            let a = _mm256_sub_epi16(v_div1, _mm256_mullo_epi16(v_div5, mul5));
            // b[i] = (c[i] * 5^1 + b[i] * 5^0) - (c[i] * 5^1)
            let b = _mm256_sub_epi16(v_div5, _mm256_mullo_epi16(v_div25, mul5));
            // c[i] = c[i] * 5^0
            let c = v_div25;

            let arr = AlignedArray{v: a}; println!("{:?}", arr.a);
            let arr = AlignedArray{v: b}; println!("{:?}", arr.a);
            let b = _mm256_slli_epi16(b, 8);
            let ab = _mm256_or_si256(a, b);
            let ab = _mm256_shuffle_epi8(ab, ab_shuffle_mask);
            let c = _mm256_shuffle_epi8(c, c_shuffle_mask);
            let abc = _mm256_or_si256(ab, c);

            let v = _mm256_permutevar8x32_epi32(abc, permute_mask);
            let v = _mm256_shuffle_epi8(lut, v);

            _mm256_storeu_si256(ptr as *mut __m256i, v);
            ptr = ptr.offset(27);
        }

        Vec::from_raw_parts(res_ptr, len, bits.len() * 27 + 5)
    }
}

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

    #[test]
    fn test_bits_to_n2_pdep() {
        assert_eq!(bits_to_n2_pdep(&vec![0b110011101110110010001010110110101101100111011101100100010101101, 0b1000101011011010110], 35),
                "ATCGNATCGNATCGNATCGNATCGNATCGNATCGN".as_bytes());
    }
}
