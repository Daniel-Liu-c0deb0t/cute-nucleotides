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

    let end_idx = if len < (bits.len() << 5) {
        bits.len() - 1
    } else {
        bits.len()
    };

    let mut res = Vec::with_capacity(len);

    unsafe {
        for i in 0..end_idx {
            let curr = *bits.get_unchecked(i);

            for j in 0..32 {
                res.push(*BITS_LUT.get_unchecked(((curr >> (j << 1)) & 0b11) as usize));
            }
        }

        if (end_idx << 5) < len {
            let curr = *bits.get_unchecked(end_idx);

            for i in 0..(len - (end_idx << 5)) {
                res.push(*BITS_LUT.get_unchecked(((curr >> (i << 1)) & 0b11) as usize));
            }
        }
    }

    res
}

union AlignedArray {
    v: __m256i,
    a: [u64; 4]
}

pub fn n_to_bits_pext(n: &[u8]) -> Vec<u64> {
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5;
    let mut res = Vec::with_capacity(end_idx + if n.len() & 31 == 0 {0} else {1});

    let ascii_mask = 0x0606060606060606; // 0b...00000110

    unsafe {
        let mut arr = AlignedArray{v: _mm256_undefined_si256()};

        for i in 0..end_idx as isize {
            // fast conversion of unaligned data to aligned
            arr.v = _mm256_loadu_si256(ptr.offset(i));

            // ascii_mask uses a special property of ATCG ASCII characters in binary
            // hide latency
            let a = _pext_u64(arr.a[0], ascii_mask);
            let b = _pext_u64(arr.a[1], ascii_mask);
            let c = _pext_u64(arr.a[2], ascii_mask);
            let d = _pext_u64(arr.a[3], ascii_mask);

            // combine low 16 bits in each 64 bit chunk
            res.push(a | (b << 16) | (c << 32) | (d << 48));
        }

        if n.len() & 31 > 0 {
            res.push(*n_to_bits_lut(&n[end_idx..]).get_unchecked(0));
        }
    }

    res
}

pub fn n_to_bits_mul(n: &[u8]) -> Vec<u64> {
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5;
    let mut res = Vec::with_capacity(end_idx + if n.len() & 31 == 0 {0} else {1});

    unsafe {
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

        let mut arr = AlignedArray{v: _mm256_undefined_si256()};

        for i in 0..end_idx as isize {
            arr.v = _mm256_loadu_si256(ptr.offset(i));
            // mask out unimportant bits
            arr.v = _mm256_and_si256(arr.v, ascii_mask);
            // multiply to pack left exactly 4 nucleotides (8 bits)
            arr.v = _mm256_mullo_epi32(arr.v, mul_mask);
            // extract last 8 bits of every 32 bit integer
            arr.v = _mm256_shuffle_epi8(arr.v, shuffle_mask);
            // combine first 32 bits from both lanes
            res.push(arr.a[0] | (arr.a[2] << 32));
        }

        if n.len() & 31 > 0 {
            res.push(*n_to_bits_lut(&n[end_idx..]).get_unchecked(0));
        }
    }

    res
}

/*pub fn bits_to_n_shuffle(b: &[u64], len: usize) -> Vec<u8> {
    return unimplemented!();

    let end_idx = len >> 6;
    let mut res = vec![];

    let nybble_mask = 0x0F0F0F0F0F0F0F0F;

    unsafe {
        for i in 0..end_idx {
            let curr = *b.get_unchecked(i);
            let v = _mm256_loadu_si256(ptr.offset(i));
            
            let a = _mm256_shuffle_epi8(nybbles1, nybble_lut);
            let b = _mm256_shuffle_epi8(nybbles2, nybble_lut);
            let c = _mm256_shuffle_epi8(nybbles3, nybble_lut);
            let d = _mm256_shuffle_epi8(nybbles4, nybble_lut);
        }
    }

    // TODO: handle leftover bytes

    res
}*/

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
                b"ATCGATCGATCGATCGATCGATCGATCGATCG");
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
    fn test_bits_to_n_pdep() {
        assert_eq!(bits_to_n_pdep(&vec![0b1101100011011000110110001101100011011000110110001101100011011000], 32),
                b"ATCGATCGATCGATCGATCGATCGATCGATCG");
    }
}
