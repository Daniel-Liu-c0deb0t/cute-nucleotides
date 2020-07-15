#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

pub fn n_to_bits_lut(n: &[u8]) -> Vec<u64> {
    unimplemented!()
}

pub fn n_to_bits_hash(n: &[u8]) -> Vec<u64> {
    unimplemented!()
}

union AlignedArray {
    v: __m256i,
    a: [u64; 4]
}

pub fn n_to_bits_pext(n: &[u8]) -> Vec<u64> {
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5;
    let mut res = Vec::with_capacity(end_idx + if n.len() & 31 == 0 {0} else {1});

    let mask = 0x0606060606060606; // 0b...00000110

    unsafe {
        let mut arr = AlignedArray{v: _mm256_undefined_si256()};

        for i in 0..end_idx as isize {
            // fast conversion of unaligned data to aligned
            arr.v = _mm256_loadu_si256(ptr.offset(i));

            // mask uses a special property of ATCG ASCII characters in binary
            // hide latency
            let a = _pext_u64(arr.a[0], mask);
            let b = _pext_u64(arr.a[1], mask);
            let c = _pext_u64(arr.a[2], mask);
            let d = _pext_u64(arr.a[3], mask);

            // combine low 16 bits in each 64 bit chunk
            res.push(a | (b << 16) | (c << 32) | (d << 48));
        }
    }

    // TODO: handle leftover bytes

    res
}

pub fn n_to_bits_mul(n: &[u8]) -> Vec<u64> {
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5;
    let mut res = Vec::with_capacity(end_idx + if n.len() & 31 == 0 {0} else {1});

    unsafe {
        let mask = _mm256_set1_epi8(0b00000110);
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
            arr.v = _mm256_and_si256(arr.v, mask);
            // multiply and truncate high bits to move left
            arr.v = _mm256_mullo_epi32(arr.v, mul_mask);
            // extract last 8 bits of every 32 bit integer
            arr.v = _mm256_shuffle_epi8(arr.v, shuffle_mask);
            // combine first 32 bits from both lanes
            res.push(arr.a[0] | (arr.a[2] << 32));
        }
    }

    // TODO: handle leftover bytes

    res
}

// A = 00, T = 10, C = 01, G = 11

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n_to_bits_pext() {
        assert_eq!(n_to_bits_pext(b"ATCGATCGATCGATCGATCGATCGATCGATCG"), vec![0b1101100011011000110110001101100011011000110110001101100011011000]);
    }

    #[test]
    fn test_n_to_bits_mul() {
        assert_eq!(n_to_bits_mul(b"ATCGATCGATCGATCGATCGATCGATCGATCG"), vec![0b1101100011011000110110001101100011011000110110001101100011011000]);
    }
}
