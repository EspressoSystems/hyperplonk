//! useful macros.

/// Takes as input a struct, and converts them to a series of bytes. All traits
/// that implement `CanonicalSerialize` can be automatically converted to bytes
/// in this manner.
#[macro_export]
macro_rules! to_bytes {
    ($x:expr) => {{
        let mut buf = ark_std::vec![];
        ark_serialize::CanonicalSerialize::serialize($x, &mut buf).map(|_| buf)
    }};
}

/// Decompose an integer into a binary vector in little endian.
#[allow(dead_code)]
pub(crate) fn bit_decompose(input: u64, num_var: usize) -> Vec<bool> {
    let mut res = Vec::with_capacity(num_var);
    let mut i = input;
    for _ in 0..num_var {
        res.push(i & 1 == 1);
        i >>= 1;
    }
    res
}

/// Project a little endian binary vector into an integer.
#[allow(dead_code)]
pub(crate) fn project(input: &[bool]) -> u64 {
    let mut res = 0;
    for &e in input.iter().rev() {
        res <<= 1;
        res += e as u64;
    }
    res
}

// Input index
// - `i := (i_0, ...i_{n-1})`,
// - `num_vars := n`
// return three elements:
// - `a:= (i_1, ..., i_{n-1}, 0)`
// - `b:= (i_1, ..., i_{n-1}, 1)`
// - `sign:= i_0`
#[inline]
pub(crate) fn get_index(i: usize, num_vars: usize) -> (usize, usize, bool) {
    let bit_sequence = bit_decompose(i as u64, num_vars);
    let x0 = project(&[bit_sequence[1..].as_ref(), [false].as_ref()].concat()) as usize;
    let x1 = project(&[bit_sequence[1..].as_ref(), [true].as_ref()].concat()) as usize;

    (x0, x1, bit_sequence[0])
}

#[cfg(test)]
mod test {
    use ark_bls12_381::Fr;
    use ark_serialize::CanonicalSerialize;
    use ark_std::{rand::RngCore, test_rng, One};

    use super::{bit_decompose, project};

    #[test]
    fn test_to_bytes() {
        let f1 = Fr::one();

        let mut bytes = ark_std::vec![];
        f1.serialize(&mut bytes).unwrap();
        assert_eq!(bytes, to_bytes!(&f1).unwrap());
    }

    #[test]
    fn test_decomposition() {
        let mut rng = test_rng();
        for _ in 0..100 {
            let t = rng.next_u64();
            let b = bit_decompose(t, 64);
            let r = project(&b);
            assert_eq!(t, r)
        }
    }
}
