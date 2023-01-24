// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

use ark_ff::PrimeField;
use ark_std::log2;

/// Decompose an integer into a binary vector in little endian.
pub fn bit_decompose(input: u64, num_var: usize) -> Vec<bool> {
    let mut res = Vec::with_capacity(num_var);
    let mut i = input;
    for _ in 0..num_var {
        res.push(i & 1 == 1);
        i >>= 1;
    }
    res
}

/// given the evaluation input `point` of the `index`-th polynomial,
/// obtain the evaluation point in the merged polynomial
pub fn gen_eval_point<F: PrimeField>(index: usize, index_len: usize, point: &[F]) -> Vec<F> {
    let index_vec: Vec<F> = bit_decompose(index as u64, index_len)
        .into_iter()
        .map(|x| F::from(x))
        .collect();
    [point, &index_vec].concat()
}

/// Return the number of variables that one need for an MLE to
/// batch the list of MLEs
#[inline]
pub fn get_batched_nv(num_var: usize, polynomials_len: usize) -> usize {
    num_var + log2(polynomials_len) as usize
}

// Input index
// - `i := (i_0, ...i_{n-1})`,
// - `num_vars := n`
// return three elements:
// - `x0 := (i_1, ..., i_{n-1}, 0)`
// - `x1 := (i_1, ..., i_{n-1}, 1)`
// - `sign := i_0`
#[inline]
pub fn get_index(i: usize, num_vars: usize) -> (usize, usize, bool) {
    let bit_sequence = bit_decompose(i as u64, num_vars);

    // the last bit comes first here because of LE encoding
    let x0 = project(&[[false].as_ref(), bit_sequence[..num_vars - 1].as_ref()].concat()) as usize;
    let x1 = project(&[[true].as_ref(), bit_sequence[..num_vars - 1].as_ref()].concat()) as usize;

    (x0, x1, bit_sequence[num_vars - 1])
}

/// Project a little endian binary vector into an integer.
#[inline]
pub(crate) fn project(input: &[bool]) -> u64 {
    let mut res = 0;
    for &e in input.iter().rev() {
        res <<= 1;
        res += e as u64;
    }
    res
}

#[cfg(test)]
mod test {
    use super::{bit_decompose, get_index, project};
    use ark_std::{rand::RngCore, test_rng};

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

    #[test]
    fn test_get_index() {
        let a = 0b1010;
        let (x0, x1, sign) = get_index(a, 4);
        assert_eq!(x0, 0b0100);
        assert_eq!(x1, 0b0101);
        assert!(sign);

        let (x0, x1, sign) = get_index(a, 5);
        assert_eq!(x0, 0b10100);
        assert_eq!(x1, 0b10101);
        assert!(!sign);

        let a = 0b1111;
        let (x0, x1, sign) = get_index(a, 4);
        assert_eq!(x0, 0b1110);
        assert_eq!(x1, 0b1111);
        assert!(sign);
    }
}
