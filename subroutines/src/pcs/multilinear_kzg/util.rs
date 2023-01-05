// Copyright (c) 2022 Espresso Systems (espressosys.com)
// This file is part of the Jellyfish library.

// You should have received a copy of the MIT License
// along with the Jellyfish library. If not, see <https://mit-license.org/>.

//! Useful utilities for KZG PCS
use ark_ff::PrimeField;
use ark_std::{end_timer, start_timer};

use crate::PCSError;



/// Evaluate eq polynomial. use the public one later
pub(crate) fn eq_eval<F: PrimeField>(x: &[F], y: &[F]) -> Result<F, PCSError> {
    if x.len() != y.len() {
        return Err(PCSError::InvalidParameters(
            "x and y have different length".to_string(),
        ));
    }
    let start = start_timer!(|| "eq_eval");
    let mut res = F::one();
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let xi_yi = xi * yi;
        res *= xi_yi + xi_yi - xi - yi + F::one();
    }
    end_timer!(start);
    Ok(res)
}
