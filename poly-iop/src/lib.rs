#![allow(dead_code)]

use std::marker::PhantomData;

use ark_ff::PrimeField;

mod errors;
mod structs;
mod sum_check;
mod transcript;
mod utils;
mod vertual_poly;
mod zero_check;

pub use vertual_poly::VirtualPolynomial;

/// Struct for PolyIOP.
pub struct PolyIOP<F: PrimeField> {
    phantom: PhantomData<F>,
}
