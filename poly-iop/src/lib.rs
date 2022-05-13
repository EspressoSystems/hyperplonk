#![allow(dead_code)]

use std::marker::PhantomData;

use ark_ff::PrimeField;

mod errors;
mod structs;
mod sum_check;
mod transcript;
mod utils;
mod virtual_poly;
mod zero_check;

pub use virtual_poly::VirtualPolynomial;

/// Struct for PolyIOP protocol.
/// It is instantiated with
/// - SumCheck protocol.
/// - ZeroCheck protocol. (WIP)
pub struct PolyIOP<F: PrimeField> {
    phantom: PhantomData<F>,
}
