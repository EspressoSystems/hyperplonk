use ark_ff::PrimeField;
use std::marker::PhantomData;

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
/// - ZeroCheck protocol.
pub struct PolyIOP<F: PrimeField> {
    phantom: PhantomData<F>,
}
