use ark_ff::PrimeField;
use std::marker::PhantomData;

mod errors;
mod structs;
mod sum_check;
mod transcript;
mod utils;
mod virtual_poly;
mod zero_check;

pub use errors::PolyIOPErrors;
pub use sum_check::SumCheck;
pub use virtual_poly::VirtualPolynomial;
pub use zero_check::{build_eq_x_r, ZeroCheck};

/// Struct for PolyIOP protocol.
/// It is instantiated with
/// - SumCheck protocol.
/// - ZeroCheck protocol.
pub struct PolyIOP<F: PrimeField> {
    phantom: PhantomData<F>,
}
