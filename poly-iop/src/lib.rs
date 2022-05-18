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

#[derive(Clone, Debug, Default, Copy)]
/// Struct for PolyIOP protocol.
/// It has an associated type `F` that defines the prime field the multi-variate
/// polynomial operates on.
///
/// An PolyIOP may be instantiated with one of the following:
/// - SumCheck protocol.
/// - ZeroCheck protocol.
/// - PermutationCheck protocol.
///
/// Those individual protocol may have similar or identical APIs.
/// The systematic way to invoke specific protocol is, for example
///     `<PolyIOP<F> as SumCheck<F>>::prove()`
pub struct PolyIOP<F: PrimeField> {
    /// Associated field
    #[doc(hidden)]
    phantom: PhantomData<F>,
}
