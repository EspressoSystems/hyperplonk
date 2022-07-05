use ark_ff::PrimeField;
use std::marker::PhantomData;

mod errors;
mod perm_check;
pub mod prelude;
mod prod_check;
mod structs;
mod sum_check;
mod utils;
mod virtual_poly;
mod zero_check;

pub use errors::PolyIOPErrors;
pub use hyperplonk::HyperPlonkPIOP;
pub use perm_check::{
    util::{identity_permutation_mle, random_permutation_mle},
    PermutationCheck,
};
pub use prod_check::ProductCheck;
pub use sum_check::SumCheck;
pub use utils::*;
pub use virtual_poly::{VPAuxInfo, VirtualPolynomial};
pub use zero_check::ZeroCheck;
use errors::PolyIOPErrors;
use virtual_poly::VirtualPolynomial;
use zero_check::ZeroCheck;

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
