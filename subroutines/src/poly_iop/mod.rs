// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

use ark_ff::PrimeField;
use std::marker::PhantomData;

mod errors;
mod perm_check;
pub mod prelude;
mod prod_check;
pub mod structs;
pub mod sum_check;
mod utils;
mod zero_check;

#[derive(Clone, Debug, Default, Copy, PartialEq, Eq)]
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
