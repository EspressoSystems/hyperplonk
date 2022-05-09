#![allow(dead_code)]

use std::marker::PhantomData;

use ark_ff::PrimeField;

mod errors;
mod poly_list;
mod structs;
mod sum_check;
mod transcript;
mod utils;

/// Struct for PolyIOP.
pub struct PolyIOP<F: PrimeField> {
    phantom: PhantomData<F>,
}
