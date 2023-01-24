// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

mod errors;
mod multilinear_polynomial;
mod univariate_polynomial;
mod util;
mod virtual_polynomial;

pub use errors::ArithErrors;
pub use multilinear_polynomial::{
    evaluate_no_par, evaluate_opt, fix_last_variables, fix_last_variables_no_par, fix_variables,
    identity_permutation, identity_permutation_mles, merge_polynomials, random_mle_list,
    random_permutation, random_permutation_mles, random_zero_mle_list, DenseMultilinearExtension,
};
pub use univariate_polynomial::{build_l, get_uni_domain};
pub use util::{bit_decompose, gen_eval_point, get_batched_nv, get_index};
pub use virtual_polynomial::{
    build_eq_x_r, build_eq_x_r_vec, eq_eval, VPAuxInfo, VirtualPolynomial,
};
