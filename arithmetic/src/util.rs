use crate::ArithErrors;
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, Radix2EvaluationDomain,
};
use ark_std::log2;

/// Decompose an integer into a binary vector in little endian.
pub(crate) fn bit_decompose(input: u64, num_var: usize) -> Vec<bool> {
    let mut res = Vec::with_capacity(num_var);
    let mut i = input;
    for _ in 0..num_var {
        res.push(i & 1 == 1);
        i >>= 1;
    }
    res
}

/// Given a list of points, build `l(points)` which is a list of univariate
/// polynomials that goes through the points
pub fn build_l<F: PrimeField>(
    num_var: usize,
    points: &[Vec<F>],
    domain: &Radix2EvaluationDomain<F>,
) -> Result<Vec<DensePolynomial<F>>, ArithErrors> {
    let prefix_len = log2(points.len()) as usize;
    let mut uni_polys = Vec::new();

    // 1.1 build the indexes and the univariate polys that go through the indexes
    let indexes: Vec<Vec<bool>> = (0..points.len())
        .map(|x| bit_decompose(x as u64, prefix_len))
        .collect();
    for i in 0..prefix_len {
        let eval: Vec<F> = indexes
            .iter()
            .map(|x| F::from(x[prefix_len - i - 1]))
            .collect();

        uni_polys.push(Evaluations::from_vec_and_domain(eval, *domain).interpolate());
    }

    // 1.2 build the actual univariate polys that go through the points
    for i in 0..num_var {
        let mut eval: Vec<F> = points.iter().map(|x| x[i]).collect();
        eval.extend_from_slice(vec![F::zero(); domain.size as usize - eval.len()].as_slice());
        uni_polys.push(Evaluations::from_vec_and_domain(eval, *domain).interpolate())
    }

    Ok(uni_polys)
}

/// Return the number of variables that one need for an MLE to
/// batch the list of MLEs
#[inline]
pub fn get_batched_nv(num_var: usize, polynomials_len: usize) -> usize {
    num_var + log2(polynomials_len) as usize
}

/// get the domain for the univariate polynomial
#[inline]
pub(crate) fn get_uni_domain<F: PrimeField>(
    uni_poly_degree: usize,
) -> Result<Radix2EvaluationDomain<F>, ArithErrors> {
    let domain = match Radix2EvaluationDomain::<F>::new(uni_poly_degree) {
        Some(p) => p,
        None => {
            return Err(ArithErrors::InvalidParameters(
                "failed to build radix 2 domain".to_string(),
            ))
        },
    };
    Ok(domain)
}
