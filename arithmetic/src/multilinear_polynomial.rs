use crate::{build_l, get_batched_nv, util::get_uni_domain, ArithErrors};
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, MultilinearExtension, Polynomial,
    Radix2EvaluationDomain,
};
use ark_std::{end_timer, log2, rand::RngCore, start_timer};
use std::rc::Rc;

pub use ark_poly::DenseMultilinearExtension;

/// Sample a random list of multilinear polynomials.
/// Returns
/// - the list of polynomials,
/// - its sum of polynomial evaluations over the boolean hypercube.
pub fn random_mle_list<F: PrimeField, R: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut R,
) -> (Vec<Rc<DenseMultilinearExtension<F>>>, F) {
    let start = start_timer!(|| "sample random mle list");
    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv))
    }
    let mut sum = F::zero();

    for _ in 0..(1 << nv) {
        let mut product = F::one();

        for e in multiplicands.iter_mut() {
            let val = F::rand(rng);
            e.push(val);
            product *= val;
        }
        sum += product;
    }

    let list = multiplicands
        .into_iter()
        .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    (list, sum)
}

// Build a randomize list of mle-s whose sum is zero.
pub fn random_zero_mle_list<F: PrimeField, R: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut R,
) -> Vec<Rc<DenseMultilinearExtension<F>>> {
    let start = start_timer!(|| "sample random zero mle list");

    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv))
    }
    for _ in 0..(1 << nv) {
        multiplicands[0].push(F::zero());
        for e in multiplicands.iter_mut().skip(1) {
            e.push(F::rand(rng));
        }
    }

    let list = multiplicands
        .into_iter()
        .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    list
}

/// Input a list of polynomials and a same list of points, built a merged
/// polynomial, and compute the evaluation of the points at each polynomial
pub fn batch_evaluate<F: PrimeField>(
    polynomials: &[Rc<DenseMultilinearExtension<F>>],
    points: &[Vec<F>],
) -> Result<(Vec<F>, Rc<DenseMultilinearExtension<F>>), ArithErrors> {
    if polynomials.len() != points.len() {
        return Err(ArithErrors::InvalidParameters(
            "polynomials and points have different sizes".to_string(),
        ));
    }
    let num_var = polynomials[0].num_vars;
    let domain = get_uni_domain::<F>(points.len())?;

    let poly_merged = merge_polynomials(polynomials)?;
    let l = build_l(num_var, points, &domain)?;
    let wl = compute_w_circ_l(&poly_merged, &l)?;
    let evals: Vec<F> = (0..points.len())
        .map(|i| wl.evaluate(&domain.element(i)))
        .collect();

    Ok((evals, poly_merged))
}

/// merge a set of polynomials. Returns an error if the
/// polynomials do not share a same number of nvs.
pub fn merge_polynomials<F: PrimeField>(
    polynomials: &[Rc<DenseMultilinearExtension<F>>],
) -> Result<Rc<DenseMultilinearExtension<F>>, ArithErrors> {
    let nv = polynomials[0].num_vars();
    for poly in polynomials.iter() {
        if nv != poly.num_vars() {
            return Err(ArithErrors::InvalidParameters(
                "num_vars do not match for polynomials".to_string(),
            ));
        }
    }

    let merged_nv = get_batched_nv(nv, polynomials.len());
    let mut scalars = vec![];
    for poly in polynomials.iter() {
        scalars.extend_from_slice(poly.to_evaluations().as_slice());
    }
    scalars.extend_from_slice(vec![F::zero(); (1 << merged_nv) - scalars.len()].as_ref());
    Ok(Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        merged_nv, scalars,
    )))
}

/// Compute W \circ l.
///
/// Given an MLE W, and a list of univariate polynomials l, generate the
/// univariate polynomial that composes W with l.
///
/// Returns an error if l's length does not matches number of variables in W.
pub(crate) fn compute_w_circ_l<F: PrimeField>(
    w: &DenseMultilinearExtension<F>,
    l: &[DensePolynomial<F>],
) -> Result<DensePolynomial<F>, ArithErrors> {
    let timer = start_timer!(|| "compute W \\circ l");

    if w.num_vars != l.len() {
        return Err(ArithErrors::InvalidParameters(format!(
            "l's length ({}) does not match num_variables ({})",
            l.len(),
            w.num_vars(),
        )));
    }

    let mut res_eval: Vec<F> = vec![];

    // TODO: consider to pass this in from caller
    // uni_degree is (product of each prefix's) + (2 * MLEs)
    // = (l.len() - (num_vars - log(l.len())) + 2) * l[0].degree
    let uni_degree = (l.len() - w.num_vars + log2(l.len()) as usize + 2) * l[0].degree();

    let domain = match Radix2EvaluationDomain::<F>::new(uni_degree) {
        Some(p) => p,
        None => {
            return Err(ArithErrors::InvalidParameters(
                "failed to build radix 2 domain".to_string(),
            ))
        },
    };
    for point in domain.elements() {
        // we reverse the order here because the coefficient vec are stored in
        // bit-reversed order
        let l_eval: Vec<F> = l.iter().rev().map(|x| x.evaluate(&point)).collect();
        res_eval.push(w.evaluate(l_eval.as_ref()).unwrap())
    }
    let evaluation = Evaluations::from_vec_and_domain(res_eval, domain);
    let res = evaluation.interpolate();

    end_timer!(timer);
    Ok(res)
}

#[cfg(test)]
mod tests {

    use super::*;
    use ark_bls12_381::Fr;
    use ark_std::test_rng;

    #[test]
    fn test_merge_poly() {
        test_merge_poly_helper::<Fr>()
    }

    fn test_merge_poly_helper<F: PrimeField>() {
        let mut rng = test_rng();
        let num_mle = 3;
        let nv = 5;

        let polynomials: Vec<Rc<DenseMultilinearExtension<F>>> = (0..num_mle)
            .map(|_| Rc::new(DenseMultilinearExtension::rand(nv, &mut rng)))
            .collect();
        let points: Vec<Vec<F>> = (0..num_mle)
            .map(|_| (0..nv).map(|_| F::rand(&mut rng)).collect())
            .collect();
        let eval = polynomials
            .iter()
            .zip(points.iter())
            .map(|(poly, point)| poly.evaluate(point).unwrap());

        let (eval_rec, _merged) = batch_evaluate(&polynomials, &points).unwrap();

        for (i, e) in eval.enumerate() {
            println!("{} {}", i, e)
        }
        for (i, e) in eval_rec.iter().enumerate() {
            println!("{} {}", i, e)
        }

        // for (a, &b) in eval.zip(eval_rec.iter()) {
        //     assert_eq!(a, b)
        // }
    }
}
