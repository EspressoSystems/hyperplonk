use crate::{
    multilinear_kzg::{
        open_internal,
        srs::MultilinearProverParam,
        util::{get_quotient_domain, get_uni_domain},
    },
    prelude::PCSError,
    univariate_kzg::{srs::UnivariateProverParam, UnivariateKzgPCS},
    PolynomialCommitmentScheme,
};
use arithmetic::DenseMultilinearExtension;
use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
    MultilinearExtension, Radix2EvaluationDomain,
};
use ark_std::{end_timer, log2, start_timer};
use std::rc::Rc;
use transcript::IOPTranscript;

/// input: k points `p`, each of dim n
/// input: k MLEs `f`, each of dim n
///
/// steps:
/// 1. define `g(y, x0,...x_{n-1}) := \sum_{i=0}^{k-1} L_i(y) f_i` which is an
/// (n+1) mle (L_i(y) is implicit)
/// 2. define `h(y) := \sum_{i=0}^{k-1} L_i(y) p_i` which are n (h(y) is
/// implicit) univariate polynomials
/// 3. define `q(y) := g(y, h1, ...hn)`
/// 4. define `q'(y) := (q(y) - \sum_{i=0}^{k-1} L_i(y) y_i)/Z_k(x)` and obtain
/// a univariate polynomial this is done via interpolation
///   4.1 y in [w..w^k] -> n evaluation: y_1...y_k
///   4.2 evaluate y at alpha...alpha^(d*mu +1 )(k-1) over the quotient domain
///   4.3 interpolate to get `q'(y)` explicit form
/// 4. commit to `q'(y)` and append commitment to transcript
///
/// 5. sample a random field element `r` from transcript
/// 6. generate `g(r, x0,...x_{n-1}) := \sum_{i=0}^{k-1} L_i(r) f_i` which is an
/// `n` dim MLE
/// 7. compute `h(r) := \sum_{i=0}^{k-1} L_i(r) p_i` which is an `n` dim point
///
/// 8. open `q'(y)` at `r` and outputs its evaluation and proof
/// 9. open `g(y, x0,...x_{n-1})` at `r, h1(r),...hn(r)` and outputs its
/// evaluations and proof  
fn open_new<E: PairingEngine>(
    uni_prover_param: &UnivariateProverParam<E::G1Affine>,
    ml_prover_param: &MultilinearProverParam<E>,
    polynomials: &[Rc<DenseMultilinearExtension<E::Fr>>],
    points: &[Vec<E::Fr>],
) -> Result<(), PCSError> {
    // 4. define `q'(y) := (q(y) - \sum_{i=0}^{k-1} L_i(y) y_i)/Z_k(x)` and obtain a
    // univariate polynomial this is done via interpolation
    //   4.1 y in [w..w^k] -> n evaluation: y_1...y_k
    //   4.2 evaluate y at alpha...alpha^(d*mu +1 )(k-1) over the quotient domain
    //   4.3 interpolate to get `q'(y)` explicit form
    let q_prime_y = compute_q_prime_y(polynomials, points)?;

    let mut transcript = IOPTranscript::<E::Fr>::new(b"ml kzg");
    for point in points {
        // println!("points proving {:?}", points);
        transcript.append_serializable_element(b"w", point)?;
    }

    let q_prime_y_commit = UnivariateKzgPCS::<E>::commit(uni_prover_param, &q_prime_y)?;
    transcript.append_serializable_element(b"q(y)", &q_prime_y_commit)?;

    // 5. sample a random field element `r` from transcript
    let r = transcript.get_and_append_challenge(b"r")?;

    // 6. generate `g(r, x0,...x_{n-1}) := \sum_{i=0}^{k-1} L_i(r) f_i` which is an
    // `n` dim MLE
    let g_y = build_g_mle(polynomials, &r)?;

    // 7. compute `h(r) := \sum_{i=0}^{k-1} L_i(r) p_i` which is an `n` dim point
    let h_r = build_h_r(points, &r)?;

    // 8. open `q(y)` at `r` and outputs its evaluation and proof
    let (q_r_proof, q_r_eval) = UnivariateKzgPCS::<E>::open(uni_prover_param, &q_prime_y, &r)?;
    println!("q'(r): {}", q_r_eval);

    // 9. open `g(y, x0,...x_{n-1})` at `r, h1(r),...hn(r)` and outputs its
    // evaluations and proof
    let (g_r_proof, g_r_eval) = open_internal(ml_prover_param, &g_y, &h_r)?;

    println!("g(r): {}", g_r_eval);

    Ok(())
}

/// 7. compute `h(r) := \sum_{i=0}^{k-1} L_i(r) p_i` which is an `n` dim point
// TODO: merge this step with step 6?
fn build_h_r<F: PrimeField>(points: &[Vec<F>], r: &F) -> Result<Vec<F>, PCSError> {
    let k = points.len();
    let nv = points[0].len();
    let domain_k: Radix2EvaluationDomain<F> = get_uni_domain(k)?;
    let l_i_r_evals = domain_k.evaluate_all_lagrange_coefficients(*r);
    let mut h_r = vec![F::zero(); nv];

    for i in 0..nv {
        for j in 0..k {
            h_r[i] += points[j][i] * l_i_r_evals[j];
        }
    }

    Ok(h_r)
}

/// 6. generate `g(r, x0,...x_{n-1}) := \sum_{i=0}^{k-1} L_i(r) f_i` which is an
/// `n` dim MLE
fn build_g_mle<F: PrimeField>(
    polynomials: &[Rc<DenseMultilinearExtension<F>>],
    r: &F,
) -> Result<Rc<DenseMultilinearExtension<F>>, PCSError> {
    let k = polynomials.len();
    let nv = polynomials[0].num_vars;
    let domain_k: Radix2EvaluationDomain<F> = get_uni_domain(k)?;
    let l_i_r_evals = domain_k.evaluate_all_lagrange_coefficients(*r);
    let mut evaluations = vec![F::zero(); 1 << nv];

    // TODO: optimization
    for i in 0..1 << nv {
        for j in 0..k {
            evaluations[i] += polynomials[j].evaluations[i] * l_i_r_evals[j];
        }
    }

    let res = Rc::new(DenseMultilinearExtension::from_evaluations_vec(
        nv,
        evaluations,
    ));

    Ok(res)
}

/// compute the degree of q(y) which is (k-1)(d*n+1)
/// - k: number of points
/// - n: number of variables
#[inline]
fn compute_q_y_degree(k: usize, n: usize) -> usize {
    (k - 1) * (n + 1)
}

/// Evaluate `q(y) := g(y, h1, ...hn)` and obtain a univariate polynomial
/// this is done via interpolation
/// Inputs:
/// - points
/// - degree: (k - 1)(n + 1)
fn compute_q_prime_y<F: PrimeField>(
    polynomials: &[Rc<DenseMultilinearExtension<F>>],
    points: &[Vec<F>],
) -> Result<DensePolynomial<F>, PCSError> {
    let k = points.len();
    let n = points[0].len();
    let timer = start_timer!(|| format!("compute q(y) for {} MLEs each of {} nv", k, m));

    let q_y_degree = compute_q_y_degree(k, n);
    // large domain for evaluating q(y) in alpha roots
    let domain: GeneralEvaluationDomain<F> = get_quotient_domain(1 << log2(q_y_degree))?;
    let domain_k: Radix2EvaluationDomain<F> = get_uni_domain(k)?;

    let eval_points = eval_h_y_at_roots(points)?;
    let mut q_prime_y_evals = vec![];
    for root in domain.elements() {
        let mut cur_value = F::zero();
        // L_i(y) evaluated at the large domain's root
        let l_i_y_evals = domain_k.evaluate_all_lagrange_coefficients(root);
        for i in 0..k {
            // f_i(h1...hn)
            // println!("{}", i);
            let mle_eval = polynomials[i].evaluate(&eval_points[i]).unwrap();
            cur_value += l_i_y_evals[i] * mle_eval;
            cur_value -= l_i_y_evals[i] * root;
        }
        let z_k_eval = domain.evaluate_vanishing_polynomial(F::multiplicative_generator() * root);
        // println!("z_k {}", z_k_eval);
        q_prime_y_evals.push(cur_value / z_k_eval)
    }

    domain.coset_fft_in_place(&mut q_prime_y_evals);
    let res = DensePolynomial {
        coeffs: q_prime_y_evals,
    };

    // Evaluations::from_vec_and_domain(q_prime_y_eval, domain).interpolate();

    end_timer!(timer);
    Ok(res)
}

/// 2. define `h(y) := \sum_{i=0}^{k-1} L_i(y) p_i` which are n (h(y) is
/// implicit) univariate polynomials
///
/// we evaluate y at omega^1 .... omega^q_y_degree and obtain q_y_degree
/// number of points, each of dimension n
///  
fn eval_h_y_at_roots<F: PrimeField>(points: &[Vec<F>]) -> Result<Vec<Vec<F>>, PCSError> {
    let k = points.len();
    let n = points[0].len();

    let q_y_degree = compute_q_y_degree(k, n);

    // large domain for evaluating q(y) in alpha roots
    let domain: Radix2EvaluationDomain<F> = get_uni_domain(q_y_degree)?;

    // small domain for interpolating L_i(y)
    let domain_k: Radix2EvaluationDomain<F> = get_uni_domain(k)?;

    let mut res = vec![];

    for root in domain.elements() {
        let mut cur_h_y = vec![F::zero(); n];
        let l_i_list = domain_k.evaluate_all_lagrange_coefficients(root);
        // println!("l_i len {}", l_i_list.len());
        // todo: parallelization
        for i in 0..k {
            // println!("{} {}", i, k);
            cur_h_y
                .iter_mut()
                .zip(points[i].iter())
                .for_each(|(h_ij, p_ij)| *h_ij += *p_ij * l_i_list[i]);
        }
        res.push(cur_h_y)
    }
    Ok(res)
}

#[cfg(test)]
mod test {
    use crate::{
        multilinear_kzg::srs::MultilinearUniversalParams, prelude::PCSError,
        univariate_kzg::srs::UnivariateUniversalParams, StructuredReferenceString,
    };
    use arithmetic::DenseMultilinearExtension;
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ec::PairingEngine;
    use ark_ff::PrimeField;
    use ark_std::{test_rng, One, Zero};
    use std::rc::Rc;

    use super::{build_g_mle, build_h_r, compute_q_prime_y, eval_h_y_at_roots, open_new};

    #[test]
    fn test_compute_q_y() -> Result<(), PCSError> {
        test_compute_q_y_helper::<Fr>()
    }

    fn test_compute_q_y_helper<F: PrimeField>() -> Result<(), PCSError> {
        // point 1 is [1, 2]
        let point1 = vec![F::from(1u64), F::from(2u64)];

        // point 2 is [3, 4]
        let point2 = vec![F::from(3u64), F::from(4u64)];

        // point 3 is [5, 6]
        let point3 = vec![F::from(5u64), F::from(6u64)];

        // W1 = 3x1x2 + 2x2 whose evaluations are
        // 0, 0 |-> 0
        // 1, 0 |-> 0
        // 0, 1 |-> 2
        // 1, 1 |-> 5
        let w_eval = vec![F::zero(), F::zero(), F::from(2u64), F::from(5u64)];
        let w1 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        // w2 = x1x2 + x1 + x2 whose evaluations are
        // 0, 0 |-> 0
        // 1, 0 |-> 1
        // 0, 1 |-> 1
        // 1, 1 |-> 2
        let w_eval = vec![F::zero(), F::one(), F::one(), F::from(2u64)];
        let w2 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        let points = &[point1.clone(), point2.clone()];
        let polynomials = &[w1.clone(), w2.clone()];

        let point = eval_h_y_at_roots(points)?;
        let q_y = compute_q_prime_y(polynomials, points);

        let r = F::from(10u64);
        let g_mle = build_g_mle(polynomials, &r)?;
        let hr = build_h_r(points, &r)?;

        Ok(())
    }

    #[test]
    fn test_commitment() -> Result<(), PCSError> {
        test_commitment_helper::<Bls12_381>()
    }

    fn test_commitment_helper<E: PairingEngine>() -> Result<(), PCSError> {
        // point 1 is [1, 2]
        let point1 = vec![E::Fr::from(1u64), E::Fr::from(2u64)];

        // point 2 is [3, 4]
        let point2 = vec![E::Fr::from(3u64), E::Fr::from(4u64)];

        // point 3 is [5, 6]
        let point3 = vec![E::Fr::from(5u64), E::Fr::from(6u64)];

        // W1 = 3x1x2 + 2x2 whose evaluations are
        // 0, 0 |-> 0
        // 1, 0 |-> 0
        // 0, 1 |-> 2
        // 1, 1 |-> 5
        let w_eval = vec![
            E::Fr::zero(),
            E::Fr::zero(),
            E::Fr::from(2u64),
            E::Fr::from(5u64),
        ];
        let w1 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        // w2 = x1x2 + x1 + x2 whose evaluations are
        // 0, 0 |-> 0
        // 1, 0 |-> 1
        // 0, 1 |-> 1
        // 1, 1 |-> 2
        let w_eval = vec![E::Fr::zero(), E::Fr::one(), E::Fr::one(), E::Fr::from(2u64)];
        let w2 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        let points = &[point1.clone(), point2.clone()];
        let polynomials = &[w1.clone(), w2.clone()];

        let mut rng = test_rng();
        let uni_params =
            UnivariateUniversalParams::<E>::gen_srs_for_testing(&mut rng, 1usize << 15)?;
        let ml_params = MultilinearUniversalParams::<E>::gen_srs_for_testing(&mut rng, 15)?;

        let (uni_ck, uni_vk) = uni_params.trim(32)?;
        let (ml_ck, ml_vk) = ml_params.trim(3)?;

        open_new(&uni_ck, &ml_ck, polynomials, points)
    }
}
