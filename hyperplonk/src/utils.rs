#![allow(dead_code)]

use arithmetic::VirtualPolynomial;
use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use pcs::PolynomialCommitmentScheme;
use std::{borrow::Borrow, rc::Rc};

use crate::{
    custom_gate::CustomizedGates, errors::HyperPlonkErrors, structs::HyperPlonkParams,
    witness::WitnessColumn,
};

use poly_iop::prelude::bit_decompose;

/// An accumulator structure that holds a polynomial and
/// its opening points
#[derive(Debug)]
pub(super) struct PcsAccumulator<E: PairingEngine, PCS: PolynomialCommitmentScheme<E>> {
    pub(crate) polynomial: Option<PCS::Polynomial>,
    pub(crate) poly_commit: Option<PCS::Commitment>,
    pub(crate) points: Vec<PCS::Point>,
}

impl<E: PairingEngine, PCS: PolynomialCommitmentScheme<E>> PcsAccumulator<E, PCS> {
    /// Create an empty accumulator.
    pub(super) fn new() -> Self {
        Self {
            polynomial: None,
            poly_commit: None,
            points: vec![],
        }
    }

    /// Initialize the polynomial; requires both the polynomial
    /// and its commitment.
    pub(super) fn init_poly(
        &mut self,
        polynomial: PCS::Polynomial,
        commitment: PCS::Commitment,
    ) -> Result<(), HyperPlonkErrors> {
        if self.polynomial.is_some() || self.poly_commit.is_some() {
            return Err(HyperPlonkErrors::InvalidProver(
                "poly already set for accumulator".to_string(),
            ));
        }

        self.polynomial = Some(polynomial);
        self.poly_commit = Some(commitment);
        Ok(())
    }

    /// Push a new evaluation point into the accumulator
    pub(super) fn insert_point(&mut self, point: &PCS::Point) {
        self.points.push(point.clone())
    }

    /// Batch open all the points over a merged polynomial.
    /// A simple wrapper of PCS::multi_open
    pub(super) fn batch_open(
        &self,
        prover_param: impl Borrow<PCS::ProverParam>,
    ) -> Result<(PCS::BatchProof, Vec<PCS::Evaluation>), HyperPlonkErrors> {
        let poly = match &self.polynomial {
            Some(p) => p,
            None => {
                return Err(HyperPlonkErrors::InvalidProver(
                    "poly is set for accumulator".to_string(),
                ))
            },
        };

        let commitment = match &self.poly_commit {
            Some(p) => p,
            None => {
                return Err(HyperPlonkErrors::InvalidProver(
                    "poly is set for accumulator".to_string(),
                ))
            },
        };
        Ok(PCS::multi_open_single_poly(
            prover_param.borrow(),
            commitment,
            poly,
            &self.points,
        )?)
    }
}

/// Build MLE from matrix of witnesses.
///
/// Given a matrix := [row1, row2, ...] where
/// row1:= (a1, a2, ...)
/// row2:= (b1, b2, ...)
/// row3:= (c1, c2, ...)
///
/// output mle(a1,b1,c1, ...), mle(a2,b2,c2, ...), ...
#[macro_export]
macro_rules! build_mle {
    ($rows:expr) => {{
        let mut res = Vec::with_capacity($rows.len());
        let num_vars = log2($rows.len()) as usize;
        let num_mles = $rows[0].0.len();

        for i in 0..num_mles {
            let mut cur_coeffs = Vec::new();
            for row in $rows.iter() {
                cur_coeffs.push(row.0[i])
            }
            res.push(Rc::new(DenseMultilinearExtension::from_evaluations_vec(
                num_vars, cur_coeffs,
            )))
        }

        Ok(res)
    }};
}

/// Sanity-check for HyperPlonk SNARK proving
pub(crate) fn prover_sanity_check<F: PrimeField>(
    params: &HyperPlonkParams,
    pub_input: &[F],
    witnesses: &[WitnessColumn<F>],
) -> Result<(), HyperPlonkErrors> {
    // public input length must be no greater than num_constraints

    if pub_input.len() > params.num_constraints {
        return Err(HyperPlonkErrors::InvalidProver(format!(
            "Public input length {} is greater than num constraits {}",
            pub_input.len(),
            params.num_pub_input
        )));
    }

    // public input length
    if pub_input.len() != params.num_pub_input {
        return Err(HyperPlonkErrors::InvalidProver(format!(
            "Public input length is not correct: got {}, expect {}",
            pub_input.len(),
            params.num_pub_input
        )));
    }
    // witnesses length
    for (i, w) in witnesses.iter().enumerate() {
        if w.0.len() != params.num_constraints {
            return Err(HyperPlonkErrors::InvalidProver(format!(
                "{}-th witness length is not correct: got {}, expect {}",
                i,
                w.0.len(),
                params.num_constraints
            )));
        }
    }
    // check public input matches witness[0]'s first 2^ell elements
    for (i, (&pi, &w)) in pub_input
        .iter()
        .zip(witnesses[0].0.iter().take(pub_input.len()))
        .enumerate()
    {
        if pi != w {
            return Err(HyperPlonkErrors::InvalidProver(format!(
                "The {:?}-th public input {:?} does not match witness[0] {:?}",
                i, pi, w
            )));
        }
    }

    Ok(())
}

/// build `f(w_0(x),...w_d(x))` where `f` is the constraint polynomial
/// i.e., `f(a, b, c) = q_l a(x) + q_r b(x) + q_m a(x)b(x) - q_o c(x)` in
/// vanilla plonk
pub(crate) fn build_f<F: PrimeField>(
    gates: &CustomizedGates,
    num_vars: usize,
    selector_mles: &[Rc<DenseMultilinearExtension<F>>],
    witness_mles: &[Rc<DenseMultilinearExtension<F>>],
) -> Result<VirtualPolynomial<F>, HyperPlonkErrors> {
    // TODO: check that selector and witness lengths match what is in
    // the gate definition

    for selector_mle in selector_mles.iter() {
        if selector_mle.num_vars != num_vars {
            return Err(HyperPlonkErrors::InvalidParameters(format!(
                "selector has different number of vars: {} vs {}",
                selector_mle.num_vars, num_vars
            )));
        }
    }

    for witness_mle in witness_mles.iter() {
        if witness_mle.num_vars != num_vars {
            return Err(HyperPlonkErrors::InvalidParameters(format!(
                "selector has different number of vars: {} vs {}",
                witness_mle.num_vars, num_vars
            )));
        }
    }

    let mut res = VirtualPolynomial::<F>::new(num_vars);

    for (coeff, selector, witnesses) in gates.gates.iter() {
        let coeff_fr = if *coeff < 0 {
            -F::from(-*coeff as u64)
        } else {
            F::from(*coeff as u64)
        };
        let mut mle_list = vec![];
        if let Some(s) = *selector {
            mle_list.push(selector_mles[s].clone())
        }
        for &witness in witnesses.iter() {
            mle_list.push(witness_mles[witness].clone())
        }
        res.add_mle_list(mle_list, coeff_fr)?;
    }

    Ok(res)
}

pub(crate) fn eval_f<F: PrimeField>(
    gates: &CustomizedGates,
    selector_evals: &[F],
    witness_evals: &[F],
) -> Result<F, HyperPlonkErrors> {
    let mut res = F::zero();
    for (coeff, selector, witnesses) in gates.gates.iter() {
        let mut cur_value = if *coeff < 0 {
            -F::from(-*coeff as u64)
        } else {
            F::from(*coeff as u64)
        };
        cur_value *= match selector {
            Some(s) => selector_evals[*s],
            None => F::one(),
        };
        for &witness in witnesses.iter() {
            cur_value *= witness_evals[witness]
        }
        res += cur_value;
    }
    Ok(res)
}

/// given the evaluation input `point` of the `index`-th polynomial,
/// obtain the evaluation point in the merged polynomial
pub(crate) fn gen_eval_point<F: PrimeField>(index: usize, index_len: usize, point: &[F]) -> Vec<F> {
    let index_vec: Vec<F> = bit_decompose(index as u64, index_len)
        .into_iter()
        .map(|x| F::from(x))
        .collect();
    [point, &index_vec].concat()
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_ff::PrimeField;
    use ark_poly::MultilinearExtension;
    #[test]
    fn test_build_gate() -> Result<(), HyperPlonkErrors> {
        test_build_gate_helper::<Fr>()
    }

    fn test_build_gate_helper<F: PrimeField>() -> Result<(), HyperPlonkErrors> {
        let num_vars = 2;

        // ql = 3x1x2 + 2x2 whose evaluations are
        // 0, 0 |-> 0
        // 0, 1 |-> 2
        // 1, 0 |-> 0
        // 1, 1 |-> 5
        let ql_eval = vec![F::zero(), F::from(2u64), F::zero(), F::from(5u64)];
        let ql = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, ql_eval));

        // W1 = x1x2 + x1 whose evaluations are
        // 0, 0 |-> 0
        // 0, 1 |-> 0
        // 1, 0 |-> 1
        // 1, 1 |-> 2
        let w_eval = vec![F::zero(), F::zero(), F::from(1u64), F::from(2u64)];
        let w1 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        // W2 = x1 + x2 whose evaluations are
        // 0, 0 |-> 0
        // 0, 1 |-> 1
        // 1, 0 |-> 1
        // 1, 1 |-> 2
        let w_eval = vec![F::zero(), F::one(), F::from(1u64), F::from(2u64)];
        let w2 = Rc::new(DenseMultilinearExtension::from_evaluations_vec(2, w_eval));

        // Example:
        //     q_L(X) * W_1(X)^5 - W_2(X)
        // is represented as
        // vec![
        //     ( 1,    Some(id_qL),    vec![id_W1, id_W1, id_W1, id_W1, id_W1]),
        //     (-1,    None,           vec![id_W2])
        // ]
        let gates = CustomizedGates {
            gates: vec![(1, Some(0), vec![0, 0, 0, 0, 0]), (-1, None, vec![1])],
        };
        let f = build_f(&gates, num_vars, &[ql.clone()], &[w1.clone(), w2.clone()])?;

        // Sanity check on build_f
        // f(0, 0) = 0
        assert_eq!(f.evaluate(&[F::zero(), F::zero()])?, F::zero());
        // f(0, 1) = 2 * 0^5 + (-1) * 1 = -1
        assert_eq!(f.evaluate(&[F::zero(), F::one()])?, -F::one());
        // f(1, 0) = 0 * 1^5 + (-1) * 1 = -1
        assert_eq!(f.evaluate(&[F::one(), F::zero()])?, -F::one());
        // f(1, 1) = 5 * 2^5 + (-1) * 2 = 158
        assert_eq!(f.evaluate(&[F::one(), F::one()])?, F::from(158u64));

        // test eval_f
        {
            let point = [F::zero(), F::zero()];
            let selector_evals = ql.evaluate(&point).unwrap();
            let witness_evals = [w1.evaluate(&point).unwrap(), w2.evaluate(&point).unwrap()];
            let eval_f = eval_f(&gates, &[selector_evals], &witness_evals)?;
            // f(0, 0) = 0
            assert_eq!(eval_f, F::zero());
        }
        {
            let point = [F::zero(), F::one()];
            let selector_evals = ql.evaluate(&point).unwrap();
            let witness_evals = [w1.evaluate(&point).unwrap(), w2.evaluate(&point).unwrap()];
            let eval_f = eval_f(&gates, &[selector_evals], &witness_evals)?;
            // f(0, 1) = 2 * 0^5 + (-1) * 1 = -1
            assert_eq!(eval_f, -F::one());
        }
        {
            let point = [F::one(), F::zero()];
            let selector_evals = ql.evaluate(&point).unwrap();
            let witness_evals = [w1.evaluate(&point).unwrap(), w2.evaluate(&point).unwrap()];
            let eval_f = eval_f(&gates, &[selector_evals], &witness_evals)?;
            // f(1, 0) = 0 * 1^5 + (-1) * 1 = -1
            assert_eq!(eval_f, -F::one());
        }
        {
            let point = [F::one(), F::one()];
            let selector_evals = ql.evaluate(&point).unwrap();
            let witness_evals = [w1.evaluate(&point).unwrap(), w2.evaluate(&point).unwrap()];
            let eval_f = eval_f(&gates, &[selector_evals], &witness_evals)?;
            // f(1, 1) = 5 * 2^5 + (-1) * 2 = 158
            assert_eq!(eval_f, F::from(158u64));
        }
        Ok(())
    }
}
