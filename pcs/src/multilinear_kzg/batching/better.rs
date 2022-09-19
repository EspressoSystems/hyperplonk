//! Aggregate and batch opening PCS, tailored to Hyperplonk.
//!
//! TODO: refactor this file and move some into hyperplonk repo

use crate::{
    multilinear_kzg::{open_internal, util::gen_eval_point},
    prelude::{
        compute_w_circ_l, fix_variables_reverse_ord, merge_polynomials, MultilinearKzgBatchProof,
        MultilinearProverParam, PCSError, UnivariateKzgPCS, UnivariateProverParam,
    },
    PolynomialCommitmentScheme,
};
use arithmetic::{build_l, get_uni_domain, DenseMultilinearExtension};
use ark_ec::PairingEngine;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, MultilinearExtension, Polynomial};
use ark_std::{end_timer, log2, start_timer, Zero};
use std::rc::Rc;
use transcript::IOPTranscript;

#[derive(Debug, Clone)]
pub struct PolyAndPoints<E: PairingEngine> {
    polynomial: Rc<DenseMultilinearExtension<E::Fr>>,
    points: Vec<Vec<E::Fr>>,
}

impl<E: PairingEngine> PolyAndPoints<E> {
    pub fn init_poly(polynomial: Rc<DenseMultilinearExtension<E::Fr>>) -> Self {
        Self {
            polynomial,
            points: vec![],
        }
    }
    pub fn init_poly_and_point(
        polynomial: Rc<DenseMultilinearExtension<E::Fr>>,
        point: &[E::Fr],
    ) -> Self {
        Self {
            polynomial,
            points: vec![point.to_vec()],
        }
    }

    pub fn insert_point(&mut self, point: &[E::Fr]) {
        self.points.push(point.to_vec())
    }
}

#[derive(Debug)]
pub struct ProverPolysAndPoints<E: PairingEngine> {
    /// the points are already padded for w_merge_nv size
    pub witness_merge_pap: PolyAndPoints<E>,
    /// the points do not require padding
    pub prod_pap: PolyAndPoints<E>,
    /// the points are not padded
    pub selectors_pap: Vec<PolyAndPoints<E>>,
}

impl<E: PairingEngine> ProverPolysAndPoints<E> {
    #[inline]
    fn get_num_points(&self) -> usize {
        self.witness_merge_pap.points.len() + self.prod_pap.points.len() + self.selectors_pap.len()
    }

    /// get the length to pad selector poly to prod(x)
    #[inline]
    fn get_selector_padding_len(&self) -> usize {
        self.prod_pap.polynomial.num_vars - self.selectors_pap[0].polynomial.num_vars
    }

    /// get selectors chunk size
    #[inline]
    fn get_selector_chunk_size(&self) -> usize {
        // each chunk stores up to 2^ {prod(x)'s dim -  q(s)'s dim} polynomials
        1 << self.get_selector_padding_len()
    }

    /// get number of selector chunks
    #[inline]
    fn get_selector_num_chunks(&self) -> usize {
        let selectors_len = self.selectors_pap.len();
        let chunk_size = self.get_selector_chunk_size();
        println!("chunk_size {}", chunk_size);
        (selectors_len + chunk_size - 1) / chunk_size
    }

    /// get the total length
    #[inline]
    fn get_total_len(&self) -> usize {
        let selector_chunks = self.get_selector_num_chunks();
        println!("selector_chunks {}", selector_chunks);
        // 1 for witness
        // 1 for prod(x)
        // selectors_pap.len() for selector
        self.prod_pap.polynomial.num_vars + log2(2 + selector_chunks) as usize
    }

    /// get the layer 1 padded length
    #[inline]
    fn get_lay_1_padded_len(&self) -> usize {
        let selector_chunks = self.get_selector_num_chunks();
        // 1 for witness
        // 1 for prod(x)
        // selectors_pap.len() for selector
        log2(2 + selector_chunks) as usize
    }

    /// merge the polynomials; following the same sequence as points
    fn merge_poly(&self) -> Result<DenseMultilinearExtension<E::Fr>, PCSError> {
        let timer = start_timer!(|| "merge polynomials");
        let nv = self.prod_pap.polynomial.num_vars;
        let mut to_be_merged = vec![];
        let w_merged_ext_one = Rc::new(DenseMultilinearExtension::from_evaluations_slice(
            nv,
            [
                self.witness_merge_pap.polynomial.evaluations.as_slice(),
                vec![E::Fr::zero(); 1 << self.witness_merge_pap.polynomial.num_vars].as_slice(),
            ]
            .concat()
            .as_ref(),
        ));
        to_be_merged.push(w_merged_ext_one);
        to_be_merged.push(self.prod_pap.polynomial.clone());

        for block_of_paps in self.selectors_pap.chunks(self.get_selector_chunk_size()) {
            let mut eval = vec![];
            for e in block_of_paps.iter() {
                eval.extend_from_slice(e.polynomial.evaluations.as_ref())
            }
            eval.resize(1 << nv, E::Fr::zero());

            to_be_merged.push(Rc::new(DenseMultilinearExtension::from_evaluations_vec(
                nv, eval,
            )))
        }

        for (i, e) in to_be_merged.iter().enumerate() {
            println!("{} {}", i, e.num_vars)
        }

        let res = merge_polynomials(&to_be_merged)?;
        end_timer!(timer);
        Ok(res)
    }

    // fn compute_r(&self)->Result<Vec<E::Fr>, HyperPlonkErrors>{

    // }

    fn build_l(&self) -> Result<Vec<DensePolynomial<E::Fr>>, PCSError> {
        let timer = start_timer!(|| "build l for poly and points");

        let lay1_to_pad = self.get_lay_1_padded_len();
        // 1. pad w_merged's points
        // Since w_merged is one dim smaller than prod(x), we will need to pad a zero
        // first
        let mut points = self
            .witness_merge_pap
            .points
            .iter()
            .map(|p| gen_eval_point(0, lay1_to_pad, &[p.as_slice(), &[E::Fr::zero()]].concat()))
            .collect::<Vec<_>>();

        // 2. add prod(x) points with
        points.extend_from_slice(
            &self
                .prod_pap
                .points
                .iter()
                .map(|p| gen_eval_point(1, lay1_to_pad, p))
                .collect::<Vec<_>>(),
        );

        // 3. pad selectors_points
        let mut current_index = 2;
        let selector_padding_len = self.get_selector_padding_len();
        for block_of_paps in self.selectors_pap.chunks(self.get_selector_chunk_size()) {
            for (index, selector) in block_of_paps.iter().enumerate() {
                // we first pad selectors to the size of prod(x)
                let point = gen_eval_point(index, selector_padding_len, &selector.points[0]);
                // then add the lay1 padding
                points.push(gen_eval_point(current_index, lay1_to_pad, &point));
            }
            current_index += 1;
        }

        println!("witness len: {}", self.witness_merge_pap.points.len());
        println!("+++++");
        println!("prod len: {}", self.prod_pap.points.len());
        println!("+++++");
        for (i, p) in points.iter().enumerate() {
            println!("{}", i);
            for s in p.iter() {
                println!("{}", s)
            }
            println!();
        }
        // 4. build l from all points
        let domain = get_uni_domain::<E::Fr>(points.len())?;
        let res = build_l(&points, &domain)?;

        end_timer!(timer);
        Ok(res)
    }

    pub fn batch_open(
        &self,
        uni_prover_param: &UnivariateProverParam<E::G1Affine>,
        ml_prover_param: &MultilinearProverParam<E>,
    ) -> Result<(MultilinearKzgBatchProof<E>, Vec<E::Fr>), PCSError> {
        let domain = get_uni_domain::<E::Fr>(self.get_num_points())?;

        // let points = self.

        // build l: a univariate that go through all points
        let uni_polys = self.build_l()?;
        println!("here");
        // build merged polynomial
        let merged_poly = self.merge_poly()?;
        println!("here");

        // build `q(x)` which is a univariate polynomial `W circ l`
        let mut transcript = IOPTranscript::new(b"ml kzg");
        let q_x = compute_w_circ_l(&merged_poly, &uni_polys, self.get_num_points())?;
        let q_x_commit = UnivariateKzgPCS::<E>::commit(uni_prover_param, &q_x)?;
        transcript.append_serializable_element(b"w", &q_x_commit)?;

        for point in self.witness_merge_pap.points.iter() {
            transcript.append_serializable_element(b"points", point)?;
        }
        for point in self.prod_pap.points.iter() {
            transcript.append_serializable_element(b"points", point)?;
        }
        for pap in self.selectors_pap.iter() {
            transcript.append_serializable_element(b"points", &pap.points[0])?;
        }
        let r = transcript.get_and_append_challenge(b"r")?;

        // 5. build q(omega^i) and their openings
        let mut q_x_opens = vec![];
        let mut q_x_evals = vec![];
        for i in 0..self.get_num_points() {
            let (q_x_open, q_x_eval) =
                UnivariateKzgPCS::<E>::open(uni_prover_param, &q_x, &domain.element(i))?;
            q_x_opens.push(q_x_open);
            q_x_evals.push(q_x_eval);
            #[cfg(feature = "extensive_sanity_checks")]
            {
                // sanity check
                let point: Vec<E::Fr> = uni_polys
                    .iter()
                    .map(|poly| poly.evaluate(&domain.element(i)))
                    .collect();
                let mle_eval = merged_poly.evaluate(&point).unwrap();
                if mle_eval != q_x_eval {
                    return Err(PCSError::InvalidProver(
                        "Q(omega) does not match W(l(omega))".to_string(),
                    ));
                }
            }
        }

        // 6. build q(r) and its opening
        let (q_x_open, q_r_value) = UnivariateKzgPCS::<E>::open(uni_prover_param, &q_x, &r)?;
        q_x_opens.push(q_x_open);
        q_x_evals.push(q_r_value);

        // 7. get a point `p := l(r)`
        let point: Vec<E::Fr> = uni_polys.iter().map(|poly| poly.evaluate(&r)).collect();

        // 9. output value that is `w` evaluated at `p` (which should match `q(r)`)
        #[cfg(feature = "extensive_sanity_checks")]
        if merged_poly.evaluate(&point).unwrap() != q_r_value {
            return Err(PCSError::InvalidProver(
                "Q(r) does not match W(l(r))".to_string(),
            ));
        }

        // 8. output an opening of `w` over point `p`
        let merged_poly_partial_eval = fix_variables_reverse_ord(
            Rc::new(merged_poly.clone()),
            point[self.prod_pap.points.len()..].as_ref(),
        );

        let (mle_opening, mle_eval) = open_internal(
            &ml_prover_param,
            &merged_poly_partial_eval,
            &point[..self.prod_pap.points.len()],
        )?;

        #[cfg(feature = "extensive_sanity_checks")]
        {
            let a = self
                .witness_merge_pap
                .polynomial
                .evaluate(&self.witness_merge_pap.points[0])
                .unwrap();
            println!("a {}", a);
            let b = merged_poly
                .evaluate(
                    [
                        self.witness_merge_pap.points[0].as_slice(),
                        vec![
                            E::Fr::zero();
                            self.get_total_len() - self.witness_merge_pap.points[0].len()
                        ]
                        .as_ref(),
                    ]
                    .concat()
                    .as_ref(),
                )
                .unwrap();
            println!("b {}", b);
        }

        println!("here asdfasd");
        todo!()
    }
}
