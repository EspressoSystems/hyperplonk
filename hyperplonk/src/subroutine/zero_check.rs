use std::rc::Rc;

use arithmetic::DenseMultilinearExtension;
use ark_ec::PairingEngine;
use ark_poly::MultilinearExtension;
use ark_std::{end_timer, start_timer};
use pcs::PolynomialCommitmentScheme;
use poly_iop::prelude::{IOPProof, PolyIOP, ZeroCheck};
use transcript::IOPTranscript;

use crate::{errors::HyperPlonkErrors, structs::HyperPlonkProvingKey, utils::build_f};

#[allow(clippy::type_complexity)]
pub(crate) fn zero_check_prover_subroutine<E, PCS>(
    pk: &HyperPlonkProvingKey<E, PCS>,
    witness_polys: &[PCS::Polynomial],
    transcript: &mut IOPTranscript<E::Fr>,
) -> Result<
    (
        IOPProof<E::Fr>,
        Vec<PCS::Proof>,
        Vec<PCS::Evaluation>,
        Vec<PCS::Proof>,
        Vec<PCS::Evaluation>,
    ),
    HyperPlonkErrors,
>
where
    E: PairingEngine,
    PCS: PolynomialCommitmentScheme<
        E,
        Polynomial = Rc<DenseMultilinearExtension<E::Fr>>,
        Point = Vec<E::Fr>,
        Evaluation = E::Fr,
    >,
{
    let step = start_timer!(|| "ZeroCheck on f");

    let fx = build_f(
        &pk.params.gate_func,
        pk.params.nv,
        &pk.selector_oracles,
        &witness_polys,
    )?;

    let zero_check_proof = <PolyIOP<E::Fr> as ZeroCheck<E::Fr>>::prove(&fx, transcript)?;

    let mut witness_zero_check_evals = vec![];
    let mut witness_zero_check_openings = vec![];
    // 4.2 open zero check proof
    // TODO: batch opening
    for wire_poly in witness_polys {
        // Open zero check proof
        let (zero_proof, zero_eval) =
            PCS::open(&pk.pcs_param, &wire_poly, &zero_check_proof.point)?;
        {
            let eval = wire_poly.evaluate(&zero_check_proof.point).ok_or_else(|| {
                HyperPlonkErrors::InvalidParameters(
                    "evaluation dimension does not match".to_string(),
                )
            })?;
            if eval != zero_eval {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }
        witness_zero_check_evals.push(zero_eval);
        witness_zero_check_openings.push(zero_proof);
    }

    // Open selector polynomial at zero_check_point
    let mut selector_oracle_openings = vec![];
    let mut selector_oracle_evals = vec![];

    // TODO: parallelization
    for selector_poly in pk.selector_oracles.iter() {
        // Open zero check proof
        // during verification, use this eval against subclaim
        let (zero_proof, zero_eval) =
            PCS::open(&pk.pcs_param, selector_poly, &zero_check_proof.point)?;

        #[cfg(feature = "extensive_sanity_checks")]
        {
            let eval = selector_poly
                .evaluate(&zero_check_proof.point)
                .ok_or_else(|| {
                    HyperPlonkErrors::InvalidParameters(
                        "evaluation dimension does not match".to_string(),
                    )
                })?;
            if eval != zero_eval {
                return Err(HyperPlonkErrors::InvalidProver(
                    "Evaluation is different from PCS opening".to_string(),
                ));
            }
        }
        selector_oracle_openings.push(zero_proof);
        selector_oracle_evals.push(zero_eval);
    }

    end_timer!(step);

    Ok((
        zero_check_proof,
        witness_zero_check_openings,
        witness_zero_check_evals,
        selector_oracle_openings,
        selector_oracle_evals,
    ))
}
