use std::time::Instant;

use ark_bls12_381::{Bls12_381, Fr};
use ark_std::test_rng;
use hyperplonk::{
    prelude::{CustomizedGates, HyperPlonkErrors, MockCircuit},
    HyperPlonkSNARK,
};
use pcs::{
    prelude::{MultilinearKzgPCS, MultilinearUniversalParams, UnivariateUniversalParams},
    PolynomialCommitmentScheme,
};
use poly_iop::PolyIOP;

fn main() -> Result<(), HyperPlonkErrors> {
    bench_vanilla_plonk()?;

    Ok(())
}

fn bench_vanilla_plonk() -> Result<(), HyperPlonkErrors> {
    let mut rng = test_rng();
    let pcs_srs = MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, 16)?;

    for nv in 1..10 {
        let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
        bench_mock_circuit_zkp_helper(nv, &vanilla_gate, &pcs_srs)?;
    }

    for nv in 1..10 {
        let vanilla_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
        bench_mock_circuit_zkp_helper(nv, &vanilla_gate, &pcs_srs)?;
    }

    Ok(())
}

fn bench_mock_circuit_zkp_helper(
    nv: usize,
    gate: &CustomizedGates,
    pcs_srs: &(
        MultilinearUniversalParams<Bls12_381>,
        UnivariateUniversalParams<Bls12_381>,
    ),
) -> Result<(), HyperPlonkErrors> {
    let repetition = if nv < 10 {
        100
    } else if nv < 20 {
        50
    } else {
        10
    };

    //==========================================================
    let start = Instant::now();
    for _ in 0..repetition {
        let circuit = MockCircuit::<Fr>::new(1 << nv, gate);
        assert!(circuit.is_satisfied());
    }
    println!(
        "mock circuit gen for {} variables: {} ns",
        nv,
        start.elapsed().as_nanos() / repetition as u128
    );

    let circuit = MockCircuit::<Fr>::new(1 << nv, gate);
    assert!(circuit.is_satisfied());
    let index = circuit.index;
    //==========================================================
    // generate pk and vks
    let start = Instant::now();
    for _ in 0..repetition {
        let (_pk, _vk) = <PolyIOP<Fr> as HyperPlonkSNARK<
            Bls12_381,
            MultilinearKzgPCS<Bls12_381>,
        >>::preprocess(&index, &pcs_srs)?;
    }
    println!(
        "key extraction for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    let (pk, vk) =
        <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::preprocess(
            &index, &pcs_srs,
        )?;
    //==========================================================
    // generate a proof
    let start = Instant::now();
    for _ in 0..repetition {
        let _proof =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::prove(
                &pk,
                &circuit.witnesses[0].coeff_ref(),
                &circuit.witnesses,
            )?;
    }
    println!(
        "proving for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    let proof = <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::prove(
        &pk,
        &circuit.witnesses[0].coeff_ref(),
        &circuit.witnesses,
    )?;
    //==========================================================
    // verify a proof
    let start = Instant::now();
    for _ in 0..repetition {
        let verify =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::verify(
                &vk,
                &circuit.witnesses[0].coeff_ref(),
                &proof,
            )?;
        assert!(verify);
    }
    println!(
        "verifying for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    Ok(())
}
