use std::{env, fs::File, time::Instant};

use ark_bls12_381::{Bls12_381, Fr};
use ark_serialize::Write;
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
use rayon::ThreadPoolBuilder;

fn main() -> Result<(), HyperPlonkErrors> {
    let args: Vec<String> = env::args().collect();
    let thread = args[1].parse().unwrap_or(12);

    ThreadPoolBuilder::new()
        .num_threads(thread)
        .build_global()
        .unwrap();
    bench_vanilla_plonk(thread)?;

    Ok(())
}

fn bench_vanilla_plonk(thread: usize) -> Result<(), HyperPlonkErrors> {
    let mut rng = test_rng();
    let pcs_srs = MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, 22)?;

    let filename = format!("vanilla {}.txt", thread);
    let mut file = File::create(filename).unwrap();
    for nv in 1..20 {
        let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
        bench_mock_circuit_zkp_helper(&mut file, nv, &vanilla_gate, &pcs_srs)?;
    }

    let filename = format!("turbo {}.txt", thread);
    let mut file = File::create(filename).unwrap();
    for nv in 1..20 {
        let vanilla_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
        bench_mock_circuit_zkp_helper(&mut file, nv, &vanilla_gate, &pcs_srs)?;
    }

    Ok(())
}

fn bench_mock_circuit_zkp_helper(
    file: &mut File,
    nv: usize,
    gate: &CustomizedGates,
    pcs_srs: &(
        MultilinearUniversalParams<Bls12_381>,
        UnivariateUniversalParams<Bls12_381>,
    ),
) -> Result<(), HyperPlonkErrors> {
    let repetition = if nv < 10 {
        10
    } else if nv < 20 {
        5
    } else {
        2
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
    let t = start.elapsed().as_micros() / repetition as u128;
    println!("proving for {} variables: {} us", nv, t);
    file.write_all(format!("{} {}\n", nv, t).as_ref()).unwrap();

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
