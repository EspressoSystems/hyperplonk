use std::{env, fs::File, io, time::Instant};

use ark_bls12_381::{Bls12_381, Fr};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Write};
use ark_std::test_rng;
use hyperplonk::{
    prelude::{CustomizedGates, HyperPlonkErrors, MockCircuit},
    HyperPlonkSNARK,
};
use rayon::ThreadPoolBuilder;
use subroutines::{
    pcs::{
        prelude::{MultilinearKzgPCS, MultilinearUniversalParams},
        PolynomialCommitmentScheme,
    },
    poly_iop::PolyIOP,
};

// For 32 G of ram we can do nv = 24
// For 128 G of ram we can do nv = 26
const SUPPORTED_SIZE: usize = 26;
const MIN_NUM_VARS: usize = 8;
const MAX_NUM_VARS: usize = SUPPORTED_SIZE;
const MIN_CUSTOM_DEGREE: usize = 1;
const MAX_CUSTOM_DEGREE: usize = 32;

fn main() -> Result<(), HyperPlonkErrors> {
    let args: Vec<String> = env::args().collect();
    let thread = args[1].parse().unwrap_or(0);
    let mut rng = test_rng();

    let pcs_srs = {
        match read_srs() {
            Ok(p) => p,
            Err(_e) => {
                let srs =
                    MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;
                write_srs(&srs);
                srs
            },
        }
    };

    ThreadPoolBuilder::new()
        .num_threads(thread)
        .build()
        .unwrap();

    bench_jellyfish_plonk(&pcs_srs, thread)?;
    bench_vanilla_plonk(&pcs_srs, thread)?;
    for degree in MIN_CUSTOM_DEGREE..MAX_CUSTOM_DEGREE {
        bench_high_degree_plonk(&pcs_srs, degree, thread)?;
    }
    Ok(())
}

fn read_srs() -> Result<MultilinearUniversalParams<Bls12_381>, io::Error> {
    let mut f = File::open("srs.params")?;
    Ok(MultilinearUniversalParams::<Bls12_381>::deserialize_unchecked(&mut f).unwrap())
}

fn write_srs(pcs_srs: &MultilinearUniversalParams<Bls12_381>) {
    let mut f = File::create("srs.params").unwrap();
    pcs_srs.serialize_uncompressed(&mut f).unwrap();
}

fn bench_jellyfish_plonk(
    pcs_srs: &MultilinearUniversalParams<Bls12_381>,
    thread: usize,
) -> Result<(), HyperPlonkErrors> {
    let filename = format!("jellyfish threads {}.txt", thread);
    let mut file = File::create(filename).unwrap();
    for nv in MIN_NUM_VARS..=MAX_NUM_VARS {
        let vanilla_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
        bench_mock_circuit_zkp_helper(&mut file, nv, &vanilla_gate, &pcs_srs)?;
    }

    Ok(())
}

fn bench_vanilla_plonk(
    pcs_srs: &MultilinearUniversalParams<Bls12_381>,
    thread: usize,
) -> Result<(), HyperPlonkErrors> {
    let filename = format!("vanilla threads {}.txt", thread);
    let mut file = File::create(filename).unwrap();
    for nv in MIN_NUM_VARS..=MAX_NUM_VARS {
        let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
        bench_mock_circuit_zkp_helper(&mut file, nv, &vanilla_gate, &pcs_srs)?;
    }

    Ok(())
}

fn bench_high_degree_plonk(
    pcs_srs: &MultilinearUniversalParams<Bls12_381>,
    degree: usize,
    thread: usize,
) -> Result<(), HyperPlonkErrors> {
    let filename = format!("high degree {} thread {}.txt", degree, thread);
    let mut file = File::create(filename).unwrap();
    for nv in MIN_NUM_VARS..=MAX_NUM_VARS {
        let vanilla_gate = CustomizedGates::mock_gate(2, degree);
        bench_mock_circuit_zkp_helper(&mut file, nv, &vanilla_gate, &pcs_srs)?;
    }

    Ok(())
}

fn bench_mock_circuit_zkp_helper(
    file: &mut File,
    nv: usize,
    gate: &CustomizedGates,
    pcs_srs: &MultilinearUniversalParams<Bls12_381>,
) -> Result<(), HyperPlonkErrors> {
    let repetition = if nv < 10 {
        5
    } else if nv < 20 {
        2
    } else {
        1
    };

    //==========================================================
    // let start = Instant::now();
    // for _ in 0..repetition {
    //     let circuit = MockCircuit::<Fr>::new(1 << nv, gate);
    //     assert!(circuit.is_satisfied());
    // }
    // println!(
    //     "mock circuit gen for {} variables: {} us",
    //     nv,
    //     start.elapsed().as_micros() / repetition as u128
    // );

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
    let mut proofs = vec![];
    for _ in 0..repetition {
        let proof =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::prove(
                &pk,
                &circuit.witnesses[0].coeff_ref(),
                &circuit.witnesses,
            )?;
        proofs.push(proof)
    }
    let t = start.elapsed().as_micros() / repetition as u128;
    println!(
        "proving for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    file.write_all(format!("{} {}\n", nv, t).as_ref()).unwrap();

    //==========================================================
    // verify a proof
    let start = Instant::now();
    for _ in 0..repetition {
        let verify =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::verify(
                &vk,
                &circuit.witnesses[0].coeff_ref(),
                &proofs[0],
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
