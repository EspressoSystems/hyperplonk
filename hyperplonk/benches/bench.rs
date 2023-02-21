// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

use std::{fs::File, io, time::Instant};

use ark_bls12_381::{Bls12_381, Fr};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Write};
use ark_std::test_rng;
use hyperplonk::{
    prelude::{CustomizedGates, HyperPlonkErrors, MockCircuit},
    HyperPlonkSNARK,
};
use subroutines::{
    pcs::{
        prelude::{MultilinearKzgPCS, MultilinearUniversalParams},
        PolynomialCommitmentScheme,
    },
    poly_iop::PolyIOP,
};

const SUPPORTED_SIZE: usize = 20;
const MIN_NUM_VARS: usize = 8;
const MAX_NUM_VARS: usize = 20;
const MIN_CUSTOM_DEGREE: usize = 1;
const MAX_CUSTOM_DEGREE: usize = 32;
const HIGH_DEGREE_TEST_NV: usize = 15;

fn main() -> Result<(), HyperPlonkErrors> {
    let thread = rayon::current_num_threads();
    println!("start benchmark with #{} threads", thread);
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
    bench_jellyfish_plonk(&pcs_srs, thread)?;
    println!();
    bench_vanilla_plonk(&pcs_srs, thread)?;
    println!();
    for degree in MIN_CUSTOM_DEGREE..=MAX_CUSTOM_DEGREE {
        bench_high_degree_plonk(&pcs_srs, degree, thread)?;
        println!();
    }
    println!();

    Ok(())
}

fn read_srs() -> Result<MultilinearUniversalParams<Bls12_381>, io::Error> {
    let mut f = File::open("srs.params")?;
    Ok(MultilinearUniversalParams::<Bls12_381>::deserialize_compressed_unchecked(&mut f).unwrap())
}

fn write_srs(pcs_srs: &MultilinearUniversalParams<Bls12_381>) {
    let mut f = File::create("srs.params").unwrap();
    pcs_srs.serialize_uncompressed(&mut f).unwrap();
}

fn bench_vanilla_plonk(
    pcs_srs: &MultilinearUniversalParams<Bls12_381>,
    thread: usize,
) -> Result<(), HyperPlonkErrors> {
    let filename = format!("vanilla threads {}.txt", thread);
    let mut file = File::create(filename).unwrap();
    for nv in MIN_NUM_VARS..=MAX_NUM_VARS {
        let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
        bench_mock_circuit_zkp_helper(&mut file, nv, &vanilla_gate, pcs_srs)?;
    }

    Ok(())
}

fn bench_jellyfish_plonk(
    pcs_srs: &MultilinearUniversalParams<Bls12_381>,
    thread: usize,
) -> Result<(), HyperPlonkErrors> {
    let filename = format!("jellyfish threads {}.txt", thread);
    let mut file = File::create(filename).unwrap();
    for nv in MIN_NUM_VARS..=MAX_NUM_VARS {
        let jf_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
        bench_mock_circuit_zkp_helper(&mut file, nv, &jf_gate, pcs_srs)?;
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
    println!("custom gate of degree {}", degree);
    let vanilla_gate = CustomizedGates::mock_gate(2, degree);
    bench_mock_circuit_zkp_helper(&mut file, HIGH_DEGREE_TEST_NV, &vanilla_gate, pcs_srs)?;

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
        >>::preprocess(&index, pcs_srs)?;
    }
    println!(
        "key extraction for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    let (pk, vk) =
        <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::preprocess(
            &index, pcs_srs,
        )?;
    //==========================================================
    // generate a proof
    let start = Instant::now();
    for _ in 0..repetition {
        let _proof =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::prove(
                &pk,
                &circuit.public_inputs,
                &circuit.witnesses,
            )?;
    }
    let t = start.elapsed().as_micros() / repetition as u128;
    println!(
        "proving for {} variables: {} us",
        nv,
        start.elapsed().as_micros() / repetition as u128
    );
    file.write_all(format!("{} {}\n", nv, t).as_ref()).unwrap();

    let proof = <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::prove(
        &pk,
        &circuit.public_inputs,
        &circuit.witnesses,
    )?;
    //==========================================================
    // verify a proof
    let start = Instant::now();
    for _ in 0..repetition {
        let verify =
            <PolyIOP<Fr> as HyperPlonkSNARK<Bls12_381, MultilinearKzgPCS<Bls12_381>>>::verify(
                &vk,
                &circuit.public_inputs,
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
