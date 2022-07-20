use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::UniformRand;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::test_rng;
use pcs::{
    prelude::{KZGMultilinearPCS, PCSErrors, PolynomialCommitmentScheme},
    StructuredReferenceString,
};
use std::time::Instant;

fn main() -> Result<(), PCSErrors> {
    bench_pcs()
}

fn bench_pcs() -> Result<(), PCSErrors> {
    let mut rng = test_rng();

    // normal polynomials
    let uni_params = KZGMultilinearPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, 18)?;

    for nv in 4..19 {
        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let poly = DenseMultilinearExtension::rand(nv, &mut rng);
        let (ml_ck, ml_vk) = uni_params.0.trim(nv)?;
        let (uni_ck, uni_vk) = uni_params.1.trim(nv)?;
        let ck = (ml_ck, uni_ck);
        let vk = (ml_vk, uni_vk);

        let point: Vec<_> = (0..nv).map(|_| Fr::rand(&mut rng)).collect();

        // commit
        let com = {
            let start = Instant::now();
            for _ in 0..repetition {
                let _commit = KZGMultilinearPCS::commit(&ck, &poly)?;
            }

            println!(
                "KZG commit for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );

            KZGMultilinearPCS::commit(&ck, &poly)?
        };

        // open
        let proof = {
            let start = Instant::now();
            for _ in 0..repetition {
                let _open = KZGMultilinearPCS::open(&ck, &poly, &point)?;
            }

            println!(
                "KZG open for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
            KZGMultilinearPCS::open(&ck, &poly, &point)?
        };
        let value = poly.evaluate(&point).unwrap();

        // verify

        {
            let start = Instant::now();
            for _ in 0..repetition {
                assert!(KZGMultilinearPCS::verify(
                    &vk, &com, &point, &value, &proof
                )?);
            }
            println!(
                "KZG verify for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }

    Ok(())
}
