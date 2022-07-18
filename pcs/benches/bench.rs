use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::UniformRand;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::test_rng;
use pcs::{
    prelude::{KZGMultilinearPC, PCSErrors, PCSScheme},
    StructuredReferenceString,
};
use std::time::Instant;

fn main() -> Result<(), PCSErrors> {
    bench_pcs()
}

fn bench_pcs() -> Result<(), PCSErrors> {
    let mut rng = test_rng();

    // normal polynomials
    let uni_params = KZGMultilinearPC::<Bls12_381>::setup(&mut rng, 18)?;

    for nv in 4..19 {
        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let poly = DenseMultilinearExtension::rand(nv, &mut rng);
        let (ck, vk) = uni_params.trim(nv)?;
        let point: Vec<_> = (0..nv).map(|_| Fr::rand(&mut rng)).collect();

        // commit
        let com = {
            let start = Instant::now();
            for _ in 0..repetition {
                let _commit = KZGMultilinearPC::commit(&ck, &poly)?;
            }

            println!(
                "KZG commit for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );

            KZGMultilinearPC::commit(&ck, &poly)?
        };

        // open
        let proof = {
            let start = Instant::now();
            for _ in 0..repetition {
                let _open = KZGMultilinearPC::open(&ck, &poly, &point)?;
            }

            println!(
                "KZG open for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
            KZGMultilinearPC::open(&ck, &poly, &point)?
        };
        let value = poly.evaluate(&point).unwrap();

        // verify

        {
            let start = Instant::now();
            for _ in 0..repetition {
                assert!(KZGMultilinearPC::verify(&vk, &com, &point, &value, &proof)?);
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
