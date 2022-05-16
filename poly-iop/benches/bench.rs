use std::time::Instant;

use ark_bls12_381::Fr;
use ark_std::test_rng;
use poly_iop::{PolyIOP, PolyIOPErrors, SumCheck, VirtualPolynomial, ZeroCheck};

fn main() -> Result<(), PolyIOPErrors> {
    bench_sum_check()?;
    println!("\n\n");
    bench_zero_check()
}

fn bench_sum_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();

    for nv in 4..25 {
        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let (poly, asserted_sum) = VirtualPolynomial::rand(nv, (2, 3), 2, &mut rng)?;
        let poly_info = poly.domain_info.clone();
        let proof = {
            let start = Instant::now();
            let mut transcript = <PolyIOP<Fr> as SumCheck<Fr>>::init_transcript();
            let proof = <PolyIOP<Fr> as SumCheck<Fr>>::prove(&poly, &mut transcript)?;

            println!(
                "sum check proving time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
            proof
        };

        {
            let start = Instant::now();
            let mut transcript = <PolyIOP<Fr> as SumCheck<Fr>>::init_transcript();
            let subclaim = <PolyIOP<Fr> as SumCheck<Fr>>::verify(
                asserted_sum,
                &proof,
                &poly_info,
                &mut transcript,
            )?;
            assert!(
                poly.evaluate(&subclaim.point).unwrap() == subclaim.expected_evaluation,
                "wrong subclaim"
            );

            println!(
                "sum check verification time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }
    Ok(())
}

fn bench_zero_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();

    for nv in 4..20 {
        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let poly = VirtualPolynomial::rand_zero(nv, (2, 3), 2, &mut rng)?;

        let poly_info = poly.domain_info.clone();
        let proof = {
            let start = Instant::now();
            let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let proof = <PolyIOP<Fr> as ZeroCheck<Fr>>::prove(&poly, &mut transcript)?;

            println!(
                "zero check proving time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
            proof
        };

        {
            let start = Instant::now();
            let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let subclaim =
                <PolyIOP<Fr> as ZeroCheck<Fr>>::verify(&proof, &poly_info, &mut transcript)?.0;
            assert!(
                poly.evaluate(&subclaim.point)? == subclaim.expected_evaluation,
                "wrong subclaim"
            );
            println!(
                "zero check verification time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }

    Ok(())
}
