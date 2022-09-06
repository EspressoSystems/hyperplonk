use arithmetic::{VPAuxInfo, VirtualPolynomial};
use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::test_rng;
use jf_primitives::pcs::{prelude::MultilinearKzgPCS, PolynomialCommitmentScheme};
use poly_iop::prelude::{
    identity_permutation_mle, PermutationCheck, PolyIOP, PolyIOPErrors, ProductCheck, SumCheck,
    ZeroCheck,
};
use std::{marker::PhantomData, rc::Rc, time::Instant};

type Kzg = MultilinearKzgPCS<Bls12_381>;

fn main() -> Result<(), PolyIOPErrors> {
    bench_permutation_check()?;
    println!("\n\n");
    bench_sum_check()?;
    println!("\n\n");
    bench_prod_check()?;
    println!("\n\n");
    bench_zero_check()
}

fn bench_sum_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();
    for degree in 2..4 {
        for nv in 4..25 {
            let repetition = if nv < 10 {
                100
            } else if nv < 20 {
                50
            } else {
                10
            };

            let (poly, asserted_sum) =
                VirtualPolynomial::rand(nv, (degree, degree + 1), 2, &mut rng)?;
            let poly_info = poly.aux_info.clone();
            let proof = {
                let start = Instant::now();
                for _ in 0..repetition {
                    let mut transcript = <PolyIOP<Fr> as SumCheck<Fr>>::init_transcript();
                    let _proof = <PolyIOP<Fr> as SumCheck<Fr>>::prove(&poly, &mut transcript)?;
                }

                println!(
                    "sum check proving time for {} variables and {} degree: {} ns",
                    nv,
                    degree,
                    start.elapsed().as_nanos() / repetition as u128
                );
                let mut transcript = <PolyIOP<Fr> as SumCheck<Fr>>::init_transcript();
                <PolyIOP<Fr> as SumCheck<Fr>>::prove(&poly, &mut transcript)?
            };

            {
                let start = Instant::now();

                for _ in 0..repetition {
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
                }
                println!(
                    "sum check verification time for {} variables and {} degree: {} ns",
                    nv,
                    degree,
                    start.elapsed().as_nanos() / repetition as u128
                );
            }

            println!("====================================");
        }
    }
    Ok(())
}

fn bench_zero_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();
    for degree in 2..4 {
        for nv in 4..20 {
            let repetition = if nv < 10 {
                100
            } else if nv < 20 {
                50
            } else {
                10
            };

            let poly = VirtualPolynomial::rand_zero(nv, (degree, degree + 1), 2, &mut rng)?;
            let poly_info = poly.aux_info.clone();
            let proof = {
                let start = Instant::now();
                let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
                transcript.append_message(b"testing", b"initializing transcript for testing")?;
                let proof = <PolyIOP<Fr> as ZeroCheck<Fr>>::prove(&poly, &mut transcript)?;

                println!(
                    "zero check proving time for {} variables and {} degree: {} ns",
                    nv,
                    degree,
                    start.elapsed().as_nanos() / repetition as u128
                );
                proof
            };

            {
                let start = Instant::now();
                let mut transcript = <PolyIOP<Fr> as ZeroCheck<Fr>>::init_transcript();
                transcript.append_message(b"testing", b"initializing transcript for testing")?;
                let zero_subclaim =
                    <PolyIOP<Fr> as ZeroCheck<Fr>>::verify(&proof, &poly_info, &mut transcript)?;
                assert!(
                    poly.evaluate(&zero_subclaim.point)? == zero_subclaim.expected_evaluation,
                    "wrong subclaim"
                );
                println!(
                    "zero check verification time for {} variables and {} degree: {} ns",
                    nv,
                    degree,
                    start.elapsed().as_nanos() / repetition as u128
                );
            }

            println!("====================================");
        }
    }
    Ok(())
}

fn bench_permutation_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();

    for nv in 4..20 {
        let srs = Kzg::gen_srs_for_testing(&mut rng, nv + 1)?;
        let (pcs_param, _) = Kzg::trim(&srs, nv + 1, Some(nv + 1))?;

        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let w = Rc::new(DenseMultilinearExtension::rand(nv, &mut rng));

        // s_perm is the identity map
        let s_perm = identity_permutation_mle(nv);

        let proof = {
            let start = Instant::now();
            let mut transcript =
                <PolyIOP<Fr> as PermutationCheck<Bls12_381, Kzg>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;

            let (proof, _q_x) = <PolyIOP<Fr> as PermutationCheck<Bls12_381, Kzg>>::prove(
                &pcs_param,
                &w,
                &w,
                &s_perm,
                &mut transcript,
            )?;

            println!(
                "permutation check proving time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
            proof
        };

        {
            let poly_info = VPAuxInfo {
                max_degree: 2,
                num_variables: nv,
                phantom: PhantomData::default(),
            };

            let start = Instant::now();
            let mut transcript =
                <PolyIOP<Fr> as PermutationCheck<Bls12_381, Kzg>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let _perm_check_sum_claim = <PolyIOP<Fr> as PermutationCheck<Bls12_381, Kzg>>::verify(
                &proof,
                &poly_info,
                &mut transcript,
            )?;
            println!(
                "permutation check verification time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }

    Ok(())
}

fn bench_prod_check() -> Result<(), PolyIOPErrors> {
    let mut rng = test_rng();

    for nv in 4..20 {
        let srs = Kzg::gen_srs_for_testing(&mut rng, nv + 1)?;
        let (pcs_param, _) = Kzg::trim(&srs, nv + 1, Some(nv + 1))?;

        let repetition = if nv < 10 {
            100
        } else if nv < 20 {
            50
        } else {
            10
        };

        let f: DenseMultilinearExtension<Fr> = DenseMultilinearExtension::rand(nv, &mut rng);
        let mut g = f.clone();
        g.evaluations.reverse();
        let f = Rc::new(f);
        let g = Rc::new(g);

        let proof = {
            let start = Instant::now();
            let mut transcript = <PolyIOP<Fr> as ProductCheck<Bls12_381, Kzg>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;

            let (proof, _prod_x) = <PolyIOP<Fr> as ProductCheck<Bls12_381, Kzg>>::prove(
                &pcs_param,
                &f,
                &g,
                &mut transcript,
            )?;

            println!(
                "product check proving time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
            proof
        };

        {
            let poly_info = VPAuxInfo {
                max_degree: 2,
                num_variables: nv,
                phantom: PhantomData::default(),
            };

            let start = Instant::now();
            let mut transcript = <PolyIOP<Fr> as ProductCheck<Bls12_381, Kzg>>::init_transcript();
            transcript.append_message(b"testing", b"initializing transcript for testing")?;
            let _perm_check_sum_claim = <PolyIOP<Fr> as ProductCheck<Bls12_381, Kzg>>::verify(
                &proof,
                &poly_info,
                &mut transcript,
            )?;
            println!(
                "product check verification time for {} variables: {} ns",
                nv,
                start.elapsed().as_nanos() / repetition as u128
            );
        }

        println!("====================================");
    }

    Ok(())
}
