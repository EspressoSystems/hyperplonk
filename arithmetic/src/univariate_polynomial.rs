// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

// TODO: remove
#![allow(dead_code)]

use crate::{bit_decompose, ArithErrors};
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, Radix2EvaluationDomain,
};
use ark_std::log2;

/// Given a list of points, build `l(points)` which is a list of univariate
/// polynomials that goes through the points; extend the dimension of the points
/// by `log(points.len())` if `with_suffix` is set.
pub fn build_l<F: PrimeField>(
    points: &[Vec<F>],
    domain: &Radix2EvaluationDomain<F>,
    with_suffix: bool,
) -> Result<Vec<DensePolynomial<F>>, ArithErrors> {
    let mut uni_polys = Vec::new();
    if with_suffix {
        // 1.1 build the indexes and the univariate polys that go through the indexes
        let prefix_len = log2(points.len()) as usize;
        let indexes: Vec<Vec<bool>> = (0..points.len())
            .map(|x| bit_decompose(x as u64, prefix_len))
            .collect();
        for i in 0..prefix_len {
            let eval: Vec<F> = indexes
                .iter()
                .map(|x| F::from(x[prefix_len - i - 1]))
                .collect();

            uni_polys.push(Evaluations::from_vec_and_domain(eval, *domain).interpolate());
        }
    }
    // 1.2 build the actual univariate polys that go through the points
    uni_polys.extend_from_slice(build_l_internal(points, domain)?.as_slice());

    Ok(uni_polys)
}

/// Given a list of points, build `l(points)` which is a list of univariate
/// polynomials that goes through the points.
pub(crate) fn build_l_internal<F: PrimeField>(
    points: &[Vec<F>],
    domain: &Radix2EvaluationDomain<F>,
) -> Result<Vec<DensePolynomial<F>>, ArithErrors> {
    let mut uni_polys = Vec::new();
    let num_var = points[0].len();
    // build the actual univariate polys that go through the points
    for i in 0..num_var {
        let mut eval: Vec<F> = points.iter().map(|x| x[i]).collect();
        eval.extend_from_slice(vec![F::zero(); domain.size as usize - eval.len()].as_slice());
        uni_polys.push(Evaluations::from_vec_and_domain(eval, *domain).interpolate())
    }
    Ok(uni_polys)
}

/// get the domain for the univariate polynomial
#[inline]
pub fn get_uni_domain<F: PrimeField>(
    uni_poly_degree: usize,
) -> Result<Radix2EvaluationDomain<F>, ArithErrors> {
    let domain = match Radix2EvaluationDomain::<F>::new(uni_poly_degree) {
        Some(p) => p,
        None => {
            return Err(ArithErrors::InvalidParameters(
                "failed to build radix 2 domain".to_string(),
            ))
        },
    };
    Ok(domain)
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_ff::{MontFp, One};
    use ark_poly::DenseUVPolynomial;

    #[test]
    fn test_build_l_with_suffix() -> Result<(), ArithErrors> {
        // point 1 is [1, 2]
        let point1 = vec![Fr::from(1u64), Fr::from(2u64)];

        // point 2 is [3, 4]
        let point2 = vec![Fr::from(3u64), Fr::from(4u64)];

        // point 3 is [5, 6]
        let point3 = vec![Fr::from(5u64), Fr::from(6u64)];

        {
            let domain = get_uni_domain::<Fr>(2)?;
            let l = build_l(&[point1.clone(), point2.clone()], &domain, true)?;

            // roots: [1, -1]
            // l0 = -1/2 * x + 1/2
            // l1 = -x + 2
            // l2 = -x + 3
            let l0 = DensePolynomial::from_coefficients_vec(vec![
                Fr::one() / Fr::from(2u64),
                -Fr::one() / Fr::from(2u64),
            ]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![Fr::from(2u64), -Fr::one()]);
            let l2 = DensePolynomial::from_coefficients_vec(vec![Fr::from(3u64), -Fr::one()]);

            assert_eq!(l0, l[0], "l0 not equal");
            assert_eq!(l1, l[1], "l1 not equal");
            assert_eq!(l2, l[2], "l2 not equal");
        }

        {
            let domain = get_uni_domain::<Fr>(3)?;
            let l = build_l::<Fr>(&[point1, point2, point3], &domain, true)?;

            // sage: q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
            // sage: P.<x> = PolynomialRing(Zmod(q))
            // sage: root1 = 1
            // sage: root2 = 0x8D51CCCE760304D0EC030002760300000001000000000000
            // sage: root3 = -1
            // sage: root4 = -root2
            // Arkwork's code is a bit wired: it also interpolate (root4, 0)
            // which returns a degree 3 polynomial, instead of degree 2

            // ========================
            // l0: [0, 0, 1]
            // ========================
            // sage: points = [(root1, 0), (root2, 0), (root3, 1), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 13108968793781547619861935127046491459422638125131909455650914674984645296128*x^3 +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888385*x^2 +
            // 13108968793781547619861935127046491459422638125131909455650914674984645296128*x +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888385
            let l0 = DensePolynomial::from_coefficients_vec(vec![
                MontFp!(
                    "39326906381344642859585805381139474378267914375395728366952744024953935888385"
                ),
                MontFp!(
                    "13108968793781547619861935127046491459422638125131909455650914674984645296128"
                ),
                MontFp!(
                    "39326906381344642859585805381139474378267914375395728366952744024953935888385"
                ),
                MontFp!(
                    "13108968793781547619861935127046491459422638125131909455650914674984645296128"
                ),
            ]);

            // ========================
            // l1: [0, 1, 0]
            // ========================
            // sage: points = [(root1, 0), (root2, 1), (root3, 0), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 866286206518413079694067382671935694567563117191340490752*x^3 +
            // 13108968793781547619861935127046491459422638125131909455650914674984645296128*x^2 +
            // 52435875175126190478581454301667552757996485117855702128036095582747240693761*x +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888385
            let l1 = DensePolynomial::from_coefficients_vec(vec![
                MontFp!(
                    "39326906381344642859585805381139474378267914375395728366952744024953935888385"
                ),
                MontFp!(
                    "52435875175126190478581454301667552757996485117855702128036095582747240693761"
                ),
                MontFp!(
                    "13108968793781547619861935127046491459422638125131909455650914674984645296128"
                ),
                MontFp!("866286206518413079694067382671935694567563117191340490752"),
            ]);

            // ========================
            // l2: [1, 3, 5]
            // ========================
            // sage: points = [(root1, 1), (root2, 3), (root3, 5), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 2598858619555239239082202148015807083702689351574021472255*x^3 +
            // 13108968793781547619861935127046491459422638125131909455650914674984645296129*x^2 +
            // 52435875175126190476848881888630726598608350352511830738900969348364559712256*x +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888387
            let l2 = DensePolynomial::from_coefficients_vec(vec![
                MontFp!(
                    "39326906381344642859585805381139474378267914375395728366952744024953935888387"
                ),
                MontFp!(
                    "52435875175126190476848881888630726598608350352511830738900969348364559712256"
                ),
                MontFp!(
                    "13108968793781547619861935127046491459422638125131909455650914674984645296129"
                ),
                MontFp!("2598858619555239239082202148015807083702689351574021472255"),
            ]);

            // ========================
            // l3: [2, 4, 6]
            // ========================
            // sage: points = [(root1, 2), (root2, 4), (root3, 6), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 3465144826073652318776269530687742778270252468765361963007*x^3 +
            // x^2 +
            // 52435875175126190475982595682112313518914282969839895044333406231173219221504*x +
            // 3
            let l3 = DensePolynomial::from_coefficients_vec(vec![
                Fr::from(3u64),
                MontFp!(
                    "52435875175126190475982595682112313518914282969839895044333406231173219221504"
                ),
                Fr::one(),
                MontFp!("3465144826073652318776269530687742778270252468765361963007"),
            ]);

            assert_eq!(l0, l[0], "l0 not equal");
            assert_eq!(l1, l[1], "l1 not equal");
            assert_eq!(l2, l[2], "l2 not equal");
            assert_eq!(l3, l[3], "l3 not equal");
        }
        Ok(())
    }

    #[test]
    fn test_build_l() -> Result<(), ArithErrors> {
        // point 1 is [1, 2]
        let point1 = vec![Fr::from(1u64), Fr::from(2u64)];

        // point 2 is [3, 4]
        let point2 = vec![Fr::from(3u64), Fr::from(4u64)];

        // point 3 is [5, 6]
        let point3 = vec![Fr::from(5u64), Fr::from(6u64)];

        {
            let domain = get_uni_domain::<Fr>(2)?;
            let l = build_l(&[point1.clone(), point2.clone()], &domain, false)?;

            // roots: [1, -1]
            // l0 = -x + 2
            // l1 = -x + 3
            let l0 = DensePolynomial::from_coefficients_vec(vec![Fr::from(2u64), -Fr::one()]);
            let l1 = DensePolynomial::from_coefficients_vec(vec![Fr::from(3u64), -Fr::one()]);

            assert_eq!(l0, l[0], "l0 not equal");
            assert_eq!(l1, l[1], "l1 not equal");
        }

        {
            let domain = get_uni_domain::<Fr>(3)?;
            let l = build_l(&[point1, point2, point3], &domain, false)?;

            // sage: q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
            // sage: P.<x> = PolynomialRing(Zmod(q))
            // sage: root1 = 1
            // sage: root2 = 0x8D51CCCE760304D0EC030002760300000001000000000000
            // sage: root3 = -1
            // sage: root4 = -root2
            // Arkwork's code is a bit wired: it also interpolate (root4, 0)
            // which returns a degree 3 polynomial, instead of degree 2

            // ========================
            // l0: [1, 3, 5]
            // ========================
            // sage: points = [(root1, 1), (root2, 3), (root3, 5), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 2598858619555239239082202148015807083702689351574021472255*x^3 +
            // 13108968793781547619861935127046491459422638125131909455650914674984645296129*x^2 +
            // 52435875175126190476848881888630726598608350352511830738900969348364559712256*x +
            // 39326906381344642859585805381139474378267914375395728366952744024953935888387
            let l0 = DensePolynomial::from_coefficients_vec(vec![
                MontFp!(
                    "39326906381344642859585805381139474378267914375395728366952744024953935888387"
                ),
                MontFp!(
                    "52435875175126190476848881888630726598608350352511830738900969348364559712256"
                ),
                MontFp!(
                    "13108968793781547619861935127046491459422638125131909455650914674984645296129"
                ),
                MontFp!("2598858619555239239082202148015807083702689351574021472255"),
            ]);

            // ========================
            // l1: [2, 4, 6]
            // ========================
            // sage: points = [(root1, 2), (root2, 4), (root3, 6), (root4, 0)]
            // sage: P.lagrange_polynomial(points)
            // 3465144826073652318776269530687742778270252468765361963007*x^3 +
            // x^2 +
            // 52435875175126190475982595682112313518914282969839895044333406231173219221504*x +
            // 3
            let l1 = DensePolynomial::from_coefficients_vec(vec![
                Fr::from(3u64),
                MontFp!(
                    "52435875175126190475982595682112313518914282969839895044333406231173219221504"
                ),
                Fr::one(),
                MontFp!("3465144826073652318776269530687742778270252468765361963007"),
            ]);

            assert_eq!(l0, l[0], "l0 not equal");
            assert_eq!(l1, l[1], "l1 not equal");
        }
        Ok(())
    }
}
