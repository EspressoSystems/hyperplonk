use crate::{errors::PolyIOPErrors, structs::DomainInfo};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{
    end_timer,
    rand::{Rng, RngCore},
    start_timer,
};
use std::{cmp::max, collections::HashMap, marker::PhantomData, ops::Add, rc::Rc};

#[rustfmt::skip]
/// A virtual polynomial is a sum of products of multilinear polynomials;
/// where the multilinear polynomials are stored via their multilinear
/// extensions:  `(coefficient, DenseMultilinearExtension)`
///
/// * Number of products n = `polynomial.products.len()`,
/// * Number of multiplicands of ith product m_i =
///   `polynomial.products[i].1.len()`,
/// * Coefficient of ith product c_i = `polynomial.products[i].0`
///
/// The resulting polynomial is
///
/// $$\sum_{i=0}^{n} c_i \cdot \prod_{j=0}^{m_i} P_{ij} $$
///
/// Example:
///  f = c0 * f0 * f1 * f2 + c1 * f3 * f4
/// where f0 ... f4 are multilinear polynomials
///
/// - flattened_ml_extensions stores the multilinear extension representation of
///   f0, f1, f2, f3 and f4
/// - products is 
///     \[ 
///         (c0, \[0, 1, 2\]), 
///         (c1, \[3, 4\]) 
///     \]
/// - raw_pointers_lookup_table maps fi to i
///
#[derive(Clone, Debug, Default, PartialEq)]
pub struct VirtualPolynomial<F: PrimeField> {
    /// Aux information about the multilinear polynomial
    pub domain_info: DomainInfo<F>,
    /// list of reference to products (as usize) of multilinear extension
    pub products: Vec<(F, Vec<usize>)>,
    /// Stores multilinear extensions in which product multiplicand can refer
    /// to.
    pub flattened_ml_extensions: Vec<Rc<DenseMultilinearExtension<F>>>,
    /// Pointers to the above poly extensions
    raw_pointers_lookup_table: HashMap<*const DenseMultilinearExtension<F>, usize>,
}

impl<F: PrimeField> Add for &VirtualPolynomial<F> {
    type Output = VirtualPolynomial<F>;
    fn add(self, other: &VirtualPolynomial<F>) -> Self::Output {
        let start = start_timer!(|| "virtual poly add");
        let mut res = self.clone();
        for products in other.products.iter() {
            let cur: Vec<Rc<DenseMultilinearExtension<F>>> = products
                .1
                .iter()
                .map(|&x| other.flattened_ml_extensions[x].clone())
                .collect();

            res.add_mle_list(cur, products.0)
                .expect("add product failed");
        }
        end_timer!(start);
        res
    }
}

impl<F: PrimeField> VirtualPolynomial<F> {
    /// Returns an empty polynomial
    pub fn new(num_variables: usize) -> Self {
        VirtualPolynomial {
            domain_info: DomainInfo {
                max_degree: 0,
                num_variables,
                phantom: PhantomData::default(),
            },
            products: Vec::new(),
            flattened_ml_extensions: Vec::new(),
            raw_pointers_lookup_table: HashMap::new(),
        }
    }

    /// Returns an new virtual polynomial from a MLE
    pub fn new_from_mle(mle: Rc<DenseMultilinearExtension<F>>, coefficient: F) -> Self {
        let mle_ptr: *const DenseMultilinearExtension<F> = Rc::as_ptr(&mle);
        let mut hm = HashMap::new();
        hm.insert(mle_ptr, 0);

        VirtualPolynomial {
            domain_info: DomainInfo {
                // The max degree is the max degree of any individual variable
                max_degree: 1,
                num_variables: mle.num_vars,
                phantom: PhantomData::default(),
            },
            products: vec![(coefficient, vec![0])],
            flattened_ml_extensions: vec![mle],
            raw_pointers_lookup_table: hm,
        }
    }

    /// Add a list of multilinear extensions that is meant to be multiplied
    /// together. The resulting polynomial will be multiplied by the scalar
    /// `coefficient`.
    pub fn add_mle_list(
        &mut self,
        mle_list: impl IntoIterator<Item = Rc<DenseMultilinearExtension<F>>>,
        coefficient: F,
    ) -> Result<(), PolyIOPErrors> {
        let mle_list: Vec<Rc<DenseMultilinearExtension<F>>> = mle_list.into_iter().collect();
        let mut indexed_product = Vec::with_capacity(mle_list.len());
        assert!(!mle_list.is_empty());
        self.domain_info.max_degree = max(self.domain_info.max_degree, mle_list.len());
        for mle in mle_list {
            if mle.num_vars != self.domain_info.num_variables {
                return Err(PolyIOPErrors::InvalidParameters(format!(
                    "product has a multiplicand with wrong number of variables {} vs {}",
                    mle.num_vars, self.domain_info.num_variables
                )));
            }

            let mle_ptr: *const DenseMultilinearExtension<F> = Rc::as_ptr(&mle);
            if let Some(index) = self.raw_pointers_lookup_table.get(&mle_ptr) {
                indexed_product.push(*index)
            } else {
                let curr_index = self.flattened_ml_extensions.len();
                self.flattened_ml_extensions.push(mle.clone());
                self.raw_pointers_lookup_table.insert(mle_ptr, curr_index);
                indexed_product.push(curr_index);
            }
        }
        self.products.push((coefficient, indexed_product));
        Ok(())
    }

    /// Multiple the current VirtualPolynomial by an MLE:
    /// - add the MLE to the MLE list
    /// - multiple each product by MLE and its coefficient
    pub fn mul_by_mle(
        &mut self,
        mle: Rc<DenseMultilinearExtension<F>>,
        coefficient: F,
    ) -> Result<(), PolyIOPErrors> {
        let start = start_timer!(|| "mul by mle");
        let mle_ptr: *const DenseMultilinearExtension<F> = Rc::as_ptr(&mle);

        let mle_index = match self.raw_pointers_lookup_table.get(&mle_ptr) {
            Some(&p) => p,
            None => {
                self.raw_pointers_lookup_table
                    .insert(mle_ptr, self.flattened_ml_extensions.len());
                self.flattened_ml_extensions.push(mle);
                self.flattened_ml_extensions.len() - 1
            },
        };

        for (prod_coef, indices) in self.products.iter_mut() {
            indices.push(mle_index);
            *prod_coef *= coefficient;
        }
        self.domain_info.max_degree += 1;
        end_timer!(start);
        Ok(())
    }

    /// Evaluate the polynomial at point `point`
    pub fn evaluate(&self, point: &[F]) -> Result<F, PolyIOPErrors> {
        let start = start_timer!(|| "evaluation");

        if self.domain_info.num_variables != point.len() {
            return Err(PolyIOPErrors::InvalidParameters(format!(
                "wrong number of variables {} vs {}",
                self.domain_info.num_variables,
                point.len()
            )));
        }

        let evals: Vec<F> = self
            .flattened_ml_extensions
            .iter()
            .map(|x| {
                x.evaluate(point).unwrap() // safe unwrap here since we have
                                           // already checked that num_var
                                           // matches
            })
            .collect();

        let res = self
            .products
            .iter()
            .map(|(c, p)| *c * p.iter().map(|&i| evals[i]).product::<F>())
            .sum();

        end_timer!(start);
        Ok(res)
    }

    /// Sample a random virtual polynomial, return the polynomial and its sum.
    pub fn rand<R: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut R,
    ) -> Result<(Self, F), PolyIOPErrors> {
        let start = start_timer!(|| "sample random virtual polynomial");

        let mut sum = F::zero();
        let mut poly = VirtualPolynomial::new(nv);
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let (product, product_sum) = random_mle_list(nv, num_multiplicands, rng);
            let coefficient = F::rand(rng);
            poly.add_mle_list(product.into_iter(), coefficient)?;
            sum += product_sum * coefficient;
        }

        end_timer!(start);
        Ok((poly, sum))
    }

    /// Sample a random virtual polynomial that evaluates to zero everywhere on
    /// the boolean hypercube.
    pub fn rand_zero<R: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut R,
    ) -> Result<Self, PolyIOPErrors> {
        let mut poly = VirtualPolynomial::new(nv);
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let product = random_zero_mle_list(nv, num_multiplicands, rng);
            let coefficient = F::rand(rng);
            poly.add_mle_list(product.into_iter(), coefficient).unwrap();
        }

        Ok(poly)
    }
}

/// Sample a random list of multilinear polynomials.
/// Returns
/// - the list of polynomials,
/// - its sum of polynomial evaluations over the boolean hypercube.
fn random_mle_list<F: PrimeField, R: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut R,
) -> (Vec<Rc<DenseMultilinearExtension<F>>>, F) {
    let start = start_timer!(|| "sample random mle list");
    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv))
    }
    let mut sum = F::zero();

    for _ in 0..(1 << nv) {
        let mut product = F::one();

        for e in multiplicands.iter_mut() {
            let val = F::rand(rng);
            e.push(val);
            product *= val;
        }
        sum += product;
    }

    let list = multiplicands
        .into_iter()
        .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    (list, sum)
}

// Build a randomize list of mle-s whose sum is zero.
pub fn random_zero_mle_list<F: PrimeField, R: RngCore>(
    nv: usize,
    degree: usize,
    rng: &mut R,
) -> Vec<Rc<DenseMultilinearExtension<F>>> {
    let start = start_timer!(|| "sample random zero mle list");

    let mut multiplicands = Vec::with_capacity(degree);
    for _ in 0..degree {
        multiplicands.push(Vec::with_capacity(1 << nv))
    }
    for _ in 0..(1 << nv) {
        multiplicands[0].push(F::zero());
        for e in multiplicands.iter_mut().skip(1) {
            e.push(F::rand(rng));
        }
    }

    let list = multiplicands
        .into_iter()
        .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
        .collect();

    end_timer!(start);
    list
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_ff::UniformRand;
    use ark_std::test_rng;

    #[test]
    fn test_virtual_polynomial_additions() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();
        for nv in 2..5 {
            for num_products in 2..5 {
                let base: Vec<Fr> = (0..nv).map(|_| Fr::rand(&mut rng)).collect();

                let (a, _a_sum) =
                    VirtualPolynomial::<Fr>::rand(nv, (2, 3), num_products, &mut rng)?;
                let (b, _b_sum) =
                    VirtualPolynomial::<Fr>::rand(nv, (2, 3), num_products, &mut rng)?;
                let c = &a + &b;

                assert_eq!(
                    a.evaluate(base.as_ref())? + b.evaluate(base.as_ref())?,
                    c.evaluate(base.as_ref())?
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_virtual_polynomial_mul_by_mle() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();
        for nv in 2..5 {
            for num_products in 2..5 {
                let base: Vec<Fr> = (0..nv).map(|_| Fr::rand(&mut rng)).collect();

                let (a, _a_sum) =
                    VirtualPolynomial::<Fr>::rand(nv, (2, 3), num_products, &mut rng)?;
                let (b, _b_sum) = random_mle_list(nv, 1, &mut rng);
                let b_mle = b[0].clone();
                let coeff = Fr::rand(&mut rng);
                let b_vp = VirtualPolynomial::new_from_mle(b_mle.clone(), coeff);

                let mut c = a.clone();

                c.mul_by_mle(b_mle, coeff)?;

                assert_eq!(
                    a.evaluate(base.as_ref())? * b_vp.evaluate(base.as_ref())?,
                    c.evaluate(base.as_ref())?
                );
            }
        }

        Ok(())
    }
}
