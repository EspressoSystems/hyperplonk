use crate::{errors::PolyIOPErrors, structs::DomainInfo};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{end_timer, start_timer};
use std::{cmp::max, collections::HashMap, marker::PhantomData, rc::Rc};

/// A virtual polynomial is a list of multilinear polynomials
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

    /// Add a list of multilinear extensions that is meant to be multiplied
    /// together. The resulting polynomial will be multiplied by the scalar
    /// `coefficient`.
    pub fn add_product(
        &mut self,
        product: impl IntoIterator<Item = Rc<DenseMultilinearExtension<F>>>,
        coefficient: F,
    ) -> Result<(), PolyIOPErrors> {
        let product: Vec<Rc<DenseMultilinearExtension<F>>> = product.into_iter().collect();
        let mut indexed_product = Vec::with_capacity(product.len());
        assert!(!product.is_empty());
        self.domain_info.max_degree = max(self.domain_info.max_degree, product.len());
        for m in product {
            if m.num_vars != self.domain_info.num_variables {
                return Err(PolyIOPErrors::InvalidParameters(format!(
                    "product has a multiplicand with wrong number of variables {} vs {}",
                    m.num_vars, self.domain_info.num_variables
                )));
            }

            let m_ptr: *const DenseMultilinearExtension<F> = Rc::as_ptr(&m);
            if let Some(index) = self.raw_pointers_lookup_table.get(&m_ptr) {
                indexed_product.push(*index)
            } else {
                let curr_index = self.flattened_ml_extensions.len();
                self.flattened_ml_extensions.push(m.clone());
                self.raw_pointers_lookup_table.insert(m_ptr, curr_index);
                indexed_product.push(curr_index);
            }
        }
        self.products.push((coefficient, indexed_product));
        Ok(())
    }

    /// Evaluate the polynomial at point `point`
    pub fn evaluate(&self, point: &[F]) -> Result<F, PolyIOPErrors> {
        let start = start_timer!(|| "begin evaluation");

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
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;
    use ark_std::rand::{Rng, RngCore};

    pub fn random_product<F: PrimeField, R: RngCore>(
        nv: usize,
        num_multiplicands: usize,
        rng: &mut R,
    ) -> (Vec<Rc<DenseMultilinearExtension<F>>>, F) {
        let mut multiplicands = Vec::with_capacity(num_multiplicands);
        for _ in 0..num_multiplicands {
            multiplicands.push(Vec::with_capacity(1 << nv))
        }
        let mut sum = F::zero();

        for _ in 0..(1 << nv) {
            let mut product = F::one();
            for i in 0..num_multiplicands {
                let val = F::rand(rng);
                multiplicands[i].push(val);
                product *= val;
            }
            sum += product;
        }

        (
            multiplicands
                .into_iter()
                .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
                .collect(),
            sum,
        )
    }

    pub(crate) fn random_list_of_products<F: PrimeField, R: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut R,
    ) -> (VirtualPolynomial<F>, F) {
        let mut sum = F::zero();
        let mut poly = VirtualPolynomial::new(nv);
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let (product, product_sum) = random_product(nv, num_multiplicands, rng);
            let coefficient = F::rand(rng);
            poly.add_product(product.into_iter(), coefficient).unwrap();
            sum += product_sum * coefficient;
        }

        (poly, sum)
    }

    pub fn random_zero_product<F: PrimeField, R: RngCore>(
        nv: usize,
        _num_multiplicands: usize,
        _rng: &mut R,
    ) -> Vec<Rc<DenseMultilinearExtension<F>>> {
        let degree = 2;
        let mut multiplicands = Vec::with_capacity(degree);
        for _ in 0..degree {
            multiplicands.push(Vec::with_capacity(1 << nv))
        }
        let mut sum = F::zero();

        for _ in 0..(1 << nv) {
            let mut product = F::one();
            for i in 0..degree {
                let val = F::zero(); // F::rand(rng);
                multiplicands[i].push(val);
                product *= val;
            }
            sum += product;
        }

        // // last nv offsets the poly to 0
        // for i in 0..num_multiplicands - 1 {
        //     multiplicands[i].push(F::one());
        // }
        // multiplicands[num_multiplicands - 1].push(-sum);

        multiplicands
            .into_iter()
            .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
            .collect()
    }

    pub(crate) fn random_zero_list_of_products<F: PrimeField, R: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut R,
    ) -> VirtualPolynomial<F> {
        let mut poly = VirtualPolynomial::new(nv);
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let product = random_zero_product(nv, num_multiplicands, rng);
            let coefficient = F::rand(rng);
            poly.add_product(product.into_iter(), coefficient).unwrap();
        }

        poly
    }
}
