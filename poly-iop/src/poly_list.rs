use crate::structs::AuxInfo;
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use std::{cmp::max, collections::HashMap, marker::PhantomData, rc::Rc};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PolynomialList<F: PrimeField> {
    /// Aux information about the multilinear system
    pub aux_info: AuxInfo<F>,
    /// list of reference to products (as usize) of multilinear extension
    pub products: Vec<(F, Vec<usize>)>,
    /// Stores multilinear extensions in which product multiplicand can refer
    /// to.
    pub flattened_ml_extensions: Vec<Rc<DenseMultilinearExtension<F>>>,
    /// Pointers to the above poly extensions
    raw_pointers_lookup_table: HashMap<*const DenseMultilinearExtension<F>, usize>,
}

impl<F: PrimeField> PolynomialList<F> {
    /// Returns an empty polynomial
    pub fn new(num_variables: usize) -> Self {
        PolynomialList {
            aux_info: AuxInfo {
                max_multiplicands: 0,
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
    ) {
        let product: Vec<Rc<DenseMultilinearExtension<F>>> = product.into_iter().collect();
        let mut indexed_product = Vec::with_capacity(product.len());
        assert!(product.len() > 0);
        self.aux_info.max_multiplicands = max(self.aux_info.max_multiplicands, product.len());
        for m in product {
            assert_eq!(
                m.num_vars, self.aux_info.num_variables,
                "product has a multiplicand with wrong number of variables"
            );
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
    }

    /// Evaluate the polynomial at point `point`
    pub fn evaluate(&self, point: &[F]) -> F {
        self.products
            .iter()
            .map(|(c, p)| {
                *c * p
                    .iter()
                    .map(|&i| self.flattened_ml_extensions[i].evaluate(point).unwrap())
                    .product::<F>()
            })
            .sum()
    }
}
