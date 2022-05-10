use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;

/// Prover State
pub struct ProverState<F: PrimeField> {
    /// sampled randomness given by the verifier
    pub challenges: Vec<F>,
    /// Stores the list of products that is meant to be added together. Each
    /// multiplicand is represented by the index in flattened_ml_extensions
    pub list_of_products: Vec<(F, Vec<usize>)>,
    /// Stores a list of multilinear extensions in which `self.list_of_products`
    /// points to
    pub flattened_ml_extensions: Vec<DenseMultilinearExtension<F>>,
    pub(crate) num_vars: usize,
    pub(crate) max_degree: usize,
    pub(crate) round: usize,
}
