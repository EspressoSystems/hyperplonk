use ark_ff::PrimeField;

/// Verifier State
pub struct VerifierState<F: PrimeField> {
    round: usize,
    nv: usize,
    max_multiplicands: usize,
    finished: bool,
    /// a list storing the univariate polynomial in evaluation form sent by the
    /// prover at each round
    polynomials_received: Vec<Vec<F>>,
    /// a list storing the randomness sampled by the verifier at each round
    challenges: Vec<F>,
}
