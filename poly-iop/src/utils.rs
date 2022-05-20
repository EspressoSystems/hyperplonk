//! useful macros.

/// Takes as input a struct, and converts them to a series of bytes. All traits
/// that implement `CanonicalSerialize` can be automatically converted to bytes
/// in this manner.
#[macro_export]
macro_rules! to_bytes {
    ($x:expr) => {{
        let mut buf = ark_std::vec![];
        ark_serialize::CanonicalSerialize::serialize($x, &mut buf).map(|_| buf)
    }};
}

#[cfg(test)]
pub(crate) fn bit_decompose(input: u64, num_var: usize) -> Vec<bool> {
    let mut res = Vec::with_capacity(num_var);
    let mut i = input;
    for _ in 0..num_var {
        res.push(i & 1 == 1);
        i >>= 1;
    }
    res
}

#[test]
fn test_to_bytes() {
    use ark_bls12_381::Fr;
    use ark_serialize::CanonicalSerialize;
    use ark_std::One;
    let f1 = Fr::one();

    let mut bytes = ark_std::vec![];
    f1.serialize(&mut bytes).unwrap();
    assert_eq!(bytes, to_bytes!(&f1).unwrap());
}
