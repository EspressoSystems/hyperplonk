// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

//! useful macros.

/// Takes as input a struct, and converts them to a series of bytes. All traits
/// that implement `CanonicalSerialize` can be automatically converted to bytes
/// in this manner.
#[macro_export]
macro_rules! to_bytes {
    ($x:expr) => {{
        let mut buf = ark_std::vec![];
        ark_serialize::CanonicalSerialize::serialize_compressed($x, &mut buf).map(|_| buf)
    }};
}

#[cfg(test)]
mod test {
    use ark_bls12_381::Fr;
    use ark_serialize::CanonicalSerialize;
    use ark_std::One;

    #[test]
    fn test_to_bytes() {
        let f1 = Fr::one();

        let mut bytes = ark_std::vec![];
        f1.serialize_compressed(&mut bytes).unwrap();
        assert_eq!(bytes, to_bytes!(&f1).unwrap());
    }
}
