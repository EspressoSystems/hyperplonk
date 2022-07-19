//! Module for PolyIOP transcript.
//! TODO(ZZ): move this module to HyperPlonk where the transcript will also be
//! useful.
//! TODO(ZZ): decide which APIs need to be public.

use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use merlin::Transcript;
use std::marker::PhantomData;

use crate::{errors::PolyIOPErrors, structs::IOPProverMessage, to_bytes, virtual_poly::VPAuxInfo};

/// An IOP transcript consists of a Merlin transcript and a flag `is_empty` to
/// indicate that if the transcript is empty.
///
/// It is associated with a prime field `F` for which challenges are generated
/// over.
///
/// The `is_empty` flag is useful in the case where a protocol is initiated by
/// the verifier, in which case the prover should start its phase by receiving a
/// `non-empty` transcript.
#[derive(Clone)]
pub struct IOPTranscript<F: PrimeField> {
    transcript: Transcript,
    is_empty: bool,
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

impl<F: PrimeField> IOPTranscript<F> {
    /// Create a new IOP transcript.
    pub fn new(label: &'static [u8]) -> Self {
        Self {
            transcript: Transcript::new(label),
            is_empty: true,
            phantom: PhantomData::default(),
        }
    }

    // Append the message to the transcript.
    pub fn append_message(
        &mut self,
        label: &'static [u8],
        msg: &[u8],
    ) -> Result<(), PolyIOPErrors> {
        self.transcript.append_message(label, msg);
        self.is_empty = false;
        Ok(())
    }

    // Append the aux information for a virtual polynomial.
    pub(crate) fn append_aux_info(&mut self, aux_info: &VPAuxInfo<F>) -> Result<(), PolyIOPErrors> {
        let message = format!(
            "max_mul {} num_var {}",
            aux_info.max_degree, aux_info.num_variables
        );
        self.append_message(b"aux info", message.as_bytes())?;

        Ok(())
    }

    // Append the message to the transcript.
    pub fn append_field_element(
        &mut self,
        label: &'static [u8],
        field_elem: &F,
    ) -> Result<(), PolyIOPErrors> {
        self.append_message(label, &to_bytes!(field_elem)?)
    }

    // Append the message to the transcript.
    pub fn append_serializable_element<S: CanonicalSerialize>(
        &mut self,
        label: &'static [u8],
        group_elem: &S,
    ) -> Result<(), PolyIOPErrors> {
        self.append_message(label, &to_bytes!(group_elem)?)
    }

    // Append a prover message to the transcript.
    pub(crate) fn append_prover_message(
        &mut self,
        prover_message: &IOPProverMessage<F>,
    ) -> Result<(), PolyIOPErrors> {
        for e in prover_message.evaluations.iter() {
            self.append_field_element(b"prover_message", e)?;
        }
        Ok(())
    }

    // Generate the challenge from the current transcript
    // and append it to the transcript.
    //
    // The output field element is statistical uniform as long
    // as the field has a size less than 2^384.
    pub fn get_and_append_challenge(&mut self, label: &'static [u8]) -> Result<F, PolyIOPErrors> {
        //  we need to reject when transcript is empty
        if self.is_empty {
            return Err(PolyIOPErrors::InvalidTranscript(
                "transcript is empty".to_string(),
            ));
        }

        let mut buf = [0u8; 64];
        self.transcript.challenge_bytes(label, &mut buf);
        let challenge = F::from_le_bytes_mod_order(&buf);
        self.transcript
            .append_message(label, &to_bytes!(&challenge)?);
        Ok(challenge)
    }

    // Generate a list of challenges from the current transcript
    // and append them to the transcript.
    //
    // The output field element are statistical uniform as long
    // as the field has a size less than 2^384.
    pub(crate) fn get_and_append_challenge_vectors(
        &mut self,
        label: &'static [u8],
        len: usize,
    ) -> Result<Vec<F>, PolyIOPErrors> {
        //  we need to reject when transcript is empty
        if self.is_empty {
            return Err(PolyIOPErrors::InvalidTranscript(
                "transcript is empty".to_string(),
            ));
        }

        let mut res = vec![];
        for _ in 0..len {
            res.push(self.get_and_append_challenge(label)?)
        }
        Ok(res)
    }
}
