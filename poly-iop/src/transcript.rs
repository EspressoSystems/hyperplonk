use std::marker::PhantomData;

use ark_ff::PrimeField;
use merlin::Transcript;

use crate::{
    errors::PolyIOPErrors,
    structs::{DomainInfo, IOPProverMessage},
    to_bytes,
};

pub struct IOPTranscript<F: PrimeField> {
    transcript: Transcript,
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

impl<F: PrimeField> IOPTranscript<F> {
    /// create a new IOP transcript
    pub(crate) fn new(label: &'static [u8]) -> Self {
        Self {
            transcript: Transcript::new(label),
            phantom: PhantomData::default(),
        }
    }

    // append the message to the transcript
    pub(crate) fn append_message(
        &mut self,
        label: &'static [u8],
        msg: &[u8],
    ) -> Result<(), PolyIOPErrors> {
        self.transcript.append_message(label, msg);

        Ok(())
    }

    pub(crate) fn append_domain_info(
        &mut self,
        domain_info: &DomainInfo<F>,
    ) -> Result<(), PolyIOPErrors> {
        let message = format!(
            "max_mul {} num_var {}",
            domain_info.max_multiplicands, domain_info.num_variables
        );
        self.append_message(b"aux info", message.as_bytes())?;

        Ok(())
    }

    // append the message to the transcript
    pub(crate) fn append_field_element(
        &mut self,
        label: &'static [u8],
        field_elem: &F,
    ) -> Result<(), PolyIOPErrors> {
        self.append_message(label, &to_bytes!(field_elem)?)
    }

    pub(crate) fn append_prover_message(
        &mut self,
        prover_message: &IOPProverMessage<F>,
    ) -> Result<(), PolyIOPErrors> {
        for e in prover_message.evaluations.iter() {
            self.append_field_element(b"prover_message", e)?;
        }
        Ok(())
    }

    // generate the challenge for the current transcript
    // and append it to the transcript
    pub(crate) fn get_and_append_challenge(
        &mut self,
        label: &'static [u8],
    ) -> Result<F, PolyIOPErrors> {
        let mut buf = [0u8; 64];
        self.transcript.challenge_bytes(label, &mut buf);
        let challenge = F::from_le_bytes_mod_order(&buf);
        self.transcript
            .append_message(label, &to_bytes!(&challenge)?);
        Ok(challenge)
    }
}
