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
    is_empty: bool,
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

impl<F: PrimeField> IOPTranscript<F> {
    /// create a new IOP transcript
    pub(crate) fn new(label: &'static [u8]) -> Self {
        Self {
            transcript: Transcript::new(label),
            is_empty: true,
            phantom: PhantomData::default(),
        }
    }

    // append the message to the transcript
    pub fn append_message(
        &mut self,
        label: &'static [u8],
        msg: &[u8],
    ) -> Result<(), PolyIOPErrors> {
        self.transcript.append_message(label, msg);
        self.is_empty = false;
        Ok(())
    }

    pub(crate) fn append_domain_info(
        &mut self,
        domain_info: &DomainInfo<F>,
    ) -> Result<(), PolyIOPErrors> {
        let message = format!(
            "max_mul {} num_var {}",
            domain_info.max_degree, domain_info.num_variables
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

    // generate a list of challenges for the current transcript
    // and append it to the transcript
    pub(crate) fn get_and_append_challenge_vectors(
        &mut self,
        label: &'static [u8],
        len: usize,
    ) -> Result<Vec<F>, PolyIOPErrors> {
        //  we need to reject when transcript is empty
        let mut res = vec![];
        for _ in 0..len {
            res.push(self.get_and_append_challenge(label)?)
        }
        Ok(res)
    }
}
