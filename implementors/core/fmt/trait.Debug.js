(function() {var implementors = {};
implementors["arithmetic"] = [{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"enum\" href=\"arithmetic/enum.ArithErrors.html\" title=\"enum arithmetic::ArithErrors\">ArithErrors</a>","synthetic":false,"types":["arithmetic::errors::ArithErrors"]},{"text":"impl&lt;F:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PrimeField&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"arithmetic/struct.VirtualPolynomial.html\" title=\"struct arithmetic::VirtualPolynomial\">VirtualPolynomial</a>&lt;F&gt;","synthetic":false,"types":["arithmetic::virtual_polynomial::VirtualPolynomial"]},{"text":"impl&lt;F:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PrimeField&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"arithmetic/struct.VPAuxInfo.html\" title=\"struct arithmetic::VPAuxInfo\">VPAuxInfo</a>&lt;F&gt;","synthetic":false,"types":["arithmetic::virtual_polynomial::VPAuxInfo"]}];
implementors["pcs"] = [{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"enum\" href=\"pcs/prelude/enum.PCSErrors.html\" title=\"enum pcs::prelude::PCSErrors\">PCSErrors</a>","synthetic":false,"types":["pcs::errors::PCSErrors"]},{"text":"impl&lt;E:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PairingEngine&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearUniversalParams.html\" title=\"struct pcs::prelude::MultilinearUniversalParams\">MultilinearUniversalParams</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G2Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,&nbsp;</span>","synthetic":false,"types":["pcs::multilinear_kzg::srs::MultilinearUniversalParams"]},{"text":"impl&lt;E:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PairingEngine&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearProverParam.html\" title=\"struct pcs::prelude::MultilinearProverParam\">MultilinearProverParam</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G1Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G1Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G2Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,&nbsp;</span>","synthetic":false,"types":["pcs::multilinear_kzg::srs::MultilinearProverParam"]},{"text":"impl&lt;E:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PairingEngine&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearVerifierParam.html\" title=\"struct pcs::prelude::MultilinearVerifierParam\">MultilinearVerifierParam</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G1Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G2Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G2Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,&nbsp;</span>","synthetic":false,"types":["pcs::multilinear_kzg::srs::MultilinearVerifierParam"]},{"text":"impl&lt;E:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PairingEngine&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.Proof.html\" title=\"struct pcs::prelude::Proof\">Proof</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G1Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,&nbsp;</span>","synthetic":false,"types":["pcs::multilinear_kzg::Proof"]},{"text":"impl&lt;E:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PairingEngine&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.BatchProof.html\" title=\"struct pcs::prelude::BatchProof\">BatchProof</a>&lt;E&gt;","synthetic":false,"types":["pcs::multilinear_kzg::BatchProof"]},{"text":"impl&lt;E:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PairingEngine&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.Commitment.html\" title=\"struct pcs::prelude::Commitment\">Commitment</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G1Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,&nbsp;</span>","synthetic":false,"types":["pcs::structs::Commitment"]},{"text":"impl&lt;E:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PairingEngine&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.UnivariateUniversalParams.html\" title=\"struct pcs::prelude::UnivariateUniversalParams\">UnivariateUniversalParams</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G1Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G2Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G2Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,&nbsp;</span>","synthetic":false,"types":["pcs::univariate_kzg::srs::UnivariateUniversalParams"]},{"text":"impl&lt;C:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + AffineCurve&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.UnivariateProverParam.html\" title=\"struct pcs::prelude::UnivariateProverParam\">UnivariateProverParam</a>&lt;C&gt;","synthetic":false,"types":["pcs::univariate_kzg::srs::UnivariateProverParam"]},{"text":"impl&lt;E:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PairingEngine&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"pcs/prelude/struct.UnivariateVerifierParam.html\" title=\"struct pcs::prelude::UnivariateVerifierParam\">UnivariateVerifierParam</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G1Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G2Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;E::G2Affine: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a>,&nbsp;</span>","synthetic":false,"types":["pcs::univariate_kzg::srs::UnivariateVerifierParam"]}];
implementors["poly_iop"] = [{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"enum\" href=\"poly_iop/prelude/enum.PolyIOPErrors.html\" title=\"enum poly_iop::prelude::PolyIOPErrors\">PolyIOPErrors</a>","synthetic":false,"types":["poly_iop::errors::PolyIOPErrors"]},{"text":"impl&lt;F:&nbsp;<a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> + PrimeField&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"poly_iop/struct.PolyIOP.html\" title=\"struct poly_iop::PolyIOP\">PolyIOP</a>&lt;F&gt;","synthetic":false,"types":["poly_iop::PolyIOP"]}];
implementors["transcript"] = [{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.62.1/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"enum\" href=\"transcript/enum.TranscriptErrors.html\" title=\"enum transcript::TranscriptErrors\">TranscriptErrors</a>","synthetic":false,"types":["transcript::errors::TranscriptErrors"]}];
if (window.register_implementors) {window.register_implementors(implementors);} else {window.pending_implementors = implementors;}})()