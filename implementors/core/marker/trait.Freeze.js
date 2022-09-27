(function() {var implementors = {};
implementors["arithmetic"] = [{"text":"impl Freeze for <a class=\"enum\" href=\"arithmetic/enum.ArithErrors.html\" title=\"enum arithmetic::ArithErrors\">ArithErrors</a>","synthetic":true,"types":["arithmetic::errors::ArithErrors"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"arithmetic/struct.VirtualPolynomial.html\" title=\"struct arithmetic::VirtualPolynomial\">VirtualPolynomial</a>&lt;F&gt;","synthetic":true,"types":["arithmetic::virtual_polynomial::VirtualPolynomial"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"arithmetic/struct.VPAuxInfo.html\" title=\"struct arithmetic::VPAuxInfo\">VPAuxInfo</a>&lt;F&gt;","synthetic":true,"types":["arithmetic::virtual_polynomial::VPAuxInfo"]}];
implementors["hyperplonk"] = [{"text":"impl Freeze for <a class=\"struct\" href=\"hyperplonk/prelude/struct.CustomizedGates.html\" title=\"struct hyperplonk::prelude::CustomizedGates\">CustomizedGates</a>","synthetic":true,"types":["hyperplonk::custom_gate::CustomizedGates"]},{"text":"impl Freeze for <a class=\"enum\" href=\"hyperplonk/prelude/enum.HyperPlonkErrors.html\" title=\"enum hyperplonk::prelude::HyperPlonkErrors\">HyperPlonkErrors</a>","synthetic":true,"types":["hyperplonk::errors::HyperPlonkErrors"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"hyperplonk/prelude/struct.MockCircuit.html\" title=\"struct hyperplonk::prelude::MockCircuit\">MockCircuit</a>&lt;F&gt;","synthetic":true,"types":["hyperplonk::mock::MockCircuit"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"hyperplonk/prelude/struct.SelectorColumn.html\" title=\"struct hyperplonk::prelude::SelectorColumn\">SelectorColumn</a>&lt;F&gt;","synthetic":true,"types":["hyperplonk::selectors::SelectorColumn"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"hyperplonk/prelude/struct.WitnessColumn.html\" title=\"struct hyperplonk::prelude::WitnessColumn\">WitnessColumn</a>&lt;F&gt;","synthetic":true,"types":["hyperplonk::witness::WitnessColumn"]}];
implementors["pcs"] = [{"text":"impl Freeze for <a class=\"enum\" href=\"pcs/prelude/enum.PCSError.html\" title=\"enum pcs::prelude::PCSError\">PCSError</a>","synthetic":true,"types":["pcs::errors::PCSError"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearUniversalParams.html\" title=\"struct pcs::prelude::MultilinearUniversalParams\">MultilinearUniversalParams</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["pcs::multilinear_kzg::srs::MultilinearUniversalParams"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearProverParam.html\" title=\"struct pcs::prelude::MultilinearProverParam\">MultilinearProverParam</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["pcs::multilinear_kzg::srs::MultilinearProverParam"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearVerifierParam.html\" title=\"struct pcs::prelude::MultilinearVerifierParam\">MultilinearVerifierParam</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["pcs::multilinear_kzg::srs::MultilinearVerifierParam"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearKzgPCS.html\" title=\"struct pcs::prelude::MultilinearKzgPCS\">MultilinearKzgPCS</a>&lt;E&gt;","synthetic":true,"types":["pcs::multilinear_kzg::MultilinearKzgPCS"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearKzgProof.html\" title=\"struct pcs::prelude::MultilinearKzgProof\">MultilinearKzgProof</a>&lt;E&gt;","synthetic":true,"types":["pcs::multilinear_kzg::MultilinearKzgProof"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.MultilinearKzgBatchProof.html\" title=\"struct pcs::prelude::MultilinearKzgBatchProof\">MultilinearKzgBatchProof</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["pcs::multilinear_kzg::MultilinearKzgBatchProof"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.Commitment.html\" title=\"struct pcs::prelude::Commitment\">Commitment</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["pcs::structs::Commitment"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.UnivariateUniversalParams.html\" title=\"struct pcs::prelude::UnivariateUniversalParams\">UnivariateUniversalParams</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["pcs::univariate_kzg::srs::UnivariateUniversalParams"]},{"text":"impl&lt;C&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.UnivariateProverParam.html\" title=\"struct pcs::prelude::UnivariateProverParam\">UnivariateProverParam</a>&lt;C&gt;","synthetic":true,"types":["pcs::univariate_kzg::srs::UnivariateProverParam"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.UnivariateVerifierParam.html\" title=\"struct pcs::prelude::UnivariateVerifierParam\">UnivariateVerifierParam</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["pcs::univariate_kzg::srs::UnivariateVerifierParam"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.UnivariateKzgPCS.html\" title=\"struct pcs::prelude::UnivariateKzgPCS\">UnivariateKzgPCS</a>&lt;E&gt;","synthetic":true,"types":["pcs::univariate_kzg::UnivariateKzgPCS"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"pcs/prelude/struct.UnivariateKzgProof.html\" title=\"struct pcs::prelude::UnivariateKzgProof\">UnivariateKzgProof</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["pcs::univariate_kzg::UnivariateKzgProof"]}];
implementors["poly_iop"] = [{"text":"impl Freeze for <a class=\"enum\" href=\"poly_iop/prelude/enum.PolyIOPErrors.html\" title=\"enum poly_iop::prelude::PolyIOPErrors\">PolyIOPErrors</a>","synthetic":true,"types":["poly_iop::errors::PolyIOPErrors"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"poly_iop/struct.PolyIOP.html\" title=\"struct poly_iop::PolyIOP\">PolyIOP</a>&lt;F&gt;","synthetic":true,"types":["poly_iop::PolyIOP"]}];
implementors["transcript"] = [{"text":"impl Freeze for <a class=\"enum\" href=\"transcript/enum.TranscriptError.html\" title=\"enum transcript::TranscriptError\">TranscriptError</a>","synthetic":true,"types":["transcript::errors::TranscriptError"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"transcript/struct.IOPTranscript.html\" title=\"struct transcript::IOPTranscript\">IOPTranscript</a>&lt;F&gt;","synthetic":true,"types":["transcript::IOPTranscript"]}];
if (window.register_implementors) {window.register_implementors(implementors);} else {window.pending_implementors = implementors;}})()