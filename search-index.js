var searchIndex = new Map(JSON.parse('[\
["arithmetic",{"t":"GFPPPFFNNNNNNNONHNNNNNNNNHHNHNNNNNNNNNNNNNNNNNNNNNNNNNNHNNNNNNHHOHHHNONNNNNNNNNNNNHHHHNHHNNNNNNNNNNNNOHNNNNONONONNNHHHHNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN","n":["ArithErrors","DenseMultilinearExtension","InvalidParameters","SerializationErrors","ShouldNotArrive","VPAuxInfo","VirtualPolynomial","add","add","add","add_assign","add_assign","add_assign","add_mle_list","aux_info","batch_check","bit_decompose","borrow","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","borrow_mut","build_eq_x_r","build_eq_x_r_vec","build_f_hat","build_l","check","clone","clone","clone","clone_into","clone_into","clone_into","default","default","default","deref","deref","deref","deref","deref_mut","deref_mut","deref_mut","deref_mut","deserialize_with_mode","drop","drop","drop","drop","eq","eq","eq","eq_eval","equivalent","equivalent","equivalent","equivalent","evaluate","evaluate","evaluate_no_par","evaluate_opt","evaluations","fix_last_variables","fix_last_variables_no_par","fix_variables","fix_variables","flattened_ml_extensions","fmt","fmt","fmt","fmt","fmt","from","from","from","from","from","from_evaluations_slice","from_evaluations_vec","gen_eval_point","get_batched_nv","get_index","get_uni_domain","hash","identity_permutation","identity_permutation_mles","index","init","init","init","init","into","into","into","into","is_zero","iter","iter_mut","max_degree","merge_polynomials","mul_by_mle","neg","new","new_from_mle","num_variables","num_vars","num_vars","print_evals","products","rand","rand","rand_zero","random_mle_list","random_permutation","random_permutation_mles","random_zero_mle_list","relabel","relabel_in_place","serialize_with_mode","serialize_with_mode","serialized_size","serialized_size","sub","sub","sub_assign","sub_assign","to_evaluations","to_owned","to_owned","to_owned","to_string","try_from","try_from","try_from","try_from","try_into","try_into","try_into","try_into","type_id","type_id","type_id","type_id","vzip","vzip","vzip","vzip","zero"],"q":[[0,"arithmetic"],[151,"arithmetic::virtual_polynomial"],[152,"ark_ff::fields::prime"],[153,"ark_poly::evaluations::multivariate::multilinear::dense"],[154,"ark_ff::fields"],[155,"arithmetic::errors"],[156,"core::result"],[157,"alloc::sync"],[158,"core::iter::traits::collect"],[159,"ark_serialize::error"],[160,"core::iter::traits::iterator"],[161,"core::marker"],[162,"alloc::vec"],[163,"ark_poly::domain::radix2"],[164,"ark_poly::polynomial::univariate::dense"],[165,"core::clone"],[166,"core::default"],[167,"ark_serialize"],[168,"std::io"],[169,"core::cmp"],[170,"core::option"],[171,"core::fmt"],[172,"core::hash"],[173,"core::slice::iter"],[174,"rand_core"],[175,"rand::rng"],[176,"alloc::string"],[177,"core::any"],[178,"arithmetic::util"],[179,"arithmetic::univariate_polynomial"],[180,"arithmetic::multilinear_polynomial"]],"i":[0,0,9,9,9,0,0,1,4,4,4,4,4,1,1,4,0,9,1,25,4,9,1,25,4,0,0,1,0,4,1,25,4,1,25,4,1,25,4,9,1,25,4,9,1,25,4,4,9,1,25,4,1,25,4,0,25,25,4,4,1,4,0,0,4,0,0,0,4,1,9,9,1,25,4,9,9,1,25,4,4,4,0,0,0,0,4,0,0,4,9,1,25,4,9,1,25,4,4,4,4,25,0,1,4,1,1,25,4,4,1,1,1,4,1,0,0,0,0,4,4,25,4,25,4,4,4,4,4,4,1,25,4,9,9,1,25,4,9,1,25,4,9,1,25,4,9,1,25,4,4],"f":"```````{{{d{{b{c}}}}{d{{b{c}}}}}ef{}}{{{h{c}}{h{c}}}{{h{c}}}j}{{{d{{h{c}}}}{d{{h{c}}}}}{}j}{{{d{l{h{c}}}}{n{c{d{{h{c}}}}}}}A`j}{{{d{l{h{c}}}}{d{{h{c}}}}}A`j}{{{d{l{h{c}}}}{h{c}}}A`j}{{{d{l{b{c}}}}ec}{{Ad{A`Ab}}}f{{Aj{}{{Af{{Ah{{h{c}}}}}}}}}}`{e{{Ad{A`Al}}}j{{An{}{{Af{{d{{h{c}}}}}}}}B`}}{{BbBd}{{Bh{Bf}}}}{{{d{c}}}{{d{e}}}{}{}}000{{{d{lc}}}{{d{le}}}{}{}}000{{{d{{Bj{c}}}}}{{Ad{{Ah{{h{c}}}}Ab}}}f}{{{d{{Bj{c}}}}}{{Ad{{Bh{c}}Ab}}}f}{{{d{{b{c}}}}{d{{Bj{c}}}}}{{Ad{{b{c}}Ab}}}f}{{{d{{Bj{{Bh{c}}}}}}{d{{Bl{c}}}}Bf}{{Ad{{Bh{{Bn{c}}}}Ab}}}f}{{{d{{h{c}}}}}{{Ad{A`Al}}}j}{{{d{{b{c}}}}}{{b{c}}}{C`f}}{{{d{{Cb{c}}}}}{{Cb{c}}}{C`f}}{{{d{{h{c}}}}}{{h{c}}}{C`j}}{{{d{c}}{d{le}}}A`{}{}}00{{}{{b{c}}}{Cdf}}{{}{{Cb{c}}}{Cdf}}{{}{{h{c}}}{Cdj}}{Bd{{d{c}}}{}}000{Bd{{d{lc}}}{}}000{{cCfCh}{{Ad{{h{e}}Al}}}Cjj}{BdA`}000{{{d{{b{c}}}}{d{{b{c}}}}}Bf{Clf}}{{{d{{Cb{c}}}}{d{{Cb{c}}}}}Bf{Clf}}{{{d{{h{c}}}}{d{{h{c}}}}}Bf{Clj}}{{{d{{Bj{c}}}}{d{{Bj{c}}}}}{{Ad{cAb}}}f}{{{d{c}}{d{e}}}Bf{}{}}000{{{d{{b{c}}}}{d{{Bj{c}}}}}{{Ad{cAb}}}f}{{{d{{h{c}}}}{d{{Bj{c}}}}}{{Cn{c}}}j}{{{d{{h{c}}}}{d{{Bj{c}}}}}cj}0`{{{d{{h{c}}}}{d{{Bj{c}}}}}{{h{c}}}f}0{{{d{{h{c}}}}{d{{Bj{c}}}}}{{h{c}}}j}0`{{{d{Ab}}{d{lD`}}}Db}0{{{d{{b{c}}}}{d{lD`}}}Db{Ddf}}{{{d{{Cb{c}}}}{d{lD`}}}Db{Ddf}}{{{d{{h{c}}}}{d{lD`}}}{{Ad{A`Df}}}j}{AlAb}{cc{}}000{{Bd{d{{Bj{c}}}}}{{h{c}}}j}{{Bd{Bh{c}}}{{h{c}}}j}{{BdBd{d{{Bj{c}}}}}{{Bh{c}}}f}{{BdBd}Bd}{{BdBd}{{n{BdBdBf}}}}{Bd{{Ad{{Bl{c}}Ab}}}f}{{{d{{h{c}}}}{d{le}}}A`{Dhj}Dj}{{BdBd}{{Bh{c}}}f}{{BdBd}{{Bh{{Ah{{h{c}}}}}}}f}{{{d{{h{c}}}}Bd}dj}{{}Bd}000{ce{}{}}000{{{d{{h{c}}}}}Bfj}{{{d{{h{c}}}}}{{Dl{c}}}j}{{{d{l{h{c}}}}}{{Dn{c}}}j}`{{{d{{Bj{{Ah{{h{c}}}}}}}}}{{Ad{{Ah{{h{c}}}}Ab}}}f}{{{d{l{b{c}}}}{Ah{{h{c}}}}c}{{Ad{A`Ab}}}f}{{{h{c}}}{}j}{Bd{{b{c}}}f}{{{d{{Ah{{h{c}}}}}}c}{{b{c}}}f}`{{{d{{h{c}}}}}Bdj}`{{{d{{b{c}}}}}A`f}`{{Bd{n{BdBd}}Bd{d{lc}}}{{Ad{{n{{b{e}}e}}Ab}}}E`f}{{Bd{d{lc}}}{{h{e}}}Ebj}{{Bd{n{BdBd}}Bd{d{lc}}}{{Ad{{b{e}}Ab}}}E`f}{{BdBd{d{lc}}}{{n{{Bh{{Ah{{h{e}}}}}}e}}}E`f}{{BdBd{d{lc}}}{{Bh{e}}}E`f}{{BdBd{d{lc}}}{{Bh{{Ah{{h{e}}}}}}}E`f}0{{{d{{h{c}}}}BdBdBd}{{h{c}}}j}{{{d{l{h{c}}}}BdBdBd}A`j}{{{d{{Cb{c}}}}eCf}{{Ad{A`Al}}}fEd}{{{d{{h{c}}}}eCf}{{Ad{A`Al}}}jEd}{{{d{{Cb{c}}}}Cf}Bdf}{{{d{{h{c}}}}Cf}Bdj}{{{d{{h{c}}}}{d{{h{c}}}}}{}j}{{{h{c}}{h{c}}}{{h{c}}}j}{{{d{l{h{c}}}}{h{c}}}A`j}{{{d{l{h{c}}}}{d{{h{c}}}}}A`j}{{{d{{h{c}}}}}{{Bh{c}}}j}{{{d{c}}}e{}{}}00{{{d{c}}}Ef{}}{c{{Ad{e}}}{}{}}0000000{{{d{c}}}Eh{}}000{ce{}{}}000{{}{{h{c}}}j}","D":"F`","p":[[5,"VirtualPolynomial",0,151],[1,"reference"],[10,"PrimeField",152],[5,"DenseMultilinearExtension",0,153],[10,"Field",154],[0,"mut"],[1,"tuple"],[1,"unit"],[6,"ArithErrors",0,155],[6,"Result",156],[17,"Item"],[5,"Arc",157],[10,"IntoIterator",158],[6,"SerializationError",159],[10,"Iterator",160],[10,"Send",161],[1,"u64"],[1,"usize"],[1,"bool"],[5,"Vec",162],[1,"slice"],[5,"Radix2EvaluationDomain",163],[5,"DensePolynomial",164],[10,"Clone",165],[5,"VPAuxInfo",0,151],[10,"Default",166],[6,"Compress",167],[6,"Validate",167],[10,"Read",168],[10,"PartialEq",169],[6,"Option",170],[5,"Formatter",171],[8,"Result",171],[10,"Debug",171],[5,"Error",171],[10,"Hash",172],[10,"Hasher",172],[5,"Iter",173],[5,"IterMut",173],[10,"RngCore",174],[10,"Rng",175],[10,"Write",168],[5,"String",176],[5,"TypeId",177]],"r":[[0,155],[1,153],[5,151],[6,151],[16,178],[25,151],[26,151],[28,179],[55,151],[62,180],[63,180],[65,180],[66,180],[67,180],[82,178],[83,178],[84,178],[85,179],[87,180],[88,180],[102,180],[115,180],[116,180],[117,180],[118,180]],"b":[[8,"impl-Add-for-DenseMultilinearExtension%3CF%3E"],[9,"impl-Add%3C%26DenseMultilinearExtension%3CF%3E%3E-for-%26DenseMultilinearExtension%3CF%3E"],[10,"impl-AddAssign%3C(F,+%26DenseMultilinearExtension%3CF%3E)%3E-for-DenseMultilinearExtension%3CF%3E"],[11,"impl-AddAssign%3C%26DenseMultilinearExtension%3CF%3E%3E-for-DenseMultilinearExtension%3CF%3E"],[12,"impl-AddAssign-for-DenseMultilinearExtension%3CF%3E"],[70,"impl-Debug-for-ArithErrors"],[71,"impl-Display-for-ArithErrors"],[125,"impl-Sub%3C%26DenseMultilinearExtension%3CF%3E%3E-for-%26DenseMultilinearExtension%3CF%3E"],[126,"impl-Sub-for-DenseMultilinearExtension%3CF%3E"],[127,"impl-SubAssign-for-DenseMultilinearExtension%3CF%3E"],[128,"impl-SubAssign%3C%26DenseMultilinearExtension%3CF%3E%3E-for-DenseMultilinearExtension%3CF%3E"]],"c":"OjAAAAAAAAA=","e":"OzAAAAEAAGcAFAAAAAAACAAFABAAAAASAAcAHAAAAB4AGQA5AAMAPgACAEIAAwBHAAUAVQAAAFcAAQBbAAMAYwAAAGkAAABtAAAAcgAAAHUAAAB3AAEAegAdAA=="}],\
["hyperplonk",{"t":"KRRRRQCMMMPFGEPPPPFPPFPPFNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNONNNNNNNNNNNNNNNNNNNONNNNNNNNNNNNNNNNNNNNNNNNNNO","n":["HyperPlonkSNARK","Index","Proof","ProvingKey","VerifyingKey","build_mle","prelude","preprocess","prove","verify","ArithmeticErrors","CustomizedGates","HyperPlonkErrors","HyperPlonkSNARK","InvalidParameters","InvalidProof","InvalidProver","InvalidVerifier","MockCircuit","PCSErrors","PolyIOPErrors","SelectorColumn","SerializationError","TranscriptError","WitnessColumn","append","append","borrow","borrow","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","clone","clone","clone","clone_into","clone_into","clone_into","coeff_ref","default","default","default","degree","deref","deref","deref","deref","deref","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","drop","drop","drop","drop","drop","eq","eq","equivalent","equivalent","equivalent","equivalent","fmt","fmt","fmt","fmt","fmt","from","from","from","from","from","from","from","from","from","from","from_selector_rows","from_witness_rows","get_nv","get_nv","index","init","init","init","init","init","into","into","into","into","into","is_satisfied","jellyfish_turbo_plonk_gate","mock_gate","new","num_selector_columns","num_selector_columns","num_variables","num_witness_columns","num_witness_columns","public_inputs","super_long_selector_gate","to_owned","to_owned","to_owned","to_string","try_from","try_from","try_from","try_from","try_from","try_into","try_into","try_into","try_into","try_into","type_id","type_id","type_id","type_id","type_id","vanilla_plonk_gate","vzip","vzip","vzip","vzip","vzip","witnesses"],"q":[[0,"hyperplonk"],[10,"hyperplonk::prelude"],[136,"hyperplonk::errors"],[137,"core::result"],[138,"hyperplonk::witness"],[139,"hyperplonk::selectors"],[140,"ark_ff::fields::prime"],[141,"hyperplonk::custom_gate"],[142,"core::clone"],[143,"core::default"],[144,"core::cmp"],[145,"core::fmt"],[146,"ark_serialize::error"],[147,"subroutines::pcs::errors"],[148,"arithmetic::errors"],[149,"subroutines::poly_iop::errors"],[150,"transcript::errors"],[151,"alloc::vec"],[152,"hyperplonk::mock"],[153,"alloc::string"],[154,"core::any"]],"i":[0,29,29,29,29,0,0,29,29,29,3,0,0,0,3,3,3,3,0,3,3,0,3,3,0,9,6,26,12,3,9,6,26,12,3,9,6,12,9,6,12,9,6,6,12,9,6,12,26,12,3,9,6,26,12,3,9,6,26,12,3,9,6,12,9,12,12,9,9,12,3,3,9,6,26,12,3,3,3,3,3,3,9,6,9,6,9,6,26,26,12,3,9,6,26,12,3,9,6,26,12,12,26,26,12,26,26,12,26,12,12,9,6,3,26,12,3,9,6,26,12,3,9,6,26,12,3,9,6,12,26,12,3,9,6,26],"f":"```````{{{b{c}}b}{{h{{d{eg}}f}}}{}{}{}}{{{b{c}}{b{j}}{b{{j{l}}}}}{{h{ef}}}{}{}}{{{b{c}}{b{j}}{b{e}}}{{h{nf}}}{}{}}```````````````{{{b{A`{Ab{c}}}}c}AdAf}{{{b{A`{l{c}}}}c}AdAf}{{{b{c}}}{{b{e}}}{}{}}0000{{{b{A`c}}}{{b{A`e}}}{}{}}0000{{{b{Ah}}}Ah}{{{b{{Ab{c}}}}}{{Ab{c}}}{AjAf}}{{{b{{l{c}}}}}{{l{c}}}{AjAf}}{{{b{c}}{b{A`e}}}Ad{}{}}00{{{b{{l{c}}}}}{{b{{j{c}}}}}Af}{{}Ah}{{}{{Ab{c}}}{AlAf}}{{}{{l{c}}}{AlAf}}{{{b{Ah}}}An}{An{{b{c}}}{}}0000{An{{b{A`c}}}{}}0000{AnAd}0000{{{b{Ah}}{b{Ah}}}n}{{{b{{Ab{c}}}}{b{{Ab{c}}}}}n{B`Af}}{{{b{c}}{b{e}}}n{}{}}000{{{b{Ah}}{b{A`Bb}}}Bd}{{{b{f}}{b{A`Bb}}}Bd}0{{{b{{Ab{c}}}}{b{A`Bb}}}Bd{BfAf}}{{{b{{l{c}}}}{b{A`Bb}}}Bd{BfAf}}{cc{}}0{Bhf}1{Bjf}{Blf}{Bnf}{C`f}55{{{b{{j{{`{c}}}}}}}{{h{{Cb{{Ab{c}}}}f}}}Af}{{{b{{j{{`{c}}}}}}}{{h{{Cb{{l{c}}}}f}}}Af}{{{b{{Ab{c}}}}}AnAf}{{{b{{l{c}}}}}AnAf}`{{}An}0000{ce{}{}}0000{{{b{{Cd{c}}}}}nAf}{{}Ah}{{AnAn}Ah}{{An{b{Ah}}}{{Cd{c}}}Af}{{{b{{Cd{c}}}}}AnAf}{{{b{Ah}}}An}110`4{{{b{c}}}e{}{}}00{{{b{c}}}Cf{}}{c{{h{e}}}{}{}}000000000{{{b{c}}}Ch{}}00008:::::`","D":"El","p":[[1,"reference"],[1,"tuple"],[6,"HyperPlonkErrors",10,136],[6,"Result",137],[1,"slice"],[5,"WitnessColumn",10,138],[1,"bool"],[0,"mut"],[5,"SelectorColumn",10,139],[1,"unit"],[10,"PrimeField",140],[5,"CustomizedGates",10,141],[10,"Clone",142],[10,"Default",143],[1,"usize"],[10,"PartialEq",144],[5,"Formatter",145],[8,"Result",145],[10,"Debug",145],[6,"SerializationError",146],[6,"PCSError",147],[6,"ArithErrors",148],[6,"PolyIOPErrors",149],[6,"TranscriptError",150],[5,"Vec",151],[5,"MockCircuit",10,152],[5,"String",153],[5,"TypeId",154],[10,"HyperPlonkSNARK",0]],"r":[[11,141],[12,136],[13,0],[18,152],[21,139],[24,138]],"b":[[70,"impl-Display-for-HyperPlonkErrors"],[71,"impl-Debug-for-HyperPlonkErrors"],[76,"impl-From%3CSerializationError%3E-for-HyperPlonkErrors"],[78,"impl-From%3CPCSError%3E-for-HyperPlonkErrors"],[79,"impl-From%3CArithErrors%3E-for-HyperPlonkErrors"],[80,"impl-From%3CPolyIOPErrors%3E-for-HyperPlonkErrors"],[81,"impl-From%3CTranscriptError%3E-for-HyperPlonkErrors"]],"c":"OjAAAAAAAAA=","e":"OzAAAAEAAFoADQACAAMABwAAAA4AAAATAAAAHAATADEAGQBNAAAATwADAFkABQBkAAAAbQAAAG8AEgCDAAUA"}],\
["subroutines",{"t":"CCQRRRRRKRRRRKRRNMMMMMNMCMMMPFFPPPPFFFFFGEPEPIFFFFFNNNNNNNNNNOONNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNONNNNNNNNNNNNNNNNNNNNNNNNNNNNOOONNNNOOOOOONNNNNNNNNNNNNNNNNNNNNNNNNNNNNOONNOOOOOONNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNFNNNNNNNNNNNNNNNNNNNNCNNNNNNNNNNNNNPFPPPPPRPKRREGKRRPPKRRRPRRKRRNNNNNNNNNNNNNNNNMNNNNNNNNNNNMMMMNNOOMMMMQNNNNNNNNMMMMNN","n":["pcs","poly_iop","to_bytes","BatchProof","Commitment","Evaluation","Point","Polynomial","PolynomialCommitmentScheme","Proof","ProverParam","ProverParam","SRS","StructuredReferenceString","VerifierParam","VerifierParam","batch_verify","commit","extract_prover_param","extract_verifier_param","gen_srs_for_testing","gen_srs_for_testing","multi_open","open","prelude","trim","trim","verify","ArithErrors","BatchProof","Commitment","InvalidParameters","InvalidProof","InvalidProver","InvalidVerifier","MultilinearKzgPCS","MultilinearKzgProof","MultilinearProverParam","MultilinearUniversalParams","MultilinearVerifierParam","PCSError","PolynomialCommitmentScheme","SerializationError","StructuredReferenceString","TranscriptError","UnivariateKzgBatchProof","UnivariateKzgPCS","UnivariateKzgProof","UnivariateProverParam","UnivariateUniversalParams","UnivariateVerifierParam","batch_check","batch_check","batch_check","batch_check","batch_check","batch_check","batch_check","batch_check","batch_check","batch_verify","beta_h","beta_h","borrow","borrow","borrow","borrow","borrow","borrow","borrow","borrow","borrow","borrow","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","check","check","check","check","check","check","check","check","check","clone","clone","clone","clone","clone","clone","clone","clone","clone","clone","clone_into","clone_into","clone_into","clone_into","clone_into","clone_into","clone_into","clone_into","clone_into","clone_into","commit","commit","default","default","default","default","default","deref","deref","deref","deref","deref","deref","deref","deref","deref","deref","deref","deref","deref","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deref_mut","deserialize_with_mode","deserialize_with_mode","deserialize_with_mode","deserialize_with_mode","deserialize_with_mode","deserialize_with_mode","deserialize_with_mode","deserialize_with_mode","deserialize_with_mode","drop","drop","drop","drop","drop","drop","drop","drop","drop","drop","drop","drop","drop","eq","eq","eq","eq","eq","eq","eq","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","equivalent","extract_prover_param","extract_prover_param","extract_verifier_param","extract_verifier_param","f_i_eval_at_point_i","fmt","fmt","fmt","fmt","fmt","fmt","fmt","fmt","fmt","fmt","fmt","fmt","from","from","from","from","from","from","from","from","from","from","from","from","from","from","from","from","g","g","g","gen_srs_for_testing","gen_srs_for_testing","gen_srs_for_testing","gen_srs_for_testing","h","h","h","h","h_mask","h_mask","hash","init","init","init","init","init","init","init","init","init","init","init","init","init","into","into","into","into","into","into","into","into","into","into","into","into","into","max_degree","multi_open","num_vars","num_vars","open","open","powers_of_g","powers_of_g","powers_of_g","proof","proofs","prover_param","serialize_with_mode","serialize_with_mode","serialize_with_mode","serialize_with_mode","serialize_with_mode","serialize_with_mode","serialize_with_mode","serialize_with_mode","serialize_with_mode","serialized_size","serialized_size","serialized_size","serialized_size","serialized_size","serialized_size","serialized_size","serialized_size","serialized_size","to_owned","to_owned","to_owned","to_owned","to_owned","to_owned","to_owned","to_owned","to_owned","to_owned","to_string","trim","trim","trim","trim","try_from","try_from","try_from","try_from","try_from","try_from","try_from","try_from","try_from","try_from","try_from","try_from","try_from","try_into","try_into","try_into","try_into","try_into","try_into","try_into","try_into","try_into","try_into","try_into","try_into","try_into","type_id","type_id","type_id","type_id","type_id","type_id","type_id","type_id","type_id","type_id","type_id","type_id","type_id","verify","verify","vzip","vzip","vzip","vzip","vzip","vzip","vzip","vzip","vzip","vzip","vzip","vzip","vzip","PolyIOP","borrow","borrow_mut","clone","clone_into","default","deref","deref_mut","drop","eq","equivalent","equivalent","extract_sum","fmt","from","init","init_transcript","init_transcript","init_transcript","init_transcript","into","prelude","prove","prove","prove","prove","to_owned","try_from","try_into","type_id","verify","verify","verify","verify","vzip","ArithmeticErrors","IOPProof","InvalidChallenge","InvalidParameters","InvalidProof","InvalidProver","InvalidVerifier","MultilinearExtension","PCSErrors","PermutationCheck","PermutationCheckSubClaim","PermutationProof","PolyIOP","PolyIOPErrors","ProductCheck","ProductCheckProof","ProductCheckSubClaim","SerializationErrors","ShouldNotArrive","SumCheck","SumCheckProof","SumCheckSubClaim","Transcript","TranscriptErrors","VPAuxInfo","VirtualPolynomial","ZeroCheck","ZeroCheckProof","ZeroCheckSubClaim","borrow","borrow","borrow_mut","borrow_mut","clone","clone_into","default","deref","deref","deref_mut","deref_mut","drop","drop","eq","equivalent","equivalent","extract_sum","fmt","fmt","fmt","from","from","from","from","from","from","init","init","init_transcript","init_transcript","init_transcript","init_transcript","into","into","point","proofs","prove","prove","prove","prove","to_bytes","to_owned","to_string","try_from","try_from","try_into","try_into","type_id","type_id","verify","verify","verify","verify","vzip","vzip"],"q":[[0,"subroutines"],[3,"subroutines::pcs"],[28,"subroutines::pcs::prelude"],[366,"subroutines::poly_iop"],[401,"subroutines::poly_iop::prelude"],[485,"transcript"],[486,"subroutines::pcs::errors"],[487,"core::result"],[488,"core::borrow"],[489,"rand::rng"],[490,"core::option"],[491,"ark_serialize::error"],[492,"ark_ec::pairing"],[493,"subroutines::pcs::multilinear_kzg::srs"],[494,"core::iter::traits::iterator"],[495,"core::marker"],[496,"subroutines::pcs::multilinear_kzg"],[497,"subroutines::pcs::structs"],[498,"subroutines::pcs::univariate_kzg::srs"],[499,"ark_ec"],[500,"subroutines::pcs::univariate_kzg"],[501,"subroutines::pcs::multilinear_kzg::batching"],[502,"core::clone"],[503,"core::default"],[504,"ark_serialize"],[505,"std::io"],[506,"core::cmp"],[507,"core::fmt"],[508,"arithmetic::errors"],[509,"transcript::errors"],[510,"core::hash"],[511,"alloc::string"],[512,"core::any"],[513,"ark_ff::fields::prime"],[514,"subroutines::poly_iop::errors"],[515,"arithmetic::virtual_polynomial"],[516,"subroutines::poly_iop::structs"],[517,"subroutines::poly_iop::perm_check"],[518,"subroutines::poly_iop::prod_check"],[519,"subroutines::poly_iop::sum_check"],[520,"subroutines::poly_iop::zero_check"]],"i":[0,0,0,34,34,34,34,34,0,34,34,11,34,0,34,11,34,34,11,11,34,11,34,34,0,34,11,34,6,0,0,6,6,6,6,0,0,0,0,0,0,0,6,0,6,0,0,0,0,0,0,20,23,24,25,26,27,29,30,31,46,27,30,46,55,6,32,20,23,24,25,26,27,29,30,31,46,55,6,32,20,23,24,25,26,27,29,30,31,20,23,24,25,26,27,29,30,31,32,20,23,24,25,26,27,29,30,31,32,20,23,24,25,26,27,29,30,31,46,55,32,26,27,29,30,46,55,6,32,20,23,24,25,26,27,29,30,31,46,55,6,32,20,23,24,25,26,27,29,30,31,20,23,24,25,26,27,29,30,31,46,55,6,32,20,23,24,25,26,27,29,30,31,32,25,26,27,29,30,31,32,32,25,25,26,26,27,27,29,29,30,30,31,31,20,27,20,27,32,6,6,32,20,23,24,25,26,27,29,30,31,46,55,6,6,6,6,32,20,23,24,25,26,27,29,30,31,23,24,30,46,55,20,27,23,24,27,30,20,24,26,46,55,6,32,20,23,24,25,26,27,29,30,31,46,55,6,32,20,23,24,25,26,27,29,30,31,27,46,23,24,46,55,23,27,29,31,25,20,20,23,24,25,26,27,29,30,31,20,23,24,25,26,27,29,30,31,32,20,23,24,25,26,27,29,30,31,6,46,55,20,27,46,55,6,32,20,23,24,25,26,27,29,30,31,46,55,6,32,20,23,24,25,26,27,29,30,31,46,55,6,32,20,23,24,25,26,27,29,30,31,46,55,46,55,6,32,20,23,24,25,26,27,29,30,31,0,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,0,50,50,50,50,50,50,50,50,50,50,50,50,50,52,0,52,52,52,52,52,56,52,0,57,57,0,0,0,58,58,52,52,0,56,56,56,52,56,56,0,59,59,52,54,52,54,54,54,54,52,54,52,54,52,54,54,54,54,56,52,52,54,52,52,52,52,52,54,52,54,57,58,56,59,52,54,54,54,57,58,56,59,0,54,52,52,54,52,54,52,54,57,58,56,59,52,54],"f":"````````````````{{{b{c}}{b{{d{e}}}}{b{{d{g}}}}{b{i}}{b{fh}}}{{n{jl}}}{}{}{}{}}{{e{b{g}}}{{n{il}}}{}{{A`{c}}}{}{}}{{{b{{Af{}{{Ab{c}}{Ad{e}}}}}}Ah}c{}{}}{{{b{{Af{}{{Ab{c}}{Ad{e}}}}}}Ah}e{}{}}{{{b{fc}}Ah}{{n{el}}}Aj{}}{{{b{fc}}Ah}{{n{{Af{}{{Ab{e}}{Ad{g}}}}l}}}Aj{}{}}{{e{b{{d{g}}}}{b{{d{i}}}}{b{{d{k}}}}{b{fh}}}{{n{ml}}}{}{{A`{c}}}{}{}{}{}}{{e{b{g}}{b{i}}}{{n{{Al{km}}l}}}{}{{A`{c}}}{}{}{}{}}`{{e{An{Ah}}{An{Ah}}}{{n{{Al{gi}}l}}}{}{{A`{c}}}{}{}}{{{b{{Af{}{{Ab{c}}{Ad{e}}}}}}Ah}{{n{{Al{ce}}l}}}{}{}}{{{b{c}}{b{e}}{b{g}}b{b{i}}}{{n{jl}}}{}{}{}{}}```````````````````````{e{{n{B`Bb}}}Bd{{Bj{}{{Bf{{b{{Bh{c}}}}}}}}Bl}}{e{{n{B`Bb}}}Bd{{Bj{}{{Bf{{b{{Bn{c}}}}}}}}Bl}}{e{{n{B`Bb}}}Bd{{Bj{}{{Bf{{b{{C`{c}}}}}}}}Bl}}{e{{n{B`Bb}}}Bd{{Bj{}{{Bf{{b{{Cb{c}}}}}}}}Bl}}{e{{n{B`Bb}}}Bd{{Bj{}{{Bf{{b{{Cd{c}}}}}}}}Bl}}{e{{n{B`Bb}}}Bd{{Bj{}{{Bf{{b{{Cf{c}}}}}}}}Bl}}{e{{n{B`Bb}}}Ch{{Bj{}{{Bf{{b{{Cj{c}}}}}}}}Bl}}{e{{n{B`Bb}}}Bd{{Bj{}{{Bf{{b{{Cl{c}}}}}}}}Bl}}{e{{n{B`Bb}}}Bd{{Bj{}{{Bf{{b{{Cn{c}}}}}}}}Bl}}{{{b{c}}{b{{d{e}}}}{b{{d{g}}}}{b{i}}{b{fh}}}{{n{jl}}}{}{}{}{}}``{{{b{c}}}{{b{e}}}{}{}}000000000000{{{b{fc}}}{{b{fe}}}{}{}}000000000000{{{b{{Bh{c}}}}}{{n{B`Bb}}}Bd}{{{b{{Bn{c}}}}}{{n{B`Bb}}}Bd}{{{b{{C`{c}}}}}{{n{B`Bb}}}Bd}{{{b{{Cb{c}}}}}{{n{B`Bb}}}Bd}{{{b{{Cd{c}}}}}{{n{B`Bb}}}Bd}{{{b{{Cf{c}}}}}{{n{B`Bb}}}Bd}{{{b{{Cj{c}}}}}{{n{B`Bb}}}Ch}{{{b{{Cl{c}}}}}{{n{B`Bb}}}Bd}{{{b{{Cn{c}}}}}{{n{B`Bb}}}Bd}{{{b{{D`{ce}}}}}{{D`{ce}}}{BdDb}{{Dd{c}}Db}}{{{b{{Bh{c}}}}}{{Bh{c}}}{DbBd}}{{{b{{Bn{c}}}}}{{Bn{c}}}{DbBd}}{{{b{{C`{c}}}}}{{C`{c}}}{DbBd}}{{{b{{Cb{c}}}}}{{Cb{c}}}{DbBd}}{{{b{{Cd{c}}}}}{{Cd{c}}}Bd}{{{b{{Cf{c}}}}}{{Cf{c}}}{DbBd}}{{{b{{Cj{c}}}}}{{Cj{c}}}{DbCh}}{{{b{{Cl{c}}}}}{{Cl{c}}}Bd}{{{b{{Cn{c}}}}}{{Cn{c}}}{DbBd}}{{{b{c}}{b{fe}}}B`{}{}}000000000{{e{b{g}}}{{n{il}}}{}{{A`{c}}}{}{}}0{{}{{D`{ce}}}{BdDf}{{Dd{c}}Df}}{{}{{Cd{c}}}Bd}{{}{{Cf{c}}}{DfBd}}{{}{{Cj{c}}}{DfCh}}{{}{{Cl{c}}}Bd}{Ah{{b{c}}}{}}000000000000{Ah{{b{fc}}}{}}000000000000{{cDhDj}{{n{{Bh{e}}Bb}}}DlBd}{{cDhDj}{{n{{Bn{e}}Bb}}}DlBd}{{cDhDj}{{n{{C`{e}}Bb}}}DlBd}{{cDhDj}{{n{{Cb{e}}Bb}}}DlBd}{{cDhDj}{{n{{Cd{e}}Bb}}}DlBd}{{cDhDj}{{n{{Cf{e}}Bb}}}DlBd}{{cDhDj}{{n{{Cj{e}}Bb}}}DlCh}{{cDhDj}{{n{{Cl{e}}Bb}}}DlBd}{{cDhDj}{{n{{Cn{e}}Bb}}}DlBd}{AhB`}000000000000{{{b{{D`{ce}}}}{b{{D`{ce}}}}}j{BdDn}{{Dd{c}}Dn}}{{{b{{Cb{c}}}}{b{{Cb{c}}}}}j{DnBd}}{{{b{{Cd{c}}}}{b{{Cd{c}}}}}jBd}{{{b{{Cf{c}}}}{b{{Cf{c}}}}}j{DnBd}}{{{b{{Cj{c}}}}{b{{Cj{c}}}}}j{DnCh}}{{{b{{Cl{c}}}}{b{{Cl{c}}}}}jBd}{{{b{{Cn{c}}}}{b{{Cn{c}}}}}j{DnBd}}{{{b{c}}{b{e}}}j{}{}}0000000000000{{{b{{Bh{c}}}}Ah}eBd{}}{{{b{{Cf{c}}}}Ah}eBd{}}10`{{{b{l}}{b{fE`}}}Eb}0{{{b{{D`{ce}}}}{b{fE`}}}Eb{BdEd}{{Dd{c}}Ed}}{{{b{{Bh{c}}}}{b{fE`}}}Eb{EdBd}}{{{b{{Bn{c}}}}{b{fE`}}}Eb{EdBd}}{{{b{{C`{c}}}}{b{fE`}}}Eb{EdBd}}{{{b{{Cb{c}}}}{b{fE`}}}Eb{EdBd}}{{{b{{Cd{c}}}}{b{fE`}}}EbBd}{{{b{{Cf{c}}}}{b{fE`}}}Eb{EdBd}}{{{b{{Cj{c}}}}{b{fE`}}}Eb{EdCh}}{{{b{{Cl{c}}}}{b{fE`}}}EbBd}{{{b{{Cn{c}}}}{b{fE`}}}Eb{EdBd}}{cc{}}0{Efl}{Ehl}2{Bbl}3333333333```{{{b{fc}}Ah}{{n{el}}}Aj{}}0{{{b{fc}}Ah}{{n{{Bh{e}}l}}}AjBd}{{{b{fc}}Ah}{{n{{Cf{e}}l}}}AjBd}``````{{{b{{Cd{c}}}}{b{fe}}}B`BdEj}{{}Ah}000000000000{ce{}{}}000000000000{{{b{{Cf{c}}}}}AhBd}{{e{b{{d{g}}}}{b{{d{i}}}}{b{{d{k}}}}{b{fh}}}{{n{{D`{m{El{m}}}}l}}}{}{{A`{c}}}{}{}{}Bd}``{{e{b{g}}{b{i}}}{{n{{Al{km}}l}}}{}{{A`{c}}}{}{}{}{}}0``````{{{b{{Bh{c}}}}eDh}{{n{B`Bb}}}BdEn}{{{b{{Bn{c}}}}eDh}{{n{B`Bb}}}BdEn}{{{b{{C`{c}}}}eDh}{{n{B`Bb}}}BdEn}{{{b{{Cb{c}}}}eDh}{{n{B`Bb}}}BdEn}{{{b{{Cd{c}}}}eDh}{{n{B`Bb}}}BdEn}{{{b{{Cf{c}}}}eDh}{{n{B`Bb}}}BdEn}{{{b{{Cj{c}}}}eDh}{{n{B`Bb}}}ChEn}{{{b{{Cl{c}}}}eDh}{{n{B`Bb}}}BdEn}{{{b{{Cn{c}}}}eDh}{{n{B`Bb}}}BdEn}{{{b{{Bh{c}}}}Dh}AhBd}{{{b{{Bn{c}}}}Dh}AhBd}{{{b{{C`{c}}}}Dh}AhBd}{{{b{{Cb{c}}}}Dh}AhBd}{{{b{{Cd{c}}}}Dh}AhBd}{{{b{{Cf{c}}}}Dh}AhBd}{{{b{{Cj{c}}}}Dh}AhCh}{{{b{{Cl{c}}}}Dh}AhBd}{{{b{{Cn{c}}}}Dh}AhBd}{{{b{c}}}e{}{}}000000000{{{b{c}}}F`{}}{{e{An{Ah}}{An{Ah}}}{{n{{Al{gi}}l}}}{}{{A`{c}}}{}{}}0{{{b{{Bh{c}}}}Ah}{{n{{Al{eg}}l}}}Bd{}{}}{{{b{{Cf{c}}}}Ah}{{n{{Al{eg}}l}}}Bd{}{}}{c{{n{e}}}{}{}}0000000000000000000000000{{{b{c}}}Fb{}}000000000000{{{b{c}}{b{e}}{b{g}}b{b{i}}}{{n{jl}}}{}{}{}{}}0{ce{}{}}000000000000`{{{b{c}}}{{b{e}}}{}{}}{{{b{fc}}}{{b{fe}}}{}{}}{{{b{{Fd{c}}}}}{{Fd{c}}}{DbFf}}{{{b{c}}{b{fe}}}B`{}{}}{{}{{Fd{c}}}{DfFf}}{Ah{{b{c}}}{}}{Ah{{b{fc}}}{}}{AhB`}{{{b{{Fd{c}}}}{b{{Fd{c}}}}}j{DnFf}}{{{b{c}}{b{e}}}j{}{}}0{{{b{c}}}e{}Ff}{{{b{{Fd{c}}}}{b{fE`}}}Eb{EdFf}}{cc{}}{{}Ah}{{}c{}}000?`{{{b{c}}{b{fe}}}{{n{gFh}}}{}{}{}}0{{b{b{{d{c}}}}{b{{d{c}}}}{b{{d{c}}}}{b{fh}}}{{n{{Al{ecc}}Fh}}}{}{}}{{b{b{{d{c}}}}{b{{d{c}}}}{b{fh}}}{{n{{Al{ecc}}Fh}}}{}{}}{{{b{c}}}e{}{}}{c{{n{e}}}{}{}}0{{{b{c}}}Fb{}}{{{b{c}}{b{e}}{b{fg}}}{{n{iFh}}}{}{}{}{}}{{c{b{e}}{b{g}}{b{fi}}}{{n{kFh}}}Ff{}{}{}{}}{{{b{c}}{b{Fj}}{b{fe}}}{{n{gFh}}}{}{}{}}2{ce{}{}}`````````````````````````````{{{b{c}}}{{b{e}}}{}{}}0{{{b{fc}}}{{b{fe}}}{}{}}0{{{b{{Fl{c}}}}}{{Fl{c}}}{DbFf}}{{{b{c}}{b{fe}}}B`{}{}}{{}{{Fl{c}}}{DfFf}}{Ah{{b{c}}}{}}0{Ah{{b{fc}}}{}}0{AhB`}0{{{b{{Fl{c}}}}{b{{Fl{c}}}}}j{DnFf}}{{{b{c}}{b{e}}}j{}{}}0{{{b{c}}}e{}Ff}{{{b{Fh}}{b{fE`}}}Eb}0{{{b{{Fl{c}}}}{b{fE`}}}Eb{EdFf}}{cc{}}{lFh}{EfFh}{EhFh}{BbFh}4{{}Ah}0{{}c{}}000{ce{}{}}0``{{b{b{{d{c}}}}{b{{d{c}}}}{b{{d{c}}}}{b{fh}}}{{n{{Al{ecc}}Fh}}}{}{}}{{b{b{{d{c}}}}{b{{d{c}}}}{b{fh}}}{{n{{Al{ecc}}Fh}}}{}{}}{{{b{c}}{b{fe}}}{{n{gFh}}}{}{}{}}0`{{{b{c}}}e{}{}}{{{b{c}}}F`{}}{c{{n{e}}}{}{}}000{{{b{c}}}Fb{}}0{{{b{c}}{b{e}}{b{fg}}}{{n{iFh}}}{}{}{}{}}{{{b{c}}{b{Fj}}{b{fe}}}{{n{gFh}}}{}{}{}}{{c{b{e}}{b{g}}{b{fi}}}{{n{kFh}}}Ff{}{}{}{}}2::","D":"ABl","p":[[1,"reference"],[1,"slice"],[0,"mut"],[5,"IOPTranscript",485],[1,"bool"],[6,"PCSError",28,486],[6,"Result",487],[10,"Borrow",488],[17,"ProverParam"],[17,"VerifierParam"],[10,"StructuredReferenceString",3],[1,"usize"],[10,"Rng",489],[1,"tuple"],[6,"Option",490],[1,"unit"],[6,"SerializationError",491],[10,"Pairing",492],[17,"Item"],[5,"MultilinearUniversalParams",28,493],[10,"Iterator",494],[10,"Send",495],[5,"MultilinearProverParam",28,493],[5,"MultilinearVerifierParam",28,493],[5,"MultilinearKzgProof",28,496],[5,"Commitment",28,497],[5,"UnivariateUniversalParams",28,498],[10,"AffineRepr",499],[5,"UnivariateProverParam",28,498],[5,"UnivariateVerifierParam",28,498],[5,"UnivariateKzgProof",28,500],[5,"BatchProof",28,501],[10,"Clone",502],[10,"PolynomialCommitmentScheme",3],[10,"Default",503],[6,"Compress",504],[6,"Validate",504],[10,"Read",505],[10,"PartialEq",506],[5,"Formatter",507],[8,"Result",507],[10,"Debug",507],[6,"ArithErrors",508],[6,"TranscriptError",509],[10,"Hasher",510],[5,"MultilinearKzgPCS",28,496],[10,"Write",505],[5,"String",511],[5,"TypeId",512],[5,"PolyIOP",366],[10,"PrimeField",513],[6,"PolyIOPErrors",401,514],[5,"VPAuxInfo",515],[5,"IOPProof",401,516],[5,"UnivariateKzgPCS",28],[10,"SumCheck",401],[10,"PermutationCheck",401],[10,"ProductCheck",401],[10,"ZeroCheck",401]],"r":[[29,501],[30,497],[35,496],[36,496],[37,493],[38,493],[39,493],[40,486],[41,3],[43,3],[45,500],[46,500],[47,500],[48,498],[49,498],[50,498],[402,516],[410,517],[413,366],[414,514],[415,518],[420,519],[427,520],[470,0]],"b":[[199,"impl-Debug-for-PCSError"],[200,"impl-Display-for-PCSError"],[213,"impl-From%3CArithErrors%3E-for-PCSError"],[214,"impl-From%3CTranscriptError%3E-for-PCSError"],[216,"impl-From%3CSerializationError%3E-for-PCSError"],[382,"impl-SumCheck%3CF%3E-for-PolyIOP%3CF%3E"],[383,"impl-ProductCheck%3CE,+PCS%3E-for-PolyIOP%3C%3CE+as+Pairing%3E::ScalarField%3E"],[384,"impl-PermutationCheck%3CE,+PCS%3E-for-PolyIOP%3C%3CE+as+Pairing%3E::ScalarField%3E"],[385,"impl-ZeroCheck%3CF%3E-for-PolyIOP%3CF%3E"],[388,"impl-SumCheck%3CF%3E-for-PolyIOP%3CF%3E"],[389,"impl-ZeroCheck%3CF%3E-for-PolyIOP%3CF%3E"],[390,"impl-PermutationCheck%3CE,+PCS%3E-for-PolyIOP%3C%3CE+as+Pairing%3E::ScalarField%3E"],[391,"impl-ProductCheck%3CE,+PCS%3E-for-PolyIOP%3C%3CE+as+Pairing%3E::ScalarField%3E"],[396,"impl-ZeroCheck%3CF%3E-for-PolyIOP%3CF%3E"],[397,"impl-SumCheck%3CF%3E-for-PolyIOP%3CF%3E"],[398,"impl-ProductCheck%3CE,+PCS%3E-for-PolyIOP%3C%3CE+as+Pairing%3E::ScalarField%3E"],[399,"impl-PermutationCheck%3CE,+PCS%3E-for-PolyIOP%3C%3CE+as+Pairing%3E::ScalarField%3E"],[447,"impl-Display-for-PolyIOPErrors"],[448,"impl-Debug-for-PolyIOPErrors"],[451,"impl-From%3CPCSError%3E-for-PolyIOPErrors"],[452,"impl-From%3CArithErrors%3E-for-PolyIOPErrors"],[453,"impl-From%3CTranscriptError%3E-for-PolyIOPErrors"],[454,"impl-From%3CSerializationError%3E-for-PolyIOPErrors"]],"c":"OjAAAAAAAAA=","e":"OzAAAAEAAE8BHQAAAAIAHgAAACoAAAAsAAAANAAIAEAANgB5AEkAyAALANYAAQDZAAAA8QANABgBHAA5ASYAYgEMAHABDAB+AQQAhAENAJkBAACcAQIAoQEBAKYBAgCqAQEArQERAMABAgDEAQMAyQEBANEBAQDYAQcA5AEBAA=="}],\
["transcript",{"t":"FPPGNNNNNNNNNNNNNNNNNNNNNNNNNNNQNNNNNNNNNN","n":["IOPTranscript","InvalidTranscript","SerializationError","TranscriptError","append_field_element","append_message","append_serializable_element","borrow","borrow","borrow_mut","borrow_mut","clone","clone_into","deref","deref","deref_mut","deref_mut","drop","drop","fmt","fmt","from","from","from","get_and_append_challenge","get_and_append_challenge_vectors","init","init","into","into","new","to_bytes","to_owned","to_string","try_from","try_from","try_into","try_into","type_id","type_id","vzip","vzip"],"q":[[0,"transcript"],[42,"transcript::errors"],[43,"core::result"],[44,"ark_ff::fields::prime"],[45,"ark_serialize"],[46,"core::clone"],[47,"core::fmt"],[48,"ark_serialize::error"],[49,"alloc::vec"],[50,"alloc::string"],[51,"core::any"]],"i":[0,7,7,0,2,2,2,7,2,7,2,2,2,7,2,7,2,7,2,7,7,7,7,2,2,2,7,2,7,2,2,0,2,7,7,2,7,2,7,2,7,2],"f":"````{{{f{b{d{c}}}}{f{{j{h}}}}{f{c}}}{{A`{ln}}}Ab}{{{f{b{d{c}}}}{f{{j{h}}}}{f{{j{h}}}}}{{A`{ln}}}Ab}{{{f{b{d{c}}}}{f{{j{h}}}}{f{e}}}{{A`{ln}}}AbAd}{{{f{c}}}{{f{e}}}{}{}}0{{{f{bc}}}{{f{be}}}{}{}}0{{{f{{d{c}}}}}{{d{c}}}{AfAb}}{{{f{c}}{f{be}}}l{}{}}{Ah{{f{c}}}{}}0{Ah{{f{bc}}}{}}0{Ahl}0{{{f{n}}{f{bAj}}}Al}0{Ann}{cc{}}0{{{f{b{d{c}}}}{f{{j{h}}}}}{{A`{cn}}}Ab}{{{f{b{d{c}}}}{f{{j{h}}}}Ah}{{A`{{B`{c}}n}}}Ab}{{}Ah}0{ce{}{}}0{{{f{{j{h}}}}}{{d{c}}}Ab}`{{{f{c}}}e{}{}}{{{f{c}}}Bb{}}{c{{A`{e}}}{}{}}000{{{f{c}}}Bd{}}055","D":"Af","p":[[0,"mut"],[5,"IOPTranscript",0],[1,"reference"],[1,"u8"],[1,"slice"],[1,"unit"],[6,"TranscriptError",0,42],[6,"Result",43],[10,"PrimeField",44],[10,"CanonicalSerialize",45],[10,"Clone",46],[1,"usize"],[5,"Formatter",47],[8,"Result",47],[6,"SerializationError",48],[5,"Vec",49],[5,"String",50],[5,"TypeId",51]],"r":[[3,42]],"b":[[19,"impl-Debug-for-TranscriptError"],[20,"impl-Display-for-TranscriptError"]],"c":"OjAAAAAAAAA=","e":"OzAAAAEAAB8AAwAFABEAGQADACEACQA="}],\
["util",{"t":"H","n":["parallelizable_slice_iter"],"q":[[0,"util"],[1,"rayon::slice"],[2,"core::marker"]],"i":[0],"f":"{{{d{{b{c}}}}}{{f{c}}}h}","D":"d","p":[[1,"slice"],[1,"reference"],[5,"Iter",1],[10,"Sync",2]],"r":[],"b":[],"c":"OjAAAAAAAAA=","e":"OjAAAAAAAAA="}]\
]'));
if (typeof exports !== 'undefined') exports.searchIndex = searchIndex;
else if (window.initSearch) window.initSearch(searchIndex);
