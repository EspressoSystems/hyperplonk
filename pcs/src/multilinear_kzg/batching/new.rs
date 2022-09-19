/// input: k points `p`, each of dim n
/// input: k MLEs `f`, each of dim n
///
///
/// steps:
/// 1. for `i in [0..k-1]`, set `L_i(x) = LI(p[i])`, where `LI` is Lagrange
/// interpolation of p[i].
/// 2. define `g(y, x0,...x_{n-1}) := \sum_{i=0}^{k-1} LI(y) f_i` which is an
/// (n+1) mle
/// 3. define `h(y) := \sum_{i=0}^{k-1} LI(y) p_i` which are n
/// univariate polynomials
/// 4. evaluate `q(y) := g(y, h1, ...hn)` and obtain a univariate polynomial
/// this is done via `compute_w_circ_l` api
/// 5. commit to `q(y)` and append commitment to transcript
///
/// 6. sample a random field element `r` from transcript
/// 7. generate `g(r, x0,...x_{n-1}) := \sum_{i=0}^{k-1} LI(r) f_i` which is an
/// `n` dim MLE
/// 8. compute `h(r) := \sum_{i=0}^{k-1} LI(r) p_i` which is an `n` dim point
///
/// 9. open `q(y)` at `r` and outputs its evaluation and proof
/// 10. open `g(y, x0,...x_{n-1})` at `r, h1(r),...hn(r)` and outputs its
/// evaluations and proof  
fn open() {}
