/// input: k points `p`, each of dim n
/// input: k MLEs `f`, each of dim n
///
///
/// steps:
/// 2. define `g(y, x0,...x_{n-1}) := \sum_{i=0}^{k-1} L_i(y) f_i` which is an
/// (n+1) mle (L_i(y) is implicit)
/// 3. define `h(y) := \sum_{i=0}^{k-1} L_i(y) p_i` which are n (h(y) is implicit)
/// univariate polynomials
/// 4. evaluate `q(y) := g(y, h1, ...hn)` and obtain a univariate polynomial
/// this is done via ~~`compute_w_circ_l` api~~ interpolation
///   4.1 y in [w..w^k] -> n evalation: y_1...y_k
///   4.2 evaluate y at alpha...alpha (d* mu +1 )(k-1) - k
///   4.3 interpolate to get `q(y)` explict form
/// 5. commit to `q(y)` and append commitment to transcript
///
/// 6. sample a random field element `r` from transcript
/// 7. generate `g(r, x0,...x_{n-1}) := \sum_{i=0}^{k-1} L_i(r) f_i` which is an
/// `n` dim MLE
/// 8. compute `h(r) := \sum_{i=0}^{k-1} L_i(r) p_i` which is an `n` dim point
///
/// 9. open `q(y)` at `r` and outputs its evaluation and proof
/// 10. open `g(y, x0,...x_{n-1})` at `r, h1(r),...hn(r)` and outputs its
/// evaluations and proof  
fn open() {}
