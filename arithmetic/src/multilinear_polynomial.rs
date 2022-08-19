use ark_ff::{Field, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{
    cfg_iter, end_timer, fmt,
    fmt::Formatter,
    ops::{Add, AddAssign, Index, Neg, Sub, SubAssign},
    rand::{Rng, RngCore},
    slice::{Iter, IterMut},
    start_timer,
    vec::Vec,
};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::rc::Rc;

use crate::{prelude::ArithErrors, util::swap_bits, MultilinearExtension};

/// Stores a multilinear polynomial in dense evaluation form.
#[derive(Clone, PartialEq, Eq, Hash, Default, CanonicalSerialize, CanonicalDeserialize)]
pub struct DenseMultilinearExtension<F: Field> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: Vec<F>,
    /// Number of variables
    pub num_vars: usize,
}

impl<F: Field> DenseMultilinearExtension<F> {
    /// Construct a new polynomial from a list of evaluations where the index
    /// represents a point in {0,1}^`num_vars` in little endian form. For
    /// example, `0b1011` represents `P(1,1,0,1)`
    pub fn from_evaluations_slice(num_vars: usize, evaluations: &[F]) -> Self {
        Self::from_evaluations_vec(num_vars, evaluations.to_vec())
    }

    /// Construct a new polynomial from a list of evaluations where the index
    /// represents a point in {0,1}^`num_vars` in little endian form. For
    /// example, `0b1011` represents `P(1,1,0,1)`
    pub fn from_evaluations_vec(num_vars: usize, evaluations: Vec<F>) -> Self {
        // assert that the number of variables matches the size of evaluations
        assert_eq!(
            evaluations.len(),
            1 << num_vars,
            "The size of evaluations should be 2^num_vars."
        );

        Self {
            num_vars,
            evaluations,
        }
    }
    /// Relabel the point inplace by switching `k` scalars from position `a` to
    /// position `b`, and from position `b` to position `a` in vector.
    ///
    /// This function turns `P(x_1,...,x_a,...,x_{a+k - 1},...,x_b,...,x_{b+k -
    /// 1},...,x_n)` to `P(x_1,...,x_b,...,x_{b+k - 1},...,x_a,...,x_{a+k -
    /// 1},...,x_n)`
    pub fn relabel_inplace(&mut self, mut a: usize, mut b: usize, k: usize) {
        // enforce order of a and b
        if a > b {
            ark_std::mem::swap(&mut a, &mut b);
        }
        assert!(
            a + k < self.num_vars && b + k < self.num_vars,
            "invalid relabel argument"
        );
        if a == b || k == 0 {
            return;
        }
        assert!(a + k <= b, "overlapped swap window is not allowed");
        for i in 0..self.evaluations.len() {
            let j = swap_bits(i, a, b, k);
            if i < j {
                self.evaluations.swap(i, j);
            }
        }
    }

    /// Returns an iterator that iterates over the evaluations over
    /// {0,1}^`num_vars`
    pub fn iter(&self) -> Iter<'_, F> {
        self.evaluations.iter()
    }

    /// Returns a mutable iterator that iterates over the evaluations over
    /// {0,1}^`num_vars`
    pub fn iter_mut(&mut self) -> IterMut<'_, F> {
        self.evaluations.iter_mut()
    }
}

impl<F: Field> MultilinearExtension<F> for DenseMultilinearExtension<F> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn evaluate(&self, point: &[F]) -> Option<F> {
        if point.len() == self.num_vars {
            Some(self.fix_variables(point)[0])
        } else {
            None
        }
    }

    fn rand<R: Rng>(num_vars: usize, rng: &mut R) -> Self {
        Self::from_evaluations_vec(
            num_vars,
            (0..(1 << num_vars)).map(|_| F::rand(rng)).collect(),
        )
    }

    fn relabel(&self, a: usize, b: usize, k: usize) -> Self {
        let mut copied = self.clone();
        copied.relabel_inplace(a, b, k);
        copied
    }

    fn fix_variables(&self, partial_point: &[F]) -> Self {
        assert!(
            partial_point.len() <= self.num_vars,
            "invalid size of partial point"
        );
        let mut poly = self.evaluations.to_vec();
        let nv = self.num_vars;
        let dim = partial_point.len();
        // evaluate single variable of partial point from left to right
        for i in 1..dim + 1 {
            let r = partial_point[i - 1];
            for b in 0..(1 << (nv - i)) {
                poly[b] = poly[b << 1] * (F::one() - r) + poly[(b << 1) + 1] * r;
            }
        }
        Self::from_evaluations_slice(nv - dim, &poly[..(1 << (nv - dim))])
    }

    fn to_evaluations(&self) -> Vec<F> {
        self.evaluations.to_vec()
    }
}

impl<F: Field> Index<usize> for DenseMultilinearExtension<F> {
    type Output = F;

    /// Returns the evaluation of the polynomial at a point represented by
    /// index.
    ///
    /// Index represents a vector in {0,1}^`num_vars` in little endian form. For
    /// example, `0b1011` represents `P(1,1,0,1)`
    ///
    /// For dense multilinear polynomial, `index` takes constant time.
    fn index(&self, index: usize) -> &Self::Output {
        &self.evaluations[index]
    }
}

impl<F: Field> Add for DenseMultilinearExtension<F> {
    type Output = DenseMultilinearExtension<F>;

    fn add(self, other: DenseMultilinearExtension<F>) -> Self {
        &self + &other
    }
}

impl<'a, 'b, F: Field> Add<&'a DenseMultilinearExtension<F>> for &'b DenseMultilinearExtension<F> {
    type Output = DenseMultilinearExtension<F>;

    fn add(self, rhs: &'a DenseMultilinearExtension<F>) -> Self::Output {
        // handle constant zero case
        if rhs.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return rhs.clone();
        }
        assert_eq!(self.num_vars, rhs.num_vars);
        let result: Vec<F> = cfg_iter!(self.evaluations)
            .zip(cfg_iter!(rhs.evaluations))
            .map(|(a, b)| *a + *b)
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result)
    }
}

impl<F: Field> AddAssign for DenseMultilinearExtension<F> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + &other;
    }
}

impl<'a, 'b, F: Field> AddAssign<&'a DenseMultilinearExtension<F>>
    for DenseMultilinearExtension<F>
{
    fn add_assign(&mut self, other: &'a DenseMultilinearExtension<F>) {
        *self = &*self + other;
    }
}

impl<'a, 'b, F: Field> AddAssign<(F, &'a DenseMultilinearExtension<F>)>
    for DenseMultilinearExtension<F>
{
    fn add_assign(&mut self, (f, other): (F, &'a DenseMultilinearExtension<F>)) {
        let other = Self {
            num_vars: other.num_vars,
            evaluations: cfg_iter!(other.evaluations).map(|x| f * x).collect(),
        };
        *self = &*self + &other;
    }
}

impl<F: Field> Neg for DenseMultilinearExtension<F> {
    type Output = DenseMultilinearExtension<F>;

    fn neg(self) -> Self::Output {
        Self::Output {
            num_vars: self.num_vars,
            evaluations: cfg_iter!(self.evaluations).map(|x| -*x).collect(),
        }
    }
}

impl<F: Field> Sub for DenseMultilinearExtension<F> {
    type Output = DenseMultilinearExtension<F>;

    fn sub(self, other: DenseMultilinearExtension<F>) -> Self {
        &self - &other
    }
}

impl<'a, 'b, F: Field> Sub<&'a DenseMultilinearExtension<F>> for &'b DenseMultilinearExtension<F> {
    type Output = DenseMultilinearExtension<F>;

    fn sub(self, rhs: &'a DenseMultilinearExtension<F>) -> Self::Output {
        self + &rhs.clone().neg()
    }
}

impl<F: Field> SubAssign for DenseMultilinearExtension<F> {
    fn sub_assign(&mut self, other: Self) {
        *self = &*self - &other;
    }
}

impl<'a, 'b, F: Field> SubAssign<&'a DenseMultilinearExtension<F>>
    for DenseMultilinearExtension<F>
{
    fn sub_assign(&mut self, other: &'a DenseMultilinearExtension<F>) {
        *self = &*self - other;
    }
}

impl<F: Field> fmt::Debug for DenseMultilinearExtension<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "DenseML(nv = {}, evaluations = [", self.num_vars)?;
        for i in 0..ark_std::cmp::min(4, self.evaluations.len()) {
            write!(f, "{} ", self.evaluations[i])?;
        }
        if self.evaluations.len() < 4 {
            write!(f, "])")?;
        } else {
            write!(f, "...])")?;
        }
        Ok(())
    }
}

impl<F: Field> Zero for DenseMultilinearExtension<F> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: vec![F::zero()],
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations[0].is_zero()
    }
}

impl<F: Field> DenseMultilinearExtension<F> {
    // Build a randomize list of mle-s whose sum is zero.
    pub fn random_zero_mle_list<R: RngCore>(
        nv: usize,
        degree: usize,
        rng: &mut R,
    ) -> Vec<Rc<DenseMultilinearExtension<F>>> {
        let start = start_timer!(|| "sample random zero mle list");

        let mut multiplicands = Vec::with_capacity(degree);
        for _ in 0..degree {
            multiplicands.push(Vec::with_capacity(1 << nv))
        }
        for _ in 0..(1 << nv) {
            multiplicands[0].push(F::zero());
            for e in multiplicands.iter_mut().skip(1) {
                e.push(F::rand(rng));
            }
        }

        let list = multiplicands
            .into_iter()
            .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
            .collect();

        end_timer!(start);
        list
    }

    /// Sample a random list of multilinear polynomials.
    /// Returns
    /// - the list of polynomials,
    /// - its sum of polynomial evaluations over the boolean hypercube.
    pub fn random_mle_list<R: RngCore>(
        nv: usize,
        degree: usize,
        rng: &mut R,
    ) -> (Vec<Rc<DenseMultilinearExtension<F>>>, F) {
        let start = start_timer!(|| "sample random mle list");
        let mut multiplicands = Vec::with_capacity(degree);
        for _ in 0..degree {
            multiplicands.push(Vec::with_capacity(1 << nv))
        }
        let mut sum = F::zero();

        for _ in 0..(1 << nv) {
            let mut product = F::one();

            for e in multiplicands.iter_mut() {
                let val = F::rand(rng);
                e.push(val);
                product *= val;
            }
            sum += product;
        }

        let list = multiplicands
            .into_iter()
            .map(|x| Rc::new(DenseMultilinearExtension::from_evaluations_vec(nv, x)))
            .collect();

        end_timer!(start);
        (list, sum)
    }

    // This function build the eq(x, r) polynomial for any given r.
    //
    // Evaluate
    //      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
    // over r, which is
    //      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
    pub fn build_eq_x_r(r: &[F]) -> Result<Rc<DenseMultilinearExtension<F>>, ArithErrors> {
        let start = start_timer!(|| "zero check build eq_x_r");

        // we build eq(x,r) from its evaluations
        // we want to evaluate eq(x,r) over x \in {0, 1}^num_vars
        // for example, with num_vars = 4, x is a binary vector of 4, then
        //  0 0 0 0 -> (1-r0)   * (1-r1)    * (1-r2)    * (1-r3)
        //  1 0 0 0 -> r0       * (1-r1)    * (1-r2)    * (1-r3)
        //  0 1 0 0 -> (1-r0)   * r1        * (1-r2)    * (1-r3)
        //  1 1 0 0 -> r0       * r1        * (1-r2)    * (1-r3)
        //  ....
        //  1 1 1 1 -> r0       * r1        * r2        * r3
        // we will need 2^num_var evaluations

        let mut eval = Vec::new();
        build_eq_x_r_helper::<F>(r, &mut eval)?;

        let mle = DenseMultilinearExtension::from_evaluations_vec(r.len(), eval);

        let res = Rc::new(mle);
        end_timer!(start);
        Ok(res)
    }
}

/// A helper function to build eq(x, r) recursively.
/// This function takes `r.len()` steps, and for each step it requires a maximum
/// `r.len()-1` multiplications.
fn build_eq_x_r_helper<F: Field>(r: &[F], buf: &mut Vec<F>) -> Result<(), ArithErrors> {
    if r.is_empty() {
        return Err(ArithErrors::InvalidParameters("r length is 0".to_string()));
    } else if r.len() == 1 {
        // initializing the buffer with [1-r_0, r_0]
        buf.push(F::one() - r[0]);
        buf.push(r[0]);
    } else {
        build_eq_x_r_helper(&r[1..], buf)?;

        // suppose at the previous step we received [b_1, ..., b_k]
        // for the current step we will need
        // if x_0 = 0:   (1-r0) * [b_1, ..., b_k]
        // if x_0 = 1:   r0 * [b_1, ..., b_k]

        let mut res = vec![];
        for &b_i in buf.iter() {
            let tmp = r[0] * b_i;
            res.push(b_i - tmp);
            res.push(tmp);
        }
        *buf = res;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use crate::{
        prelude::{DenseMultilinearExtension, MultilinearExtension},
        util::bit_decompose,
    };
    use ark_bls12_381::Fr;
    use ark_ff::{Field, Zero};
    use ark_std::{end_timer, ops::Neg, start_timer, test_rng, vec::Vec, UniformRand};

    /// utility: evaluate multilinear extension (in form of data array) at a
    /// random point
    fn evaluate_data_array<F: Field>(data: &[F], point: &[F]) -> F {
        if data.len() != (1 << point.len()) {
            panic!("Data size mismatch with number of variables. ")
        }

        let nv = point.len();
        let mut a = data.to_vec();

        for i in 1..nv + 1 {
            let r = point[i - 1];
            for b in 0..(1 << (nv - i)) {
                a[b] = a[b << 1] * (F::one() - r) + a[(b << 1) + 1] * r;
            }
        }
        a[0]
    }

    #[test]
    fn evaluate_at_a_point() {
        let mut rng = test_rng();
        let poly = DenseMultilinearExtension::rand(10, &mut rng);
        for _ in 0..10 {
            let point: Vec<_> = (0..10).map(|_| Fr::rand(&mut rng)).collect();
            assert_eq!(
                evaluate_data_array(&poly.evaluations, &point),
                poly.evaluate(&point).unwrap()
            )
        }
    }

    #[test]
    fn relabel_polynomial() {
        let mut rng = test_rng();
        for _ in 0..20 {
            let mut poly = DenseMultilinearExtension::rand(10, &mut rng);
            let mut point: Vec<_> = (0..10).map(|_| Fr::rand(&mut rng)).collect();

            let expected = poly.evaluate(&point);

            poly.relabel_inplace(2, 2, 1); // should have no effect
            assert_eq!(expected, poly.evaluate(&point));

            poly.relabel_inplace(3, 4, 1); // should switch 3 and 4
            point.swap(3, 4);
            assert_eq!(expected, poly.evaluate(&point));

            poly.relabel_inplace(7, 5, 1);
            point.swap(7, 5);
            assert_eq!(expected, poly.evaluate(&point));

            poly.relabel_inplace(2, 5, 3);
            point.swap(2, 5);
            point.swap(3, 6);
            point.swap(4, 7);
            assert_eq!(expected, poly.evaluate(&point));

            poly.relabel_inplace(7, 0, 2);
            point.swap(0, 7);
            point.swap(1, 8);
            assert_eq!(expected, poly.evaluate(&point));
        }
    }

    #[test]
    fn arithmetic() {
        const NV: usize = 10;
        let mut rng = test_rng();
        for _ in 0..20 {
            let point: Vec<_> = (0..NV).map(|_| Fr::rand(&mut rng)).collect();
            let poly1 = DenseMultilinearExtension::rand(NV, &mut rng);
            let poly2 = DenseMultilinearExtension::rand(NV, &mut rng);
            let v1 = poly1.evaluate(&point).unwrap();
            let v2 = poly2.evaluate(&point).unwrap();
            // test add
            assert_eq!((&poly1 + &poly2).evaluate(&point).unwrap(), v1 + v2);
            // test sub
            assert_eq!((&poly1 - &poly2).evaluate(&point).unwrap(), v1 - v2);
            // test negate
            assert_eq!(poly1.clone().neg().evaluate(&point).unwrap(), -v1);
            // test add assign
            {
                let mut poly1 = poly1.clone();
                poly1 += &poly2;
                assert_eq!(poly1.evaluate(&point).unwrap(), v1 + v2)
            }
            // test sub assign
            {
                let mut poly1 = poly1.clone();
                poly1 -= &poly2;
                assert_eq!(poly1.evaluate(&point).unwrap(), v1 - v2)
            }
            // test add assign with scalar
            {
                let mut poly1 = poly1.clone();
                let scalar = Fr::rand(&mut rng);
                poly1 += (scalar, &poly2);
                assert_eq!(poly1.evaluate(&point).unwrap(), v1 + scalar * v2)
            }
            // test additive identity
            {
                assert_eq!(&poly1 + &DenseMultilinearExtension::zero(), poly1);
                assert_eq!(&DenseMultilinearExtension::zero() + &poly1, poly1);
                {
                    let mut poly1_cloned = poly1.clone();
                    poly1_cloned += &DenseMultilinearExtension::zero();
                    assert_eq!(&poly1_cloned, &poly1);
                    let mut zero = DenseMultilinearExtension::zero();
                    let scalar = Fr::rand(&mut rng);
                    zero += (scalar, &poly1);
                    assert_eq!(zero.evaluate(&point).unwrap(), scalar * v1);
                }
            }
        }
    }

    #[test]
    fn test_eq_xr() {
        let mut rng = test_rng();
        for nv in 4..10 {
            let r: Vec<Fr> = (0..nv).map(|_| Fr::rand(&mut rng)).collect();
            let eq_x_r = DenseMultilinearExtension::build_eq_x_r(r.as_ref()).unwrap();
            let eq_x_r2 = build_eq_x_r_for_test(r.as_ref());
            assert_eq!(eq_x_r, eq_x_r2);
        }
    }

    /// Naive method to build eq(x, r).
    /// Only used for testing purpose.
    // Evaluate
    //      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
    // over r, which is
    //      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
    fn build_eq_x_r_for_test<F: Field>(r: &[F]) -> Rc<DenseMultilinearExtension<F>> {
        let start = start_timer!(|| "zero check naive build eq_x_r");

        // we build eq(x,r) from its evaluations
        // we want to evaluate eq(x,r) over x \in {0, 1}^num_vars
        // for example, with num_vars = 4, x is a binary vector of 4, then
        //  0 0 0 0 -> (1-r0)   * (1-r1)    * (1-r2)    * (1-r3)
        //  1 0 0 0 -> r0       * (1-r1)    * (1-r2)    * (1-r3)
        //  0 1 0 0 -> (1-r0)   * r1        * (1-r2)    * (1-r3)
        //  1 1 0 0 -> r0       * r1        * (1-r2)    * (1-r3)
        //  ....
        //  1 1 1 1 -> r0       * r1        * r2        * r3
        // we will need 2^num_var evaluations

        // First, we build array for {1 - r_i}
        let one_minus_r: Vec<F> = r.iter().map(|ri| F::one() - ri).collect();

        let num_var = r.len();
        let mut eval = vec![];

        for i in 0..1 << num_var {
            let mut current_eval = F::one();
            let bit_sequence = bit_decompose(i, num_var);

            for (&bit, (ri, one_minus_ri)) in
                bit_sequence.iter().zip(r.iter().zip(one_minus_r.iter()))
            {
                current_eval *= if bit { *ri } else { *one_minus_ri };
            }
            eval.push(current_eval);
        }

        let mle = DenseMultilinearExtension::from_evaluations_vec(num_var, eval);

        let res = Rc::new(mle);
        end_timer!(start);
        res
    }
}
