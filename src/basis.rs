use ndarray::{Array1, Array2, array};
use crate::integrals;
use crate::math::linalg;


#[derive(Default)]
pub struct RadialBasisFunctions{
    /* Array of radial basis functions g_nl with the structure
            
            (g_00, g_01, ...)   <- n=0
    g   =   (g_10, g_11, ...)
            ( .     .     . )

              ^
              l=0

    and primitives (φ) with the same structure.
    */

    pub g:    Array2<RadialBasisFunction>,
    pub phi:  Array2<PrimitiveBasisFunction>,
}


impl RadialBasisFunctions{

    pub fn n_max(&self) -> usize{
        // n_max is just the number of rows
        self.phi.nrows()
    }


    pub fn l_max(&self) -> usize{
        // l_max is just the number of columns
        self.phi.ncols() - 1
    }


    pub fn construct(&mut self,
                     n_max: usize,
                     l_max: usize,
                     r_cut: f64){
        /* Generate a set of radial basis functions given a maximum value of 
        the angular (l) and radial (n) 'quantum' numbers


        Arguments:
            n_max (int): [1, 10]
            l_max (int): [0, 10]
            r_cut (float): Cut-off for the RBFs
        */
        if n_max < 2{
            panic!("n_max must be at least 2")
        }

        let mut g_flat = Vec::with_capacity(n_max*l_max);
        let mut phi_flat = Vec::with_capacity(n_max*l_max);

        for n in 1..=n_max{
            for l in 0..=l_max{

                // g_nl
                g_flat.push(RadialBasisFunction{beta: Default::default()});

                
                // φ_nl
                let mut phi: PrimitiveBasisFunction = Default::default();
                phi.alpha = PrimitiveBasisFunction::alpha(n, l, r_cut, n_max);
                phi_flat.push(phi);
            }// l
        }// n

        // Reshape into matricies
        self.g = Array2::from_shape_vec((n_max, l_max+1), g_flat).unwrap();
        self.phi = Array2::from_shape_vec((n_max, l_max+1), phi_flat).unwrap();


        self.orthogonalise();
    }


    fn overlap(&self, 
               n: usize, 
               n_prime: usize, 
               l: usize) -> f64{
        /* Compute the overlap of two basis functions with the same angular with
        the same angular number (l) 

                    ∞
        S_nn'l  =   ∫ dr r^2 g_nl(r) g_n'l(r)
                    0 

        for a properly orthogonalised set then for n=n' -> c and n≠n' -> 0 where
        c is a constant. And g_nl(r) is a sum of primitives

               n_max
        g_nl =   Σ   β_nli φ_il
                 i
        */
        let mut s_nnl = 0.0_f64;

        for i in 0..self.n_max(){
            for j in 0..self.n_max(){

                s_nnl += self.g[[n-1, l]].beta[i]
                         * self.g[[n_prime-1, l]].beta[j]
                         * self.phi_overlap(i, j, l);
            }
        }
        return s_nnl
    }


    fn orthogonalise(&mut self){
        /* Orthogonalise the basis functions by setting the β elements
        of each g_nl, such that the overlap integral is zero for two 
        different g_nl functions. Uses Löwdin orthogonalisation 

        β = S^(-1/2)

        */
        let n_max = self.n_max();

        for l in 0..=self.l_max(){

            // Overlap matrix for a particular l
            let mut s = Array2::<f64>::zeros((n_max, n_max));

            // which is symmetric, thus only needs to be looped over the
            // unique elements
            for i in 0..n_max{
                for j in i..n_max{

                    s[[i, j]] = self.phi_overlap(i, j, l);
                    s[[j, i]] = s[[i, j]].clone();
                }
            }

            let beta = linalg::inverse_sqrt_db(&s);
            
            // Set the list of betas for each RBF (g_nl) which are the sum
            // of primitives
            for n in 0..n_max{
                self.g[[n, l]].beta =  beta.column(n).to_owned();
            }// n

        }// l
    }


    fn phi_overlap(&self, 
                   n: usize, 
                   n_prime: usize, 
                   l: usize) -> f64{
        /* Compute the overlap of two primitive radial basis functions with
        the same angular number (l)

                    ∞
        S_nn'l  =   ∫ dr r^2 φ_nl(r) φ_n'l(r)
                    0 

        where the gaussian primitive (φ_nl) has the form

        φ_nl(r) = r^l exp(-α_nl r^2)
        */

        integrals::gaussian_2m((l as i32) + 1, 
                               self.phi[[n, l]].alpha + self.phi[[n_prime, l]].alpha)
    }
}


#[derive(Default)]
pub struct RadialBasisFunction{   
    /*
               n_max
        g_nl =   Σ   β_nil φ_il
                 i
    */
    pub beta: Array1<f64>,       // {β_nil} over i
}


#[derive(Default)]
pub struct PrimitiveBasisFunction{   // φ_nl(r) = r^l exp(-α_nl r^2)
    pub alpha: f64,
}


impl PrimitiveBasisFunction{

    fn alpha(n: usize, l: usize, r_cut: f64, n_max: usize) -> f64{
        /*
        For a radial basis function of the form

        ϕ_nl(r) = r^l e^(-α_nl r^2)   -> α_nl = -ln(ϕ/r^l)/r^2

        the decay parameters are chosen for each l so that the functions decay 
        to a value 1E-3 at

        r_nl = 1 + n(r_cut - 1)/n_max

        {n} = [1, 2, ..., n_max]
        {l} = [0, 1, ..., l_max]
        */

        let r = 1_f64 + (n as f64 - 1_f64) * (r_cut - 1_f64) / (n_max as f64);
        -(0.001_f64 / r.powi(l as i32)).log(std::f64::consts::E) / (r*r)
    }
}



/*
      /$$$$$$$$ /$$$$$$$$  /$$$$$$  /$$$$$$$$ /$$$$$$
     |__  $$__/| $$_____/ /$$__  $$|__  $$__//$$__  $$
        | $$   | $$      | $$  \__/   | $$  | $$  \__/
        | $$   | $$$$$   |  $$$$$$    | $$  |  $$$$$$
        | $$   | $$__/    \____  $$   | $$   \____  $$
        | $$   | $$       /$$  \ $$   | $$   /$$  \ $$
        | $$   | $$$$$$$$|  $$$$$$/   | $$  |  $$$$$$/
        |__/   |________/ \______/    |__/   \______/
*/
#[cfg(test)]
mod tests{

    use super::*;

    fn is_close(x: f64, y: f64, atol: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        (x - y).abs() <= atol
    }


    fn is_far(x: f64, y: f64, atol: f64) -> bool{
        // Are two numbers far away to within an absolute tolerance? 
        (x - y).abs() > atol
    }


    #[test]
    fn test_primitive_overlap(){

        let phi_00 = PrimitiveBasisFunction{alpha: 1.0};
        let phi_10 = PrimitiveBasisFunction{alpha: 1.3};

        let mut rbfs: RadialBasisFunctions = Default::default();

        // *************************** l = 0 ******************************
        rbfs.phi = array![[phi_00], 
                          [phi_10]];

        // Check that the alpha parameters have been set correctly in the array 
        assert!(is_close(rbfs.phi[[0, 0]].alpha, 1.0, 1E-8));
        assert!(is_close(rbfs.phi[[1, 0]].alpha, 1.3, 1E-8));


        let args = integrals::IntegrandArguments{integer_args: vec![0],
                                                 float_args:   vec![1.0, 1.3]};

        fn integrand(r: f64, args: &integrals::IntegrandArguments) -> f64{
            // NOTE: the unneeded moves are intended for readability 
            let l = args.integer_args[0];

            let alpha_nl = args.float_args[0];
            let alpha_nprimel = args.float_args[1];

            r.powi(2*l + 2) * (-(alpha_nl + alpha_nprimel) * r.powi(2)).exp()
        }
        
        let analytic_integral = rbfs.phi_overlap(0, 1, 0);
        let numeric_integral = integrals::trapz(&integrand, 0.0, 10.0, &args, 1E-4);

        assert!(is_close(analytic_integral, numeric_integral, 1E-5));


        // *************************** l = 2 ******************************
        let phi_00 = PrimitiveBasisFunction{alpha: 1.0};
        let phi_10 = PrimitiveBasisFunction{alpha: 1.3};
        let phi_01 = PrimitiveBasisFunction{alpha: 1.0};
        let phi_02 = PrimitiveBasisFunction{alpha: 1.8};
        let phi_11 = PrimitiveBasisFunction{alpha: 0.6};
        let phi_12 = PrimitiveBasisFunction{alpha: 1.9};

        rbfs.phi = array![[phi_00, phi_01, phi_02], 
                          [phi_10, phi_11, phi_12]];

        let args = integrals::IntegrandArguments{integer_args: vec![2],
                                                 float_args:   vec![1.8, 1.9]};

        let analytic_integral = rbfs.phi_overlap(0, 1, 2);
        let numeric_integral = integrals::trapz(&integrand, 0.0, 10.0, &args, 1E-4);

        assert!(is_close(analytic_integral, numeric_integral, 1E-5));

    }


    #[test]
    fn test_radial_overlap(){

        let mut rbfs: RadialBasisFunctions = Default::default();

        // *************************** l = 0 ******************************
        let phi_00 = PrimitiveBasisFunction{alpha: 1.0};
        let phi_10 = PrimitiveBasisFunction{alpha: 1.3};

        rbfs.phi = array![[phi_00], 
                          [phi_10]];

        assert_eq!(rbfs.phi.ncols(), 1);
        assert_eq!(rbfs.phi.nrows(), 2);

        let g_00 = RadialBasisFunction{beta: array![2.0, 1.0]};
        let g_10 = RadialBasisFunction{beta: array![1.0, 2.0]};
        rbfs.g = array![[g_00], 
                        [g_10]];

        assert_eq!(rbfs.g.ncols(), 1);
        assert_eq!(rbfs.g.nrows(), 2);

        for n in 0..2{
            for l in 0..1{
                assert!(is_far(rbfs.phi[[n, l]].alpha, 0.0, 1E-3));
            }
        }

        // Without orthogonalisation the overlap is non-zero: n != n'
        assert!(is_far(rbfs.overlap(1, 2, 0), 0.0, 1E-8));

        rbfs.orthogonalise();

        // For n != n' the overlap should be ~0
        assert!(is_close(rbfs.overlap(1, 2, 0), 0.0, 1E-8));
    }
    

    #[test]
    fn test_simple_construction(){
        // Ensure that the constructed RBFS are orthogonal

        let mut rbfs: RadialBasisFunctions = Default::default();

        // n_max = 2, l_max = 0
        rbfs.construct(2, 0, 3.0);

        assert!(is_close(rbfs.overlap(1, 2, 0), 0.0, 1E-8));

        assert!(is_far(rbfs.overlap(1, 1, 0), 0.0, 1E-8));
        assert!(is_far(rbfs.overlap(2, 2, 0), 0.0, 1E-8));

    }


    #[test]
    fn test_large_n_lmax_construction(){
        // Ensure that the constructed RBFS are orthogonal

        let mut rbfs: RadialBasisFunctions = Default::default();

        let n_max = 8;
        rbfs.construct(n_max, 0, 3.0);

        // Should be orthogonal to all n' not equal to 10 
        for n in 1..n_max{
            //println!("<g_{}| g_{}>_l=0 = {}", 
            //         n_max, n, rbfs.overlap(n_max, n, 0));
                     
            assert!(is_close(rbfs.overlap(n_max, n, 0), 0.0, 1E-8));
        }
    }


    #[test]
    fn test_almost_orthogonal(){
        // For large n_max the orthgonalisation isn't quite perfect..

        let mut rbfs: RadialBasisFunctions = Default::default();

        rbfs.construct(10, 0, 3.0);
        assert!(is_close(rbfs.overlap(10, 1, 0), 0.0, 1E-4));
    }


    #[test]
    #[should_panic]
    fn test_n_max_less_than_2(){
        let mut rbfs: RadialBasisFunctions = Default::default();

        rbfs.construct(1, 0, 3.0);
    }
}
