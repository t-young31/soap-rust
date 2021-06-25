use ndarray::{Array1, Array2, array};
use crate::integrals;


#[derive(Default)]
pub struct RadialBasisFunctions{
    /* Array of radial basis functions g_nl with the structure
            
            (g_00, g_01, ...)   <- n=0
    g   =   (g_10, g_11, ...)
            ( .     .     . )

              ^
              l=0

    and primitives (φ) with the same structure. Beta is a 3-index
    tensor containing the orthogonalising weights (nn'l)
    */

    pub g:    Array2<RadialBasisFunction>,
    pub phi:  Array2<PrimitiveBasisFunction>,
}


impl RadialBasisFunctions{

    fn n_max(&self) -> usize{
        // n_max is just the number of rows
        self.phi.shape()[0]
    }

    fn l_max(&self) -> usize{
        // l_max is just the number of columns
        self.phi.shape()[1]
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

                s_nnl += self.g[[n, l]].beta[i] 
                     * self.g[[n_prime, l]].beta[j]
                     * self.phi_overlap(i, j, l);
            }
        }
        return s_nnl
    }


    fn orthogonalise(&self){
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


            /*

            Why are there no working linear algebra libraries?!

            */

            // let diag_elems = s.eigenvalues(); // .sqrt().inv();

            // let (lamdas, u) = a.clone().eigh(UPLO::Upper).unwrap();

            // println!("{}", lamdas);
        
        }




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


    fn phi_value(&self, 
                 n: usize,
                 l: usize,
                 r: f64) -> f64{

        /*
         φ_nl(r) = r^l exp(-α_nl r^2)

        Returns:
            (float): Value at a particular distance, r
        */
        r.powi(l as i32) * (-self.phi[[n, l]].alpha * r.powi(2)).exp()
    }

}


#[derive(Default)]
pub struct RadialBasisFunction{   // g_nl
    pub beta: Array1<f64>,       // {β_nn'l} over n'
}


#[derive(Default)]
pub struct PrimitiveBasisFunction{   // φ_nl
    alpha: f64,
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

    /*
    TODO: fix


    #[test]
    fn test_radial_overlap(){

        let mut rbfs: RadialBasisFunctions = Default::default();

        // *************************** l = 0 ******************************
        let phi_00 = PrimitiveBasisFunction{alpha: 1.0};
        let phi_10 = PrimitiveBasisFunction{alpha: 1.3};

        rbfs.phi = array![[phi_00], 
                          [phi_10]];

        println!("{}", rbfs.n_max());

        // Without orthogonalisation the overlap is non-zero: n != n'
        assert!(is_far(rbfs.overlap(0, 1, 0), 0.0, 1E-8));

        return
        rbfs.orthogonalise();

        // For n != n' the overlap should be ~0
        assert!(is_close(rbfs.overlap(0, 1, 0), 0.0, 1E-8));
        assert!(is_far(rbfs.overlap(0, 1, 0), 0.0, 1E-8));

    }
    */

}
