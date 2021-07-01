/*
Smooth overlap of atomic positions. Theory from 

Bartok et al, 2013, 10.1103/PhysRevB.87.184115
*/
use crate::geometry::SphericalPolarCoordinate;
use crate::basis::RadialBasisFunctions;
use crate::math::sphr::sphr_harm;


fn c_nlm(sphr_coords: Vec<SphericalPolarCoordinate>,
         rbfs:        RadialBasisFunctions,
         n:           usize,
         l:           i32,
         m:           i32,
         sigma_at:    f64) -> f64{
    /*
    Calculate an element of the SOAP vector as a c_nlm coefficent
    considering the environment around an atom given as their 
    spherical coordinates

    Arguments:
        sphr_coords: Spherical coordinates
        
        rbfs: Radial basis functions

        n (int): >0

        l (int): >= 0

        m (int): [-l, .., 0, .., +l]

        sigma_at (float): Spead of the Gaussians in the SOAP. 0.5 Å
                          is usually a good value
    
    Returns:
        (float): c_nlm
    */

    let sqrt_eight_pi_cubed = 15.7496099457224197442_f64;
    let a = 1_f64 / (2_f64 * sigma_at * sigma_at);    
    let mut sum_i = 0_f64;

    for coord in sphr_coords.iter(){   // Σ_ι
        
        let mut sum = 0_f64;
        
        for n_prime in 0..rbfs.n_max(){
            sum += rbfs.g[[n-1, l as usize]].beta[n_prime] 
                   * ((2_f64*a*coord.r).powi(l)
                      / (2_f64*(rbfs.phi[[n_prime, l as usize]].alpha + a)).powf(l as f64+1.5))
                    * ((a * coord.r).powi(2) / (rbfs.phi[[n_prime, l as usize]].alpha + a)).exp();
        }// n'
        sum_i += sphr_harm(coord, l, m) * (-a*coord.r*coord.r).exp() * sum;
    }// i
    
    sqrt_eight_pi_cubed * sum_i
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

    fn is_very_close(x: f64, y: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        println!("\nleft = {}\nright = {}", x, y);
        (x - y).abs() <= 1E-8
    }
    

    #[test]
    fn test_c_100(){
        // Simple SOAP coefficent

        let coord = SphericalPolarCoordinate{r: 1.0,
                                             theta: 1.5707963267948966, // π/2 
                                             phi: 1.0471975511965976};  // π/3

        let mut rbfs: RadialBasisFunctions = Default::default();
         // n_max = 2, l_max = 0, r_cut = 2.0
         rbfs.construct(2, 0, 2.0);
        
        assert!(is_very_close(c_nlm(vec![coord], rbfs, 1, 0, 0, 0.5),
                              -0.07133018467726725  // Python implementation
                              ));
    }


}
