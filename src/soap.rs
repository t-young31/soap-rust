/*
Smooth overlap of atomic positions. Theory from 

Bartok et al, 2013, 10.1103/PhysRevB.87.184115
*/
use crate::geometry::SphericalPolarCoordinate;
use crate::structure::Structure;
use crate::basis::RadialBasisFunctions;
use crate::math::sphr::sphr_harm;


fn c_nlm(sphr_coords: &Vec<SphericalPolarCoordinate>,
         rbfs:        &RadialBasisFunctions,
         n:           usize,
         l:           usize,
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
    let l_f64 = l as f64;
    let l_i32 = l as i32;
    
    let mut sum_i = 0_f64;

    for coord in sphr_coords.iter(){   // Σ_ι
        
        let mut sum = 0_f64;
        
        for n_prime in 0..rbfs.n_max(){

            sum += rbfs.g[[n-1, l as usize]].beta[n_prime] 
                   * ((2_f64*a*coord.r).powi(l_i32)
                      / (2_f64*(rbfs.phi[[n_prime, l]].alpha + a)).powf(l_f64+1.5))
                    * ((a * coord.r).powi(2) / (rbfs.phi[[n_prime, l]].alpha + a)).exp();
        }// n'
        sum_i += sphr_harm(coord, l_i32, m) * (-a*coord.r*coord.r).exp() * sum;
    
    }// i
    
    sqrt_eight_pi_cubed * sum_i
}


fn power_spectrum(structure:      &Structure,
                  atom_idx:       usize,
                  nbr_element:    &str,
                  n_max:          usize,
                  l_max:          usize,
                  r_cut:          f64,
                  sigma_at:       f64) -> Vec<f64>{
    /*
    Compute the vector comprising the elements of the power spectrum
    p_nn'l, where n, n'∈[1, n_max] and l ∈ [0, l_max]. n_max must be at
    least 2.

    Arguments:
        structure: strucutre (box/molecule/system) defining the geometry

        atom_idx (int): Index of the atom to expand about

        nbr_element (str): String of the element about which to expand
                           the neighbour density e.g. "C"
        
        n_max (int): Maximum radial expansion basis (>1)

        l_max (int): Maximum angular expansion basis (≥0)
    
        r_cut (float): Cut-off distance used in the construction of the
                       radial basis functions

        sigma_at (float): Width of the Gaussians used in the SOAP

    Returns:
        (vector(float)): Concatenated elements
    */
   
    let sphr_nbrs = structure.sphr_neighbours(atom_idx, 
                                              nbr_element,
                                              r_cut);
 
    let mut rbfs: RadialBasisFunctions = Default::default();
    rbfs.construct(n_max, l_max, r_cut);

    
    let mut p = Vec::<f64>::with_capacity(sphr_nbrs.len() * n_max * (l_max + 1));

    for n in 1..=n_max{
        for n_prime in n..=n_max{
            for l in 0..=l_max{

                let l_i32 = l as i32;
                let mut sum_m = 0_f64;
            
                for m_idx in 0..=2*l{   // 2l + 1 terms

                    let m = (m_idx as i32) - l_i32;
                    let _c_nlm = c_nlm(&sphr_nbrs, &rbfs, n, l, m, sigma_at);

                    sum_m += c_nlm(&sphr_nbrs, &rbfs, n, l, m, sigma_at) 
                             * c_nlm(&sphr_nbrs, &rbfs, n_prime, l, m, sigma_at);

                    println!("n = {}, l = {}, m = {}", n, l, m); 

                }// m

                // Append the normalised component of the power spectrum
                p.push(sum_m / ((2*l + 1) as f64).sqrt())        

            }// l
        }// n'
    }// n
   
    p
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
        is_close(x, y, 1E-8)
    }
    

    fn is_close(x: f64, y: f64, atol: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        println!("\nleft = {}\nright = {}", x, y);
        (x - y).abs() <= atol
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
        
        assert!(is_very_close(c_nlm(&vec![coord], &rbfs, 1, 0, 0, 0.5),
                              -0.07133018467726725  // Python implementation
                              ));
    }


    #[test]
    fn test_c_200(){
        // Simple SOAP coefficent

        let coord = SphericalPolarCoordinate{r: 1.0,
                                             theta: 1.5707963267948966, // π/2 
                                             phi: 1.0471975511965976};  // π/3

        let mut rbfs: RadialBasisFunctions = Default::default();
         // n_max = 2, l_max = 0, r_cut = 2.0
        rbfs.construct(2, 0, 2.0);
        
        assert!(is_very_close(c_nlm(&vec![coord], &rbfs, 2, 0, 0, 0.5),
                              0.28493141517019677  // Python implementation
                              ));
    }


    #[test]
    fn test_c_211(){
        // Simple SOAP coefficent

        let coord = SphericalPolarCoordinate{r: 1.0,
                                             theta: 1.5707963267948966, // π/2 
                                             phi: 1.0471975511965976};  // π/3

        let mut rbfs: RadialBasisFunctions = Default::default();
         // n_max = 2, l_max = 1, r_cut = 2.0
        rbfs.construct(2, 1, 2.0);
        
        assert!(is_very_close(c_nlm(&vec![coord], &rbfs, 2, 1, 1, 0.5),
                              0.1795563639871126  // Python implementation
                              ));
    }


    #[test]
    fn test_power_spectrum_methane(){
        // Test the SOAP vector generated for a methane strucutre 
        // (just RDKit generated)
        
	    std::fs::write("methane.xyz", 
                       "5\n\n\
                        C     0.00000   0.00000   0.00000\n\
                        H    -0.65860  -0.85220  -0.30120\n\
                        H    -0.45940   0.97110  -0.28590\n\
                        H     0.08440  -0.02940   1.10060\n\
                        H     1.02910  -0.10990  -0.41250\n")
                       .expect("Failed to write methane.xyz!");
  
 
        let p = power_spectrum(&Structure::from("methane.xyz"),
                               0,          // Atom index to expand about
                               "H",        // Neighbour denisty
                               2,          // n_max
                               0,          // l_max
                               2.0_f64,    // r_cut
                               0.5_f64);   // σ_atom        
 
        assert_eq!(p.len(), 3); // only (n=1, n'=1), (n=1, n'=2), (n=2, n'=2)
                                // using only the unique pairs
     
        let p_expected = vec![0.07695406_f64, -0.24412914_f64, 0.7744755_f64];
        
        for i in 0..p.len(){
            assert!(is_close(p[i], p_expected[i], 1E-5));
        }

        std::fs::remove_file("methane.xyz").expect("Could not remove file!");
    }


    #[test]
    fn test_power_spectrum_methane_trns(){
        // Test the SOAP vector is translationally invariant
        
	    std::fs::write("methane_trns.xyz", 
                       "5\n\n\
                        C     1.00000   0.00000   0.00000\n\
                        H     0.34140  -0.85220  -0.30120\n\
                        H     0.54060   0.97110  -0.28590\n\
                        H     1.08440  -0.02940   1.10060\n\
                        H     2.02910  -0.10990  -0.41250\n")
                       .expect("Failed to write methane.xyz!");
  
 
        let p = power_spectrum(&Structure::from("methane_trns.xyz"),
                               0, "H", 2, 0, 2.0_f64, 0.5_f64);        
 
        let p_expected = vec![0.07695406_f64, -0.24412914_f64, 0.7744755_f64];
        
        for i in 0..p.len(){
            assert!(is_close(p[i], p_expected[i], 1E-5));
        }

        std::fs::remove_file("methane_trns.xyz").expect("Could not remove file!");
    }

    #[test]
    fn test_power_spectrum_methane_rot(){
        // Test the SOAP vector is rotationally invariant
        
	    std::fs::write("methane_rot.xyz",
                       "5\n\n\
                        C     0.00000   0.00000   0.00000\n\
                        H    -0.82213065 -0.46572569 -0.598265 \n\
                        H    -0.41587856  0.67019903  0.78339049 \n\
                        H     0.54036726 -0.82165307  0.50219273\n\
                        H     0.73410205  0.52779665 -0.65100379\n")
                       .expect("Failed to write methane.xyz!");
  
 
        let p = power_spectrum(&Structure::from("methane_rot.xyz"),
                               0, "H", 2, 0, 2.0_f64, 0.5_f64);        
 
        let p_expected = vec![0.07695406_f64, -0.24412914_f64, 0.7744755_f64];
        
        for i in 0..p.len(){
            assert!(is_close(p[i], p_expected[i], 1E-5));
        }

        std::fs::remove_file("methane_rot.xyz").expect("Could not remove file!");
    }

    #[test]
    fn test_power_spectrum_methane_disp(){
        // Test the SOAP vector isn't the same for all input coordinates
        
	    std::fs::write("methane_disp.xyz",
                       "5\n\n\
                        C     0.00000   0.00000   0.00000\n\
                        H    -0.52213065 -0.46572569 -0.598265 \n\
                        H    -0.713487856  0.67019903  0.78339049 \n\
                        H     0.84036726 -0.82165307  0.50219273\n\
                        H     0.13410205  0.52779665 -0.65100379\n")
                       .expect("Failed to write methane.xyz!");
  
 
        let p = power_spectrum(&Structure::from("methane_disp.xyz"),
                               0, "H", 2, 0, 2.0_f64, 0.5_f64);        
 
        let p_expected = vec![0.07695406_f64, -0.24412914_f64, 0.7744755_f64];
        
        for i in 0..p.len(){
            assert!(! is_close(p[i], p_expected[i], 1E-5));
        }

        std::fs::remove_file("methane_disp.xyz").expect("Could not remove file!");
    }


}
