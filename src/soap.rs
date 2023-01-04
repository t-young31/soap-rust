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


pub fn power_spectrum(
    structure: &Structure,
    atom_idx: usize,
    species: &Vec<String>,
    n_max: usize,
    l_max: usize,
    r_cut: f64,
    sigma_at: f64,
) -> Vec<f64> {
    /*
    Compute the vector comprising the elements of the power spectrum
    p_nn'l, where n, n'∈[1, n_max] and l ∈ [0, l_max]. n_max must be at
    least 2.

    Arguments:
        structure: strucutre (box/molecule/system) defining the geometry

        atom_idx (int): Index of the atom to expand about

        species (Vec<String>): The symbols of the elements to consider. e.g. H

        n_max (int): Maximum radial expansion basis (>1)

        l_max (int): Maximum angular expansion basis (≥0)

        r_cut (float): Cut-off distance used in the construction of the
                       radial basis functions

        sigma_at (float): Width of the Gaussians used in the SOAP

    Returns:
        (vector(float)): Concatenated elements
    */

    let mut rbfs: RadialBasisFunctions = Default::default();
    rbfs.construct(n_max, l_max, r_cut);

    let mut c = vec![];

    for atom_number in species {
        let mut cz = vec![];
        let sphr_nbrs = structure.sphr_neighbours(atom_idx, atom_number, r_cut);
        for n in 1..=n_max {
            let mut c_n = vec![];

            for l in 0..=l_max {
                let mut c_nl = vec![];
                let l_i32 = l as i32;

                for m_idx in 0..=2 * l {
                    // 2l + 1 terms

                    let m = (m_idx as i32) - l_i32;

                    c_nl.push(c_nlm(&sphr_nbrs, &rbfs, n, l, m, sigma_at));
                } // m
                c_n.push(c_nl);
            } // l

            cz.push(c_n);
        } // n
        c.push(cz); // z
    }

    let mut p = Vec::<f64>::with_capacity(species.len() * n_max * n_max * (l_max + 1));

    for zi in 0..species.len() {
        for zj in 0..species.len() {
            for n in 1..=n_max {
                for n_prime in n..=n_max {
                    for l in 0..=l_max {
                        if zi <= zj {
                            let mut sum_m = 0_f64;

                            for m_idx in 0..=2 * l {
                                // 2l + 1 terms

                                sum_m += c[zi][n - 1][l][m_idx] * c[zj][n_prime - 1][l][m_idx];
                                // println!("n = {}, l = {}, m = {}", n, l, (m_idx as i32)-l as i32);
                            } // m

                            p.push(sum_m / ((2 * l + 1) as f64).sqrt())
                        }
                    } // l
                } // n'
            } // n
        }
    }

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
    use crate::geometry::CartesianCoordinate;
    use rand::thread_rng;
    use rand::distributions::{Distribution, Uniform};
 

    fn is_very_close(x: f64, y: f64) -> bool{
        is_close(x, y, 1E-8)
    }
    

    fn is_close(x: f64, y: f64, atol: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        println!("\nleft = {}\nright = {}", x, y);
        (x - y).abs() <= atol
    }

    pub fn power_spectrum_old(structure:      &Structure,
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

        let mut c = vec![];

        for n in 1..=n_max{
            let mut c_n = vec![];

            for l in 0..=l_max{
                let mut c_nl = vec![];
                let l_i32 = l as i32;

                for m_idx in 0..=2*l{   // 2l + 1 terms

                    let m = (m_idx as i32) - l_i32;
                    c_nl.push(c_nlm(&sphr_nbrs, &rbfs, n, l, m, sigma_at));

                }// m
                c_n.push(c_nl);

            }// l

            c.push(c_n);
        }// n


        let mut p = Vec::<f64>::with_capacity(n_max * n_max * (l_max + 1));

        for n in 1..=n_max{
            for n_prime in n..=n_max{
            for l in 0..=l_max{

                let mut sum_m = 0_f64;

                for m_idx in 0..=2*l{   // 2l + 1 terms
                    sum_m += c[n-1][l][m_idx] * c[n_prime-1][l][m_idx]; 
                    //println!("n = {}, l = {}, m = {}", n, l, (m_idx as i32)-l as i32); 
                }// m

                    p.push(sum_m / ((2*l + 1) as f64).sqrt())        

                }// l
            }// n'
        }// n

        p
    }

    fn write_methane_xyz(i: usize){
        std::fs::write(format!("methane{}.xyz", i), 
                       "5\n\n\
                        C     0.00000   0.00000   0.00000\n\
                        H    -0.65860  -0.85220  -0.30120\n\
                        H    -0.45940   0.97110  -0.28590\n\
                        H     0.08440  -0.02940   1.10060\n\
                        H     1.02910  -0.10990  -0.41250\n")
                       .expect("Failed to write methane.xyz!");
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
        
	    write_methane_xyz(0); 
 
        let p = power_spectrum_old(&Structure::from("methane0.xyz"),
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

        std::fs::remove_file("methane0.xyz").expect("Could not remove file!");
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
  
 
        let p = power_spectrum_old(&Structure::from("methane_trns.xyz"),
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
  
 
        let p = power_spectrum_old(&Structure::from("methane_rot.xyz"),
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
  
 
        let p = power_spectrum_old(&Structure::from("methane_disp.xyz"),
                               0, "H", 2, 0, 2.0_f64, 0.5_f64);        
 
        let p_expected = vec![0.07695406_f64, -0.24412914_f64, 0.7744755_f64];
        
        for i in 0..p.len(){
            assert!(! is_close(p[i], p_expected[i], 1E-5));
        }

        std::fs::remove_file("methane_disp.xyz").expect("Could not remove file!");
    }

    fn cnlm_integrand(coord:       CartesianCoordinate,
                      cart_nbrs:   &Vec<CartesianCoordinate>,
                      rbfs:        &RadialBasisFunctions, 
                      n:           usize,
                      l:           usize,
                      m:           i32,
                      a:           f64) -> f64 {
        /*

        c_nlm = ∭  g_n(r) Y_lm(θ, φ) Σ exp(-1/2σ^2 |r - R_i|^2)  dr
               all                   i

        where i sums over neighbours and R_i is the position of the neighbour and
        a = 1 / 2σ^2
        */
        
        let mut g_n = 0_f64;
        let sphr_coord = coord.to_polar();
        let l_i32 = l as i32;

        for n_prime in 0..rbfs.n_max(){
            g_n += rbfs.g[[n-1, l]].beta[n_prime]       // Β_nn'l r^l e^(-αr^2)
                   * sphr_coord.r.powi(l_i32)
                   * (-rbfs.phi[[n_prime, l]].alpha * sphr_coord.r * sphr_coord.r).exp();
        }
        
        let y_lm = sphr_harm(&sphr_coord, l_i32, m);

        let mut sum = 0_f64;
        for i in 0..cart_nbrs.len(){
            sum += (-a * ((coord.x - cart_nbrs[i].x).powi(2)
                          + (coord.y - cart_nbrs[i].y).powi(2)
                          + (coord.z - cart_nbrs[i].z).powi(2))
                   ).exp();
        }

        g_n * y_lm * sum
    }


    fn mc_integral(cart_nbrs:   &Vec<CartesianCoordinate>,
                   rbfs:        &RadialBasisFunctions, 
                   n:           usize,
                   l:           usize,
                   m:           i32,
                   a:           f64) -> f64{
        /* Integrate the triple integral in 3D
        */
        let mut rng = thread_rng();

        let length = 4_f64;
        let v = length.powi(3);       // Uniform around (0, 0, 0)  x in [-l/2, l/2]
        let uniform = Uniform::<f64>::new_inclusive(-length/2_f64, length/2_f64);
        
        let n_points = 100000_usize;
        let mut sum = 0_f64;

        for _ in 0..n_points{
            let coord = CartesianCoordinate{x: uniform.sample(&mut rng),
                                            y: uniform.sample(&mut rng),
                                            z: uniform.sample(&mut rng)};

            sum += cnlm_integrand(coord, cart_nbrs, rbfs, n, l, m, a);
        }

        v * sum / (n_points as f64)
     }

    #[test]
    fn test_numerical_c_n00(){
   
        write_methane_xyz(1);
        let methane = Structure::from("methane1.xyz");
 
        let mut rbfs: RadialBasisFunctions = Default::default();
        rbfs.construct(2, 0, 2.0);
    
        let a = 1_f64 / (2_f64 * 0.5_f64.powi(2));   // 1/2σ^2, σ = 0.5 Å


        for n in vec![1, 2]{
            let analytic_c_100 = c_nlm(&methane.sphr_neighbours(0, "H", 2.0),
                                       &rbfs,
                                       n,
                                       0,
                                       0,
                                       0.5_f64);

            // As the C atom is at the origin the neighbours r - R_i are just the 
            // coordinates
            let numerical_c_100 = mc_integral(&(methane.coordinates[1..].to_vec()),
                                              &rbfs, n, 0, 0, a);

            // Monte Carlo integration is slow, so use a sloppy tolerance
            assert!(is_close(numerical_c_100, analytic_c_100, 0.1));
        }

        std::fs::remove_file("methane1.xyz").expect("Could not remove file!");
    }


    #[test]
    fn test_numerical_c_11m(){
   
        write_methane_xyz(2);
        let methane = Structure::from("methane2.xyz");
 
        let mut rbfs: RadialBasisFunctions = Default::default();
        rbfs.construct(2, 1, 2.0);
    
        let numerical_c11m_s = vec![-0.19955, -0.23227, -0.010524];
        
        for m in vec![-1, 0, 1]{
            let analytic_c_11m = c_nlm(&methane.sphr_neighbours(0, "H", 2.0),
                                       &rbfs, 1, 1, m, 0.5_f64);

            /* Uncomment for true MC evaluation

            let a = 1_f64 / (2_f64 * 0.5_f64.powi(2));   // 1/2σ^2, σ = 0.5 Å
            let numerical_c_11m = mc_integral(&(methane.coordinates[1..].to_vec()),
                                             &rbfs, 1, 1, m, a);
            println!("m = {} I = {}", m, numerical_c_11m);
            */

            // Monte Carlo is slow, so use pre-computed values
            let numerical_c_11m = numerical_c11m_s[(m+1) as usize];

            assert!(is_close(numerical_c_11m, analytic_c_11m, 0.1));
        }

        std::fs::remove_file("methane2.xyz").expect("Could not remove file!");
    }

    
    fn normalised(p: &Vec<f64>) -> Vec<f64>{
        // Normalise: p -> p / √(p.p)
        let mut norm = 0_f64;

        for i in 0..p.len(){
            norm += p[i]*p[i];
        }
        let sqrt_norm = norm.sqrt();

        let mut p_normed = Vec::<f64>::with_capacity(p.len());
        for i in 0..p.len(){
            p_normed.push(p[i].clone()/sqrt_norm);
        }
         
        p_normed
    }


    #[test]
    fn test_soap_kernel(){
        // Test that the power spectrum can be used as a similariry measure
    

        write_methane_xyz(3);
        std::fs::write("methane_long.xyz", 
                       "5\n\n\
                        C    -3.13849        0.66522        0.00000\n\
                        H     1.41112        0.71579        0.00000\n\
                        H    -3.49515       -0.25764       -0.40745\n\
                        H    -3.49515        1.47951       -0.59550\n\
                        H    -3.49516        0.77379        1.00295\n")
                       .expect("Failed to write methane.xyz!");

        std::fs::write("methane_short.xyz", 
                       "5\n\n\
                        C    -3.13849        0.66522        0.00000\n\
                        H    -1.87886        0.67786        0.00000\n\
                        H    -3.49515       -0.25764       -0.40745\n\
                        H    -3.49515        1.47951       -0.59550\n\
                        H    -3.49516        0.77379        1.00295\n")
                       .expect("Failed to write methane.xyz!");
        
        let p_0 = normalised(&power_spectrum_old(&Structure::from("methane3.xyz"),
                                           0, "H", 6, 6, 6.0_f64, 0.5_f64)); 

        let p_1 = normalised(&power_spectrum_old(&Structure::from("methane_short.xyz"),
                                              0, "H", 6, 6, 6.0_f64, 0.5_f64));         
        
        let p_2 = normalised(&power_spectrum_old(&Structure::from("methane_long.xyz"),
                                    0, "H", 6, 6, 6.0_f64, 0.5_f64)); 

        let mut k_01 = 0_f64;
        let mut k_02 = 0_f64;

        for i in 0..p_0.len(){
            k_01 += p_0[i] * p_1[i];
            k_02 += p_0[i] * p_2[i];
        }
         
        let zeta = 4_i32;

        // The short displacement of the H atom should be a strucutre more similar
        // than the larger displacement, thus have a larger similariry (k value)
        assert!(k_01.powi(zeta) > k_02.powi(zeta));
        
        std::fs::remove_file("methane3.xyz").expect("Could not remove file!");
        std::fs::remove_file("methane_short.xyz").expect("Could not remove file!");
        std::fs::remove_file("methane_long.xyz").expect("Could not remove file!");
    }


}
