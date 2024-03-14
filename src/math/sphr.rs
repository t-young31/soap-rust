use std::f64::consts::PI;

use rgsl::legendre::associated_polynomials::legendre_sphPlm; 
use crate::geometry::SphericalPolarCoordinate;



pub fn sphr_harm(coord: &SphericalPolarCoordinate,
                 l:     i32,
                 m:     i32) -> f64{
    /*

    */
    let x = coord.theta.cos();
    let m_f64 = m as f64;
    let l_f64 = l as f64;
    
    if m == 0 {
        return ((2f64 * l_f64 + 1f64) / (4f64 * PI)).powf(0.5) 
                * rgsl::legendre::associated_polynomials::legendre_Plm(l, m, x)  
    }
    
    let factor = std::f64::consts::SQRT_2 * (-1_f64).powi(m % 2);

    if m < 0 {           // |m|
        return factor * (-m_f64 * coord.phi).sin() * legendre_sphPlm(l, -m, x) 
    }

    else { // m > 0
        return factor * (m_f64 * coord.phi).cos() * legendre_sphPlm(l, m, x)
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

    fn is_very_close(x: f64, y: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        (x - y).abs() <= 1E-8
    }


    #[test]
    fn test_sphr11(){

        let coord = SphericalPolarCoordinate{r: 1.0, 
                                             theta: 0.7853981633974483, // π/4 
                                             phi: 1.0471975511965976};  // π/3
        assert!(is_very_close(sphr_harm(&coord, 1, 1), 
                              0.1727470747356678));

        assert!(is_very_close(sphr_harm(&coord, 1, -1),
                              0.2992067103010745));
    }



    
}
