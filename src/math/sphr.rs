use rgsl::legendre::associated_polynomials::legendre_sphPlm; 
use crate::geometry::SphericalPolarCoordinate;



pub fn sphr_harm(coord: &SphericalPolarCoordinate,
                 l:     i32,
                 m:     i32) -> f64{
    /*

    */
    if m == 0 {
        return 0.2820947917738781434740_f64   // 1/2 √(1/π)
    }
    
    let x = coord.theta.cos();
    let m_float = m as f64;
    let factor = std::f64::consts::SQRT_2 * (-1_f64).powi(m % 2);

    if m < 0 {           // |m|
        return factor * (-m_float * coord.phi).sin() * legendre_sphPlm(l, -m, x) 
    }

    else { // m > 0
        return factor * (m_float * coord.phi).cos() * legendre_sphPlm(l, m, x)
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
