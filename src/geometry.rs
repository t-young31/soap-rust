#[derive(Default, Clone)]
pub struct CartesianCoordinate{
    pub x: f64,
    pub y: f64,
    pub z: f64,
}


impl CartesianCoordinate{

    fn sq_distance(&self) -> f64{
        // Square distance from the origin (r^2 = x^2 + y^2 + z^2)
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    pub fn to_polar(&self) -> SphericalPolarCoordinate{
        // Convert this Cartesian coordinate into spherical polars

        let r: f64 = self.sq_distance().sqrt();

        SphericalPolarCoordinate{r:     r,
                                 phi: self.y.atan2(self.x),
                                 theta:   (self.z / r).acos()}
    }

    pub fn shift_then_to_polar(&self, 
                               origin: &CartesianCoordinate) -> SphericalPolarCoordinate{
        /* Convert this coordinate to spherical polars with a defined
        origin (defaults to (0, 0, 0))
        */
        CartesianCoordinate{x: &self.x - origin.x,
                            y: &self.y - origin.y,
                            z: &self.z - origin.z}.to_polar()

    }

}

#[derive(Clone)]
pub struct SphericalPolarCoordinate{
    /* Spherical polar coordinate, with notation taken from 
    https://mathworld.wolfram.com/SphericalCoordinates.html
        
    
              ^  φ
              |---- /
            z |   / 
              | / 
              ----------> y
               \
                \> x 
    
    r: Distance 
    θ: Azimuthal angle range = [0, 2π]
    φ: Polar angle range = [0, π]

    */
    pub r:     f64,
    pub theta: f64,
    pub phi:   f64,
}

impl SphericalPolarCoordinate{
    pub fn to_cartesian(&self) -> CartesianCoordinate{
        // Convert this spherical polar coordinate into Cartesian
        CartesianCoordinate{x: self.r * self.theta.sin() * self.phi.cos(),
                            y: self.r * self.theta.sin() * self.phi.sin(),
                            z: self.r * self.theta.cos()}
    }
    pub fn distance(&self, other: &SphericalPolarCoordinate) -> f64{
        // Distance between two points in spherical polars
        ((self.r*other.r).powi(2) + 
         (self.r*other.r).powi(2) - 
         2.0*self.r*other.r*(self.theta*other.theta).cos()).sqrt()
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

    // Inherit all the functions in the module above this one
    use super::*;
    use std::f64::consts::*;

    fn is_close(x: f64, y: f64, atol: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        (x - y).abs() <= atol
    }

    #[test]
    fn test_geometry_init(){
        // Check initialisation of coordinate
        assert_eq!(CartesianCoordinate{x: 0.0, y: 0.0, z: 0.0}.x, 0.0);
        assert_eq!(CartesianCoordinate{x: 0.0, y: 0.0, z: 0.0}.y, 0.0);
        assert_eq!(CartesianCoordinate{x: 0.0, y: 0.0, z: 0.0}.x, 0.0);
    }

    #[test]
    fn test_geometry_cart_to_spherical(){
        /* Check that the conversion from Cartesian (3D) to spherical
            polars coordinates is close to what it should be
        */

        // Point aligned along the z-axis
        let coordinate = CartesianCoordinate{x: 0.0, y: 0.0, z: 1.0 }.to_polar();
        
        assert!(is_close(coordinate.r, 1.0, 1E-8));
        assert!(is_close(coordinate.phi, 0.0, 1E-8));
        assert!(is_close(coordinate.theta, 0.0, 1E-8));

        // non-unity distance
        let coordinate = CartesianCoordinate{x: 1.0, y: 1.0, z: 0.0 }.to_polar();
        assert!(is_close(coordinate.r, 1.4, 0.1));

        // non-zero angles
        let coordinate = CartesianCoordinate{x: 1.0, y: 0.0, z: 0.0 }.to_polar();
        assert!(is_close(coordinate.r, 1.0, 1E-8));
        assert!(is_close(coordinate.phi, 0.0, 1E-8));
        assert!(is_close(coordinate.theta, FRAC_PI_2, 1E-8));

        // non-positive quadrant
        let coordinate = CartesianCoordinate{x: -1.0, y: 0.0, z: 0.0}.to_polar();
        assert!(is_close(coordinate.r, 1.0, 1E-8));
        assert!(is_close(coordinate.phi.abs(), PI, 1E-8));
        assert!(is_close(coordinate.theta, FRAC_PI_2, 1E-8));
    }

}
