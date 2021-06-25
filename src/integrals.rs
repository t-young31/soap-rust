use crate::math::misc::factorial_f;
use crate::math::msbf;


#[derive(Default)]
pub struct IntegrandArguments{
    pub integer_args: Vec<i32>,
    pub float_args: Vec<f64>,
}


pub fn trapz(integrand: &dyn Fn(f64, &IntegrandArguments) -> f64,
             a:         f64,
             b:         f64,
             args:      &IntegrandArguments,
             dx:        f64) -> f64{
    /* Integrate a function (the integrand) from a->b using the trapezium rule

    Arguments:
        integrand (function):

        a (float): Lower bound of the integral

        b (float): Upper bound of the integral

        args (IntegrandArguments): Arguments 

    -------------------------------------------------------------------
        Example:
                        1
                        ∫ x^2 dx   = 1/3
                        0  

        Define a function and arguments as:

        fn x_sq(x: f64) -> f64{
            x*x
        }

        let args = IntegrandArguments{integer_args: vec![], float_args:   vec![]};

        println!("{}", trapz(x_sq, 0.0, 1.0, &args))  ----> 0.3333..
    -------------------------------------------------------------------
    */
    let n = ((b - a).abs() / dx).floor() as u32;              // Number of points

    let mut integral: f64 = integrand(a, args) / 2_f64;       // First: f(a)/2
    let mut x: f64 = a;

    for _ in 1..n {                                           // Intermediate: f(x_n)
        x += dx;
        integral += integrand(x, args);

    } 
    integral += integrand(b, args) / 2_f64;                   // First: f(b)/2

    integral * dx
}


pub fn gaussian_2m(m: i32, k: f64) -> f64{
    /*
    Calculate a gaussian integral of the form below, from wikipedia:

        ∞                             ((2m - 1)!!
        ∫ dr r^(2m) e^(-kr^2) =  -----------------------
        0                       (k^m 2^(m+1)))(pi / k)^1/2

    Arguments:
        m (int): Half of the r exponent

        k (float): Gaussian exponent

    Returns:
        (float): Value of the integral
    */
    let d_factorial = factorial_f(2_i32*m) / (2_f64.powi(m) * factorial_f(m));

    d_factorial * (std::f64::consts::PI / k).sqrt() / (k.powi(m) * 2_f64.powi(m+1))
}


pub fn xn_2_gaussian_msbf(n: i32, a: f64, b: f64) -> f64{
    /*
    Calculate the value of

    ∞
     ∫ dx x^(n+2) e^(-ax^2) i_n(bx)
    0

    where i_n is the modified spherical Bessel function (MSBF) of the first kind
    */

    (std::f64::consts::PI/2_f64).sqrt()
    * (b.powi(n) / (2_f64 * a).powf((n as f64) + 1.5))
    * (b*b/(4_f64 * a)).exp()
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
    use std::f64::consts::*;

    fn is_close(x: f64, y: f64, atol: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        (x - y).abs() <= atol
    }

    #[test]
    fn test_gaussian_2m(){

        // Wolfram: integral 0 to +oo exp(-x^2) dx
        assert!(is_close(gaussian_2m(0, 1.0), PI.sqrt()/2_f64, 1E-8));

        // Wolfram: integral 0 to +oo x^2 exp(-x^2) dx
        assert!(is_close(gaussian_2m(1, 1.0), PI.sqrt()/4_f64, 1E-8));

        // Wolfram: integral 0 to +oo x^6*exp(-3*x^2) dx
        assert!(is_close(gaussian_2m(3, 3.0), 5_f64*(PI/3_f64).sqrt()/144_f64, 1E-8));
    }


    #[test]
    fn test_trapz_x_sq(){
        // Test that the trapezium integrator can integrate x^2

        fn x_sq(x: f64, _args: &IntegrandArguments) -> f64{
            x * x
        }

        let args: IntegrandArguments = Default::default();

        assert!(is_close(trapz(&x_sq, 0.0, 1.0, &args, 1E-6), 1_f64/3_f64, 1E-4));
    }


    #[test]
    fn test_trapz_gaussian(){
        // Test the integration of a Gaussian
        // Wolfram: integral from -oo to +oo exp(-a*(x+b)^2) 

        fn gaussian(x: f64, args: &IntegrandArguments) -> f64{
            let a = args.float_args[0];
            let b = args.float_args[1];

            (-a*(x + b).powi(2)).exp()    // exp(-a*(x+b)^2) 
        }

        let args = IntegrandArguments{integer_args: vec![],
                                      float_args:   vec![2.0, -3.0]};
        
        let integral = trapz(&gaussian, 0.0, 6.0, &args, 1E-4);
        let analytic_integral = (PI / args.float_args[0]).sqrt();

        assert!(is_close(integral, analytic_integral, 1E-4));
    }


    #[test]
    fn test_trapz_gaussian_x_sq(){
        // Test the integration of a polynomial multiplied by a Gaussian
        // Wolfram: integral from 0 to +oo x^2 exp(-x^2) 

        fn x_sq_gaussian(x: f64, args: &IntegrandArguments) -> f64{
            let n = args.integer_args[0];

            x.powi(n) * (-x.powi(2)).exp()
        }

        let args = IntegrandArguments{integer_args: vec![2],
                                      float_args:   vec![]};
        
        let integral = trapz(&x_sq_gaussian, 0.0, 6.0, &args, 1E-4);
        let analytic_integral = PI.sqrt() / 4_f64;

        assert!(is_close(integral, analytic_integral, 1E-6));
    }


    #[test]
    fn test_xn_2_gaussian_msbf(){
        // Test xn_2_gaussian_msbf using trapezium integration

        fn integrand(x: f64, args: &IntegrandArguments) -> f64{
            // x^(n+2) e^(-ax^2) i_n(bx)
            let n = args.integer_args[0];
            let a = args.float_args[0];
            let b = args.float_args[1];

            x.powi(n+2) * (-a*(x*x)).exp() * msbf::i_n(n, b*x)
        }

        // Check multiple values of n for the spherical bessel function
        for n in 0..4{
            let args = IntegrandArguments{integer_args: vec![n],
                                          float_args:   vec![1.5, 2.2]};

            let analytic_integral = xn_2_gaussian_msbf(n, 1.5, 2.2);
            let numeric_integral = trapz(&integrand, 0.0, 5.0, &args, 1E-4);

            assert!(is_close(analytic_integral, numeric_integral, 1E-6));
        }

    }


}