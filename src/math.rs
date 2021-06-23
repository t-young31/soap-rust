extern crate rgsl;   // GNU Scientific library


pub fn factorial(n: i32) -> i32{
    // Compute n! for n as an integer
    (1..=n).product()
}


pub fn factorial_f(n: i32) -> f64{
    // Compute n! for n as an floating point number
    factorial(n) as f64
}


pub fn msbf_first(n: i32, x: f64) -> f64{
    /* Modified spherical Bessel function of the first kind. See:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_in.html#scipy.special.spherical_in

    i_n = √(π/2x) I_(n+1/2)(x)

    Due to the way it's defined in GSL and wrapped the scaled value
    needs to be unscaled...
    */

    rgsl::bessel::il_scaled(n, x) / (-x.abs()).exp()
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

    #[test]
    fn test_factorial(){
        assert_eq!(factorial(0), 1);
        assert_eq!(factorial(1), 1);
        assert_eq!(factorial(3), 6);
    }

} // Tests

