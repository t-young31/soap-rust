extern crate rgsl;   // GNU Scientific library


pub fn i_n(n: i32, x: f64) -> f64{
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

    fn assert_is_close(x: f64, y: f64, atol: f64){
        // Are two numbers close to within an absolute tolerance? 
        assert!((x - y).abs() <= atol, "{} != {}", x, y)
    }


    #[test]
    fn test_in(){
        /* Generated using scipy with 

        from scipy.special import spherical_in

        for n in range(10):
            print(spherical_in(n=n, z=1.1))
        */

        let scipy_in = vec![1.2142249728401608,
                             0.4129941645291778,
                             0.08787725139694845,
                             0.013552112724866856,
                             0.0016365340568866291,
                             0.00016228862306716564,
                             1.364782621497245E-05,
                             9.961314356731056E-07,
                             6.421572852101202E-08,
                             3.7065403483745103E-09];
        
        for n in 0..10{
            assert_is_close(i_n(n, 1.1), scipy_in[n as usize], 1E-8);
        }

    }

} // Tests





