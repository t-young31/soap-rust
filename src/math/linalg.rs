/* Linear algebra routines, based on:

[1] https://en.wikipedia.org/wiki/Gaussian_elimination
[2]

*/
use ndarray::{Array2, array};



fn inverse_gauss_jordan(m: Array2::<f64>) -> Array2::<f64>{
    /* Use the Gauss Jordan method to find the inverse of 
    a square matrix m 


    Arguments:
        m (ndarray): Matrix

    Returns:
        (ndarray): m^-1
    */
    let (n_rows, n_cols) = (m.shape()[0], m.shape()[1]);

    if n_rows != n_cols{
        panic!("Expecting a square matrix, but {} != {}", n_cols, n_cols);
    }

    let mut tmp_m = Array2::<f64>::zeros((n_rows, 2 * n_cols));

    /* Initial augmentation on the RHS diagonal
            
            (0  0  0  0 ..)         (0  0  1  0 ..)
    tmp_m = (0  0  0  0 ..)   -->   (0  0  0  1 ..)
            (.  .  .  .   )         (.  .  .  . ..)
    */
    for i in 0..n_rows{
        tmp_m[[i, i+n_rows]] = 1.0;
    }

    // Set the LHS as the matrix to invert
    for i in 0..n_rows{

        if is_close(m[[i, i]], 0.0, 1E-8){
            panic!("Cannot invert a matrix with a zero on the diagonal");
        }

        for j in 0..n_cols{
            tmp_m[[i, j]] = m[[i, j]].clone();
        }// j
    }// i

    // And apply the Gauss-Jordan method O(n^3)
    for i in 0..n_rows{
        for j in 0..n_cols{

            if i == j{
                continue;
            }

            let ratio = tmp_m[[j, i]] / tmp_m[[i, i]];

            for k in 0..2*n_cols{
                tmp_m[[j, k]] -= ratio * tmp_m[[i, k]];
            }// k
        }// j
    }// i

    for i in 0..n_rows{
        let tmp_m_ii = tmp_m[[i, i]].clone();

        for j in 0..2*n_cols{
            tmp_m[[i, j]] = tmp_m[[i, j]] / tmp_m_ii;
        }// j
    }// i

    // Finally set the elements of the inverse matrix
    let mut inv = Array2::<f64>::zeros((n_rows, n_cols));

    for i in 0..n_rows{
        for j in n_cols..2*n_cols{
            inv[[i, j-n_cols]] = tmp_m[[i, j]];
        }// j
    }// i

    inv
}


fn is_close(x: f64, y: f64, atol: f64) -> bool{
    // Are two numbers close to within an absolute tolerance? 
    (x - y).abs() <= atol
}


fn f_norm(m: Array2::<f64>) -> f64{
    /*
    Calculate the Frobenius norm of a matrix, defined as per
    https://mathworld.wolfram.com/FrobeniusNorm.html

    ||m|| = √(Σ Σ |a_ij|^2)
              i j
    */
    let mut norm_sq = 0.0_f64;

    for i in 0..m.shape()[0]{
        for j in 0..m.shape()[1]{
            
            norm_sq += m[[i, j]].powi(2);
        }
    }
    norm_sq.sqrt()
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

    fn matrix_is_close(a: Array2::<f64>, b: Array2::<f64>) -> bool{
        // Are two matricies close to within an tolerance on the Frobenius norm?
        f_norm(a - b) < 1E-8
    }


    #[test]
    #[should_panic]
    fn test_non_sq_matrix(){
        // Cannot invert a non-square matrix using Gauss Jordan
        let _ = inverse_gauss_jordan(array![[0., 0.]]);
    }


    #[test]
    #[should_panic]
    fn test_diagonal_zeros(){
        // Cannot invert a zero matrix using Gauss Jordan
        let _ = inverse_gauss_jordan(array![[0., 0.],
                                            [0., 0.]]);

        // or one with a zero diagonal element
        let _ = inverse_gauss_jordan(array![[0., 1.],
                                            [2., 3.]]);
    }


    #[test]
    fn test_frob_norm(){
        // Test ||m|| of some simple matrices

        let m = array![[1.0, 0.0],
                         [0.0, 1.0]];
        assert!(is_close(f_norm(m), 2.0_f64.sqrt(), 1E-8));

        let m = array![[3.0, 4.0],
                      [0.0, 0.0]];
        assert!(is_close(f_norm(m), 5.0, 1E-8));

        // they need not be square
        let m = array![[3.0, 4.0]];
        assert!(is_close(f_norm(m), 5.0, 1E-8));
    }


    #[test]
    fn test_inv_diagonal_gj(){
        // Test the inverse of some diagonal matrices, which should just be
        // the inverse of their diagonal elements

        let m = array![[1.0, 0.0],
                       [0.0, 1.0]];

        assert!(matrix_is_close(inverse_gauss_jordan(m),
                                 array![[1.0, 0.0],
                                        [0.0, 1.0]]));

        // Check a selection of deterministic diagonal matrices
        for n in 1..4{

            for dim in 1..10{
                let mut m = Array2::<f64>::zeros((dim, dim));
                let mut expected_inv =  Array2::<f64>::zeros((dim, dim));

                for i in 0..dim{
                    m[[i, i]] = n.clone() as f64;
                    expected_inv[[i, i]] = 1_f64 / (n.clone() as f64);
                }

                assert!(matrix_is_close(inverse_gauss_jordan(m),
                                        expected_inv));
            }

        }// n

    }


    #[test]
    fn test_non_diagonal_gj(){
        // Simple examples calculated using numpy

        let m = array![[1.0, 2.0],
                       [3.0, 4.0]];

        assert!(matrix_is_close(inverse_gauss_jordan(m),
                array![[-2.0, 1.0],
                       [1.5, -0.5]]));

        // Random normal 4x4
        let m = array![[0.0887191769245145, 1.103703023166563, -0.9660151847773749, 0.3489494627314538], 
                       [2.775694781834193, 0.23875614627475789, 0.10486181505063653, -0.1816924313770859], 
                       [2.2548197390169356, -0.5776152654289637, -0.7700473235450745, -0.3230021800759852], 
                       [0.6224882824513317, 0.42422899069884973, 2.3384581631868406, 0.13358277318643444]];

        let expected_inv = array![[0.39393344045306977, -0.4937778002938748, 0.9027323710244478, 0.48214280057426767], 
                                  [-0.8373094224518616, 2.9939031039729396, -3.2266771294796954, -1.5426805038913658], 
                                  [-0.22642457170714925, 0.10540465336540064, -0.1941260062687446, 0.26544468507994606], 
                                  [4.787117030455704, -9.052169980358372, 9.438830866943361, 5.491653160900342]];

        assert!(matrix_is_close(inverse_gauss_jordan(m), expected_inv));
    }

}


