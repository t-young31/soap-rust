extern crate clap;
use crate::structure::Structure;
use clap::Parser;
mod basis;
mod geometry;
mod integrals;
mod math;
mod soap;
mod structure;

#[derive(Parser)]
#[command(name = "soap-rust")]
#[command(author = "Tom Young")]
#[command(version = "1.1.0")]
#[command(about = "Smooth overlap of atomic positions (SOAP) \
in rust", long_about = None)]
struct Cli {
    #[arg(
        long,
        short,
        help = "The .xyz to calculate a SOAP for. e.g. methane.xyz"
    )]
    filename: String,
    #[arg(
        long,
        short,
        help = "The atom index about which the neighbour\
    density will be expanded. e.g. 0"
    )]
    center: usize,
    #[arg(
        long,
        short,
        help = "The symbols of the elements to consider. e.g. H"
    )]
    species: Vec<String>,
    #[arg(long, default_value_t = 6)]
    nmax: usize,
    #[arg(long, default_value_t = 6)]
    lmax: usize,
    #[arg(
        long,
        short,
        default_value_t = 6.,
        help = "Cut off distance for the SOAP in Å"
    )]
    rcut: f64,
    #[arg(
        long,
        default_value_t = 0.5,
        help = "Smoothness of the neighbour denisty, in Å"
    )]
    sigma: f64,
}

fn main() {
    let cli = Cli::parse();

    let p = soap::power_spectrum(
        &Structure::from(&cli.filename),
        cli.center,
        &cli.species,
        cli.nmax,
        cli.lmax,
        cli.rcut,
        cli.sigma,
    );
    println!("{:?}", p);
}

#[cfg(test)]
mod tests {
    use crate::{soap, structure::Structure};
    use ndarray::array;

    #[test]
    fn test_rotate(){
        // for Test

        std::fs::write(
                format!("nacl.xyz"),
                    "8\n\n\
                        Na       0.00000000       2.82010000       0.00000000\n\
                        Cl       0.00000000       2.82010000       2.82010000\n\
                        Na       0.00000000       0.00000000       2.82010000\n\
                        Cl       0.00000000       0.00000000       0.00000000\n\
                        Na       2.82010000       2.82010000       2.82010000\n\
                        Cl       2.82010000       2.82010000       0.00000000\n\
                        Na       2.82010000       0.00000000       0.00000000\n\
                        Cl       2.82010000       0.00000000       2.82010000\n",
            )
            .expect("Failed to write nacl.xyz!");
        
        std::fs::write( // rotate with 90 degree with axis 0,1,1 
                format!("nacl_r.xyz"),
                    "8\n\n\
                    Na      -1.99411183       1.41005000       1.41005000\n\
                    Cl       0.00000000       2.82010000       2.82010000\n\
                    Na       1.99411183       1.41005000       1.41005000\n\
                    Cl       0.00000000       0.00000000       0.00000000\n\
                    Na       0.00000000       4.81421183       0.82598817\n\
                    Cl      -1.99411183       3.40416183      -0.58406183\n\
                    Na       0.00000000       1.99411183      -1.99411183\n\
                    Cl       1.99411183       3.40416183      -0.58406183\n",
            )
            .expect("Failed to write nacl_r.xyz!");

        let p = soap::power_spectrum(
            &Structure::from("./nacl.xyz"),
            0,
            &vec!["Cl".to_owned(),"Na".to_owned()],
            6,
            6,
            3.,
            0.5,
        );
        let p2 = soap::power_spectrum(
            &Structure::from("./nacl_r.xyz"),
            0,
            &vec!["Cl".to_owned(),"Na".to_owned()],
            6,
            6,
            3.,
            0.5,
        );
        let sum1:f64 = p.iter().map(|x|x.abs()).sum();
        let sum2:f64 = p2.iter().map(|x|x.abs()).sum();
        let sum:f64 = p.iter().zip(p2.iter()).map(|(x1,x2)|(x1-x2).abs()).sum();

        std::fs::remove_file("nacl.xyz").expect("Could not remove file!");
        std::fs::remove_file("nacl_r.xyz").expect("Could not remove file!");
        assert!(sum / sum1 < 1E-6);
        assert!(sum / sum2 < 1E-6);

    }


    #[test]
    fn test_rotate_2(){
        // A more complicated test, for rotating the convenrional cell of a zeolite structure called ABW
        let xyz_path = "ABW.xyz";
        std::fs::write( 
            format!("{}", xyz_path),
                "36\n\n\
                 O    4.936500    1.313500    3.089671\n\
                 O    4.936500    3.940500    5.680329\n\
                 O    0.000000    3.940500    7.474671\n\
                 O    9.873000    3.940500    7.474671\n\
                 O    0.000000    1.313500    1.295329\n\
                 O    9.873000    1.313500    1.295329\n\
                 O    3.060630    0.000000    4.385000\n\
                 O    3.060630    5.254000    4.385000\n\
                 O    6.812370    0.000000    4.385000\n\
                 O    6.812370    5.254000    4.385000\n\
                 O    6.812370    2.627000    4.385000\n\
                 O    3.060630    2.627000    4.385000\n\
                 O    7.997130    2.627000    0.000000\n\
                 O    7.997130    2.627000    8.770000\n\
                 O    1.875870    2.627000    0.000000\n\
                 O    1.875870    2.627000    8.770000\n\
                 O    1.875870    0.000000    0.000000\n\
                 O    1.875870    0.000000    8.770000\n\
                 O    1.875870    5.254000    0.000000\n\
                 O    1.875870    5.254000    8.770000\n\
                 O    7.997130    0.000000    0.000000\n\
                 O    7.997130    0.000000    8.770000\n\
                 O    7.997130    5.254000    0.000000\n\
                 O    7.997130    5.254000    8.770000\n\
                 O    2.468250    1.313500    2.192500\n\
                 O    7.404750    3.940500    6.577500\n\
                 O    7.404750    1.313500    2.192500\n\
                 O    2.468250    3.940500    6.577500\n\
                Si    3.384464    1.313500    3.514139\n\
                Si    6.488536    3.940500    5.255861\n\
                Si    6.488536    1.313500    3.514139\n\
                Si    3.384464    3.940500    5.255861\n\
                Si    8.320964    3.940500    7.899139\n\
                Si    1.552036    1.313500    0.870861\n\
                Si    1.552036    3.940500    7.899139\n\
                Si    8.320964    1.313500    0.870861\n\
                ",
        )
        .expect("Failed to write ABW.xyz!");
        let structure = Structure::from(xyz_path);
        let orth_matrix = array![
        [1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0]
        ]; // rotate 90 degree along x axis
        let structure_r = structure.rotate_coords(&orth_matrix);
        let p = soap::power_spectrum(
            &structure,
            0,
            &vec!["O".to_owned(),"Si".to_owned()],
            6,
            6,
            3_f64,
            0.5,
        );
        let p2 = soap::power_spectrum(
            &structure_r,
            0,
            &vec!["O".to_owned(),"Si".to_owned()],
            6,
            6,
            3_f64,
            0.5,
        );
        let sum1:f64 = p.iter().map(|x|x.abs()).sum();
        let sum2:f64 = p2.iter().map(|x|x.abs()).sum();
        let sum:f64 = p.iter().zip(p2.iter()).map(|(x1,x2)|(x1-x2).abs()).sum();

        std::fs::remove_file(xyz_path).expect("Could not remove file!");
        println!("{:?}", sum1);
        println!("{:?}", sum2);
        println!("{:?}", sum);
        assert!(sum / sum1 < 1E-6);
        assert!(sum / sum2 < 1E-6);
    }
}
