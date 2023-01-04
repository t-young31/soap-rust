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
#[command(version = "1.0.1")]
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
            6.,
            0.5,
        );
        let p2 = soap::power_spectrum(
            &Structure::from("./nacl_r.xyz"),
            0,
            &vec!["Cl".to_owned(),"Na".to_owned()],
            6,
            6,
            6.,
            0.5,
        );
        let sum1:f64 = p.iter().map(|x|x.abs()).sum();
        let sum2:f64 = p2.iter().map(|x|x.abs()).sum();
        let sum:f64 = p.iter().zip(p2.iter()).map(|(x1,x2)|(x1-x2).abs()).sum();
        println!("sum of soap before ratation p1: {:?}",sum1);
        println!("sum of soap after ratation p2: {:?}",sum2);
        println!("sum of (p1 - p2): {:?}",sum);

        std::fs::remove_file("nacl.xyz").expect("Could not remove file!");
        std::fs::remove_file("nacl_r.xyz").expect("Could not remove file!");
    }

}