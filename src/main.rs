extern crate clap;
use clap::{Arg, App, value_t};
use crate::structure::Structure;
mod integrals;
mod math;
mod geometry;
mod basis;
mod structure;
mod soap;


fn main() {

 let args = App::new("soap-rust")
                        .version("1.0.0a")
                        .author("Tom Young")
                        .about("Smooth overlap of atomic positions (SOAP) \
                                in rust")
                        .arg(Arg::with_name("xyz_filename")
                             .short("file")
                             .long("xyz_filename")
                             .required(true)
                             .help("The .xyz to calculate a SOAP for. e.g. methane.xyz")
                             .index(1))
                        .arg(Arg::with_name("atom_idx")
                              .short("i")
                              .long("atom_idx")
                              .required(true)
                              .help("The atom index about which the neighbour\
                                     density will be expanded. e.g. 0")
                              .index(2))
                        .arg(Arg::with_name("nbr")
                              .short("e")
                              .long("neighbour_element")
                              .required(true)
                              .help("The symbol of the elements neihgbouring\
                                     that defined by atom_idx to consider. e.g. H")
                              .index(3))
                        .arg(Arg::with_name("n_max")
                              .short("n")
                              .long("n_max")
                              .required(false)
                              .default_value("6"))
                        .arg(Arg::with_name("l_max")
                              .short("l")
                              .long("l_max")
                              .required(false)
                              .default_value("6"))
                        .arg(Arg::with_name("r_cut")
                              .short("r")
                              .long("r_cut")
                              .required(false)
                              .default_value("2")
                              .help("Cut off distance for the SOAP in Å"))
                        .arg(Arg::with_name("sigma_at")
                              .short("s")
                              .long("sigma_at")
                              .required(false)
                              .default_value("0.5")
                              .help("Smoothness of the neihgbour denisty, in Å"))
                        .get_matches();

   let p = soap::power_spectrum(&Structure::from(args.value_of("xyz_filename").unwrap()),
                                value_t!(args, "atom_idx", usize).unwrap_or_else(|e| e.exit()),
                                args.value_of("nbr").unwrap(),                
                                value_t!(args, "n_max", usize).unwrap_or_else(|e| e.exit()),  
                                value_t!(args, "l_max", usize).unwrap_or_else(|e| e.exit()), 
                                value_t!(args, "r_cut", f64).unwrap_or_else(|e| e.exit()), 
                                value_t!(args, "sigma_at", f64).unwrap_or_else(|e| e.exit()));    

    println!("{:?}", p);
}
