/*
Generate a structure from .xyz files
*/
use crate::geometry::CartesianCoordinate;
use std::fs::File;
use std::io::prelude::*;


#[derive(Default)]
struct Structure{

    pub coordinates: Vec<CartesianCoordinate>,
    pub atomic_symbols: Vec<String>,
}


impl Structure{

    pub fn from(xyz_filename: String) -> Self{
        /* 
        Construct a strucutre from a standard .xyz file
        */
        let mut structure: Structure = Default::default();

        let file = File::open(xyz_filename.clone()).expect("File not found");
        let reader = std::io::BufReader::new(file);

        let mut line_n = 0_i32;

        for line in reader.lines() {

            let _line = &line.unwrap();
            let items: Vec<&str> = _line.split(" ").collect();

            if line_n == 0{
                let n_atoms: i32 = items[0].parse().unwrap();
                println!("Had {} atoms", n_atoms);

                line_n += 1;
                continue;
            }

            if line_n == 1{
                // Should be the title line

                line_n += 1;
                continue;
            }

            
            if items.len() == 0{
                // Assume the whole file has been parsed
                break;
            }


            structure.atomic_symbols.push(items[0].clone().to_string());
            
            let coord = CartesianCoordinate{x: items[1].clone().parse().unwrap(),
                                            y: items[2].clone().parse().unwrap(),
                                            z: items[3].clone().parse().unwrap()};

            structure.coordinates.push(coord);

            line_n += 1;
        }
        
        structure
    }


    pub fn n_atoms(&self) -> usize{
        self.coordinates.len()
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
    use std::env;
    use std::path::Path;


    fn is_very_close(x: f64, y: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        (x - y).abs() <= 1E-10
    }


    #[test]
    fn test_simple_construction(){

        let path = env!("CARGO_MANIFEST_DIR").to_owned() + "/tests/h2.xyz";
        assert!(std::path::Path::is_file(Path::new(&path)));

        let strucutre = Structure::from(path);

        assert_eq!(strucutre.n_atoms(), 2);

        assert!(is_very_close(strucutre.coordinates[0].x, 0.0));
        assert!(is_very_close(strucutre.coordinates[0].y, 0.0));
        assert!(is_very_close(strucutre.coordinates[0].z, 0.0));

        assert!(is_very_close(strucutre.coordinates[1].x, 1.0));
        assert!(is_very_close(strucutre.coordinates[1].y, 0.0));
        assert!(is_very_close(strucutre.coordinates[1].z, 0.0));
    }

}

