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
        Construct a strucutre from a standard .xyz file with the form

	2
 	Title line
	H  0.0  0.0  0.0
	H  0.0  0.0  0.0
   
	where there may be any amount of whitespace between the 
        coordinates and the first item in the file should be the number
        of total atoms in the file.
        */
        let mut structure: Structure = Default::default();

        let file = File::open(xyz_filename).expect("File not found");
        let reader = std::io::BufReader::new(file);

        let mut line_n = 0_i32;
        let mut n_atoms: usize = 0;

        for line in reader.lines() {

            let _line = &line.unwrap();
            let items: Vec<&str> = _line.split(" ").collect();

            if line_n == 0{
                n_atoms = items[0].parse().unwrap();
                println!("Had {} atoms", n_atoms);

                line_n += 1;
                continue;
            }

            if line_n == 1{
                // Should be the title line

                line_n += 1;
                continue;
            }

            
            if items.len() == 0 || (items.len() == 1 && items[0] == ""){
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

        if structure.coordinates.len() != n_atoms{
	    panic!("Number of declared atoms not equal that found")
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
    use std::path::Path;


    fn is_very_close(x: f64, y: f64) -> bool{
        // Are two numbers close to within an absolute tolerance? 
        (x - y).abs() <= 1E-10
    }


    fn delete_when_exists(filename: &str){
       // Pole for a file existing and quit if it does not exist after
       // 0.5 s

       for _ in 0..50{
           std::thread::sleep(std::time::Duration::from_millis(10));
           if Path::new(filename).exists(){
               std::fs::remove_file(filename).expect("Could not remove file!");
               break; 
           } 
       }
    }


    #[test]
    fn test_simple_construction(){

	std::fs::write("h2.xyz", 
                       "2\n\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0")
                       .expect("Failed to write the test file!");

        let strucutre = Structure::from("h2.xyz".to_string());

        assert_eq!(strucutre.n_atoms(), 2);

        assert!(is_very_close(strucutre.coordinates[0].x, 0.0));
        assert!(is_very_close(strucutre.coordinates[0].y, 0.0));
        assert!(is_very_close(strucutre.coordinates[0].z, 0.0));

        assert!(is_very_close(strucutre.coordinates[1].x, 1.0));
        assert!(is_very_close(strucutre.coordinates[1].y, 0.0));
        assert!(is_very_close(strucutre.coordinates[1].z, 0.0));
    
	std::fs::remove_file("h2.xyz").expect("Could not remove file!"); 

    }

    
    
    #[test]
    fn test_correct_format_with_blank_lines(){

	std::fs::write("h2_blank_lines.xyz", 
                       "2\n\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0\n\n")
                       .expect("Failed to write the test file!");

        let strucutre = Structure::from("h2_blank_lines.xyz".to_string());

        assert_eq!(strucutre.n_atoms(), 2);
	std::fs::remove_file("h2_blank_lines.xyz").expect("Could not remove file!"); 

    }


    #[test]
    #[should_panic]
    fn test_incorrect_format_n_atoms(){

	std::fs::write("h2_broken.xyz", 
                       "1\n\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0")
                       .expect("Failed to write the test file!");

        let _ = Structure::from("h2_broken.xyz".to_string());
    }
   
    
    #[test]
    #[should_panic]
    fn test_incorrect_format_not_enough_coords(){

	std::fs::write("h2_broken2.xyz", 
                       "2\n\nH 0.0 0.0\nH 1.0 0.0 0.0")
                       .expect("Failed to write the test file!");

        let _ = Structure::from("h2_broken2.xyz".to_string());
    }


   #[test]
   fn cleanup_test(){

       delete_when_exists("h2_broken.xyz");
       delete_when_exists("h2_broken2.xyz");
   } 

}

