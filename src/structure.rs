/*
Generate a structure from .xyz files
*/
use crate::geometry::CartesianCoordinate;
use crate::geometry::SphericalPolarCoordinate;
use std::fs::File;
use std::io::prelude::*;


#[derive(Default)]
pub struct Structure{

    pub coordinates: Vec<CartesianCoordinate>,
    pub atomic_symbols: Vec<String>,
}


impl Structure{

    pub fn from(xyz_filename: &str) -> Self{
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

        let file = File::open(xyz_filename.to_string()).expect("File not found");
        let reader = std::io::BufReader::new(file);

        let mut line_n = 0_i32;
        let mut n_atoms: usize = 0;

        for line in reader.lines() {

            let _line = &line.unwrap();
            let items: Vec<&str> = _line.split_whitespace().collect();

            if line_n == 0{
                n_atoms = items[0].parse().expect("Failed to read n_atoms");
                // println!("Had {} atoms", n_atoms);

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
            
            // NOTE: Clone is delibrerate for array positioning
            let coord = CartesianCoordinate{x: items[1].clone().parse().expect("Failed to read x"),
                                            y: items[2].clone().parse().expect("Failed to read y"),
                                            z: items[3].clone().parse().expect("Failed to read z"),};

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


    pub fn sphr_neighbours(&self,
                           atom_idx:    usize, 
                           nbr_element: &str,
                           r_cut:       f64) -> Vec<SphericalPolarCoordinate>{
        /*
        Generate a vector of spherical coordinates around a single atom

        Arguments:
            atom_idx (int):

            nbr_element (str):

            r_cut (float): Cut-off distance (Å)
        
        Returns:
            (vector(coord)):
        */
        if atom_idx > self.n_atoms(){
            panic!("Cannot expand the neighbours about an atom not the structure");
        }
        
        let mut coords = Vec::<SphericalPolarCoordinate>::with_capacity(self.n_atoms());
        let origin = &self.coordinates[atom_idx];        

        for i in 0..self.n_atoms(){

            if (i == atom_idx) | (self.atomic_symbols[i] != nbr_element){
                // Only consider the neighbours of the correct element
                continue;
            }
            let sphr_coord = self.coordinates[i].shift_then_to_polar(origin);
            
            if sphr_coord.r < r_cut{
                // Only add neighbours within the cut-off distance
                coords.push(sphr_coord);
            }
        }
        
        coords
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

    
    fn write_h2(){
        // Write a correct xyz file for dihydrogen
	    std::fs::write("h2.xyz", 
                       "2\n\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0")
                       .expect("Failed to write the test file!");
    }
    

    #[test]
    fn test_simple_construction(){

        write_h2();
        let strucutre = Structure::from("h2.xyz");

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

        let strucutre = Structure::from("h2_blank_lines.xyz");

        assert_eq!(strucutre.n_atoms(), 2);
	    std::fs::remove_file("h2_blank_lines.xyz").expect("Could not remove file!"); 

    }


    #[test]
    #[should_panic]
    fn test_incorrect_format_n_atoms(){

	    std::fs::write("h2_broken.xyz", 
                       "1\n\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0")
                       .expect("Failed to write the test file!");

        let _ = Structure::from("h2_broken.xyz");
    }
   
    
    #[test]
    #[should_panic]
    fn test_incorrect_format_not_enough_coords(){

	    std::fs::write("h2_broken2.xyz", 
                       "2\n\nH 0.0 0.0\nH 1.0 0.0 0.0")
                       .expect("Failed to write the test file!");

        let _ = Structure::from("h2_broken2.xyz");
    }


   #[test]
   fn cleanup_test(){

       delete_when_exists("h2_broken.xyz");
       delete_when_exists("h2_broken2.xyz");
   } 


   #[test]
   #[should_panic]
   fn sphr_neighbours_invalid_idx(){
        // Check that neighbours cannot be defined for a atom not in the set
        write_h2();
        let structure = Structure::from("h2.xyz");
        let _ = structure.sphr_neighbours(3, "H", 2.9);
        std::fs::remove_file("h2.xyz").expect("Could not remove file!"); 
    }  


   #[test]
   fn sphr_neighbours_h2(){
        // But can be defined for a valid index, where (1, 0, 0) is the
        // coordinate of the other H atom    
    
        write_h2();
        let structure = Structure::from("h2.xyz");
        let nbr_coords = structure.sphr_neighbours(0, "H", 2.0);
        
        assert_eq!(nbr_coords.len(), 1);
        let coord = &nbr_coords[0];

        assert!(is_very_close(coord.r, 1_f64));
        assert!(is_very_close(coord.theta, 0_f64));
        assert!(is_very_close(coord.phi, std::f64::consts::FRAC_PI_2));
    
        std::fs::remove_file("h2.xyz").expect("Could not remove file!"); 

    }  


   #[test]
   fn sphr_neighbours_h2_c_atoms(){
        // In H2 the hydrogen atom has no carbon neighbours        

        write_h2();
        let structure = Structure::from("h2.xyz");
        let nbr_coords = structure.sphr_neighbours(0, "C", 2.0);

        assert_eq!(nbr_coords.len(), 0);        

        std::fs::remove_file("h2.xyz").expect("Could not remove file!"); 

    }  
 
}

