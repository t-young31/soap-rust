## soap-rust
Smooth overlap of atomic positions (SOAP)[1] evaluated in rust using GTO radial basis functions and real spherical harmonic angular functions.[2]

[1] *A. P. Bartók et al., Phys. Rev. B, 87, 184115.*  
[2] *L. Himanen et al., Comput. Phys. Commun., 247, 106949.*

### Installation

To install from source first [install rust](https://www.rust-lang.org/tools/install), clone the repository, then build:

```bash
git clone https://github.com/t-young31/soap-rust.git
cd soap-rust
cargo build --release
```

> **_NOTE:_**  Requires linking against [GSL](https://www.gnu.org/software/gsl/).

### Usage

Execute the binary built in the *targets/release/* directory. For example to generate the 'SOAP vector' for the H-atom density in a methane molecule:

```bash
target/release/soap methane.xyz 0 H
```

where the first argument is the .xyz structure file, the second the atom index of the origin and the third the element comprising the neighbours. See the help for all options:

```bash
target/release/soap --help
```

### Testing
Rust provides a awesome testing platform so the tests can be run simply with:

```bash
cargo test
```

### Disclaimer

Absolutely no guarantee of correctness. Please fork this repository and submit a pull request if there are egregious errors!
