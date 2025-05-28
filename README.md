# bader-rs (0.4.8)
![build](https://github.com/kerrigoon/bader-rs/workflows/build/badge.svg?branch=master)
[![Latest Version](https://img.shields.io/crates/v/bader.svg)](https://crates.io/crates/bader)
[![Documentation](https://docs.rs/bader/badge.svg)](https://docs.rs/bader/)
[![DOI](https://zenodo.org/badge/292534636.svg)](https://zenodo.org/badge/latestdoi/292534636)
[![MSRV: rustc 1.60.0+](https://img.shields.io/badge/MSRV-rustc_1.60.0+-lightgray.svg)](https://blog.rust-lang.org/2022/04/07/Rust-1.60.0/)

An incredibly fast, multi-threaded, Bader charge partitioning tool. Based on methods presented in [Yu Min  and Trinkle Dallas R. 2011  J. Che.m Phys. 134 064111] with adaptions for multi-threading and increased speed.
## Installation
### Pre-built Binary
There are pre-built 64bit binaries for Linux, Mac and Windows provided with the source code for the latest [release].
### Cargo
If these binaries don't cover your OS the easiest way to install is via [cargo].
```sh
$ cargo install bader
```
### From Source
To check out the lastest features not in the binaries yet you can compile from source. To do this run the following, which will create the ./target/release/bca executable.
```sh
$ git clone https://github.com/adam-kerrigan/bader-rs
$ cd bader-rs
$ cargo build --verbose --release
```
From here you can either move or link the binary to folder in your path.
```sh
$ mv ./target/release/bca ~/bin
```
### Minimum Supported Rust Version (MSRV)
This crate is guaranteed to compile on stable Rust 1.60.0 and up.
## Usage
The program takes a charge density file as input and performs Bader analysis of the data. Currently it supports density in [VASP] or [cube] formats. It is recommended to run VASP calculations with [LAECHG] = .TRUE. to print the core density and self-consistent valence density. These can then be passed as reference files to the program using the -r, --reference flag where they will be summed.
```sh
$ bca CHGCAR -r AECCAR0 -r AECCAR2
```
VASP charge density files containing spin densities will output the the partitioned spin also. To achieve this for cube files requires using the --spin flag to pass a second file to treat as the spin density.
```sh
$ bca charge-density.cube -s spin-density.cube
```
For a detailed list of usage options run
```sh
$ bca --help
```
## Output
The program outputs the Atomic Charge File (ACF.dat) which contians the charge (and spin) information for each atom.
## License
MIT

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

[release]: <https://github.com/adam-kerrigan/bader-rs/releases/latest>
[VASP]: <https://www.vasp.at/>
[cube]: <https://gaussian.com/>
[LAECHG]: <https://www.vasp.at/wiki/index.php/LAECHG>
[Yu Min  and Trinkle Dallas R. 2011  J. Che.m Phys. 134 064111]: <https://doi.org/10.1063/1.3553716>
[cargo]: <https://doc.rust-lang.org/cargo/getting-started/installation.html>
