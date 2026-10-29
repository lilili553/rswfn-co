# rswfn

[![crates.io](https://img.shields.io/crates/v/rswfn.svg)](https://crates.io/crates/rswfn)
[![docs.rs](https://docs.rs/rswfn/badge.svg)](https://docs.rs/rswfn)
[![license](https://img.shields.io/crates/l/rswfn.svg)](#license)

**rswfn** (Rust Wavefunction) is a wavefunction analysis library and
command-line tool for quantum chemistry, written in Rust. It reads
electronic-structure data produced by common quantum chemistry programs and
performs a wide range of wavefunction analyses — population analysis, bonding
analysis, aromaticity, orbital analysis, and real-space (grid) function
analysis.

It can be used both as a **library** and as a **CLI / TUI** application.

## Features

- **Readers** for Gaussian (`.gjf`, `.log`, `.out`, `.fch`, `.fchk`), ORCA
  (`.inp`, `.out`), Molden (`.molden`), cube (`.cub`), `.mol`, `.sdf`, `.xyz`,
  and a custom `.dmi` / `.dmo` format that carries per-atom local coordinate
  frames.
- **Atomic properties**: Mulliken / Löwdin / Hirshfeld charges, spin
  population, π-electron population (the POCV and MOCV orbital-transformation
  methods), Fukui functions and dual descriptors, and free valence.
- **Bond properties**: Mayer and Wiberg bond orders, per-orbital bond-order
  contributions, σ / π / δ / φ decomposition, and π bond orders.
- **Molecular & orbital properties**: π-based aromaticity (PISD), multipole
  moments, orbital composition, and per-orbital π-electron counts.
- **Real-space / grid analysis**: electron density and its derivatives,
  sign(λ₂)ρ, RDG, IRI, ELF / LOL, kinetic-energy densities, electrostatic
  potential, information entropy, Fukui-function grids, and topology search
  for critical points.
- **Integrals** over Gaussian basis sets: overlap, kinetic, nuclear-attraction
  and two-electron repulsion.
- **Writers** for cube (`.cub`), `.xyz`, `.gjf`, `.fch`, `.mol` / `.dmo` and
  `.obj`.

## Installation

### Command-line tool

```sh
cargo install rswfn
```

### As a library

```sh
cargo add rswfn
```

## Command-line usage

Running `rswfn` with no arguments launches the interactive TUI. With arguments
it runs as a CLI:

```sh
# Mulliken / Löwdin / Hirshfeld charges and spin population
rswfn atomprop charge mulliken C6H6.fch
rswfn atomprop charge lowdin C6H6.fch
rswfn atomprop charge hirshfeld C6H6.fch
rswfn atomprop charge spin C6H6.fch

# π-electron population (POCV / MOCV)
rswfn atomprop pi-pop C6H6.fch --method pocv

# Fukui functions and dual descriptor (needs the N, N+1 and N-1 wavefunctions)
rswfn atomprop fukui phenol.out phenol_N+1.out phenol_N-1.out

# Mayer / Wiberg bond orders
rswfn bondprop order mayer C6H6.fch

# σ/π/δ/φ decomposition of the bond order between atoms 1 and 2
rswfn bondprop decom 1 2 C6H6.fch

# Aromaticity (PISD)
rswfn moleprop aromacity C6H6.fch

# Composition of molecular orbital 21
rswfn orbtprop composition 21 C6H6.fch

# Export a real-space function on a grid to a cube file
rswfn gridprop density C6H6.fch rho.cub --func density --step 0.2
rswfn gridprop potential C6H6.fch esp.cub

# Extract coordinates, convert formats, extract supporting information
rswfn coord C6H6.fch coords.xyz
rswfn trans C6H6.fch C6H6.gjf
rswfn si scan1.log scan2.log -o si.txt
```

Run `rswfn --help`, or `rswfn <subcommand> --help`, for the full set of options.

## Library usage

```rust
use rswfn::atomprop::charge::Calculator;
use rswfn::base::Mole;

fn main() {
    // Read a wavefunction from a file; the format is detected from the extension.
    let mole = Mole::from_file("C6H6.fch");

    // Mulliken charges
    let charges = Calculator::new(&mole).mulliken();
    println!("{charges:?}");
}
```

## Feature flags

Both features are **disabled by default** and link against locally built native
libraries, so they are not usable from the published crate:

- `use-cint` — use the external `libcint` library for integrals.
- `use-func` — use the external `libxc` library for exchange–correlation
  functionals.

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in this crate by you, as defined in the Apache-2.0 license, shall
be dual licensed as above, without any additional terms or conditions.

## 中文简介

`rswfn` 是一个用 Rust 编写的量子化学**波函数分析库与命令行工具**，可以读取
Gaussian、ORCA、Molden 等程序的输出文件，进行各类波函数分析，包括：原子电荷
布居（Mulliken / Löwdin / Hirshfeld）、π 电子布居（POCV / MOCV）、福井函数与
双描述符、自由价、键级（Mayer / Wiberg）及其分解、芳香性、轨道成分分析，以及
实空间函数（密度、静电势、ELF、RDG 等）的格点计算与导出。既可以作为库使用，
也可以直接作为命令行 / TUI 程序运行。
