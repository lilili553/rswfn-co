mod base;
mod datas;
mod editor;
mod funcs;
mod integ;
mod loger;
mod march;
mod maths;
mod qmcalc;
mod reader;
mod tools;
mod utils;
mod writer;

mod atomprop;
mod bondprop;
mod gridprop;
mod moleprop;
mod orbtprop;
mod tests;

mod cli;
mod tui;

use std;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 {
        cli::main();
    } else {
        tui::Main::new();
    }
}
