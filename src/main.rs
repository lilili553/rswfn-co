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
mod cmd;

use std;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    // println!("欢迎使用rswfn: 基于rust的波函数分析工具");
    // tests::test();
    if args.len() > 1 {
        println!("args:{:?}", args);
        cmd::main();
    } else {
        cli::Main::new();
    }
}
