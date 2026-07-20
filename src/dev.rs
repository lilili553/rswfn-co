use std::collections::HashMap;
use std::f64::consts::PI;

use log::Level;
use ndarray::{Array2, s};
use ndarray::{Axis, Shape};
use ndarray_linalg::Inverse;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;

use crate::base::{Mole, Stm, Vec3};
use crate::bondprop::direction;
use crate::datas::consts::BOHR;
use crate::datas::{self, btrans, mocv};
use crate::datas::{consts, roots};
use crate::editor::Editor;
use crate::gridprop;
use crate::gridprop::{density, wfnfunc};
use crate::integ::{self, get_2eid};
use crate::march;
use crate::maths;
use crate::moleprop;
use crate::orbtprop;
use crate::reader::{LogReader, MdeReader, Reader, XyzReader};
use crate::tools;
use crate::utils::{self, comp_matrix};
use crate::writer;
use crate::writer::{FchWriter, GjfWriter, Writer};
use crate::{atomprop, qmcalc};
use crate::{bondprop, editor};
use crate::{funcs, loger};

/// 测试读取 log 文件并解析分子信息
pub fn test_reader() {
    let path = r"E:\code\articles\MOCV\moles\MCl2\Agcl2.log";
    let reader = LogReader::new(path);
    println!("multip: {}", reader.multip);
    println!("charge: {}", reader.charge);
    let mole = Mole::from_file(path);
    println!("{}", mole);
}

pub fn test_order() {
    let path = r"E:\code\articles\MOCV\moles\sxf-scan-fchk\fchk-nmy\fchk\dui-Ca-lianben.log_scan_1.fchk";
    let mole = Mole::from_file(path);
    let caler = bondprop::order::Calculator::new(&mole);
    let omat = caler.mayer();
    println!("omat: {:?}", omat.shape());
}

pub fn test_charge() {
    let path = String::from(r"E:\code\rswfn\moles\C6H6.fch");
    let mole = Mole::from_file(&path);
    let caler = atomprop::charge::Calculator::new(&mole);
    let res = caler.hirshfeld();
    println!("res:{:?}", res);
}

pub fn test_topology() {
    let path = String::from(r"E:\code\rswfn\moles\C6H6.fch");
    let mole = Mole::from_file(&path);
    let caler = gridprop::topolog::Calculator::new(&mole);
    let res = caler.analysis("mol_rho_dm");
    for (pos, val, cpt) in res {
        let [px, py, pz] = pos;
        println!(
            "pos: ({:>10.4}, {:>10.4}, {:>10.4}), val: {}, cpt: {}",
            px, py, pz, val, cpt
        );
    }
}

/// 运行所有测试的入口函数
pub fn test() {
    loger::init_log();
    let t0 = std::time::Instant::now();
    // test_order();
    test_topology();
    println!("消耗时间 time: {:?}", t0.elapsed());
}
