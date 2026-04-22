#![allow(non_snake_case)]
mod cub;
mod dmo;
mod fch;
mod gjf;
mod inp;
mod log;
mod mde;
mod mol;
mod non;
mod out;
mod sdf;
mod xyz;

use crate::base::Atoms;
use crate::base::Basis;
use crate::base::Coefs;
use crate::base::Mole;

pub use cub::CubReader;
pub use dmo::DmoReader;
pub use fch::FchReader;
pub use gjf::GjfReader;
pub use inp::InpReader;
pub use log::LogReader;
pub use mde::MdeReader;
pub use mol::MolReader;
pub use non::NonReader;
pub use out::OutReader;
pub use sdf::SdfReader;
pub use xyz::XyzReader;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub struct Marks {
    pub marks: HashMap<String, usize>,
}

impl Marks {
    pub fn new() -> Marks {
        Self { marks: HashMap::new() }
    }

    pub fn add(&mut self, marker: String, number: usize) {
        self.marks.insert(marker, number);
    }
}

pub trait Reader {
    fn get_atoms(&self) -> Atoms;
    fn get_basis(&self) -> Basis;
    fn get_coefs(&self) -> Coefs;
}

pub fn read_lines(path: &str) -> Vec<String> {
    let file = File::open(path).expect(format!("无法打开文件:{}\n", path).as_ref());
    let read = BufReader::new(file);
    let mut texts = Vec::new();
    for line in read.lines() {
        let line = line.expect("line read error");
        texts.push(line);
    }
    texts
}

// 根据电荷和自旋多重度计算α和β电子数
pub fn charge_number(charge: isize, multip: usize, total: usize) -> [usize; 2] {
    let N = total as isize;
    let C = charge as isize;
    let S = multip as isize;
    let nela = (N - C + S - 1) / 2;
    let nelb = (N - C - S + 1) / 2;
    [nela as usize, nelb as usize]
}

#[derive(Debug)]
pub enum FileType {
    Gjf, // 高斯输入文件
    Log, // 高斯输出文件
    Fch, // 高斯检查点文件
    Mol, // 分子文件
    Sdf,
    Xyz,
    Mde, // molden文件
    Dmo, // 自定义的包含原子局部坐标系的文件
    Inp, // Orca输入文件
    Cub, // 格点文件
}

pub fn get_fileType(path: &str) -> FileType {
    if path.ends_with(".gjf") {
        return FileType::Gjf;
    }
    if path.ends_with(".log") {
        return FileType::Log;
    }
    if path.ends_with(".out") {
        return FileType::Log;
    }
    if path.ends_with(".fch") || path.ends_with(".fchk") {
        return FileType::Fch;
    }
    if path.ends_with(".mol") {
        return FileType::Mol;
    }
    if path.ends_with(".sdf") {
        return FileType::Sdf;
    }
    if path.ends_with(".xyz") {
        return FileType::Xyz;
    }
    if path.ends_with(".molden") {
        return FileType::Mde;
    }
    if path.ends_with(".dmi") || path.ends_with(".dmo") {
        return FileType::Dmo;
    }
    if path.ends_with(".inp") {
        return FileType::Inp;
    }
    if path.ends_with(".cub") {
        return FileType::Cub;
    }
    panic!("未知文件类型:{path}");
}
