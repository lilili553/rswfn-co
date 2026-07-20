use std::process::exit;
// use std::char;

use regex::Regex;

use crate::base::{Atoms, Basis, Coefs, Stm};
use crate::datas::consts::BOHR;
use crate::datas::elems::sym2atm;
use crate::reader::charge_number;
use crate::reader::{Reader, read_lines};

/// 高斯输入文件(.gjf)读取器
pub struct GjfReader {
    path: String,
    lines: Vec<String>,
}

impl GjfReader {
    /// 创建GjfReader实例
    pub fn new(path: &str) -> Self {
        let lines = read_lines(path);
        Self {
            path: path.to_string(),
            lines,
        }
    }

    /// 读取几何结构
    pub fn read_geome(&self) -> (Vec<usize>, Vec<[f64; 3]>) {
        let mut atms = vec![];
        let mut xyzs = vec![];
        let s1 = Regex::new(r" ([A-Za-z\d]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)").expect("正则表达式生成失败");
        for line in self.lines.iter() {
            if s1.is_match(line) {
                let caps = s1.captures(line).expect("捕获失败");
                let sym = caps[1].to_string();
                let atm = sym2atm(&sym);
                let x = caps[2].parse::<f64>().expect("解析浮点数失败") * BOHR;
                let y = caps[3].parse::<f64>().expect("解析浮点数失败") * BOHR;
                let z = caps[4].parse::<f64>().expect("解析浮点数失败") * BOHR;
                atms.push(atm);
                xyzs.push([x, y, z]);
            }
        }
        (atms, xyzs)
    }

    /// 读取电荷和自旋多重度
    pub fn read_charge(&self) -> (isize, usize) {
        let s0 = Regex::new(r"^ *(-?\d+) +(-?\d+)").expect("正则表达式生成失败");
        let mut C = 0; // 电荷
        let mut S = 0; // 自旋
        let mut matched = false;
        for line in self.lines.iter() {
            if s0.is_match(line) {
                let caps = s0.captures(line).expect("捕获失败");
                C = caps[1].parse::<isize>().expect("转换失败");
                S = caps[2].parse::<usize>().expect("转换失败");
                matched = true;
                break;
            }
        }
        if !matched {
            log::error!("无法读取电荷和自旋");
            exit(0)
        }
        let (atms, _) = self.read_geome();
        (C, S)
    }

    /// 读取作业关键词
    pub fn read_job(&self) -> String {
        let mut job = String::new();
        for line in self.lines.iter() {
            if line.starts_with("#") {
                job.push_str(&line.replace("#", "").trim().to_string());
            }
        }
        job
    }
}

impl Reader for GjfReader {
    fn get_atoms(&self) -> Atoms {
        let (atms, xyzs) = self.read_geome();
        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
        let mut geome = Atoms::new(&atms, &xyzs, &crgs);
        let N = atms.iter().map(|atm| *atm as usize).sum::<usize>();
        let (C, S) = self.read_charge();

        geome.neles = charge_number(C, S, N);
        geome
    }
    fn get_basis(&self) -> Option<Basis> {
        None
    }

    fn get_coefs(&self) -> Option<Coefs> {
        None
    }
}
