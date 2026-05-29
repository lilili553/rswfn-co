#![allow(non_snake_case)]
use regex::Regex;

use crate::datas::consts::BOHR;
use crate::datas::elems;
use crate::datas::temps;
use crate::reader::LogReader;

const ENGS_STR: &str = r#"
Zero-point correction                        {E0}
Thermal correction to Energy                 {E1}
Thermal correction to Enthalpy               {E2}
Thermal correction to Gibbs Free Energy      {E3}
Sum of electronic and zero-point Energies    {E4}
Sum of electronic and thermal Energies       {E5}
Sum of electronic and thermal Enthalpies     {E6}
Sum of electronic and thermal Free Energies  {E7}"#;

pub struct Tool {}

impl Tool {
    pub fn new() -> Self {
        Self {}
    }
    pub fn get_SI(&self, paths: &Vec<String>) -> String {
        // 提取支持信息

        let mut text = String::new();
        for path in paths.iter() {
            let mut temp = temps::SI.to_string();
            temp = temp.replace("{name}", path);
            let reader = LogReader::new(path);
            let (atms, xyzs) = reader.read_geome(None);
            // 能量信息
            let mut engsStr = ENGS_STR.to_string();
            match reader.read_energys() {
                Some(engs) => {
                    for (i, eng) in engs.iter().enumerate() {
                        engsStr = engsStr.replace(&format!("{{E{i}}}"), &format!("{:>16.8}", eng));
                    }
                    temp = temp.replace("{engs}", &engsStr);
                }
                None => {
                    // text += &format!("{}:获取能量失败！\n", path);
                    temp = temp.replace("{engs}", &format!("{}:获取能量失败！\n", path));
                }
            };
            // 频率信息
            let mut freqsStr = String::new();
            match reader.read_freqs() {
                Some(freqs) => {
                    for (i, freq) in freqs.iter().enumerate() {
                        freqsStr += &format!("{:20.2}", freq);
                        if (i + 1) % 3 == 0 && (i + 1) != freqs.len() {
                            freqsStr += "\n";
                        }
                    }
                    freqsStr = freqsStr.trim_end_matches("\n").to_string();
                    temp = temp.replace("{freqs}", &freqsStr);
                }
                None => {
                    temp = temp.replace("{freqs}", &format!("{}:获取频率失败！\n", path));
                }
            };
            // 结构信息
            let mut geomeStr = String::new();

            for i in 0..atms.len() {
                let atm = atms[i];
                let xyz = xyzs[i];
                let sym = elems::atm2sym(atm);
                geomeStr += &format!("{:<12}{:>16.8}{:>16.8}{:>16.8}\n", sym, xyz[0] / BOHR, xyz[1] / BOHR, xyz[2] / BOHR);
            }
            geomeStr = geomeStr.trim_end_matches("\n").to_string();
            temp = temp.replace("{geome}", &geomeStr);

            text += &format!("{}\n", temp);
        }
        text
    }

    /// 获取每一个优化步骤的结构和能量
    pub fn split_opt(&self, path: String) -> Vec<(Vec<usize>, Vec<[f64; 3]>, f64)> {
        let reader = LogReader::new(&path);
        let lines = &reader.lines;
        let mut geomes = Vec::new();
        let mut atms = Vec::new();
        let mut xyzs = Vec::new();
        let mut eng = 0.0;
        let reg_xyz = Regex::new(r" +\d+ +(\d+) +\d+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)").expect("reg new err");
        let reg_eng = Regex::new(r" SCF Done:.*=\s*(-?\d+\.\d+)").expect("reg new err");
        for line in lines {
            if line == "" {
                continue;
            }
            if line.trim() == "Input orientation:" {
                atms.clear();
                xyzs.clear();
                continue;
            }
            if line.trim() == "Standard orientation:" {
                atms.clear();
                xyzs.clear();
                continue;
            }
            if reg_xyz.is_match(line) {
                let caps = reg_xyz.captures(line).expect("cap err");
                let atm = caps[1].parse::<usize>().expect("use parse err");
                let x = caps[2].parse::<f64>().expect("f64 parse error");
                let y = caps[3].parse::<f64>().expect("f64 parse error");
                let z = caps[4].parse::<f64>().expect("f64 parse error");
                atms.push(atm);
                xyzs.push([x, y, z]);
                continue;
            }
            if line.starts_with(" SCF Done:") {
                let caps = reg_eng.captures(line).expect("cap err");
                eng = caps[1].parse::<f64>().expect("parse f64 err");
            }

            if line.starts_with(" Predicted change in Energy=") {
                // let geome = Geome::new(&atms, &xyzs);
                geomes.push((atms.clone(), xyzs.clone(), eng));
                continue;
            }
            if line.starts_with(" Optimization completed.") {
                break;
            }
        }
        geomes
    }
}
