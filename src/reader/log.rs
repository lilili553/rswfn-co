// 读取log文件的脚本

use crate::base::{Atoms, Basis, BasisData, Coefs};
use crate::datas::consts::BOHR;
use crate::maths;
use crate::reader::{self, Reader, read_lines};
use ndarray::{Array2, Axis, concatenate};
use regex::Regex;
use std::collections::HashMap;
use std::process::exit;

// 一个振动频率
#[derive(Debug, Clone)]
pub struct Vibe {
    pub freq: f64,
    pub redm: f64,
    pub frcc: f64,
    pub irin: f64,
    pub vibs: Vec<[f64; 3]>,
}

impl Vibe {
    // 获取震动幅度最大的原子及方向
    pub fn get_max(&self) -> (usize, [f64; 3], f64) {
        let mut max_idx = 0;
        let mut max_vib = [0.0; 3];
        let mut max_len = 0.0;
        for (i, vib) in self.vibs.iter().enumerate() {
            let li = maths::vec_len(vib);
            if li > max_len {
                max_len = li;
                max_idx = i;
                max_vib = vib.clone();
            }
        }
        (max_idx, max_vib, max_len)
    }
}

pub struct LogReader {
    path: String,
    pub marks: HashMap<String, usize>,
    pub lines: Vec<String>,
    nato: usize,
    nobt: usize,
    pub charge: isize,
    pub multip: usize,
}

impl LogReader {
    // 搜索标题行，找到每块数据的位置
    pub fn new(path: &str) -> LogReader {
        let mut marks: HashMap<String, usize> = HashMap::new();
        let lines = read_lines(path);
        let mut nato = 0;
        let mut nobt = 0;
        let mut charge = 0;
        let mut multip = 0;
        for (l, line) in lines.iter().enumerate() {
            if line.len() == 0 {
                continue;
            }
            if line == "                          Input orientation:                          " {
                *marks.entry("geome".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == "                         Standard orientation:                         " {
                *marks.entry("geome".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == " AO basis set in the form of general basis input (Overlap normalization):" {
                *marks.entry("basis".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == " Effective core potential in the form of general ECP input:" {
                *marks.entry("ecp".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == "     Molecular Orbital Coefficients:" {
                *marks.entry("coefs".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == "     Alpha Molecular Orbital Coefficients:" {
                *marks.entry("coefs_a".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == "     Beta Molecular Orbital Coefficients:" {
                *marks.entry("coefs_b".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line.len() > 2 {
                if line.starts_with(" #") {
                    marks.entry("job".to_string()).or_insert(l as usize); // 只在第一次遇到的时候设置
                }
            }
            if line.contains(" NBasis=") {
                nato = line[8..14].trim().parse::<usize>().expect("NBasis转换失败");
                continue;
            }
            if line.contains(" NBsUse=") {
                nobt = line[8..14].trim().parse::<usize>().expect("NBasis转换失败");
                continue;
            }
            if line.contains(" Charge =") {
                charge = line[9..12].trim().parse::<isize>().expect("电荷解析失败");
                multip = line[27..].trim().parse::<usize>().expect("自旋解析失败");
                continue;
            }
            if line == "     Density Matrix:" {
                *marks.entry("density".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == " *** Overlap *** " {
                *marks.entry("overlap".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == " *** Kinetic Energy *** " {
                *marks.entry("kinetic".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == " ***** Potential Energy ***** " {
                *marks.entry("Potential".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == " ******    1 imaginary frequencies (negative Signs) ****** " {
                *marks.entry("vibrations".to_string()).or_insert(l as usize) = l as usize;
            }
        }
        LogReader {
            path: path.to_string(),
            marks: marks,
            lines: lines,
            nato: nato,
            nobt: nobt,
            charge: charge,
            multip: multip,
        }
    }

    pub fn read_job(&self) -> String {
        let title_number = match self.marks.get("job") {
            Some(n) => *n,
            None => {
                return "".to_string();
            }
        };
        let mut job = "".to_string();
        let re = Regex::new(r"^ -+$").expect("正则表达式构建失败");
        for line in self.lines.iter().skip(title_number as usize) {
            if re.is_match(line) {
                break;
            }
            job.push_str(&line.replace("#", "").trim());
        }
        println!("job={}", job);
        job
    }

    // 读取系数矩阵
    pub fn read_coefs(&self, title_number: usize) -> (Vec<usize>, Vec<usize>, Vec<String>, Vec<f64>, Vec<bool>, Array2<f64>) {
        let job = self.read_job();
        if !job.contains("pop=full") {
            log::error!("请加入关键词pop=full");
            exit(0);
        }
        let nato = self.nato;
        let nobt = self.nobt;

        let nblock = nobt / 5 + if nobt % 5 == 0 { 0 } else { 1 }; // 块的数量
        let mut cmat = Array2::<f64>::zeros((nato as usize, nobt as usize));
        let mut ato_atms: Vec<usize> = Vec::with_capacity(nato as usize);
        let mut ato_shls: Vec<usize> = Vec::with_capacity(nato as usize);
        let mut ato_syms: Vec<String> = Vec::with_capacity(nato as usize);
        let mut obt_engs: Vec<f64> = Vec::with_capacity(nobt as usize);
        let mut obt_occs: Vec<bool> = Vec::with_capacity(nobt as usize);
        // 先获取所有的行
        // let mut lines = Vec::new();
        // 彻底解决列不同问题
        let line = self.lines.get(title_number as usize + 4).expect("获取mark行失败");
        let nnum = line.len() - if nobt > 5 { 50 } else { nobt as usize * 10 };
        let line = line[0..nnum].to_owned();
        let mut ups = Vec::new();
        let mut los = Vec::new();
        let line: Vec<char> = line.chars().collect();
        for i in 1..line.len() {
            if line[i - 1] == ' ' && line[i] != ' ' {
                ups.push(i)
            }
            if line[i - 1] != ' ' && line[i] == ' ' {
                los.push(i)
            }
        }
        assert_eq!(ups.len(), 4);
        assert_eq!(los.len(), 4);
        let cols = [ups[1], ups[2], ups[3] - 1, ups[3] + 1, ups[3] + 1, nnum];
        // println!("cols={:?}", cols);
        let mut iatm = 0;
        for i in 0..nblock {
            let loc = title_number + 1 + i * (nato + 3);

            // 原子信息
            if i == 0 {
                for j in 0..nato {
                    let line = self.lines.get((loc + j + 3) as usize).expect("获取文件行失败");
                    let atom_part = line[cols[0]..cols[1]].trim();
                    if !atom_part.is_empty() {
                        iatm = atom_part.parse::<usize>().expect("原子索引解析失败");
                    }
                    ato_atms.push(iatm - 1);
                    let shl = line[cols[2]..cols[3]].trim().parse::<usize>().unwrap_or_else(|_| {
                        let msg = format!("{},shl解析失败:{}|{:?}", loc + j + 3, line, &line[cols[2]..cols[3]]);
                        log::error!("{}", msg);
                        panic!("{}", msg);
                    });
                    let sym = line[cols[4]..cols[5]].trim().to_string();
                    ato_shls.push(shl - 1);
                    ato_syms.push(sym);
                }
            }

            // 轨道占据
            let occs_line = self.lines.get(loc as usize + 1).expect("获取文件行失败");
            let occs_part = &occs_line[cols[5]..];
            // println!("occ_part={}", occs_part);
            let occs: Vec<bool> = occs_part
                .chars()
                .collect::<Vec<_>>()
                .chunks(10)
                .map(|chunk| {
                    let s: String = chunk.iter().collect();
                    // println!("s={}", s);
                    s.trim().ends_with('O')
                })
                .collect();
            obt_occs.extend_from_slice(&occs);

            // 轨道能量
            let engs_line = self.lines.get(loc as usize + 2).expect("获取文件行失败");
            let engs_part = &engs_line[cols[5]..];
            let engs: Vec<f64> = engs_part
                .char_indices()
                .step_by(10)
                .map(|(i, _)| {
                    engs_part[i..i + 10].trim().parse::<f64>().unwrap_or_else(|_| {
                        log::error!("轨道能量浮点数解析失败,loc:{},{:?}", loc + 1, engs_line);
                        exit(0);
                    })
                })
                .collect();
            obt_engs.extend_from_slice(&engs);

            // 轨道系数
            for j in 0..nobt {
                let coef_line = self.lines.get((loc + j + 3) as usize).expect("获取文件行失败");
                let coef_part = &coef_line[cols[5]..];
                let coefs: Vec<f64> = coef_part
                    .chars()
                    .collect::<Vec<_>>()
                    .chunks(10)
                    .map(|chunk| {
                        let s: String = chunk.iter().collect();
                        s.trim().parse::<f64>().expect("浮点数解析失败")
                    })
                    .collect();
                for (k, coef) in coefs.iter().enumerate() {
                    // cmat[(j as usize, i as usize * 5 + k)] = *coef;
                    if let Some(val) = cmat.get_mut((j as usize, i as usize * 5 + k)) {
                        *val = *coef;
                    } else {
                        panic!("分子轨道系数访问越界, 可能是分子轨道数量和原子轨道数量不一致")
                    }
                }
            }
        }
        (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat)
    }

    pub fn read_geome(&self, irow: Option<usize>) -> (Vec<usize>, Vec<[f64; 3]>) {
        let reg = Regex::new(r" +\d+ +(\d+) +\d +(-?\d+.\d{6}) +(-?\d+.\d{6}) +(-?\d+.\d{6})").expect("正则表达式构建失败");
        let mut atms = Vec::new();
        let mut xyzs = Vec::new();
        let irow = match irow {
            Some(i) => i,
            None => match self.marks.get("geome") {
                Some(number) => *number as usize,
                None => {
                    log::error!("没有读取到几何结构信息");
                    panic!("没有读取到几何结构信息");
                }
            },
        };
        for line in self.lines.iter().skip(irow + 5) {
            if let Some(caps) = reg.captures(&line) {
                let atm: usize = caps[1].parse().expect("整数解析失败");
                let x: f64 = caps[2].parse::<f64>().expect("浮点数解析失败") / BOHR;
                let y: f64 = caps[3].parse::<f64>().expect("浮点数解析失败") / BOHR;
                let z: f64 = caps[4].parse::<f64>().expect("浮点数解析失败") / BOHR;
                atms.push(atm);
                xyzs.push([x, y, z]);
            } else {
                break;
            }
        }

        (atms, xyzs)
    }

    pub fn read_geomes(&self) -> Vec<Atoms> {
        let mut geomes = Vec::new();
        let job = self.read_job();
        let mut mark = "Standard orientation:";
        if job.contains("nosymm") {
            mark = "Input orientation:";
        }
        if job == "" {
            return geomes;
        } else {
            if job.contains("opt") {
                for (l, line) in self.lines.iter().enumerate() {
                    if line.contains(mark) {
                        let (atms, xyzs) = self.read_geome(Some(l));
                        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
                        let geome = Atoms::new(&atms, &xyzs, &crgs, None);
                        geomes.push(geome);
                    }
                }
                return geomes;
            }
            return geomes;
        }
    }

    pub fn read_energys(&self) -> Option<Vec<f64>> {
        let marks = vec![
            "Zero-point correction",
            "Thermal correction to Energy",
            "Thermal correction to Enthalpy",
            "Thermal correction to Gibbs Free Energy",
            "Sum of electronic and zero-point Energies",
            "Sum of electronic and thermal Energies",
            "Sum of electronic and thermal Enthalpies",
            "Sum of electronic and thermal Free Energies",
        ];
        let mut engs = vec![0.0; 8];
        let mut yes = false;
        for line in self.lines.iter() {
            for (i, mark) in marks.iter().enumerate() {
                if line.contains(mark) {
                    yes = true;
                    if i < 4 {
                        let val: f64 = line[45..58].trim().parse().expect(format!("解析失败: {}", &line[45..58]).as_ref());
                        engs[i] = val;
                    } else {
                        let val = line[45..65].trim().parse().expect(format!("解析失败: {}", &line[45..65]).as_ref());
                        engs[i] = val;
                    }
                }
            }
        }
        if yes {
            return Some(engs);
        } else {
            return None;
        }
    }

    pub fn read_freqs(&self) -> Option<Vec<f64>> {
        let mut freqs = Vec::new();
        for line in self.lines.iter() {
            if line.contains("Frequencies --") {
                let vals: Vec<&str> = line[15..].trim().split_whitespace().collect();
                let vals = vals.iter().map(|val| val.parse::<f64>().expect(format!("解析失败: {}", val).as_ref()));
                freqs.extend(vals);
            }
        }
        if freqs.len() == 0 {
            return None;
        } else {
            return Some(freqs);
        }
    }

    // 读取对称矩阵
    pub fn read_symMat(&self, mut row: usize, loc: usize) -> Array2<f64> {
        let nato = self.nato as isize;
        let mut nrow = nato;
        let nobt = self.nobt;
        let mut matrix = Array2::zeros((nato as usize, nato as usize));
        let mut col = 0;
        let nblock = nobt / 5 + if nobt % 5 == 0 { 0 } else { 1 };
        for _i in 0..nblock {
            row += 1;
            for j in 0..nrow as usize {
                let line = self.lines.get(row).expect("获取行失败");
                let nums = line[loc..].trim().split_whitespace().map(|s| s.replace("D", "E").parse::<f64>().expect(&format!("解析失败:{},{}", row, line[loc..].to_string()))).collect::<Vec<f64>>();
                for (k, num) in nums.iter().enumerate() {
                    let irow = (nato - nrow) as usize + j;
                    matrix[(irow, col + k)] = *num;
                    matrix[(col + k, irow)] = *num;
                }
                row += 1;
            }
            col += 5;
            nrow -= 5;
        }
        matrix
    }

    pub fn read_dmat(&self) -> Option<Array2<f64>> {
        let row = match self.marks.get("density") {
            None => return None,
            Some(idx) => idx + 1,
        } as usize;
        let dmat = self.read_symMat(row, 21);
        Some(dmat)
    }

    pub fn read_smat(&self) -> Option<Array2<f64>> {
        let row = match self.marks.get("overlap") {
            None => return None,
            Some(idx) => idx + 1,
        } as usize;
        let dmat = self.read_symMat(row, 7);
        Some(dmat)
    }

    pub fn read_tmat(&self) -> Option<Array2<f64>> {
        let row = match self.marks.get("kinetic") {
            None => return None,
            Some(idx) => idx + 1,
        } as usize;
        let dmat = self.read_symMat(row, 7);
        Some(dmat)
    }

    pub fn read_vmat(&self) -> Option<Array2<f64>> {
        let row = match self.marks.get("Potential") {
            None => return None,
            Some(idx) => idx + 1,
        } as usize;
        let dmat = self.read_symMat(row, 7) * -1.0;
        Some(dmat)
    }

    // 读取赝势基组每个原子的有效核电荷数
    pub fn read_ecps(&self) -> Option<HashMap<usize, usize>> {
        let row = match self.marks.get("ecp") {
            None => return None,
            Some(idx) => idx + 1,
        } as usize;
        let mut crgs = HashMap::new();
        let reg1 = Regex::new(r" [A-Za-z]+(\d+)-ECP +\d +(\d+)").expect("正则表达式构建失败");
        for line in self.lines.iter().skip(row) {
            if reg1.is_match(line) {
                let caps = reg1.captures(line).expect("正则匹配失败");
                let atm = caps[1].parse::<usize>().expect("解析失败");
                let ele = caps[2].parse::<usize>().expect("解析失败");
                crgs.insert(atm, ele);
            }
            if line == " " {
                break;
            }
        }
        Some(crgs)
    }

    pub fn read_vibrations(&self) -> Option<Vec<Vibe>> {
        let mut vibrations: Vec<Vibe> = vec![];
        let irow = match self.marks.get("vibrations") {
            Some(idx) => idx + 7,
            None => return None,
        };
        let natm = self.read_geome(None).0.len();
        let ncol = natm * 3 - 6;
        let sblock = natm + 7; // 块的大小
        let nblock = ncol / 3 + if ncol % 3 == 0 { 0 } else { 1 }; // 块的数量
        let mut freqs = vec![];
        let mut redms = vec![];
        let mut frccs = vec![];
        let mut irins = vec![];
        let mut vibrs = vec![vec![]; ncol];
        for ib in 0..nblock {
            let start = irow + ib * sblock;
            let freql: Vec<f64> = self.lines[start + 2][15..].trim().split_whitespace().map(|v| v.parse().unwrap()).collect();
            let redml: Vec<f64> = self.lines[start + 3][15..].trim().split_whitespace().map(|v| v.parse().unwrap()).collect();
            let frccl: Vec<f64> = self.lines[start + 4][15..].trim().split_whitespace().map(|v| v.parse().unwrap()).collect();
            let irinl: Vec<f64> = self.lines[start + 5][15..].trim().split_whitespace().map(|v| v.parse().unwrap()).collect();
            freqs.extend(freql);
            redms.extend(redml);
            frccs.extend(frccl);
            irins.extend(irinl);
            for l in start + 7..start + 7 + natm {
                let line = &self.lines[l][10..];
                let vals: Vec<f64> = line.trim().split_whitespace().map(|v| v.parse().unwrap()).collect();
                for j in 0..vals.len() / 3 {
                    let dx = vals[j * 3 + 0];
                    let dy = vals[j * 3 + 1];
                    let dz = vals[j * 3 + 2];
                    vibrs[ib * 3 + j].push([dx, dy, dz]);
                }
            }
        }
        for i in 0..ncol {
            let vibe = Vibe { freq: freqs[i], redm: redms[i], frcc: frccs[i], irin: irins[i], vibs: vibrs[i].clone() };
            vibrations.push(vibe);
        }
        Some(vibrations)
    }
}

impl Reader for LogReader {
    // 读取结构
    fn get_atoms(&self) -> Atoms {
        let (atms, xyzs) = self.read_geome(None);

        let mut crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();

        let ecps = self.read_ecps();
        match ecps {
            None => {}
            Some(ecps) => {
                for (idx, ele) in ecps.iter() {
                    let idx = *idx as usize - 1;
                    crgs[idx] = (atms[idx] - *ele) as f64;
                }
            }
        };
        let mut atoms = Atoms::new(&atms, &xyzs, &crgs, None);
        let C = self.charge;
        let S = self.multip as isize;
        let mut N = 0;
        let (atms, _) = self.read_geome(None);
        for atm in atms.iter() {
            N += *atm as isize;
        }
        let neles = reader::charge_number(C, S as usize, N as usize);
        atoms.neles = neles;
        atoms
    }

    fn get_basis(&self) -> Basis {
        let s1 = Regex::new(r"^ +(\d+) +\d+$").expect("reg new error");
        let s2 = Regex::new(r" ([SPDFG]+) +(\d+) \d.\d{2} +\d.\d{12}").expect("reg new error");
        let s3 = Regex::new(r"^ +(( +-?\d.\d{10}D[+-]\d{2}){2,3})").expect("reg new error");
        let mut ang_map: HashMap<char, usize> = HashMap::new();
        ang_map.insert('S', 0);
        ang_map.insert('P', 1);
        ang_map.insert('D', 2);
        ang_map.insert('F', 3);
        ang_map.insert('G', 4);
        let mut data: Vec<BasisData> = Vec::new();
        let mut angs: Vec<usize> = Vec::new();
        let mut shl: usize = 0;
        let mut atm: usize = 0;
        let mut basis = Basis::blank();
        let title_number = match self.marks.get("basis") {
            Some(number) => *number,
            None => {
                log::error!("没有找到基组信息");
                return basis;
            }
        };
        for line in self.lines.iter().skip(title_number as usize + 1) {
            if s1.is_match(&line) {
                let caps = s1.captures(&line).expect("cap error");
                atm = caps[1].parse::<usize>().expect("usize解析失败");
            } else if s2.is_match(&line) {
                let caps = s2.captures(&line).expect("cap error");
                let shl_sym: String = caps[1].to_string();
                angs = shl_sym.chars().map(|s| *ang_map.get(&s).expect("未知角动量符号")).collect();
                shl += 1;
            } else if s3.is_match(&line) {
                let caps = s3.captures(&line).expect("cap error");
                let nums: Vec<f64> = caps[1]
                    .split_whitespace()
                    .map(|each| {
                        let coef = each.replace("D", "e");
                        coef.parse::<f64>().expect("f64解析失败")
                    })
                    .collect();
                if angs.len() == 1 {
                    let alp = nums[0];
                    let coe = nums[1];
                    data.push(BasisData::new(atm - 1, shl - 1, angs[0], alp, coe, (0, 0, 0)));
                }
                if angs.len() == 2 {
                    let alp = nums[0];
                    let coe1 = nums[1];
                    let coe2 = nums[2];
                    data.push(BasisData::new(atm - 1, shl - 1, angs[0], alp, coe1, (0, 0, 0)));
                    data.push(BasisData::new(atm - 1, shl - 1, angs[1], alp, coe2, (0, 0, 0)));
                }
            } else if line == " ****" {
                shl = 0;
                continue;
            } else {
                break;
            }
        }
        let mut idxs: Vec<usize> = (0..data.len()).collect();
        idxs.sort_by_key(|idx| (data[*idx].atm, data[*idx].shl, data[*idx].ang));
        data = idxs.iter().map(|idx| data[*idx].clone()).collect();
        basis.build(data);
        basis
    }

    fn get_coefs(&self) -> Coefs {
        let coefs = if self.marks.contains_key("coefs") {
            let title_number = *self.marks.get("coefs").expect("titles get coefs error");
            let (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat) = self.read_coefs(title_number);
            let obt_occs = obt_occs
                .iter()
                .map(|x| {
                    if *x {
                        return 2;
                    } else {
                        return 0;
                    }
                })
                .collect();
            Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat)
        } else if self.marks.contains_key("coefs_a") && self.marks.contains_key("coefs_b") {
            let title_number_a = *self.marks.get("coefs_a").expect("titles get coefs_a error");
            let title_number_b = *self.marks.get("coefs_b").expect("titles get coefs_b error");
            let (ato_atms, ato_shls, ato_syms, mut obt_engs, mut obt_occs, cmat_a) = self.read_coefs(title_number_a);
            let coefs_b = self.read_coefs(title_number_b);
            obt_engs.extend_from_slice(&coefs_b.3);
            obt_occs.extend_from_slice(&coefs_b.4);
            let cmat = concatenate(Axis(1), &[cmat_a.view(), coefs_b.5.view()]).expect("矩阵拼接失败");
            let obt_occs = obt_occs
                .iter()
                .map(|x| {
                    if *x {
                        return 1;
                    } else {
                        return 0;
                    }
                })
                .collect();
            Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat)
        } else {
            log::error!("没有找到分子轨道");
            Coefs::blank()
        };
        coefs
    }
}
