use ndarray::{Array2, s};

use crate::base::Atoms;
use crate::base::Mole;
use crate::datas::FACS;
use crate::datas::btrans;
use crate::datas::btrans::ang2syms;
use crate::datas::btrans::sym2lmn;
use crate::integ::RsInt;
use core::fmt;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::OnceLock;

#[derive(Debug, Clone)]
pub struct Ato {
    pub atm: usize,
    pub shl: usize,
    pub ang: usize,
    pub sym: String,
    pub loc: usize,
    pub form: String,
    pub alps: Vec<f64>,
    pub coes: Vec<f64>,
    coes_norm: OnceLock<Vec<f64>>,
    lmn: OnceLock<[usize; 3]>,
}

impl Ato {
    pub fn new(
        atm: usize,
        shl: usize,
        ang: usize,
        sym: String,
        loc: usize,
        form: String,
        alps: &Vec<f64>,
        coes: &Vec<f64>,
    ) -> Self {
        let alps = alps.clone();
        let coes = coes.clone();
        let coes_norm = OnceLock::new();
        let lmn = OnceLock::new();
        Self {
            atm,
            shl,
            ang,
            sym,
            loc,
            form,
            alps,
            coes,
            coes_norm,
            lmn,
        }
    }
    pub fn coes_norm(&self) -> &Vec<f64> {
        self.coes_norm.get_or_init(|| {
            let [l, m, n] = self.lmn();
            let mut coes = self.coes.clone();
            for i in 0..coes.len() {
                let fac = FACS[*l] * FACS[*m] * FACS[*n];
                let alp = self.alps[i];
                let coe = self.coes[i];
                let nm = (2.0 * alp / PI).powf(0.75) * f64::sqrt((4.0 * alp).powi(self.ang as i32) / fac);
                coes[i] = nm * coe;
            }
            coes
        })
    }

    pub fn lmn(&self) -> &[usize; 3] {
        self.lmn.get_or_init(|| sym2lmn(&self.sym))
    }
}

#[derive(Debug, Clone)]
pub struct Atos {
    pub raws: Vec<Ato>,
    pub sphs: Vec<Ato>,
    pub cars: Vec<Ato>,
}

impl Atos {
    pub fn len(&self, form: &str) -> usize {
        match form {
            "raw" => self.raws.len(),
            "sph" => self.sphs.len(),
            "car" => self.cars.len(),
            _ => panic!("form error"),
        }
    }

    pub fn atms(&self, form: &str) -> Vec<usize> {
        match form {
            "raw" => self.raws.iter().map(|ato| ato.atm).collect(),
            "sph" => self.sphs.iter().map(|ato| ato.atm).collect(),
            "car" => self.cars.iter().map(|ato| ato.atm).collect(),
            _ => panic!("form error"),
        }
    }

    pub fn shls(&self, form: &str) -> Vec<usize> {
        match form {
            "raw" => self.raws.iter().map(|ato| ato.shl).collect(),
            "sph" => self.sphs.iter().map(|ato| ato.shl).collect(),
            "car" => self.cars.iter().map(|ato| ato.shl).collect(),
            _ => panic!("form error"),
        }
    }

    pub fn syms(&self, form: &str) -> Vec<String> {
        match form {
            "raw" => self.raws.iter().map(|ato| ato.sym.clone()).collect(),
            "sph" => self.sphs.iter().map(|ato| ato.sym.clone()).collect(),
            "car" => self.cars.iter().map(|ato| ato.sym.clone()).collect(),
            _ => panic!("form error"),
        }
    }

    pub fn locs(&self, form: &str) -> Vec<usize> {
        match form {
            "raw" => self.raws.iter().map(|ato| ato.loc).collect(),
            "sph" => self.sphs.iter().map(|ato| ato.loc).collect(),
            "car" => self.cars.iter().map(|ato| ato.loc).collect(),
            _ => panic!("form error"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Shl {
    pub atm: usize,
    pub shl: usize,
    pub ang: usize,
    pub form: String, // 0 mix 1 sph 2 car
    pub alps: Vec<f64>,
    pub coes: Vec<f64>,
}

impl Shl {
    pub fn new(atm: usize, shl: usize, ang: usize, alps: Vec<f64>, coes: Vec<f64>) -> Self {
        Self {
            atm,
            shl,
            ang,
            form: "car".to_string(),
            alps,
            coes,
        }
    }
    pub fn coes_norm(&self) -> Vec<f64> {
        let mut coes = self.coes.clone();
        for i in 0..coes.len() {
            let alp = self.alps[i];
            let coe = self.coes[i];
            let nm = (2.0 * alp / PI).powf(0.75) * f64::sqrt((4.0 * alp).powi(self.ang as i32));
            coes[i] = nm * coe;
        }
        coes
    }

    pub fn loc(&self) -> usize {
        let alp_max = self.alps.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        if *alp_max > 10.0 { 0 } else { 1 }
    }
}

#[derive(Debug, Clone)]
pub struct Shls {
    pub shls: Vec<Shl>,
    pub uplo: OnceLock<[Vec<[usize; 2]>; 3]>,
}

const ANG_NUM_SPH: [usize; 6] = [1, 3, 5, 7, 9, 11];
const ANG_NUM_CAR: [usize; 6] = [1, 3, 6, 10, 15, 21];

impl Shls {
    pub fn new(shls: Vec<Shl>) -> Shls {
        Shls {
            shls,
            uplo: OnceLock::new(),
        }
    }

    pub fn get(&self, i: usize) -> &Shl {
        &self.shls[i]
    }

    pub fn len(&self) -> usize {
        self.shls.len()
    }
    pub fn uplo(&self, form: &str) -> &Vec<[usize; 2]> {
        let uplo = self.uplo.get_or_init(|| {
            let mut uplo_raw = vec![];
            let mut uplo_sph = vec![];
            let mut uplo_car = vec![];
            let mut idx_raw = 0;
            let mut idx_sph = 0;
            let mut idx_car = 0;
            for shl in self.shls.iter() {
                let ang = shl.ang;
                let num_raw = match shl.form.as_str() {
                    "car" => ANG_NUM_CAR[ang],
                    "sph" => ANG_NUM_SPH[ang],
                    _ => panic!("form error"),
                };
                let num_sph = ANG_NUM_SPH[ang];
                let num_car = ANG_NUM_CAR[ang];
                uplo_raw.push([idx_raw, idx_raw + num_raw]);
                uplo_sph.push([idx_sph, idx_sph + num_sph]);
                uplo_car.push([idx_car, idx_car + num_car]);
                idx_raw += num_raw;
                idx_sph += num_sph;
                idx_car += num_car;
            }
            [uplo_raw, uplo_sph, uplo_car]
        });
        match form {
            "raw" => &uplo[0],
            "sph" => &uplo[1],
            "car" => &uplo[2],
            _ => panic!("form error"),
        }
    }

    pub fn nato(&self, form: &str) -> usize {
        let uplo = self.uplo(form);
        let mut nato = 0;
        for [up, lo] in uplo.iter() {
            nato += lo - up;
        }
        nato
    }

    // 在 basis.rs 的 impl Shls 中
    pub fn atm_uplo(&self, form: &str) -> Vec<[usize; 2]> {
        let uplo = self.uplo(form);
        let mut atm_uplo = Vec::new();
        let mut prev_atm = None;
        for (i, shl) in self.shls.iter().enumerate() {
            let [start, end] = uplo[i];
            if Some(shl.atm) != prev_atm {
                atm_uplo.push([start, end]);
                prev_atm = Some(shl.atm);
            } else {
                atm_uplo.last_mut().unwrap()[1] = end;
            }
        }
        atm_uplo
    }
}

#[derive(Debug, Clone)]
pub struct Bas {
    pub atm: usize,
    pub shl: usize,
    pub ang: usize,
    pub alp: f64,
    pub coe: f64,
    // pub nlm: (usize, usize, isize), //n,l,m 兼容真实原子轨道，先这么写，后面应该要搞成枚举类型
}

impl Bas {
    pub fn new(atm: usize, shl: usize, ang: usize, alp: f64, coe: f64) -> Self {
        Self {
            atm,
            shl,
            ang,
            alp,
            coe,
        }
    }
}
#[derive(Debug, Clone)]
pub struct Bass {
    bass: Vec<Bas>,
}

#[derive(Debug)]
pub struct Basis {
    pub shls: Shls, // atm,shl,ang,alp,coe,sph
    pub atos: OnceLock<Atos>,
    pub bass: OnceLock<Bass>,
    pub smat: OnceLock<[Array2<f64>; 3]>, // 存储球谐和笛卡尔的重叠积分
                                          // pub mole: *const Mole,
}

impl Basis {
    pub fn new(shls: Vec<Shl>) -> Basis {
        let shls = Shls::new(shls);
        let atos = OnceLock::new();
        let bass = OnceLock::new();
        let smat = OnceLock::new();
        // let mole = std::ptr::null();
        Basis { shls, atos, bass, smat }
    }

    pub fn build(&mut self, shls: Vec<Shl>) {
        self.shls = Shls::new(shls);
    }

    pub fn is_none(&self) -> bool {
        self.shls.len() == 0
    }

    pub fn get_by_atm(&self, atm: usize) -> Vec<Shl> {
        let mut shls = vec![];
        let nshl = self.shls.len();
        for i in 0..nshl {
            let shl = self.shls.get(i);
            if shl.atm == atm {
                shls.push(shl.clone());
            }
        }
        shls
    }

    pub fn ato_keys(&self) -> Vec<String> {
        let mut keys = vec![];
        let nshl = self.shls.len();
        for i in 0..nshl {
            let shl = self.shls.get(i);
            let syms = ang2syms(shl.ang, &shl.form);
            for sym in syms.iter() {
                let key = format!("{},{},{}", shl.atm, shl.shl, sym);
                if !keys.contains(&key) {
                    keys.push(key);
                }
            }
        }
        keys
    }

    pub fn atos(&self) -> &Atos {
        let shls = &self.shls;
        let nshl = shls.len();
        self.atos.get_or_init(|| {
            let mut raws = vec![];
            let mut sphs = vec![];
            let mut cars = vec![];
            for i in 0..nshl {
                let atm = shls.get(i).atm;
                let shl = shls.get(i).shl;
                let ang = shls.get(i).ang;
                let loc = shls.get(i).loc();
                let form = &shls.get(i).form;
                let alps = &shls.get(i).alps;
                let coes = &shls.get(i).coes;
                let syms = ang2syms(ang, form);
                for sym in syms {
                    let ato = Ato::new(atm, shl, ang, sym.to_owned(), loc, form.to_owned(), alps, coes);
                    raws.push(ato);
                }
                let syms = ang2syms(ang, "sph");
                for sym in syms {
                    let ato = Ato::new(atm, shl, ang, sym.to_owned(), loc, form.to_owned(), alps, coes);
                    sphs.push(ato);
                }
                let syms = ang2syms(ang, "car");
                for sym in syms {
                    let ato = Ato::new(atm, shl, ang, sym.to_owned(), loc, form.to_owned(), alps, coes);
                    cars.push(ato);
                }
            }
            Atos { raws, sphs, cars }
        })
    }

    pub fn nbas(&self) -> usize {
        let mut nbas = 0;
        let nshl = self.shls.len();
        for i in 0..nshl {
            let shl = self.shls.get(i);
            nbas += match shl.form.as_str() {
                "sph" => ANG_NUM_SPH[shl.ang],
                "car" => ANG_NUM_CAR[shl.ang],
                _ => 0,
            };
            nbas += self.shls.get(i).alps.len();
        }
        nbas
    }

    pub fn nshl(&self) -> usize {
        self.shls.len()
    }

    // 判断两个基组信息是否相等
    pub fn equal(&self, other: &Basis) -> bool {
        let datas_0 = &self.shls;
        let datas_1 = &other.shls;
        if datas_0.len() != datas_1.len() {
            log::error!("两个基组信息长度不相等");
            return false;
        }
        for i in 0..datas_0.len() {
            let each_1 = datas_0.get(i);
            let each_2 = datas_1.get(i);
            if each_1.atm != each_2.atm {
                log::error!("{},atm 不相等:{}!={}", i, each_1.atm, each_2.atm);
                return false;
            }
            if each_1.shl != each_2.shl {
                log::error!("{},shl 不相等:{}!={}", i, each_1.shl, each_2.shl);
                return false;
            }
            if each_1.ang != each_2.ang {
                log::error!("{},ang 不相等:{}!={}", i, each_1.ang, each_2.ang);
                return false;
            }
        }
        return true;
    }

    pub fn get_smat(&self, form: &str, atoms: &Atoms) -> &Array2<f64> {
        let smats = self.smat.get_or_init(|| {
            let int_caler = RsInt::new(atoms, self);
            let smat_car = int_caler.smat("car");
            let smat_sph = int_caler.mat1e_sph(&smat_car);
            let smat_raw = int_caler.mat1e_raw(&smat_car);
            [smat_raw, smat_sph, smat_car]
        });
        match form {
            "raw" => &smats[0],
            "sph" => &smats[1],
            "car" => &smats[2],
            _ => panic!("Invalid form"),
        }
    }

    pub fn atm_uplo(&self, form: &str) -> Vec<[usize; 2]> {
        self.shls.atm_uplo(form)
    }
}

impl Clone for Basis {
    fn clone(&self) -> Self {
        Self {
            shls: self.shls.clone(),
            atos: self.atos.clone(),
            bass: self.bass.clone(),
            smat: self.smat.clone(),
        }
    }
}

impl fmt::Display for Basis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{:<5}{:<5}{:<5}{:<5}{:>10}", "idx", "atm", "shl", "ang", "form")?;
        let nshl = self.shls.len();
        for i in 0..nshl {
            let shl = self.shls.get(i);
            let alps = &shl.alps;
            let coes = &shl.coes;
            writeln!(
                f,
                "{:<5}{:<5}{:<5}{:<5}{:>10}({})",
                i,
                shl.atm,
                shl.shl,
                shl.ang,
                shl.form,
                alps.len()
            )?;
            for i in 0..alps.len() {
                writeln!(f, "     {:>12.6E}  {:>12.6E}", alps[i], coes[i])?;
            }
        }
        Ok(())
    }
}

// 将笛卡尔重叠矩阵转为原始重叠矩阵
// pub fn smat_raw(smat_car: &Array2<f64>, basis: &Basis) -> Array2<f64> {
//     let shls = &basis.shls;
//     let nshl = shls.len();
//     let uplo_raw = shls.uplo("raw");
//     let uplo_car = shls.uplo("car");
//     let nato_raw = shls.nato("raw");
//     let mut smat_raw = Array2::zeros((nato_raw, nato_raw));
//     for i in 0..nshl {
//         let ang_i = shls.get(i).ang;
//         let [up_car_i, lo_car_i] = uplo_car[i];
//         let [up_raw_i, lo_raw_i] = uplo_raw[i];
//         let n_car_i = lo_car_i - up_car_i;
//         let n_sph_i = lo_raw_i - up_raw_i;
//         let tmati = btrans::car2sph_mat(ang_i);
//         for j in 0..nshl {
//             let angj = shls.get(j).ang;
//             let [up_car_j, lo_car_j] = uplo_car[j];
//             let [up_raw_j, lo_raw_j] = uplo_raw[j];
//             let n_car_j = lo_car_j - up_car_j;
//             let n_sph_j = lo_raw_j - up_raw_j;
//             let tmatj = btrans::sph2car_mat(angj);
//             let car_part = smat_car.slice(s![up_car_i..lo_car_i, up_car_j..lo_car_j]);
//             let raw_part: Array2<f64> = if n_car_i == n_sph_i && n_car_j == n_sph_j { car_part.to_owned() } else { tmati.dot(&car_part).dot(&tmatj) };
//             smat_raw.slice_mut(s![up_raw_i..lo_raw_i, up_raw_j..lo_raw_j]).assign(&raw_part);
//         }
//     }
//     smat_raw
// }

unsafe impl Send for Basis {}
unsafe impl Sync for Basis {}
