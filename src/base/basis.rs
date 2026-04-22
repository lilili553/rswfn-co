use ndarray::Array2;

use crate::base::mole::Mole;
use crate::datas::FACS;
use crate::datas::btrans;
use crate::datas::btrans::ang2syms;
use crate::datas::btrans::sym2lmn;
use crate::integ::Int_mat;
use core::fmt;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::OnceLock;

#[derive(Debug, Clone)]
pub struct Ato {
    pub atm: usize,
    pub shl: usize,
    pub ang: usize,
    pub sym: String,
    pub alps: Vec<f64>,
    pub coes: Vec<f64>,
    coes_norm: OnceLock<Vec<f64>>,
    lmn: OnceLock<[usize; 3]>,
}

impl Ato {
    pub fn new(atm: usize, shl: usize, ang: usize, sym: String, alps: &Vec<f64>, coes: &Vec<f64>) -> Self {
        let alps = alps.clone();
        let coes = coes.clone();
        let coes_norm = OnceLock::new();
        let lmn = OnceLock::new();
        Self { atm, shl, ang, sym, alps, coes, coes_norm, lmn }
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
pub struct MapAto {
    pub sphs: Vec<Ato>,
    pub cars: Vec<Ato>,
}

impl MapAto {
    // fn new(atms: Vec<usize>, xyzs: Vec<[f64; 3]>, lmns: Vec<[usize; 3]>, syms: Vec<String>, alps: Vec<Vec<f64>>, coes: Vec<Vec<f64>>, lyes: Vec<usize>, shls: Vec<usize>, nlms: Vec<(usize, usize, isize)>) -> Self {
    //     Self { atms, xyzs, lmns, syms, alps, coes, lyes, shls, nlms }
    // }

    pub fn len(&self, form: &str) -> usize {
        match form {
            "sph" => self.sphs.len(),
            "car" => self.cars.len(),
            _ => panic!("form error"),
        }
    }

    pub fn atms(&self, form: &str) -> Vec<usize> {
        match form {
            "sph" => self.sphs.iter().map(|ato| ato.atm).collect(),
            "car" => self.cars.iter().map(|ato| ato.atm).collect(),
            _ => panic!("form error"),
        }
    }

    pub fn shls(&self, form: &str) -> Vec<usize> {
        match form {
            "sph" => self.sphs.iter().map(|ato| ato.shl).collect(),
            "car" => self.cars.iter().map(|ato| ato.shl).collect(),
            _ => panic!("form error"),
        }
    }

    pub fn syms(&self, form: &str) -> Vec<String> {
        match form {
            "sph" => self.sphs.iter().map(|ato| ato.sym.clone()).collect(),
            "car" => self.cars.iter().map(|ato| ato.sym.clone()).collect(),
            _ => panic!("form error"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Shl {
    pub atm: usize,
    pub shl: usize,
    pub ang: usize,
    pub alps: Vec<f64>,
    pub coes: Vec<f64>,
}

impl Shl {
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
}

#[derive(Debug, Clone)]
pub struct MapShl {
    pub shls: Vec<Shl>,
    // pub cars: Vec<Shl>,
    // pub atms: Vec<usize>,
    // pub shls: Vec<usize>,
    // pub angs: Vec<usize>,
    // pub alps: Vec<Vec<f64>>,
    // pub coes: Vec<Vec<f64>>,
    pub uplo: OnceLock<[Vec<[usize; 2]>; 2]>,
}

const ANG_NUM_SPH: [usize; 6] = [1, 3, 5, 7, 9, 11];
const ANG_NUM_CAR: [usize; 6] = [1, 3, 6, 10, 15, 21];

impl MapShl {
    // pub fn new(atms: Vec<usize>, shls: Vec<usize>, angs: Vec<usize>, alps: Vec<Vec<f64>>, coes: Vec<Vec<f64>>) -> Self {
    //     let uplo = OnceLock::new();
    //     Self { atms, shls, angs, alps, coes, uplo }
    // }
    pub fn len(&self) -> usize {
        self.shls.len()
    }
    pub fn uplo(&self, form: &str) -> &Vec<[usize; 2]> {
        let uplo = self.uplo.get_or_init(|| {
            let mut uplo_sph = vec![];
            let mut uplo_car = vec![];
            let mut idx_sph = 0;
            let mut idx_car = 0;
            for shl in self.shls.iter() {
                let ang = shl.ang;
                let num_sph = ANG_NUM_SPH[ang];
                let num_car = ANG_NUM_CAR[ang];
                uplo_sph.push([idx_sph, idx_sph + num_sph]);
                uplo_car.push([idx_car, idx_car + num_car]);
                idx_sph += num_sph;
                idx_car += num_car;
            }
            [uplo_sph, uplo_car]
        });
        match form {
            "sph" => &uplo[0],
            "car" => &uplo[1],
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
}

#[derive(Debug, Clone)]
pub struct BasisData {
    pub atm: usize,
    pub shl: usize,
    pub ang: usize,
    pub alp: f64,
    pub coe: f64,
    pub nlm: (usize, usize, isize), //n,l,m 兼容真实原子轨道，先这么写，后面应该要搞成枚举类型
}

impl BasisData {
    pub fn new(atm: usize, shl: usize, ang: usize, alp: f64, coe: f64, nlm: (usize, usize, isize)) -> Self {
        Self { atm, shl, ang, alp, coe, nlm }
    }
}

#[derive(Debug, Clone)]
pub struct Basis {
    pub data: Vec<BasisData>, // atm,shl,ang,alp,coe,sph
    map_ato: OnceLock<MapAto>,
    map_shl: OnceLock<MapShl>,
    pub smat: OnceLock<[Array2<f64>; 2]>, // 存储球谐和笛卡尔的重叠积分
    pub otp: String,                      // 轨道类型 gto or sto
}

impl Basis {
    pub fn new(data: Vec<BasisData>, otp: Option<String>) -> Basis {
        let otp = match otp {
            Some(otp) => otp,
            None => "gto".to_string(),
        };
        let map_ato = OnceLock::new();
        let map_shl = OnceLock::new();
        let smat = OnceLock::new();
        Basis { data, map_ato, map_shl, smat, otp }
    }

    pub fn blank() -> Basis {
        let data = vec![];
        let map_ato = OnceLock::new();
        let map_shl = OnceLock::new();
        let smat = OnceLock::new();
        Basis { data, map_ato, map_shl, smat, otp: "gto".to_string() }
    }

    pub fn build(&mut self, data: Vec<BasisData>) {
        self.data = data;
    }

    pub fn is_none(&self) -> bool {
        self.data.len() == 0
    }

    pub fn map_bas(&self, mole: &Mole) -> (Vec<usize>, Vec<f64>, Vec<f64>, Vec<[usize; 3]>) {
        let mut atos = Vec::new();
        let mut alps = Vec::new();
        let mut coes = Vec::new();
        let mut lmns = Vec::new();
        // let mole = mole.upg.expect("未初始化分子");
        let ato_keys = mole.coefs().ato_keys("car");
        for each in self.data.iter() {
            for (_i, lmn) in btrans::ang2lmns(each.ang).iter().enumerate() {
                let sym = btrans::lmn2sym(lmn);
                let key = format!("{},{},{}", each.atm, each.shl, sym);
                let ato = match ato_keys.iter().position(|item| item == &key) {
                    Some(ato) => ato,
                    None => panic!("未找到索引:{key}"),
                };

                atos.push(ato as usize);
                alps.push(each.alp);
                coes.push(each.coe);
                lmns.push(*lmn);
            }
        }
        let mut idxs: Vec<usize> = (0..atos.len()).collect();
        idxs.sort_by_key(|idx| atos[*idx]);

        atos = idxs.iter().map(|idx| atos[*idx]).collect();
        alps = idxs.iter().map(|idx| alps[*idx]).collect();
        coes = idxs.iter().map(|idx| coes[*idx]).collect();
        lmns = idxs.iter().map(|idx| lmns[*idx]).collect();
        (atos, alps, coes, lmns)
    }

    pub fn get_by_atm(&self, atm: usize) -> Vec<BasisData> {
        let mut data = vec![];
        for each in self.data.iter() {
            if each.atm == atm {
                data.push(each.clone());
            }
        }
        data
    }
    pub fn ato_keys(&self) -> Vec<String> {
        let mut keys = vec![];
        for each in self.data.iter() {
            let syms = match each.ang {
                0 => vec!["S"],
                1 => vec!["PX", "PY", "PZ"],
                2 => btrans::CAR_D_SYMS.to_vec(),
                3 => btrans::CAR_F_SYMS.to_vec(),
                4 => btrans::CAR_G_SYMS.to_vec(),
                5 => btrans::CAR_H_SYMS.to_vec(),
                _ => panic!("不支持的角动量"),
            };
            for sym in syms.iter() {
                let key = format!("{},{},{}", each.atm, each.shl, sym);
                if !keys.contains(&key) {
                    keys.push(key);
                }
            }
        }
        keys
    }

    pub fn map_ato(&self) -> &MapAto {
        let shls = &self.map_shl().shls;
        let nshl = shls.len();
        self.map_ato.get_or_init(|| {
            let mut sphs = vec![];
            let mut cars = vec![];
            for i in 0..nshl {
                let atm = shls[i].atm;
                let shl = shls[i].shl;
                let ang = shls[i].ang;
                let alps = &shls[i].alps;
                let coes = &shls[i].coes;
                let syms = ang2syms(ang, "sph");
                for sym in syms {
                    let ato = Ato::new(atm, shl, ang, sym.to_owned(), alps, coes);
                    sphs.push(ato);
                }
                let syms = ang2syms(ang, "car");
                for sym in syms {
                    let ato = Ato::new(atm, shl, ang, sym.to_owned(), alps, coes);
                    cars.push(ato);
                }
            }
            MapAto { sphs, cars }
        })
    }

    pub fn map_shl(&self) -> &MapShl {
        self.map_shl.get_or_init(|| {
            let mut nshl = HashMap::new();
            let mut keys = vec![];
            let mut shls = vec![];
            for data in self.data.iter() {
                let atm = data.atm;
                let shl = data.shl;
                let ang = data.ang;
                let key = (atm, shl, ang);
                let alp = data.alp;
                let coe = data.coe;
                *nshl.entry(key).or_insert(0) += 1;
                if !keys.contains(&key) {
                    keys.push(key);
                    let alps = vec![alp];
                    let coes = vec![coe];
                    shls.push(Shl { atm, shl, ang, alps, coes });
                } else {
                    let num = shls.len();
                    shls[num - 1].alps.push(alp);
                    shls[num - 1].coes.push(coe);
                }
            }
            let mut nums: Vec<usize> = vec![];
            let mut uplo = vec![];
            let mut idx = 0;
            for key in keys.iter() {
                let num = nshl.get(key).unwrap();
                nums.push(*num);
                uplo.push([idx, idx + *num]);
                idx += *num;
            }
            let uplo = OnceLock::new();
            MapShl { shls, uplo }
        })
    }

    pub fn nbas(&self) -> usize {
        self.data.len()
    }
    pub fn nshl(&self) -> usize {
        let mut shls = Vec::new();
        for each in self.data.iter() {
            let key = (each.atm, each.shl);
            if !shls.contains(&key) {
                shls.push(key);
            }
        }
        shls.len()
    }

    // 判断两个基组信息是否相等
    pub fn equal(&self, other: &Basis) -> bool {
        let datas_0 = &self.data;
        let datas_1 = &other.data;
        if datas_0.len() != datas_1.len() {
            log::error!("两个基组信息长度不相等");
            return false;
        }
        for i in 0..datas_0.len() {
            let each_1 = &datas_0[i];
            let each_2 = &datas_1[i];
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
            if (each_1.alp - each_2.alp).abs() > 1e-3 {
                log::error!("{},alp 不相等:{}!={}", i, each_1.alp, each_2.alp);
                return false;
            }
            if (each_1.coe - each_2.coe).abs() > 1e-3 {
                log::error!("{},coe 不相等:{}!={}", i, each_1.coe, each_2.coe);
                return false;
            }
        }
        return true;
    }

    pub fn get_smat(&self, form: &str, mole: &Mole) -> &Array2<f64> {
        let smats = self.smat.get_or_init(|| {
            let int_caler = Int_mat::new(mole);
            let smat_car = int_caler.smat("car");
            let smat_sph = int_caler.mat1e_sph(&smat_car);
            [smat_sph, smat_car]
        });
        match form {
            "sph" => &smats[0],
            "car" => &smats[1],
            _ => panic!("Invalid form"),
        }
    }
}

impl fmt::Display for Basis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{:>5}{:>5}{:>5}{:>5}{:>14}{:>14}{:>14}", "idx", "atm", "shl", "ang", "alp", "coe", "nlm")?;
        for (i, each) in self.data.iter().enumerate() {
            let (n, l, m) = each.nlm;
            writeln!(f, "{:>5}{:>5}{:>5}{:>5}{:>14.8}{:>14.8}{:>8}{:2}{:>2}{:>2}", i, each.atm, each.shl, each.ang, each.alp, each.coe, "", n, l, m)?;
        }
        Ok(())
    }
}
