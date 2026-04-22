#![allow(non_snake_case)]
use std::collections::HashMap;

use ndarray::Array2;
use ndarray::{Array1, s};
use rayon::prelude::*;

use crate::atomprop;
use crate::datas::rdens;
use crate::gridprop;
use crate::integ;

use crate::base::{Mole, Stm, Vec3};
use crate::maths;
use crate::moleprop::newmole;
use crate::orbtprop;
use crate::orbtprop::obtmat::Mocv;
// use crate::utils;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    pub fn spin(&self, ctype: &str) -> Vec<[f64; 3]> {
        let natm = self.mole.atoms().len();
        if !self.mole.is_open() {
            log::warn!("非开壳层分子无法计算自旋");
            return vec![[0.0; 3]; natm];
        }
        let mole_caler = newmole::Calculator::new(self.mole);
        let (mole_a, mole_b) = mole_caler.mole_spin();

        let caler_a = Calculator::new(&mole_a);
        let caler_b = Calculator::new(&mole_b);

        let (eles_a, eles_b) = (caler_a.eles(ctype), caler_b.eles(ctype));

        let mut spin = vec![[0.0; 3]; natm];
        for i in 0..natm {
            spin[i][0] = eles_a[i];
            spin[i][1] = eles_b[i];
            spin[i][2] = eles_a[i] - eles_b[i];
        }
        spin
    }

    pub fn mulliken(&self) -> Vec<f64> {
        let dmat = self.mole.coefs().get_dmat("sph");
        let smat = self.mole.basis().get_smat("sph", self.mole); // 重叠矩阵
        let uplo = self.mole.atoms().ato_uplo("sph", self.mole);
        mulliken(dmat, smat, &uplo)
    }

    // 计算lowdin 电荷，返回电子数量
    pub fn lowdin(&self) -> Vec<f64> {
        let dmat = self.mole.coefs().get_dmat("sph");
        let smat = self.mole.basis().get_smat("sph", self.mole);
        let uplo = self.mole.atoms().ato_uplo("sph", self.mole);
        lowdin(dmat, smat, &uplo, None)
    }

    pub fn hirshfeld(&self) -> Vec<f64> {
        let gridClaer = gridprop::dftgrid::Calculator::new(self.mole);
        let densCaler = gridprop::density::Calculator::new(self.mole);
        let atoms = self.mole.atoms();
        let natm = atoms.len();
        let mut eles = vec![0.0; natm]; // 原子的电子数
        // 计算原子与格点之间的距离
        for i in 0..natm {
            let (atm_grids, atm_weits) = gridClaer.atm_grids(i, 75, 434);

            let ngrid = atm_grids.len();
            let mut pro_dens: Array1<f64> = Array1::zeros(ngrid);
            let mut sum_dens: Array1<f64> = Array1::zeros(ngrid);
            let atms = self.mole.atoms().idxs();
            let mol_dens: Vec<f64> = (0..ngrid).into_par_iter().map(|g| densCaler.mol_rho_dm(&atm_grids[g], 0).0).collect();
            let mol_dens = Array1::from_vec(mol_dens);

            for j in 0..natm {
                let atom_j = atoms.get(j as usize);

                let fdens: Vec<f64> = (0..ngrid)
                    .into_par_iter()
                    .map(|g| {
                        let dist = maths::calc_dist(&atm_grids[g], &atom_j.xyz);
                        rdens::get_radDens_v1(atom_j.nuc, dist)
                    })
                    .collect();

                let fdens = Array1::from_vec(fdens.to_vec());

                if i == j {
                    pro_dens = fdens.clone();
                }
                sum_dens += &fdens;
            }
            let atom_i = atoms.get(i as usize);
            eles[i] = (0..ngrid)
                .into_par_iter()
                .map(|g| {
                    let dist = maths::calc_dist(&atm_grids[g], &atom_i.xyz);
                    if dist > 10.0 {
                        return 0.0;
                    } else if sum_dens[g].abs() < 1e-10 {
                        return 0.0;
                    } else {
                        let weit = pro_dens[g] / sum_dens[g];
                        return weit * mol_dens[g] * atm_weits[g];
                    }
                })
                .collect::<Vec<f64>>()
                .iter()
                .sum::<f64>();

            log::info!("原子{:>3}  电子数:{:>10.4}", i, eles[i]);
        }
        eles
    }

    pub fn eles(&self, ctype: &str) -> Vec<f64> {
        match ctype {
            "mulliken" => self.mulliken(),
            "lowdin" => self.lowdin(),
            "hirshfeld" => self.hirshfeld(),
            _ => panic!("未知的计算方法"),
        }
    }

    // 分子轨道内原子电子布局
    pub fn eles_obt(&self, ctype: &str, obt: usize) -> Vec<f64> {
        let smat = self.mole.basis().get_smat("sph", self.mole);
        let old_cmat = self.mole.coefs().get_cmat("sph");
        let nato = old_cmat.nrows();
        let nobt = old_cmat.ncols();
        let mut new_cmat = Array2::zeros((nato, nobt));
        new_cmat.slice_mut(s![.., obt]).assign(&old_cmat.slice(s![.., obt]));
        let dmat = old_cmat.dot(&new_cmat.t());
        let uplo = self.mole.atoms().ato_uplo("sph", self.mole);
        mulliken(&dmat, smat, &uplo)
    }

    // 计算经过pocv方法改造过后的分子轨道的电子布居
    pub fn pocv(&self, dirs: &HashMap<usize, Vec3>, keep_other_atm: bool, keep_other_sym: bool, ctype: &str) -> Vec<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_pocv = obt_caler.pocv(dirs, keep_other_atm, keep_other_sym);

        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.set_cmat("sph", &cmat_pocv);
        let ele_caler = Calculator::new(&mole_pocv);
        match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        }
    }

    pub fn mocv(
        &self,
        mocvs: &HashMap<usize, Mocv>, // 每个原子的保留信息
        ctype: &str,
    ) -> Vec<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_mocv = obt_caler.mocv(mocvs);

        let mut mole_mocv = self.mole.clone();
        mole_mocv.coefs.set_cmat("sph", &cmat_mocv);
        let ele_caler = Calculator::new(&mole_mocv);
        match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        }
    }

    // 多局部坐标系的MOCV算法
    pub fn mocv_mult(&self, atm: usize) {
        let nebs = self.mole.atoms().nebs(atm);
    }

    pub fn pi_pocv(&self, ctype: &str) -> (HashMap<usize, Vec3>, Vec<f64>) {
        let mut obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        obt_caler.atms = self.atms;
        let (dirs, cmat_pocv) = obt_caler.pi_pocv();

        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.set_cmat("sph", &cmat_pocv);
        let ele_caler = Calculator::new(&mole_pocv);
        let eles = match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        };
        (dirs, eles)
    }

    pub fn pi_mocv(&self, ctype: &str) -> (HashMap<usize, Mocv>, Vec<f64>) {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let (mocvs, cmat_mocv) = obt_caler.pi_mocv();
        let mut mole_mocv = self.mole.clone();
        mole_mocv.coefs.set_cmat("sph", &cmat_mocv);
        let ele_caler = Calculator::new(&mole_mocv);
        let eles = match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        };
        (mocvs, eles)
    }

    // /// 计算三个垂直与键轴的
    // pub fn pi_mocv_split(&self, atm: usize) -> f64 {
    //     let mut dir_caler = atomprop::direction::Calculator::new(self.mole);
    //     let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
    //     let atoms = self.mole.atoms();
    //     let nebs = atoms.nebs(atm);
    //     let obts = self.mole.coefs().obt_idxs(true);
    //     let smat = self.mole.coefs().get_smat("car", self.mole);
    //     let nato = smat.nrows();
    //     let mut dmat = Array2::<f64>::zeros((nato, nato));
    //     let uplo = self.mole.atoms().ato_uplo("car", self.mole);
    //     let [up, lo] = uplo[atm];
    //     for idx in nebs.iter() {
    //         // let dir_caler.LCS(atm, *idx);
    //         let stm = dir_caler.LCS(atm, Some(*idx));
    //         if let Some(stm) = stm {
    //             let mocv = Mocv::new(stm, [0], [0, 0, 1], [0, 0, 1, 0, 1, 0]);
    //             let mut mocvs = HashMap::with_capacity(2);
    //             mocvs.insert(atm, mocv.clone());
    //             mocvs.insert(*idx, mocv);
    //             let cmat_mocv = obt_caler.mocv(&mocvs);
    //             let mut dmat_mocv = maths::calc_dmat(&cmat_mocv, obts);
    //             dmat_mocv.slice_mut(s![up..lo, up..lo]).fill(0.0);
    //             dmat = dmat + dmat_mocv;
    //             // obt_caler.mocv()
    //         }
    //     }
    //     // let smat = smat - Array2::eye(nato); // 不考虑自己和自己的重叠
    //     let uplo = &self.mole.atoms().ato_uplo("car", self.mole);
    //     let omat = mulliken(&dmat, &smat, uplo);
    //     omat[atm]
    // }

    pub fn pi_ldao(&self) -> (HashMap<usize, Stm>, Vec<f64>) {
        let mut stms = HashMap::new();
        let mut dir_caler = atomprop::direction::Calculator::new(self.mole);
        let natm = self.mole.atoms().len();
        for atm in 0..natm {
            if let Some(stm) = dir_caler.LCS(atm, None) {
                stms.insert(atm, stm);
            } else {
                log::warn!("原子{}没有局部坐标系", atm);
            }
        }
        let cmat_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat = cmat_caler.ldao(&stms); // 局部坐标系的分子轨道
        let idxs = self.mole.coefs().obt_idxs(true);
        let dmat = maths::calc_dmat(&cmat, idxs, self.mole.is_open());
        let smat = self.mole.basis().get_smat("car", self.mole);
        let syms = self.mole.coefs().ato_syms("car");
        let atms = self.mole.coefs().ato_atms("car");
        let atoms = self.mole.atoms();
        let mut atos = vec![false; syms.len()];
        for (i, sym) in syms.iter().enumerate() {
            if sym == "PZ" {
                atos[i] = true;
            }
            if sym == "PY" && atoms.is_line(atms[i]) {
                atos[i] = true;
            }
            if sym == "XZ" || sym == "YZ" {
                atos[i] = true;
            }
        }
        let uplo = self.mole.atoms().ato_uplo("sph", self.mole);
        let eles = lowdin(&dmat, &smat, &uplo, Some(&atos));
        (stms, eles)
    }
}

pub fn mulliken(dmat: &Array2<f64>, smat: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Vec<f64> {
    assert_eq!(dmat.shape(), smat.shape());
    let qmat = dmat * smat;
    let nato = dmat.nrows(); // 基函数数量
    let natm = uplo.len();
    let mut eles = vec![0.0; natm];
    for a in 0..natm {
        let [u, l] = uplo[a];
        for i in u..l {
            for j in 0..nato {
                eles[a] += qmat[(i, j)];
            }
        }
    }
    eles
}

/// 计算lowdin 电荷，返回电子数量
/// cmat：分子轨道系数矩阵
/// smat：重叠矩阵
/// nums：每个原子的轨道数量
/// idxs：占据轨道的索引
/// atos：原子轨道是否保留
pub fn lowdin(dmat: &Array2<f64>, smat: &Array2<f64>, uplo: &Vec<[usize; 2]>, atos: Option<&Vec<bool>>) -> Vec<f64> {
    let nato = dmat.nrows();
    let atos = match atos {
        Some(atos) => atos,
        None => &vec![true; nato],
    };
    let smat_half = maths::matrix_phalf(&smat); // 重叠矩阵的平方根

    let qmat = smat_half.dot(dmat).dot(&smat_half);
    let natm = uplo.len();
    let mut eles = vec![0.0; natm];
    for a in 0..natm {
        let [ua, la] = uplo[a];
        for i in ua..la {
            if atos[i] == true {
                eles[a] += qmat[(i, i)];
            }
        }
    }
    eles
}
