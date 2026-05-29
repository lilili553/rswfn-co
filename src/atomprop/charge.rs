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
// use crate::orbtprop::obtmat::Mocv;
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
        let atoms = self.mole.atoms();
        let basis = self.mole.basis();
        let coefs = self.mole.coefs();
        let dmat = coefs.get_dmat("raw", basis);
        let smat = basis.get_smat("raw", atoms); // 重叠矩阵
        let uplo = coefs.atm_uplo("raw", basis);
        mulliken(dmat, smat, &uplo)
    }

    // 计算lowdin 电荷，返回电子数量
    pub fn lowdin(&self) -> Vec<f64> {
        let dmat = self.mole.coefs().get_dmat("raw", self.mole.basis());
        let smat = self.mole.basis().get_smat("raw", self.mole.atoms());
        let uplo = self.mole.coefs().atm_uplo("raw", self.mole.basis());
        lowdin(dmat, smat, &uplo, None)
    }

    pub fn hirshfeld(&self) -> Vec<f64> {
        let gridClaer = gridprop::dftgrid::Calculator::new(self.mole);
        let densCaler = gridprop::density::Calculator::new(self.mole);
        let atoms = self.mole.atoms();
        let natm = atoms.len();
        let mut eles = vec![0.0; natm]; // 原子的电子数
        let (mole_grids, mole_weits) = gridClaer.mol_grids(75, 434);
        let ngrid = mole_grids.len();
        let mole_dens: Vec<f64> = mole_grids
            .par_iter()
            .map(|grid| densCaler.mol_rho_dm(grid, 0).val0)
            .collect();
        let pro_mole_dens: Vec<f64> = mole_grids.par_iter().map(|grid| densCaler.pro_rho(grid)).collect();
        let mut pro_atom_dens = vec![0.0; ngrid];
        let mut atom_weits = vec![0.0; ngrid];

        // 计算原子与格点之间的距离
        for i in 0..natm {
            let nuc = atoms.get(i as usize).nuc;
            let [ax, ay, az] = atoms.get(i as usize).xyz;
            pro_atom_dens = mole_grids
                .par_iter()
                .map(|grid| {
                    let [gx, gy, gz] = grid;
                    let [dx, dy, dz] = [gx - ax, gy - ay, gz - az];
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                    rdens::get_rad_dens_v1(nuc, dist)
                })
                .collect();
            atom_weits = (0..ngrid)
                .map(|g| {
                    if pro_mole_dens[g] > 1e-12 {
                        pro_atom_dens[g] / pro_mole_dens[g]
                    } else {
                        0.0
                    }
                })
                .collect();
            for g in 0..ngrid {
                eles[i] += mole_dens[g] * mole_weits[g] * atom_weits[g];
            }

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
        let smat = self.mole.basis().get_smat("raw", self.mole.atoms());
        let old_cmat = self.mole.coefs().get_cmat("raw", self.mole.basis());
        let nato = old_cmat.nrows();
        let nobt = old_cmat.ncols();
        let mut new_cmat = Array2::zeros((nato, nobt));
        let nocc = if self.mole.is_open() { 1.0 } else { 2.0 };
        new_cmat.slice_mut(s![.., obt]).assign(&old_cmat.slice(s![.., obt]));
        let dmat = old_cmat.dot(&new_cmat.t()) * nocc;
        let uplo = self.mole.coefs().atm_uplo("raw", self.mole.basis());
        mulliken(&dmat, smat, &uplo)
    }

    // 计算经过pocv方法改造过后的分子轨道的电子布居
    pub fn pocv(
        &self,
        dirs: &HashMap<usize, Vec3>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        ctype: &str,
    ) -> Vec<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_pocv = obt_caler.pocv(dirs, keep_other_atm, keep_other_sym);

        let mut mole_pocv = self.mole.clone();
        let basis = mole_pocv.basis.as_ref().unwrap();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, basis);
        let ele_caler = Calculator::new(&mole_pocv);
        match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        }
    }

    pub fn mocv(&self, stms: &HashMap<usize, Stm>, atos: &Vec<usize>, ctype: &str) -> Vec<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_mocv = obt_caler.mocv(stms, atos);

        let mut mole_mocv = self.mole.clone();
        let basis = mole_mocv.basis.as_ref().unwrap();
        mole_mocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_mocv, basis);
        let ele_caler = Calculator::new(&mole_mocv);
        match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        }
    }

    pub fn pi_pocv(&self, ctype: &str) -> (HashMap<usize, Vec3>, Vec<f64>) {
        let mut obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        obt_caler.atms = self.atms;
        let (dirs, cmat_pocv) = obt_caler.pi_pocv();

        let mut mole_pocv = self.mole.clone();
        let basis = mole_pocv.basis.as_ref().unwrap();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, basis);
        let ele_caler = Calculator::new(&mole_pocv);
        let eles = match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        };
        (dirs, eles)
    }

    pub fn pi_mocv(&self, ctype: &str) -> (HashMap<usize, Stm>, Vec<usize>, Vec<f64>) {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let (stms, atos, cmat_mocv) = obt_caler.pi_mocv(0);
        let mut mole_mocv = self.mole.clone();
        let basis = mole_mocv.basis.as_ref().unwrap();
        mole_mocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_mocv, basis);
        let caler = Calculator::new(&mole_mocv);
        let eles = match ctype {
            "mulliken" => caler.mulliken(),
            "lowdin" => caler.lowdin(),
            "hirshfeld" => caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        };
        (stms, atos, eles)
    }

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
        let occs = self.mole.coefs().obt_occs.clone();
        let dmat = maths::calc_dmat(&cmat, &occs);
        let smat = self.mole.basis().get_smat("car", self.mole.atoms());
        let syms = self.mole.basis().atos().syms("car");
        let atms = self.mole.basis().atos().atms("car");
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
        let uplo = self.mole.coefs().atm_uplo("raw", self.mole.basis());
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

pub fn hirshfeld(
    atom_nucs: Vec<usize>,
    atom_xyzs: Vec<[f64; 3]>,
    mole_grids: &Vec<[f64; 3]>,
    mole_weits: &Vec<f64>,
    mole_dens: &Vec<f64>,
    pro_mole_dens: &Vec<f64>,
) -> Vec<f64> {
    let natm = atom_nucs.len();
    let mut eles = vec![0.0; natm]; // 原子的电子数
    let ngrid = mole_grids.len();
    let mut pro_atom_dens = vec![0.0; ngrid];
    let mut atom_weits = vec![0.0; ngrid];

    // 计算原子与格点之间的距离
    for i in 0..natm {
        let nuc = atom_nucs[i];
        let [ax, ay, az] = atom_xyzs[i];
        pro_atom_dens = mole_grids
            .par_iter()
            .map(|grid| {
                let [gx, gy, gz] = grid;
                let [dx, dy, dz] = [gx - ax, gy - ay, gz - az];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                rdens::get_rad_dens_v1(nuc, dist)
            })
            .collect();
        atom_weits = (0..ngrid)
            .map(|g| {
                if pro_mole_dens[g] > 1e-12 {
                    pro_atom_dens[g] / pro_mole_dens[g]
                } else {
                    0.0
                }
            })
            .collect();
        for g in 0..ngrid {
            eles[i] += mole_dens[g] * mole_weits[g] * atom_weits[g];
        }

        log::info!("原子{:>3}  电子数:{:>10.4}", i, eles[i]);
    }
    eles
}

fn check(res: Vec<f64>, ans: Vec<f64>) {
    assert_eq!(res.len(), ans.len());
    for i in 0..res.len() {
        assert!(
            (res[i] - ans[i]).abs() < 1e-4,
            "res[{}] = {}, ans[{}] = {}",
            i,
            res[i],
            i,
            ans[i]
        );
    }
}

#[test]
fn test_charge() {
    let path = String::from(r"E:\code\rswfn\moles\C6H6.fch");
    let mole = Mole::from_file(&path);

    let caler = Calculator::new(&mole);
    let res = caler.mulliken();
    let ans = vec![
        6.12846023, 6.12846024, 6.12846023, 6.12846023, 6.12846024, 6.12846023, 0.87153977, 0.87153977, 0.87153977,
        0.87153977, 0.87153977, 0.87153977,
    ];
    check(res, ans);
    let res = caler.lowdin();
    let ans = vec![
        6.15647202, 6.15982449, 6.15982449, 6.15647202, 6.15982449, 6.15982449, 0.84289460, 0.84049220, 0.84049220,
        0.84289460, 0.84049220, 0.84049220,
    ];
    check(res, ans);
    let mut res = caler.hirshfeld();
    let res = (0..res.len()).map(|atm| mole.atoms.get(atm).crg - res[atm]).collect();
    let ans = vec![
        -0.04268004,
        -0.04261213,
        -0.04261213,
        -0.04268004,
        -0.04261213,
        -0.04261213,
        0.04262433,
        0.04263998,
        0.04263999,
        0.04262433,
        0.04263998,
        0.04263999,
    ];
    check(res, ans);
    let (dirs, vals) = caler.pi_pocv("mulliken");
    let ans = vec![
        0.98514, 0.98514, 0.98514, 0.98514, 0.98514, 0.98514, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    println!("pi_pocv:{:?}", vals);
    check(vals, ans);
    let (dirs, atos, vals) = caler.pi_mocv("mulliken");
    let ans = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    println!("pi_mocv:{:?}", vals);
    check(vals, ans);
    // 测试自旋
    let path = String::from(r"E:\code\rswfn\moles\C6H6_N.fch");
    let mole = Mole::from_file(&path);
    let caler = Calculator::new(&mole);
    let res = caler.spin("mulliken").iter().map(|row| row[2]).collect();
    let ans = vec![
        -0.11890, 0.32872, 0.32872, -0.11890, 0.32872, 0.32872, 0.00392, -0.02123, -0.02123, 0.00392, -0.02123,
        -0.02123,
    ];
    check(res, ans);
}
