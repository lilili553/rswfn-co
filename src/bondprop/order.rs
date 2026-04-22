use ndarray::{Array2, s};
use std::collections::HashMap;

use crate::base::{Mole, Vec3};
use crate::datas::consts::BOHR;
use crate::orbtprop::obtmat::Mocv;
use crate::{maths, orbtprop};

pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    pub fn order(&self, otype: &str) -> Array2<f64> {
        let smat = self.mole.basis().get_smat("sph", self.mole);
        let dmat = self.mole.coefs().get_dmat("sph");
        assert_eq!(smat.shape(), dmat.shape());
        let uplo = self.mole.atoms().ato_uplo("sph", self.mole);
        match otype {
            "mayer" => mayer(dmat, smat, &uplo),
            "lowdin" => lowdin(dmat, smat, &uplo),
            "wiberg" => wiberg(dmat, &uplo),
            _ => panic!("未知键级类型:{}", otype),
        }
    }

    pub fn mayer_obt(&self, obt: usize) -> Array2<f64> {
        let smat = self.mole.basis().get_smat("sph", self.mole);
        let old_cmat = self.mole.coefs().get_cmat("sph");
        let nato = old_cmat.nrows();
        let nobt = old_cmat.ncols();
        let mut new_cmat = Array2::zeros((nato, nobt));
        new_cmat.slice_mut(s![.., obt]).assign(&old_cmat.slice(s![.., obt]));
        let dmat = old_cmat.dot(&new_cmat.t());
        let uplo = self.mole.atoms().ato_uplo("sph", self.mole);
        mayer(&dmat, smat, &uplo)
    }

    pub fn pocv(&self, dirs: &HashMap<usize, Vec3>, keep_other_atm: bool, keep_other_sym: bool, otype: &str) -> Array2<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_pocv = obt_caler.pocv(dirs, keep_other_atm, keep_other_sym);

        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.set_cmat("sph", &cmat_pocv);
        let caler = Calculator::new(&mole_pocv);
        caler.order(otype)
    }

    pub fn mocv(
        &self,
        mocvs: &HashMap<usize, Mocv>, // 每个原子的保留信息
        otype: &str,
    ) -> Array2<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_mocv = obt_caler.mocv(mocvs);

        let mut mole_mocv = self.mole.clone();
        mole_mocv.coefs.set_cmat("sph", &cmat_mocv);
        let caler = Calculator::new(&mole_mocv);
        caler.order(otype)
    }

    pub fn pi_pocv(&self) -> (HashMap<usize, Vec3>, Array2<f64>) {
        let mut obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        obt_caler.atms = self.atms;
        let (dirs, cmat_pocv) = obt_caler.pi_pocv();

        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.set_cmat("sph", &cmat_pocv);

        let caler = Calculator::new(&mole_pocv);
        let mut omat = caler.order("mayer");
        omat = omat.abs().sqrt(); // 取绝对值再开根号
        (dirs, omat)
    }

    pub fn pi_mocv(&self) -> (HashMap<usize, Mocv>, Array2<f64>) {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let (mocvs, cmat_mocv) = obt_caler.pi_mocv();
        let mut mole_mocv = self.mole.clone();
        mole_mocv.coefs.set_cmat("sph", &cmat_mocv);
        let caler = Calculator::new(&mole_mocv);
        let mut omat = caler.order("mayer");
        omat = omat.abs().sqrt(); // 取绝对值再开根号
        (mocvs, omat)
    }

    /// 计算束缚键级
    pub fn bound(&self, atm: usize, dir: &Vec3, otype: &str) -> (Vec<usize>, Vec<f64>) {
        log::info!("计算束缚键级");
        let mut dir_map = HashMap::new();
        dir_map.insert(atm, dir.clone());

        let omat = self.pocv(&dir_map, true, true, otype); // 获得键级矩阵
        let mut orders = Vec::new();
        // let nebs = self.mole.atoms().nebs(atm);
        // 使用距离小于5.5埃的键进行计算
        let mut nebs = Vec::new();
        for idx in self.mole.atoms().idxs() {
            let dist = self.mole.atoms().dist(atm, *idx);
            if dist < 5.5 * BOHR && *idx != atm {
                nebs.push(*idx);
            }
        }
        for idx in nebs.iter() {
            let order = omat[(atm as usize, *idx as usize)];
            log::info!("order.bound:{atm}-{idx}:{order}");
            orders.push(order);
        }

        (nebs, orders)
    }
}

pub fn mayer(dmat: &Array2<f64>, smat: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Array2<f64> {
    let natm = uplo.len();
    let psmat = dmat.dot(smat);
    let bomat = &psmat * &psmat.t();
    let mut omat = Array2::zeros((natm, natm)); // 键级矩阵

    for i in 0..natm {
        let [ui, li] = uplo[i];
        for j in 0..natm {
            let [uj, lj] = uplo[j];
            let order = bomat.slice(s![ui..li, uj..lj]).sum();
            omat[(i, j)] = order;
        }
    }
    omat
}

pub fn wiberg(dmat: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Array2<f64> {
    let p2mat = dmat * dmat;
    let natm = uplo.len();
    let mut omat = Array2::zeros((natm, natm));
    for i in 0..natm {
        let [ui, li] = uplo[i];
        for j in 0..natm {
            let [uj, lj] = uplo[j];
            let order = p2mat.slice(s![ui..li, uj..lj]).sum();
            omat[(i, j)] = order;
        }
    }
    omat
}

pub fn lowdin(dmat: &Array2<f64>, smat: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Array2<f64> {
    let smat_sqrt = maths::matrix_phalf(smat);
    let pmat_lowd = smat_sqrt.dot(dmat).dot(&smat_sqrt);
    let natm = uplo.len();
    let mut omat = Array2::zeros((natm, natm));
    for i in 0..natm {
        let [ui, li] = uplo[i];
        for j in 0..natm {
            let [uj, lj] = uplo[j];
            let order = pmat_lowd.slice(s![ui..li, uj..lj]).sum();
            omat[(i, j)] = order;
        }
    }
    omat
}
