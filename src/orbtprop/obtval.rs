use std::collections::HashMap;

use ndarray::Axis;

use crate::atomprop;
use crate::base::Mole;
use crate::integ;
use crate::maths;
use crate::orbtprop;
use crate::utils;
pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    // 计算分子轨道的pi电子数
    pub fn pi_eles_pocv(&self) -> Vec<f64> {
        let mut caler = orbtprop::obtmat::Calculator::new(self.mole);
        caler.vobt = true;
        caler.atms = self.atms;
        let (_, cmat) = caler.pi_pocv();
        let smat = self.mole.basis().get_smat("raw", self.mole.atoms());
        let mut eles = Vec::new();
        let nobt = cmat.ncols();
        let k = if self.mole.is_open() { 1.0 } else { 2.0 };
        for i in 0..nobt {
            let cvec = cmat.column(i).insert_axis(Axis(1));
            let dmat = cvec.dot(&cvec.t());
            let qmat = dmat * smat;
            eles.push(qmat.sum() * k);
        }
        eles
    }

    // 计算分子轨道的pi电子数
    pub fn pi_eles_ldao(&self) -> Vec<f64> {
        let mut stms = HashMap::new();
        let mut dir_caler = atomprop::direction::Calculator::new(self.mole);
        let natm = self.mole.atoms().len();
        for atm in 0..natm {
            if let Some(stm) = dir_caler.LCS(atm, None) {
                stms.insert(atm, stm);
            }
        }
        let cmat_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let smat_caler = integ::RsInt::new(self.mole.atoms(), self.mole.basis());
        let cmat = cmat_caler.ldao(&stms); // 局部坐标系的分子轨道
        let smat = self.mole.basis().get_smat("raw", self.mole.atoms());
        let smat_half = maths::matrix_phalf(&smat);
        let k = if self.mole.is_open() { 1.0 } else { 2.0 };
        let mut eles = Vec::new();
        let syms = self.mole.basis().atos().syms("raw");
        let mut ipzs = vec![false; syms.len()];
        for (i, sym) in syms.iter().enumerate() {
            if sym.to_string().as_str() == "PZ" {
                ipzs[i] = true;
            }
        }
        let nobt = cmat.ncols();
        let ctype = "lowdin";
        for obt in 0..nobt {
            let cvec = cmat.column(obt).insert_axis(Axis(1));
            let dmat = cvec.dot(&cvec.t());
            if ctype == "mulliken" {
                let mut qmat = dmat * smat; // 元素乘积
                for (ato, ipz) in ipzs.iter().enumerate() {
                    if *ipz == false {
                        qmat.row_mut(ato).fill(0.0);
                    }
                }
                eles.push(qmat.sum() * k);
            } else if ctype == "lowdin" {
                let mut qmat = smat_half.dot(&dmat).dot(&smat_half);
                for (ato, ipz) in ipzs.iter().enumerate() {
                    if *ipz == false {
                        qmat[(ato, ato)] = 0.0;
                    }
                }
                eles.push(qmat.diag().sum() * k);
            } else {
                panic!("未知方法")
            }
        }
        eles
    }

    pub fn eles(&self) -> Vec<f64> {
        // let mut caler = orbtprop::obtmat::Calculator::new(self.mole);
        // caler.vobt = true;
        // let atms = self.mole.atoms().idxs();
        // let (dirs, cmat) = caler.pi_pocv(atms);
        let cmat = self.mole.coefs().get_cmat("sph", self.mole.basis());
        let smat = self.mole.basis().get_smat("sph", self.mole.atoms());
        let mut eles = Vec::new();
        let nobt = cmat.ncols();
        let k = if self.mole.is_open() { 1.0 } else { 2.0 };
        for i in 0..nobt {
            let cvec = cmat.column(i).insert_axis(Axis(1));
            let dmat = cvec.dot(&cvec.t());
            let qmat = dmat * smat;
            eles.push(qmat.sum() * k);
        }
        eles
    }
}
