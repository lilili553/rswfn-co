use ndarray::{Array1, Array2, s};

use crate::base::Mole;
use crate::moleprop;
use crate::orbtprop;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Calculator {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    pub fn mole_spin(&self) -> (Mole, Mole) {
        let cmat = self.mole.coefs().get_cmat("raw", self.mole.basis());
        let [nato, nobt] = cmat.shape().try_into().expect("获取矩阵形状失败");
        let occs = &self.mole.coefs().obt_occs;
        if self.mole.is_open() {
            let mut mole_a = self.mole.clone();
            let mut mole_b = self.mole.clone();

            let mut cmat_a = cmat.clone();
            let mut cmat_b = cmat.clone();

            cmat_a.slice_mut(s![.., nobt / 2..]).fill(0.0);
            cmat_b.slice_mut(s![.., ..nobt / 2]).fill(0.0);

            mole_a.coefs.as_mut().unwrap().set_cmat("raw", &cmat_a, self.mole.basis());
            mole_b.coefs.as_mut().unwrap().set_cmat("raw", &cmat_b, self.mole.basis());
            (mole_a, mole_b)
        } else {
            let occs_a: Vec<f64> = occs.iter().map(|occ| if *occ > 0.0 { 1.0 } else { 0.0 }).collect();
            let occs_b: Vec<f64> = occs.iter().map(|occ| if *occ > 0.0 { 1.0 } else { 0.0 }).collect();
            let mut mole_a = self.mole.clone();
            let mut mole_b = self.mole.clone();
            mole_a.coefs.as_mut().unwrap().obt_occs = occs_a;
            mole_b.coefs.as_mut().unwrap().obt_occs = occs_b;
            (mole_a, mole_b)
        }
    }

    pub fn mole_pi(&self) -> Mole {
        let mut mole_caler = orbtprop::obtmat::Calculator::new(self.mole);
        mole_caler.atms = self.atms;
        let (_dirs, cmat_pocv) = mole_caler.pi_pocv();
        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, self.mole.basis());
        mole_pocv
    }

    /// 生成局部坐标系的分子
    pub fn mole_lcao(&self) {
        // 将原始的系数矩阵转为局部坐标系的系数矩阵
    }
}
