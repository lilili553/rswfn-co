// 半经验方法，直接忽略双电子积分先

use std::collections::HashMap;

use crate::base::{Coefs, Mole};
use crate::integ::RsInt;
use crate::maths;
use crate::{atomprop, utils};
use ndarray::{Array1, Array2};

pub struct RSE<'m> {
    pub mole: &'m Mole,
}

impl<'m> RSE<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        Self { mole }
    }

    fn calc_ele_eng(
        &self,
        dmat: &Array2<f64>, // 重叠矩阵
        tmat: &Array2<f64>, // 势能矩阵
        vmat: &Array2<f64>, // 动能矩阵
        fmat: &Array2<f64>, // fork矩阵
    ) -> f64 {
        (dmat * (tmat + vmat + fmat)).sum() * 0.5
    }

    // 实现RSE方法，忽略双电子积分，直接对角化核心哈密顿矩阵
    pub fn calc_rse(&self, lcs: bool) -> (f64, Coefs) {
        let atos = self.mole.basis().atos();
        // let nato = atos.len();
        let mut dir_caler = atomprop::direction::Calculator::new(self.mole);
        let atoms = self.mole.atoms();

        let ints = RsInt::new(self.mole.atoms(), self.mole.basis());
        let mut stms = HashMap::new();
        let [smat, tmat, vmat] = if lcs {
            for atm in 0..atoms.len() {
                let nebs = atoms.nebs(atm);
                let stm = match dir_caler.LCS(atm, Some(nebs[0])) {
                    Some(stm) => stm,
                    None => dir_caler.RLCS(atm),
                };
                stms.insert(atm, stm);
            }
            ints.mat1e_lcs(&stms, [1, 1, 1])
        } else {
            ints.mat1e([1, 1, 1], "car")
        };
        let hmat = tmat + vmat;

        let xmat = maths::matrix_nhalf(&smat);
        let hmat_ = xmat.dot(&hmat).dot(&xmat);
        let (engs, cmat_) = maths::eigh(&hmat_);
        let cmat = xmat.dot(&cmat_);
        let [nela, nelb] = self.mole.atoms.neles;
        let mut eng = 0.0;
        for i in 0..nela {
            eng += engs[i];
        }
        let mut mole_c = self.mole.clone();
        if stms.len() != 0 {
            let natm = mole_c.atoms().len();
            for i in 0..natm {
                mole_c.atoms.atoms[i].stm = Some(stms[&i].clone());
            }
        }
        let nobt = cmat.ncols();
        let ato_atms = atos.atms("sph");
        let ato_shls = atos.shls("sph");
        let ato_syms = atos.syms("sph");
        let obt_engs = engs.to_vec();
        let obt_occs = [vec![1.0; nela], vec![0.0; nobt - nela]].concat();
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat);
        (eng, coefs)
    }
}
