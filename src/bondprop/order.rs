use ndarray::{Array2, Axis, s};
use std::collections::HashMap;

use crate::base::{Mole, Stm, Vec3};
use crate::integ::RsInt;
use crate::maths::calc_dmat;
// use crate::orbtprop::obtmat::Mocv;
use crate::orbtprop::obtmat::MocvDest;
use crate::{atomprop, maths, orbtprop};

/// 键级计算器，支持 Mayer 和 Wiberg 键级分析
pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    /// 计算 Mayer 键级矩阵
    pub fn mayer(&self) -> Array2<f64> {
        self.order("mayer")
    }

    /// 计算 Wiberg 键级矩阵（Löwdin 正交化基组下）
    pub fn wiberg(&self) -> Array2<f64> {
        self.order("wiberg")
    }

    /// 计算指定类型的键级矩阵
    ///
    /// # 参数
    /// * `otype` - 键级类型，"mayer" 或 "wiberg"
    pub fn order(&self, otype: &str) -> Array2<f64> {
        let smat = self.mole.basis().get_smat("raw", self.mole.atoms());
        let dmat = self.mole.coefs().get_dmat("raw", self.mole.basis());
        assert_eq!(smat.shape(), dmat.shape());
        let uplo = self.mole.basis().atm_uplo("raw");
        let omat_ato = match otype {
            "mayer" => mayer(dmat, smat),
            "wiberg" => wiberg(dmat, smat),
            _ => panic!("未知键级类型:{}", otype),
        };
        ato2atm(&omat_ato, &uplo)
    }

    /// 计算某个分子轨道对键级的贡献
    ///
    /// 对于占据轨道，删去该轨道后键级减少表示正贡献；
    /// 对于空轨道，添加该轨道后键级增加表示正贡献。
    pub fn order_obt(&self, otype: &str, iobt: usize) -> Array2<f64> {
        let smat = self.mole.basis().get_smat("raw", self.mole.atoms());
        let old_cmat = self.mole.coefs().get_cmat("raw", self.mole.basis());
        let old_dmat = self.mole.coefs().get_dmat("raw", self.mole.basis());
        let nato = old_cmat.nrows();
        let nobt = old_cmat.ncols();
        let obt_occs = self.mole.coefs().obt_occs.clone();
        let new_dmat = if obt_occs[iobt] != 0.0 {
            let mut new_cmat = old_cmat.clone();
            new_cmat.slice_mut(s![.., iobt]).fill(0.0);
            maths::calc_dmat(&new_cmat, &obt_occs)
        } else {
            let mut obt_occs = obt_occs.clone();
            obt_occs[iobt] = if self.mole.is_open() { 1.0 } else { 2.0 };
            maths::calc_dmat(&old_cmat, &obt_occs)
        };
        let old_omat = order(&old_dmat, smat, otype);
        let new_omat = order(&new_dmat, smat, otype);
        let omat_ato = if obt_occs[iobt] != 0.0 {
            old_omat - new_omat
        } else {
            new_omat - old_omat
        };
        let uplo = &self.mole.basis().atm_uplo("raw");
        ato2atm(&omat_ato, uplo)
    }

    /// 使用 POCV（投影轨道系数向量）方法计算键级
    pub fn pocv(
        &self,
        dirs: &HashMap<usize, Vec3>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        otype: &str,
    ) -> Result<Array2<f64>, String> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_pocv = obt_caler.pocv(dirs, keep_other_atm, keep_other_sym)?;

        let mut mole_pocv = self.mole.clone();
        let basis = mole_pocv.basis.as_ref().unwrap();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, basis);
        let caler = Calculator::new(&mole_pocv);
        Ok(caler.order(otype))
    }

    /// 使用 MOCV（分子轨道系数向量）方法计算键级
    pub fn mocv(&self, stms: &HashMap<usize, Stm>, atos: &Vec<usize>, otype: &str) -> Array2<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_mocv = obt_caler.mocv(stms, atos);
        let mut mole_mocv = self.mole.clone();
        let basis = mole_mocv.basis.as_ref().unwrap();
        mole_mocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_mocv, basis);
        let caler = Calculator::new(&mole_mocv);
        caler.order(otype)
    }

    /// 使用 POCV 方法计算 π 键级
    pub fn pi_pocv(&self, otype: &str) -> Result<(HashMap<usize, Vec3>, Array2<f64>), String> {
        let mut obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        obt_caler.atms = self.atms;
        let (dirs, cmat_pocv) = obt_caler.pi_pocv()?;

        let mut mole_pocv = self.mole.clone();
        let basis = mole_pocv.basis.as_ref().unwrap();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, basis);

        let caler = Calculator::new(&mole_pocv);
        let mut omat = caler.order(otype);
        omat = omat.abs().sqrt();
        Ok((dirs, omat))
    }

    /// 使用 MOCV 方法计算 π 键级
    pub fn pi_mocv(&self, otype: &str) -> (HashMap<usize, Stm>, Vec<usize>, Array2<f64>) {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let dest = MocvDest::PiProperty;
        let (stms, atos, cmat_mocv) = obt_caler.pi_mocv(&dest);
        let mut mole_mocv = self.mole.clone();
        let basis = mole_mocv.basis.as_ref().unwrap();
        mole_mocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_mocv, basis);
        let caler = Calculator::new(&mole_mocv);
        let mut omat = caler.order(otype);
        omat = omat.abs().sqrt(); // 取绝对值再开根号
        (stms, atos, omat)
    }

    /// 键级分解 获得σ,pi和δ键级
    pub fn decom(&self, atms: [usize; 2]) -> Vec<f64> {
        let keep_ato_syms_plan = vec![vec!["S", "PZ", "D 0"], vec!["PX", "XZ", "D-1"], vec!["XY"]];
        let keep_ato_syms_line = vec![
            vec!["S", "PZ", "D 0"],
            vec!["PX", "PY", "XZ", "YZ", "D+2", "D-2"],
            vec!["XY", "D+2", "D-2"],
        ];
        let dest = MocvDest::OrderDecom;
        let cmat_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let mut dir_caler = atomprop::direction::Calculator::new(self.mole);
        let natm = self.mole.atoms.len();
        let mut stms = HashMap::new();
        for i in 0..natm {
            let stm = if i == atms[0] || i == atms[1] {
                if i == atms[0] {
                    dir_caler.LCS(i, Some(atms[1])).unwrap()
                } else {
                    dir_caler.LCS(i, Some(atms[0])).unwrap()
                }
            } else {
                dir_caler.LCS(i, None).unwrap()
            };
            stms.insert(i, stm);
        }
        let int_caler = RsInt::new(self.mole.atoms(), self.mole.basis());
        let [smat, _, _] = int_caler.mat1e_lcs(&stms, [1, 0, 0]);
        let smat = int_caler.mat1e_raw(&smat);
        let ato_syms = self.mole.basis().atos().syms("raw");
        let nato = ato_syms.len();
        let obt_occs = &self.mole.coefs().obt_occs;
        let cmat_mocv = cmat_caler.mocv(&stms, &(0..nato).collect());
        let dmat_mocv = calc_dmat(&cmat_mocv, obt_occs);
        let uplo = self.mole.basis().atm_uplo("raw");
        let mut vals = vec![];
        let keep_ato_syms = if self.mole.atoms().is_line(atms[0]) {
            &keep_ato_syms_line
        } else {
            &keep_ato_syms_plan
        };
        for keep_syms in keep_ato_syms {
            let mut atos = vec![];
            for (idx, sym) in ato_syms.iter().enumerate() {
                if keep_syms.contains(&sym.as_str()) {
                    atos.push(idx);
                }
            }
            log::info!("atos:{:?}", atos);
            let omat = mayer(&dmat_mocv, &smat);
            vals.push(omat[(atms[0], atms[1])]);
        }
        vals
    }

    /// 计算束缚键级
    pub fn bound(&self, atm: usize, dir: &Vec3, nebs: &Vec<usize>, otype: &str) -> Result<Vec<f64>, String> {
        let mut info = format!("计算束缚键级:{} {} {}\n", atm + 1, dir, otype);
        let mut dirs = HashMap::new();
        dirs.insert(atm, dir.clone());
        let omat = self.pocv(&dirs, true, true, otype)?; // 获得键级矩阵
        let mut orders = Vec::new();
        for idx in nebs.iter() {
            let order = omat[(atm as usize, *idx as usize)];
            info.push_str(&format!("order.bound:{}->{}:{:>10.4}\n", atm + 1, idx + 1, order));
            orders.push(order);
        }
        log::info!("{}", info);
        Ok(orders)
    }
}

/// 根据指定类型计算键级矩阵
pub fn order(dmat: &Array2<f64>, smat: &Array2<f64>, otype: &str) -> Array2<f64> {
    match otype {
        "mayer" => mayer(dmat, smat),
        "wiberg" => wiberg(dmat, smat),
        _ => mayer(dmat, smat),
    }
}

/// Mayer 键级分析
pub fn mayer(dmat: &Array2<f64>, smat: &Array2<f64>) -> Array2<f64> {
    let psmat = dmat.dot(smat);
    &psmat * &psmat.t()
}

/// Wiberg bond order analysis in Löwdin orthogonalized basis
pub fn wiberg(dmat: &Array2<f64>, smat: &Array2<f64>) -> Array2<f64> {
    let smat_sqrt = maths::matrix_phalf(smat);
    let pmat_lowd = smat_sqrt.dot(dmat).dot(&smat_sqrt);
    &pmat_lowd * &pmat_lowd
}

/// 将原子轨道（AO）键级矩阵约化为原子（ATM）键级矩阵
pub fn ato2atm(omat_ato: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Array2<f64> {
    let natm = uplo.len();
    let mut omat_atm = Array2::zeros((natm, natm));
    for i in 0..natm {
        let [ui, li] = uplo[i];
        for j in 0..natm {
            let [uj, lj] = uplo[j];
            let order = omat_ato.slice(s![ui..li, uj..lj]).sum();
            omat_atm[(i, j)] = order;
        }
    }
    omat_atm
}
