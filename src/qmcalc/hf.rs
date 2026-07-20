use std::collections::HashMap;

use crate::base::Coefs;
use crate::datas::btrans::lmn2sym;
// 实现HF方法，根据结构和基组信息，计算得到波函数
use crate::integ::{Mat2e, RsInt};
use crate::{atomprop, integ, maths};
use crate::{base::Mole, utils};
use ndarray::{Array1, Array2, Array4, s};

/// 闭壳层Hartree-Fock计算器
pub struct RHF<'m> {
    /// 分子结构引用
    pub mole: &'m Mole,
    /// 是否使用局部坐标系
    pub use_lstm: bool,
    /// 积分形式("sph"或"car")
    pub int_form: String,
}

impl<'m> RHF<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        let use_lstm = false;
        let int_form = "sph".to_string();
        Self {
            mole,
            use_lstm,
            int_form,
        }
    }

    // 计算闭壳层密度矩阵
    fn calc_rdmat(
        &self,
        nele: usize,
        cmat: &Array2<f64>, // 系数矩阵
    ) -> Array2<f64> {
        // let mut dmat:ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>=Array::<f64, _>::zeros((nmat,nmat));
        let coef = cmat.slice(s![.., ..nele / 2]).to_owned();
        let dmat = coef.dot(&coef.t());
        return dmat * 2.0;
    }

    fn calc_rfork(
        // 计算fork矩阵
        &self,
        nato: usize,
        tmat: &Array2<f64>, // 动能矩阵
        vmat: &Array2<f64>, // 势能矩阵
        dmat: &Array2<f64>, // 密度矩阵
        emat: &Mat2e,       // 双电子积分矩阵
    ) -> Array2<f64> {
        // 计算Fork矩阵

        let mut jmat = Array2::zeros((nato, nato)); // 库伦积分
        let mut kmat = Array2::zeros((nato, nato)); // 交换积分
        let mut jsum: f64;
        let mut ksum: f64;
        for i in 0..nato {
            for j in 0..nato {
                jsum = 0.0;
                ksum = 0.0;
                for k in 0..nato {
                    for l in 0..nato {
                        jsum += emat.get(i, j, k, l) * dmat[(k, l)];
                        ksum += emat.get(i, k, j, l) * dmat[(k, l)];
                    }
                }
                jmat[(i, j)] = jsum;
                kmat[(i, j)] = ksum;
            }
        }
        return tmat + vmat + jmat - kmat * 0.5;
    }

    fn calc_ele_eng(
        &self,
        dmat: &Array2<f64>, // 密度矩阵
        tmat: &Array2<f64>, // 动能矩阵
        vmat: &Array2<f64>, // 势能矩阵
        fmat: &Array2<f64>, // fork矩阵
    ) -> f64 {
        (dmat * (tmat + vmat + fmat)).sum() * 0.5
    }

    fn calc_nuc_eng(&self) -> f64 {
        let natm = self.mole.atoms.len();
        let mut eng = 0.0;
        let nucs = self.mole.atoms.nucs();
        for i in 0..natm - 1 {
            for j in i + 1..natm {
                let dist = self.mole.atoms.dist(i, j);
                eng += (nucs[i] * nucs[j]) as f64 / dist;
            }
        }
        eng
    }

    /// 实现RHF自洽场计算
    ///
    /// 通过迭代求解Roothaan方程,收敛后返回电子能量和分子轨道系数。
    ///
    /// # 参数
    /// - `emat`: 可选的双电子积分,不提供时内部自动计算
    pub fn run(&self, emat: Option<Mat2e>) -> (f64, Coefs) {
        println!("开始计算");
        // let map_ato = self.mole.basis().map_ato();
        // let nato = map_ato.len("car");
        let ints = RsInt::new(self.mole.atoms(), self.mole.basis());
        let (mat1e, emat) = if self.use_lstm {
            println!("使用局部坐标系");
            let mut dir_caler = atomprop::direction::Calculator::new(self.mole);
            let mut stms = HashMap::new();
            let natm = self.mole.atoms.len();
            for atm in 0..natm {
                let stm = dir_caler.LCS(atm, None).unwrap();
                stms.insert(atm, stm);
            }

            let mat1e = ints.mat1e_lcs(&stms, [1, 1, 1]);
            println!("单电子积分计算完成");
            let mat2e = match emat {
                Some(emat) => emat,
                None => ints.mat2e_lcs(&stms),
            };
            println!("双电子积分计算完成");

            (mat1e, mat2e)
        } else {
            println!("使用全局坐标系");
            let mat1e = ints.mat1e([1, 1, 1], &self.int_form);
            println!("单电子积分计算完成");
            let mat2e = match emat {
                Some(emat) => emat,
                None => ints.mat2e(&self.int_form),
            };
            println!("双电子积分计算完成");
            (mat1e, mat2e)
        };
        let [smat, tmat, vmat] = mat1e;
        let nato = smat.shape()[0];
        println!("单电子积分计算完成");
        // println!("双电子积分计算完成");
        let xmat = maths::matrix_nhalf(&smat);
        let mut eng = f64::INFINITY;
        let mut engs = Array1::zeros(nato);
        let [nela, nelb] = self.mole.atoms.neles; // 分子总电子数
        println!("电子数:{},{}", nela, nelb);

        let mut cmat = {
            let hmat = &tmat + &vmat; // 核心哈密顿矩阵
            let fmat_ = xmat.dot(&hmat).dot(&xmat);
            let (_engs, cmat_) = maths::eigh(&fmat_);
            xmat.dot(&cmat_)
        }; // 系数矩阵初猜
        let mut dmat = self.calc_rdmat(nela + nelb, &cmat); // 计算密度矩阵
        let nuc_eng = self.calc_nuc_eng();
        // let mut cmat = xmat.clone(); // 开壳层有两个系数矩阵
        for i in 0..100 {
            dmat = 0.5 * dmat + 0.5 * self.calc_rdmat(nela + nelb, &cmat); // 计算密度矩阵
            let fmat = self.calc_rfork(nato, &tmat, &vmat, &dmat, &emat);
            let fmat_ = xmat.dot(&fmat).dot(&xmat);
            let (_engs, cmat_) = maths::eigh(&fmat_);
            engs = _engs;
            // println!("engs:{:?}", engs.to_vec());
            cmat = xmat.dot(&cmat_);
            check_nele(&dmat, &smat);
            let _eng = self.calc_ele_eng(&dmat, &tmat, &vmat, &fmat);
            println!("epoch:{:>2},电子能量:{:>14.8},分子能量:{}", i + 1, _eng, _eng + nuc_eng);
            if (_eng - eng).abs() < 1e-8 {
                break;
            } else {
                eng = _eng;
            }
        }
        let atos = &self.mole.basis().atos().sphs;
        let ato_atms = atos.iter().map(|ato| ato.atm).collect();
        let ato_shls = atos.iter().map(|ato| ato.shl).collect();
        let ato_syms = atos.iter().map(|ato| ato.sym.clone()).collect();
        let obt_engs = engs.to_vec();
        let nobt = cmat.ncols();
        let obt_occs = [vec![1.0; nela], vec![0.0; nobt - nela]].concat();
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat);
        // println!("coefs:\n{}", coefs);
        (eng, coefs)
    }
}

/// 开壳层(自旋非限制)Hartree-Fock计算器
pub struct UHF<'m> {
    /// 分子结构引用
    pub mole: &'m Mole,
}

impl<'m> UHF<'m> {
    /// 创建新的UHF计算器
    pub fn new(mole: &'m Mole) -> Self {
        Self { mole }
    }
    /// 开壳层密度矩阵计算,对α和β电子分别计算
    pub fn calc_udmat(
        &self,
        nele: [usize; 2],        // alpha 和 beta 电子数
        cmat: &[Array2<f64>; 2], // 两个系数矩阵
    ) -> [Array2<f64>; 2] {
        let coef_a = cmat[0].slice(s![.., ..nele[0]]).to_owned();
        let coef_b = cmat[1].slice(s![.., ..nele[1]]).to_owned();
        let dmat_a = coef_a.dot(&coef_a.t());
        let dmat_b = coef_b.dot(&coef_b.t());
        [dmat_a, dmat_b]
    }

    fn calc_ufork(
        // 计算fork矩阵
        &self,
        nmat: usize,
        tmat: &Array2<f64>,      // 动能矩阵
        vmat: &Array2<f64>,      // 势能矩阵
        dmat: &[Array2<f64>; 2], // 密度矩阵
        emat: &Mat2e,            // 双电子积分矩阵
    ) -> [Array2<f64>; 2] {
        // 计算Fork矩阵

        let fmat_a = Array2::zeros((nmat, nmat));
        let fmat_b = Array2::zeros((nmat, nmat));
        let mut fmat = [fmat_a, fmat_b];
        let hmat = tmat + vmat;
        let dmat_sum = &dmat[0] + &dmat[1];
        let mut jsum;
        let mut csum;
        for s in 0..2usize {
            for i in 0..nmat {
                for j in 0..nmat {
                    jsum = 0.0; // 库伦矩阵元
                    csum = 0.0; // 交换矩阵元
                    for k in 0..nmat {
                        for l in 0..nmat {
                            jsum += emat.get(i, j, k, l) * dmat_sum[(k, l)];
                            csum += emat.get(i, k, j, l) * dmat[s][(k, l)];
                        }
                    }
                    fmat[s][(i, j)] = hmat[(i, j)] + jsum - csum;
                }
            }
        }

        return fmat;
    }

    // 开壳层的能量
    fn calc_ele_eng(
        &self,
        dmat: &[Array2<f64>; 2], // 密度矩阵
        tmat: &Array2<f64>,      // 势能矩阵
        vmat: &Array2<f64>,      // 动能矩阵
        fmat: &[Array2<f64>; 2], // fork矩阵
    ) -> [f64; 2] {
        let dmata = &dmat[0];
        let dmatb = &dmat[1];
        let fmata = &fmat[0];
        let fmatb = &fmat[1];
        let hcore = &(tmat + vmat);
        let engas = (hcore + fmata) * dmata;
        let engbs = (hcore + fmatb) * dmatb;
        return [engas.sum() * 0.5, engbs.sum() * 0.5];
    }

    /// 计算开壳层HF的自洽场
    ///
    /// 迭代求解Pople-Nesbet方程,返回α和β电子的能量、轨道能量和系数矩阵。
    pub fn calc_uhf(&self) -> ([f64; 2], [Array1<f64>; 2], [Array2<f64>; 2]) {
        let atos = &self.mole.basis().atos().sphs;
        let nato = atos.len();
        let cmat_a = Array2::zeros((nato, nato));
        let cmat_b = Array2::zeros((nato, nato));
        let mut cmat = [cmat_a, cmat_b];

        let dmat_a = Array2::zeros((nato, nato));
        let dmat_b = Array2::zeros((nato, nato));
        let mut dmat = [dmat_a, dmat_b];

        let rs_int = integ::RsInt::new(self.mole.atoms(), self.mole.basis());
        let [smat, tmat, vmat] = rs_int.mat1e([1, 1, 1], "car");
        let emat = rs_int.mat2e("car");

        let xmat = maths::matrix_nhalf(&smat);
        let mut eng = [0.0; 2];
        let mut obt_engs = [Array1::zeros(nato), Array1::zeros(nato)];
        let mut beta = 0.5; // 新旧密度矩阵的混合系数
        let mut old_eng = 100.0;
        let [nela, nelb] = self.mole.atoms.neles;
        log::info!("UHF计算");
        for i in 0..50 {
            let new_dmat = self.calc_udmat([nela, nelb], &cmat); // 计算密度矩阵
            for c in 0..2 {
                dmat[c] = &dmat[c] * (1.0 - beta) + &new_dmat[c] * beta;
            }
            let fmat = self.calc_ufork(nato, &tmat, &vmat, &dmat, &emat);
            for c in 0..2usize {
                let slice = &fmat[c];
                let fmat_ = xmat.dot(slice).dot(&xmat);
                let (_engs, _cmat) = maths::eigh(&fmat_);
                obt_engs[c] = _engs;
                cmat[c] = xmat.dot(&_cmat);
            }

            eng = self.calc_ele_eng(&dmat, &tmat, &vmat, &fmat);

            if (eng[0] + eng[1] - old_eng).abs() < 1e-2 && beta != 1.0 {
                beta = 1.0;
            }
            if i > 10 {
                beta = 1.0;
            }

            log::info!(
                "epoch:{i:>2},电子能量:{:>10.4} + {:>10.4} = {:>14.8}",
                eng[0],
                eng[1],
                eng[0] + eng[1]
            );
            if (eng[0] + eng[1] - old_eng).abs() < 1e-6 {
                break;
            }
            old_eng = eng[0] + eng[1];
        }
        (eng, obt_engs, cmat)
    }
}

fn check_nele(dmat: &Array2<f64>, smat: &Array2<f64>) {
    let ds = dmat * smat;
    let nele = ds.sum();
    println!("当前电子数:{}", nele);
}
