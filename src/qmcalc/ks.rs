use core::f64;

use crate::base::Mole;
use crate::funcs::RsFun;
use crate::gridprop;
use crate::gridprop::wfnfunc::ato_wfn;
use crate::integ::{Mat2e, RsInt};
use crate::maths::{self, calc_dmat};
use crate::{funcs, utils};
use ndarray::Array2;

pub struct RKS<'m> {
    pub mole: &'m Mole,
    pub grids: Vec<[f64; 3]>,
    pub weits: Vec<f64>,
    pub wfn0s: Vec<Vec<f64>>, // 所有原子轨道在所有格点处的波函数
    pub wfn1s: Vec<Vec<[f64; 3]>>,
    pub wfn2s: Vec<Vec<[[f64; 3]; 3]>>,
}

impl<'m> RKS<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        let grid_caler = gridprop::dftgrid::Calculator::new(mole);
        let (grids, weits) = grid_caler.mol_grids(100, 74);
        let atos = &mole.basis().atos().cars;
        let atms = mole.atoms.idxs();
        let ngrid = grids.len();
        let atoms = &mole.atoms;
        let mut wfn0s = Vec::with_capacity(ngrid);
        let mut wfn1s = Vec::with_capacity(ngrid);
        let mut wfn2s = Vec::with_capacity(ngrid);
        for i in 0..ngrid {
            let grid = &grids[i];
            let (wfn0, wfn1, wfn2) = ato_wfn(atos, grid, atms, atoms, 2);
            wfn0s.push(wfn0);
            wfn1s.push(wfn1);
            wfn2s.push(wfn2);
        }
        Self { mole, grids, weits, wfn0s, wfn1s, wfn2s }
    }

    // 计算交换相关矩阵
    fn calc_xcmat(
        &self,
        vrhos: &Vec<f64>, // 能量对电子密度的梯度 de/dρ
        rho1s: &Vec<[f64; 3]>,
        vsigmas: &Vec<f64>, // 能量对sigma的梯度 de/dσ
    ) -> Array2<f64> {
        let nato = self.wfn0s[0].len(); // 矩阵大小/原子轨道数量
        let npos = self.weits.len(); // 格点数量
        let mut xcmat = Array2::zeros((nato, nato));
        for p in 0..npos {
            let weit = self.weits[p];
            let vrho = vrhos[p];
            let vsigma = vsigmas[p];
            let rho1 = rho1s[p];
            let wfn0 = &self.wfn0s[p];
            let wfn1 = &self.wfn1s[p];
            for i in 0..nato {
                let wfn0_i = wfn0[i];
                let wfn1_i = wfn1[i];
                for j in 0..nato {
                    let wfn0_j = wfn0[j];
                    let wfn1_j = wfn1[j];
                    let val_lda = wfn0_i * wfn0_j * vrho * weit;
                    let term =  // 中间项
                    (wfn0_i*wfn1_j[0]+wfn0_j*wfn1_i[0]) * rho1[0]+ // x方向
                    (wfn0_i*wfn1_j[1]+wfn0_j*wfn1_i[1]) * rho1[1]+ // y方向
                    (wfn0_i*wfn1_j[2]+wfn0_j*wfn1_i[2]) * rho1[2]; // z方向
                    let val_gga = 2.0 * vsigma * term * weit;

                    xcmat[(i, j)] += val_lda + val_gga;
                }
            }
        }
        return xcmat;
    }

    fn calc_jkmat(&self, dmat: &Array2<f64>, emat: &Mat2e) -> [Array2<f64>; 2] {
        let nato = dmat.shape()[0];
        let mut jmat = Array2::zeros((nato, nato)); // 库伦积分
        let mut kmat = Array2::zeros((nato, nato)); // 交换积分
        for i in 0..nato {
            for j in 0..nato {
                let mut jsum = 0.0;
                let mut ksum = 0.0;
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
        kmat *= -0.5;
        [jmat, kmat]
    }

    fn calc_ele_eng(&self, dmat: &Array2<f64>, hmat: &Array2<f64>, jmat: &Array2<f64>, kmat: &Array2<f64>, xcmat: &Array2<f64>, zk: &Vec<f64>, rhos: &Vec<f64>) -> f64 {
        let hcore_eng = (dmat * hmat).sum();
        let exact_x_eng = (dmat * (jmat + kmat * 0.2)).sum() * 0.5;
        let mut xc_eng = 0.0;
        for i in 0..zk.len() {
            xc_eng += zk[i] * rhos[i] * self.weits[i];
        }
        // xc_eng *= 0.8;
        println!("计算电子能量:{:>10.4},{:>10.4},{:>10.4}", hcore_eng, exact_x_eng, xc_eng);
        hcore_eng + exact_x_eng + xc_eng
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

    pub fn run(&self) {
        // 准备电子积分

        let int_caler = RsInt::new(self.mole.atoms(), self.mole.basis());
        let [smat, tmat, vmat] = int_caler.mat1e([1, 1, 1], "car");
        let hmat = &tmat + &vmat; //核心哈密顿矩阵
        println!("核心哈密顿矩阵:");
        utils::show_matrix(&hmat);
        let xmat = maths::matrix_nhalf(&smat);
        let emat = int_caler.mat2e("car");
        let nato = smat.shape()[0];

        let [nela, nelb] = self.mole.atoms.neles;
        let nobt = (nela + nelb) / 2;
        let occs: Vec<f64> = (0..nobt).map(|_| 2.0).collect();
        let atms = self.mole.atoms.idxs();
        let ato_atms = &self.mole.basis().atos().atms("car");
        let mut dmat = {
            let fmat_ = xmat.dot(&hmat).dot(&xmat);
            let (_engs, cmat_) = maths::eigh(&fmat_);
            let cmat = xmat.dot(&cmat_);
            calc_dmat(&cmat, &occs)
        }; // 密度矩阵初猜
        check_nele(&dmat, &smat);
        let ngrid = self.grids.len();
        let mut rho0s = vec![0.0; ngrid];
        let mut rho1s = vec![[0.0; 3]; ngrid];
        let mut rho2s = vec![[[0.0; 3]; 3]; ngrid];
        let mut sigmas = vec![0.0; ngrid];
        let mut ele_eng = f64::MAX;
        let rsfun = RsFun::new();
        let nuc_eng = self.calc_nuc_eng();
        let mut engs = vec![0.0; nato];
        for epoch in 0..30 {
            // 计算所有格点的电子密度
            let mut nele = 0.0;
            for g in 0..ngrid {
                let (rho0, rho1, rho2) = gridprop::density::mol_rho_dm(ato_atms, &dmat, &self.wfn0s[g], &self.wfn1s[g], &self.wfn2s[g], atms, 2);
                nele += rho0 * self.weits[g];
                rho0s[g] = rho0;
                rho1s[g] = rho1;
                rho2s[g] = rho2;
                let [dx, dy, dz] = &rho1s[g];
                sigmas[g] = dx.powi(2) + dy.powi(2) + dz.powi(2);
            }
            println!("电子密度电子数:{:>10.4}", nele); // 通过这里，电子密度计算应该没错
            let (zk, vrho, vsigma) = rsfun.run(&rho0s, &sigmas);
            let xcmat = self.calc_xcmat(&vrho, &rho1s, &vsigma); // 交换关联矩阵
            let [jmat, kmat] = self.calc_jkmat(&dmat, &emat); // fock 矩阵
            let fmat = &hmat + &jmat + &kmat * 0.20 + &xcmat;
            let fmat_ = xmat.dot(&fmat).dot(&xmat);
            let (_engs, cmat_) = maths::eigh(&fmat_);
            engs = _engs.flatten().to_vec();
            let cmat = xmat.dot(&cmat_);
            dmat = 0.5 * dmat + 0.5 * calc_dmat(&cmat, &occs); // 更新密度矩阵
            check_nele(&dmat, &smat);
            let eng = self.calc_ele_eng(&dmat, &hmat, &jmat, &kmat, &xcmat, &zk, &rho0s);
            let del_eng = eng - ele_eng;
            ele_eng = eng;
            println!("{:>2},电子能量:{:>10.4},分子能量:{:>10.4}", epoch, ele_eng, ele_eng + nuc_eng);
            if epoch > 0 && del_eng.abs() < 1e-5 {
                break;
            }
        }
        println!("轨道能量:{:?}", engs);
    }
}

fn check_nele(dmat: &Array2<f64>, smat: &Array2<f64>) {
    let ds = dmat * smat;
    let nele = ds.sum();
    println!("密度矩阵电子数:{}", nele);
}
