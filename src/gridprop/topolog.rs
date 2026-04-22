#![allow(non_snake_case)]

/**
 * 拓扑分析，我理解就是寻找分子中的一些点
*/
use crate::base::Mole;
use crate::gridprop::density;
use ndarray::Array2;

use ndarray_linalg::Inverse;
use ndarray_linalg::{Determinant, EigValsh};
use rayon::prelude::*;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    // 提供一个初猜点，从该初猜点寻找临界点
    pub fn AIM(&self, init_grid: &[f64; 3]) -> ([f64; 3], f64, i32) {
        let rho_caler = density::Calculator::new(self.mole);
        let mut grid = init_grid.clone();
        let mut Hessian = Array2::zeros((3, 3));
        let mut rho = 0.0;
        for _i in 0..100 {
            // 最多迭代100次
            let (rho0, rho1, rho2) = rho_caler.mol_rho_dm(&grid, 2);
            rho = rho0;
            Hessian = Array2::from_shape_vec((3, 3), rho2.iter().cloned().flatten().collect()).expect("矩阵生成失败");
            // Hessian=Hessian.t().to_owned();
            let det = Hessian.det().expect("行列式计算失败");

            if det.abs() < 1e-10 {
                break;
            }
            let [dx, dy, dz] = rho1;
            if (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt() < 1e-6 {
                break;
            }
            let Hessinv = Hessian.inv().expect("逆矩阵计算失败");
            let gradent = Array2::from_shape_vec((3, 1), rho1.to_vec()).expect("矩阵生成失败");
            let delta = Hessinv.dot(&gradent).flatten().to_vec();
            let mx = delta[0];
            let my = delta[1];
            let mz = delta[2];
            grid[0] -= mx;
            grid[1] -= my;
            grid[2] -= mz;
        }
        let eig_vals = Hessian.eigvalsh(ndarray_linalg::UPLO::Lower).expect("求本征值失败");
        let mut n_pos = 0;
        let mut n_neg = 0;
        for val in eig_vals.iter() {
            if *val > 1e-8 {
                n_pos += 1;
            }
            if *val < -1e-8 {
                n_neg += 1;
            }
        }
        let cpt = n_pos - n_neg; // 临界点类型
        (grid, rho, cpt)
    }

    // 获取分子中的所有临界点
    pub fn AIM_mole(&self) -> Vec<([f64; 3], i32)> {
        let atoms = self.mole.atoms();
        let xyzs = atoms.xyzs();
        let mut grids = Vec::new();
        // 每个原子的中心
        for xyz in xyzs.iter() {
            grids.push(xyz.clone());
        }
        let natm = xyzs.len();
        for i in 0..natm - 1 {
            let [xi, yi, zi] = xyzs[i];
            for j in i + 1..natm {
                let [xj, yj, zj] = xyzs[j];
                if !atoms.is_bond(i, j) {
                    continue;
                }
                let cent = [(xi + xj) / 2.0, (yi + yj) / 2.0, (zi + zj) / 2.0];
                grids.push(cent);
            }
        }

        let ngrid = grids.len();
        let aims: Vec<([f64; 3], f64, i32)> = (0..ngrid).into_par_iter().map(|g| self.AIM(&grids[g])).collect();
        aims.iter().filter(|aim| aim.1 > 1e-3).map(|aim| (aim.0, aim.2)).collect()
    }
}
