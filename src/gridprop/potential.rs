use std::sync::OnceLock;

use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::base::Mole;
use crate::gridprop::{density, dftgrid};
use crate::maths::calc_dist;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub rhos: OnceLock<(Vec<[f64; 3]>, Vec<f64>, Vec<f64>)>,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole, rhos: OnceLock::new() }
    }

    pub fn nuc_potential(
        &self,
        grid: &[f64; 3], // 请求点（要计算的点）的坐标
    ) -> f64 {
        // 计算每个格点的核势能
        let natm = self.mole.atoms().len();
        let nucs = self.mole.atoms().nucs();
        let xyzs = self.mole.atoms().xyzs();
        let mut val = 0.0;
        for a in 0..natm {
            let dist = calc_dist(grid, &xyzs[a]); // 计算距离
            if dist < 1e-6 {
                continue;
            }
            val += nucs[a] as f64 / dist;
        }
        val
    }

    pub fn ele_potential(
        &self,
        grid: &[f64; 3], // query position 请求点（要计算的点）的坐标
    ) -> f64 {
        let (grids, weits, dens) = self.rhos.get_or_init(|| {
            let rho_caler = density::Calculator::new(self.mole);
            let xyz_caler = dftgrid::Calculator::new(self.mole);
            let (grids, weits) = xyz_caler.mol_grids(5, 74);
            let atms = self.mole.atoms().idxs();
            let dens: Vec<f64> = grids.par_iter().map(|grid| rho_caler.mol_rho_dm(grid, 0).0).collect();
            (grids, weits, dens)
        });

        // 计算每个格点的电子势能
        let ngrid = grids.len();
        let mut val = 0.0;
        for g in 0..ngrid {
            let mut dist = calc_dist(grid, &grids[g]);
            if dist < 1e-6 {
                dist = 1e-6;
            }
            val += weits[g] * dens[g] / dist;
        }
        val
    }
}
