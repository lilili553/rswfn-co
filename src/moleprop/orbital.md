#![allow(non_snake_case)]

use nalgebra::Vector3;
use ndarray::{Array2, s};
use ndarray_linalg::Inverse;
use std::collections::HashMap;
use std::fmt;
use std::hash::Hash;

use crate::atomprop;
use crate::base::Mole;
use crate::datas;
use crate::datas::mmocv;
use crate::moleprop;
use crate::{maths, utils};

pub struct Calculator<'m> {
pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
pub fn new(mole: &'m Mole) -> Self {
Calculator { mole }
}

    pub fn pocv(
        &self,
        dirs: &HashMap<usize, Vector3<f64>>, // 每个原子对应的方向向量
        keep_other_atm: bool,              // 其它原子轨道系数是否保留
        keep_other_sym: bool,              //是否保留其它类型的轨道
    ) -> Array2<f64> {
        let atoms = self.mole.atoms();
        let uls = atoms.ato_uls("raw", self.mole);
        let mut atms: Vec<usize> = dirs.keys().cloned().collect();
        atms.sort_by(|x, y| x.cmp(y));
        let cmat_raw = self.mole.coefs().get_cmat("raw");
        let [nato, nobt] = cmat_raw.shape().try_into().expect("获取矩阵形状失败");
        // 初始化投影的分子轨道
        let mut cmat_pocv: Array2<f64> = Array2::zeros([nato, nobt]);

        let ato_atms = self.mole.coefs().ato_atms("raw");
        let ato_syms = self.mole.coefs().ato_syms("raw");

        let obts = self.mole.coefs().occ_idxs(true); // 占据轨道的索引
        for obt in obts.iter() {
            for ato in 0..nato {
                let atm = &ato_atms[ato];
                let sym = &ato_syms[ato];
                if atms.contains(atm) && sym.contains("P") {
                    cmat_pocv[(ato, *obt as usize)] = cmat_raw[(ato, *obt as usize)]
                }
                if !atms.contains(atm) && keep_other_atm {
                    cmat_pocv[(ato, *obt as usize)] = cmat_raw[(ato, *obt as usize)]
                }
                if !sym.contains("P") && keep_other_sym {
                    cmat_pocv[(ato, *obt as usize)] = cmat_raw[(ato, *obt as usize)]
                }
            }
        }

        let mut pidxs = Vec::new(); // p轨道的索引
        let ato_syms = self.mole.coefs().ato_syms("raw");

        assert_eq!(uls.len(), atoms.len(), "原子再轨道内索引的数量");

        for atm in atms.iter() {
            let [u, l] = uls[*atm as usize];
            for ato in u..l {
                let ato_sym = &ato_syms[ato];
                if ato_sym.contains("P") {
                    pidxs.push([*atm as usize, ato]);
                }
            }
        }
        assert_eq!(pidxs.len() % 3, 0, "原子P轨道的数量必须是3的倍数");
        for i in (0..pidxs.len()).step_by(3) {
            let [atm, ato] = pidxs[i];
            let dir = dirs.get(&(atm as usize)).expect("dir获取atm失败");
            let dx = dir[0];
            let dy = dir[1];
            let dz = dir[2];
            for obt in obts.iter() {
                let px = &cmat_pocv[(ato + 0, *obt as usize)];
                let py = &cmat_pocv[(ato + 1, *obt as usize)];
                let pz = &cmat_pocv[(ato + 2, *obt as usize)];
                let [tx, ty, tz] = maths::projection([*px, *py, *pz], [dx, dy, dz]);
                cmat_pocv[(ato + 0, *obt as usize)] = tx;
                cmat_pocv[(ato + 1, *obt as usize)] = ty;
                cmat_pocv[(ato + 2, *obt as usize)] = tz;
            }
        }
        cmat_pocv
    }

    pub fn pi_pocv(&self, atms: Option<&Vec<usize>>) -> (HashMap<usize, Vector3<f64>>, Array2<f64>) {
        let dir_caler = atomprop::direction::Calculator::new(self.mole);
        let mut dirs = HashMap::new();
        let atoms = self.mole.atoms();

        let atms = match atms {
            Some(atms) => atms,
            None => &atoms.idxs(),
        };
        // dbg!("atms:{:?}", atms);
        for atm in atms {
            if atoms.get(*atm).sym == "H" {
                continue;
            }
            match dir_caler.normal_vector(*atm) {
                Some(dir) => {
                    dirs.insert(*atm, dir);
                }
                None => {}
            }
        }
        let cmat = self.pocv(&dirs, false, false);
        (dirs, cmat)
    }

    pub fn mocv(&self, mocvs: &HashMap<usize, mocv>) -> Array2<f64> {
        let ato_syms = self.mole.coefs().ato_syms("car");
        let ato_atms = self.mole.coefs().ato_atms("car");
        let cmat_car = self.mole.coefs().get_cmat("car");
        let [nrow, ncol] = cmat_car.shape().try_into().expect("获取矩阵形状失败");
        let mut cmat_mocv = Array2::zeros([nrow, ncol]);
        let mut sidx = Vec::new();
        let mut pidx = Vec::new();
        let mut didx = Vec::new();
        let mut fidx = Vec::new();
        // 获取每种类型轨道的索引
        for (i, sym) in ato_syms.iter().enumerate() {
            if datas::btrans::MIX_S_SYMS.contains(&sym.as_str()) {
                sidx.push(i);
                continue;
            }
            if datas::btrans::MIX_P_SYMS.contains(&sym.as_str()) {
                pidx.push(i);
                continue;
            }
            if datas::btrans::CAR_D_SYMS.contains(&sym.as_str()) {
                didx.push(i);
                continue;
            }
            if datas::btrans::CAR_F_SYMS.contains(&sym.as_str()) {
                fidx.push(i);
                continue;
            }
        }
        let atms: Vec<usize> = mocvs.keys().copied().collect();
        // 对每种类型分别处理
        for ato in sidx.iter() {
            // S轨道的索引
            let atm = ato_atms[*ato];
            if !atms.contains(&atm) {
                continue;
            }
            if mocvs.get(&atm).expect("mocvs获取atm失败").skeep[0] == 1 {
                // 如果保留S轨道
                let rcoefs: Array2<f64> = cmat_car.slice(s![*ato..*ato + 0, ..]).to_owned();
                cmat_mocv.slice_mut(s![*ato, ..]).assign(&rcoefs);
            }
        }
        for ato in pidx.iter().step_by(3) {
            // P轨道的索引
            let atm = ato_atms[*ato];
            if !atms.contains(&atm) {
                continue;
            }
            let pkeep: Vec<f64> = mocvs
                .get(&atm)
                .expect("mocvs获取atm失败")
                .pkeep
                .iter()
                .map(|x| *x as f64)
                .collect();
            let pkeep = Array2::from_shape_vec((3, 1), pkeep).expect("生成矩阵失败");
            let rcoefs: Array2<f64> = cmat_car.slice(s![*ato..*ato + 3, ..]).to_owned();
            let T = &mocvs.get(&atm).expect("mocvs获取atm失败").base;
            let tcoefs = mocv(T, &rcoefs, &pkeep);
            cmat_mocv.slice_mut(s![*ato..*ato + 3, ..]).assign(&tcoefs);
        }
        for ato in didx.iter().step_by(6) {
            // D轨道的索引
            let atm = ato_atms[*ato];
            if !atms.contains(&atm) {
                continue;
            }
            let dkeep: Vec<f64> = mocvs
                .get(&atm)
                .expect("mocvs获取atm失败")
                .dkeep
                .iter()
                .map(|x| *x as f64)
                .collect();
            let dkeep = Array2::from_shape_vec((6, 1), dkeep).expect("生成矩阵失败");
            let rcoefs: Array2<f64> = cmat_car.slice(s![*ato..*ato + 6, ..]).to_owned();
            let T = &mocvs.get(&atm).expect("mocvs获取atm失败").base;
            let tcoefs = mocv(T, &rcoefs, &dkeep);
            cmat_mocv.slice_mut(s![*ato..*ato + 6, ..]).assign(&tcoefs);
        }
        cmat_mocv
    }

    pub fn pi_mocv(&self) -> (HashMap<usize, mocv>, Array2<f64>) {
        let dir_caler = atomprop::direction::Calculator::new(self.mole);
        let mut mocvs = HashMap::new();
        let atoms = self.mole.atoms();
        for i in 0..atoms.len() {
            if atoms.get(i as usize).atm == 1 {
                continue;
            }
            // match dir_caler.local_coord_system(i as usize, None) {
            //     Some(bse) => {
            //         let mut mocv = mocv {
            //             base: bse,
            //             skeep: [0],
            //             pkeep: [0, 0, 1],
            //             dkeep: [0, 0, 0, 1, 1, 1],
            //         };
            //         if atoms.is_linear(i as usize) {
            //             mocv.pkeep = [0, 1, 1]
            //         }
            //         mocvs.insert(i as usize, mocv);
            //     }
            //     None => {}
            // }
            let stm = dir_caler.rand_lcs(i as usize);
            let mut mocv = mocv {
                base: stm,
                skeep: [0],
                pkeep: [0, 0, 1],
                dkeep: [0, 0, 0, 1, 1, 1],
            };
            if atoms.is_linear(i as usize) {
                mocv.pkeep = [0, 1, 1]
            }
            mocvs.insert(i as usize, mocv);
        }

        let cmat = self.mocv(&mocvs);
        (mocvs, cmat)
    }

    // 轨道定域化
    pub fn LMO(&self) {
        let mut cmat = self.mole.coefs().get_cmat("raw").to_owned();
        let smat = self.mole.coefs().get_smat("raw", self.mole).to_owned();
        let natm = self.mole.atoms().len();
        let uplo = self.mole.atoms().ato_uls("raw", self.mole);
        // let [nato, nobt] = cmat.shape().try_into().expect("获取矩阵形状失败");
        let obts = self.mole.coefs().occ_idxs(true);
        let nobt = obts.len();
        let mut pval_old = 0.0;
        for _cyl in 0..30 {
            // 最多循环30次
            for i in 0..nobt - 1 {
                let oi = obts[i] as usize;
                let cmi = cmat
                    .slice(s![.., oi])
                    .insert_axis(ndarray::Axis(1))
                    .to_owned(); // 第i个分子轨道系数
                for j in i + 1..nobt {
                    let oj = obts[j] as usize;
                    let cmj = cmat
                        .slice(s![.., oj])
                        .insert_axis(ndarray::Axis(1))
                        .to_owned(); // 第j个分子轨道系数
                    let dij: Array2<f64> = cmi.dot(&cmj.t()); // i和j的密度矩阵
                    let dii: Array2<f64> = cmi.dot(&cmi.t()); // i的密度矩阵
                    let djj: Array2<f64> = cmj.dot(&cmj.t()); // j的密度矩阵
                    let qij = dij * &smat; // 电子布局矩阵
                    let qii = dii * &smat;
                    let qjj = djj * &smat;

                    let mut Aval = 0.0;
                    let mut Bval = 0.0;
                    for a in 0..uplo.len() {
                        let [u, l] = uplo[a];
                        let Qij = qij.slice(s![u..l, u..l]).sum(); // 第a个原子的轨道i和j上的电子数量
                        let Qii = qii.slice(s![u..l, u..l]).sum(); // 第a个原子的轨道i上的电子数量
                        let Qjj = qjj.slice(s![u..l, u..l]).sum(); // 第a个原子的轨道j上的电子数量
                        Aval += Qii.powi(2) - Qjj.powi(2);
                        Bval += 2.0 * Qij * (Qii - Qjj);
                    }
                    if Aval.powi(2) + Bval.powi(2) < 1e-12 {
                        continue;
                    }
                    let gamma = 0.25 * f64::atan2(Bval, -Aval); // 旋转角度
                    let cmi = cmat.slice(s![.., oi]);
                    let cmj = cmat.slice(s![.., oj]);
                    let veci = f64::cos(gamma) * &cmi + f64::sin(gamma) * &cmj; // 旋转两个分子轨道
                    let vecj = f64::cos(gamma) * &cmj - f64::sin(gamma) * &cmi;
                    cmat.slice_mut(s![.., oi]).assign(&veci);
                    cmat.slice_mut(s![.., oj]).assign(&vecj);
                }
            }

            let mut pval_new = 0.0;
            for imo in obts.iter() {
                let imo = *imo as usize;
                let cmi = cmat
                    .slice(s![.., imo])
                    .insert_axis(ndarray::Axis(1))
                    .to_owned();
                let dmi: Array2<f64> = cmi.dot(&cmi.t());
                let qmi = &dmi * &smat;
                for iatm in 0..natm {
                    let [u, l] = uplo[iatm];
                    let Qii = qmi.slice(s![u..l, u..l]).sum();
                    pval_new += Qii.powi(2); // 一个分子轨道内所有原子的电子数的平方和
                }
            }
            if (pval_old - pval_new).abs() < 1e-5 {
                break;
            } else {
                pval_old = pval_new;
            }
        }
    }

}

pub fn mocv(T: &Array2<f64>, rcoefs: &Array2<f64>, keeps: &Array2<f64>) -> Array2<f64> {
let tmat = match rcoefs.shape()[0] {
1 => mmocv::tmat\*s(T),
3 => mmocv::tmat_p(T),
6 => mmocv::tmat_d(T),
10 => mmocv::tmat_f(T),
15 => mmocv::tmat_g(T),

- => panic!("Invalid orbital type"),
  };
  let Mr = tmat.inv().expect("矩阵求逆失败");
  let mut tcoefs = Mr.dot(rcoefs); //系数映射到局部坐标系
  tcoefs = &tcoefs \* keeps; // 保留某个维度
  let Mi = Mr.inv().expect("矩阵求逆失败");
  let fcoefs = Mi.dot(&tcoefs); // 再映射回全局坐标系
  fcoefs
  }

#[derive(Debug, Clone)]
pub struct mocv {
pub base: Array2<f64>,
pub skeep: [usize; 1],
pub pkeep: [usize; 3],
pub dkeep: [usize; 6],
}

impl fmt::Display for mocv {
fn fmt(&self, f: &mut fmt::Formatter<'\_>) -> fmt::Result {
let ds = ["dx", "dy", "dz"];
write!(f, "\n")?;
for j in 0..3 {
write!(f, "{}:", ds[j])?;
for i in 0..3 {
write!(f, "{:>10.4}", self.base[(i, j)])?;
}
write!(f, "\n")?;
}
writeln!(f, "s:{:>3}", self.skeep[0])?;

        write!(f, "p:")?;
        for i in 0..3 {
            write!(f, "{:>3}", self.pkeep[i])?;
        }
        write!(f, "\n")?;

        write!(f, "d:")?;
        for i in 0..6 {
            write!(f, "{:>3}", self.dkeep[i])?;
        }
        write!(f, "\n")?;

        Ok(())
    }

}
