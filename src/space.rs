#![allow(non_snake_case)]

use crate::maths::calc_dist;

/// 生成三维空间格点坐标
pub fn get_grids(nx: usize, ny: usize, nz: usize) -> Vec<[f64; 3]> {
    let mut grids = Vec::with_capacity(nx * ny * nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                grids.push([i as f64, j as f64, k as f64]);
            }
        }
    }
    grids
}

/// 计算原子到分子的权重分配（原子权重乘以Becke权重）
pub fn a2m_weits(
    iatm: usize,
    atm_grids: &Vec<[f64; 3]>, // 原子格点坐标
    atm_weits: &Vec<f64>,      // 原子格点权重
    atm_pos: &Vec<[f64; 3]>,   // 原子坐标
    atm_rad: &Vec<f64>,        // 原子半径
) -> Vec<f64> {
    let ngrid = atm_grids.len();
    let natm = atm_pos.len();
    let mut a2m_weits = Vec::with_capacity(ngrid);
    let mut su: Vec<Vec<f64>>;
    for g in 0..ngrid {
        let gp = atm_grids[g];
        su = vec![vec![1.0; natm]; natm];
        for i in 0..natm {
            let pi = atm_pos[i];
            let ri = calc_dist(&pi, &gp);
            for j in 0..natm {
                if j == i {
                    continue;
                };
                let pj = atm_pos[j];
                let rj = calc_dist(&pj, &gp);
                let dij = calc_dist(&pi, &pj);
                let miu_ij = (ri - rj) / dij;
                let chi = atm_rad[i] / atm_rad[j]; // 两原子半径之比
                let mut nu_ij;
                if (chi - 1.0).abs() < 1e-6 {
                    nu_ij = miu_ij;
                } else {
                    let u_ij = (chi - 1.0) / (chi + 1.0);
                    let mut a_ij = u_ij / (u_ij.powi(2) - 1.0);
                    if a_ij > 0.5 {
                        a_ij = 0.5
                    };
                    if a_ij < -0.5 {
                        a_ij = -0.5
                    };
                    nu_ij = miu_ij + a_ij * (1.0 - miu_ij.powi(2));
                }
                nu_ij = 1.5 * nu_ij - 0.5 * nu_ij.powi(3);
                nu_ij = 1.5 * nu_ij - 0.5 * nu_ij.powi(3);
                nu_ij = 1.5 * nu_ij - 0.5 * nu_ij.powi(3);

                su[j][i] = 0.5 * (1.0 - nu_ij);
            }
        }
        // 计算权重乘积
        let mut wt = vec![1.0; natm];
        for i in 0..natm {
            for j in 0..natm {
                wt[j] *= su[i][j];
            }
        }

        // 归一化并保存结果
        let sum_wt: f64 = wt.iter().sum();
        let rat = if sum_wt.abs() > 1e-12 { wt[iatm] / sum_wt } else { 0.0 };
        a2m_weits.push(atm_weits[g] * rat);
    }

    a2m_weits
}

/// 计算原子（片段）间的静电相互作用能
/// 包括电子-电子排斥、核-电子吸引和核-核排斥三项
pub fn calc_EIEBA(
    fragA: &Vec<[f64; 4]>,  // 片段A的原子(核电荷数及坐标)
    fragB: &Vec<[f64; 4]>,  // 片段B的原子(核电荷数及坐标)
    rhoA: &Vec<f64>,        // 片段A的电子密度
    rhoB: &Vec<f64>,        //片段B的电子密度
    gridsA: &Vec<[f64; 3]>, // 格点坐标
    gridsB: &Vec<[f64; 3]>, // 格点坐标
    weitsA: &Vec<f64>,      // 格点权重
    weitsB: &Vec<f64>,      // 格点权重
) -> (f64, f64, f64) {
    let mut eng_ee = 0.0; // 电子-电子排斥
    let mut eng_en = 0.0; // 核-电子吸引
    let mut eng_nn = 0.0;

    for i in 0..rhoA.len() {
        // 电子-电子排斥
        for j in 0..rhoB.len() {
            let dist_ee = calc_dist(&gridsA[i], &gridsB[j]); // 两个格点之间的距离
            if dist_ee < 1e-6 {
                continue;
            }
            eng_ee += rhoA[i] * weitsA[i] * rhoB[j] * weitsB[j] / dist_ee
        }
    }

    // A核-B电子 吸引
    for g in 0..rhoB.len() {
        for a in 0..fragA.len() {
            let [c, x, y, z] = fragA[a];
            let dist_ne = calc_dist(&[x, y, z], &gridsB[g]);
            if dist_ne < 1e-6 {
                continue;
            }
            eng_en -= rhoB[g] * weitsB[g] * c / dist_ne;
        }
    }

    // B核-A电子 吸引
    for g in 0..rhoA.len() {
        for a in 0..fragB.len() {
            let [c, x, y, z] = fragB[a];
            let dist_ne = calc_dist(&[x, y, z], &gridsA[g]);
            if dist_ne < 1e-6 {
                continue;
            }
            eng_en -= rhoA[g] * weitsA[g] * c / dist_ne;
        }
    }

    // A核-B核 排斥
    for i in 0..fragA.len() {
        let [ca, xa, ya, za] = fragA[i];
        for j in 0..fragB.len() {
            let [cb, xb, yb, zb] = fragB[j];
            let dist_ij = calc_dist(&[xa, ya, za], &[xb, yb, zb]);
            eng_nn += ca * cb / dist_ij;
        }
    }

    (eng_ee, eng_en, eng_nn)
}
