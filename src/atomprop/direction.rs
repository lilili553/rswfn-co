#![allow(non_snake_case)]

use core::f64;
use std::sync::Mutex;
use std::{collections::HashMap, f64::consts::PI, vec};

use ndarray::Array2;

use crate::base::{Mole, Stm, Vec3};
use crate::gridprop;
use crate::maths;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub rand_vz: Vec3,
    normals: Mutex<HashMap<usize, Option<Vec3>>>, // 记录每个原子的法向量
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        // let rand_vz = maths::rand_vect();
        let rand_vz = Vec3::rand();
        Calculator { mole: mole, rand_vz: rand_vz, normals: Mutex::new(HashMap::new()) }
    }

    pub fn normal_vector(&self, atm: usize) -> Option<Vec3> {
        let atoms = self.mole.atoms();
        let nebs = atoms.nebs(atm);
        let mut vz = None;

        match self.get_normal(atm) {
            Some(normal) => vz = Some(normal),
            // 如果没法直接计算的话则迭代搜索
            None => {
                if nebs.len() == 4 {
                    vz = None;
                } else if nebs.len() == 5 {
                    vz = None;
                } else {
                    let mut searchd = vec![atm]; // 已经搜索过的
                    let mut searchs = atoms.nebs(atm);
                    while searchs.len() > 0 {
                        let idx = searchs.pop().expect("弹出列表项失败");
                        searchd.push(idx);
                        let normal = self.get_normal(idx);
                        match normal {
                            Some(normal) => vz = Some(normal),
                            None => {
                                for neb in atoms.nebs(idx) {
                                    if searchd.contains(&neb) {
                                        continue;
                                    };
                                    if searchs.contains(&neb) {
                                        continue;
                                    };
                                    searchs.push(neb);
                                }
                            }
                        }
                    }
                }
            }
        };
        // 保证只有一个原子的时候法向量垂直于键轴
        vz = match vz {
            Some(vz) => {
                if nebs.len() == 1 {
                    let vx = atoms.vec(nebs[0], atm).normalize();
                    let vz = if vx.angle(&vz).abs() < 0.1 { Vec3::rand() } else { vz };
                    let vy = vz.cross(&vx).normalize();
                    Some(vx.cross(&vy).normalize())
                } else {
                    Some(vz)
                }
            }
            None => None,
        };
        // assert!()
        vz
    }

    pub fn get_normal(&self, atm: usize) -> Option<Vec3> {
        let mut normals = self.normals.lock().expect("锁定失败");
        normals
            .entry(atm)
            .or_insert_with(|| {
                let atoms = self.mole.atoms();
                let nebs = atoms.nebs(atm);
                if nebs.len() == 1 {
                    return None;
                } else if nebs.len() == 2 {
                    let ai = nebs[0];
                    let aj = nebs[1];
                    let vi = atoms.vec(atm, ai).normalize();
                    let vj = atoms.vec(atm, aj).normalize();
                    let angle = vi.angle(&vj);
                    if (angle - 1.0).abs() < 0.01 {
                        return None;
                    } else {
                        let mut normal = vi.cross(&vj).normalize();
                        // let [dx, dy, dz] = self.rand_vz.clone().into();
                        // let angle = maths::vec_angle(&normal.into(), &[dx, dy, dz]);
                        let angle = normal.angle(&self.rand_vz);
                        if angle.abs() > 0.5 {
                            normal = normal * -1.0;
                        }
                        return Some(normal);
                    }
                } else if nebs.len() == 3 {
                    let pi: Vec3 = atoms.get(nebs[0]).vec();
                    let pj: Vec3 = atoms.get(nebs[1]).vec();
                    let pk: Vec3 = atoms.get(nebs[2]).vec();
                    let vij: Vec3 = (&pj - &pi).normalize();
                    let vik: Vec3 = (&pk - &pi).normalize();
                    let mut normal = vij.cross(&vik).normalize();
                    let [A, B, C] = normal.into_arr();
                    let [cx, cy, cz] = ((pi + pj + pk) / 3.0).into_arr();
                    let D = -A * cx - B * cy - C * cz;
                    let [ax, ay, az] = atoms.get(atm).xyz;
                    let dist = f64::abs(A * ax + B * ay + C * az + D) / f64::sqrt(A.powi(2) + B.powi(2) + C.powi(2)); // 点到平面的距离
                    if dist < 0.1 {
                        // 是一个平面
                        let rand_vect = self.rand_vz.clone();
                        // let angle = maths::vec_angle(&normal.into(), &rand_vect.into());
                        let angle = normal.angle(&rand_vect);
                        if angle > 0.5 {
                            normal *= -1.0;
                        }
                        return Some(normal);
                    } else {
                        let vect = atoms.get(atm).vec() - Vec3::new(cx, cy, cz);
                        let angle = normal.angle(&vect) / PI;
                        if angle > 0.5 {
                            normal *= -1.0;
                        }
                        return Some(normal);
                    }
                } else if nebs.len() == 6 {
                    // 碳硼烷的特例，计算每一个键和其它原子的最小排斥方向
                    let mut min_val = f64::MIN;
                    let mut min_idx = nebs[0];
                    for i in 0..nebs.len() {
                        let mut val = f64::MAX;
                        for j in 0..nebs.len() {
                            if i == j {
                                continue;
                            }
                            let angle = atoms.ang(nebs[i], atm, nebs[j]);
                            if angle < val {
                                val = angle;
                            }
                        }
                        if val > min_val {
                            min_val = val;
                            min_idx = nebs[i];
                        }
                    }
                    let vect = atoms.vec(atm, min_idx);
                    Some(vect.normalize())
                } else {
                    return None;
                }
            })
            .clone()
    }

    pub fn dmat_fix(&self) -> Array2<f64> {
        let mut dirs = vec![];
        let atoms = self.mole.atoms();
        let natm = atoms.len();
        for atm in 0..natm {
            let normal = self.get_normal(atm as usize);
            dirs.push(normal);
        }
        let mut dmat_fix = Array2::<f64>::ones((natm, natm));
        for ai in 0..natm {
            for aj in 0..natm {
                if !atoms.is_bond(ai as usize, aj as usize) {
                    continue;
                }
                if atoms.nebs(ai as usize).len() < 5 && atoms.nebs(aj as usize).len() < 5 {
                    continue;
                }
                let c = match (dirs[ai].clone(), dirs[aj].clone()) {
                    (Some(di), Some(dj)) => {
                        let ang = di.angle(&dj);
                        ang.cos()
                    }
                    _ => 1.0,
                };
                dmat_fix[(ai, aj)] = c;
            }
        }
        dmat_fix
    }

    /// 计算原子所在的局部坐标系 local coord system
    pub fn LCS(&mut self, atm: usize, neb: Option<usize>) -> Option<Stm> {
        let atoms = self.mole.atoms();
        let atom = atoms.get(atm);
        if !atom.stm.non {
            return Some(atom.stm.clone());
        }
        let nebs = atoms.nebs(atm);

        let normal: Option<Vec3> = match self.normal_vector(atm) {
            Some(normal) => Some(normal),
            None => {
                let idx1;
                let idx2;
                if nebs.len() == 1 {
                    idx1 = nebs[0];
                    idx2 = atm;
                } else if nebs.len() == 2 {
                    idx1 = nebs[0];
                    idx2 = nebs[1];
                } else {
                    return None;
                }

                let vx = (atoms.get(idx1).vec() - atoms.get(idx2).vec()).normalize();
                let vr = Vec3::rand();
                Some(vx.cross(&vr).normalize())
            }
        };
        let stm = match normal {
            None => None,
            Some(normal) => {
                let vz = normal.normalize(); // 法向量
                let [vx, vy, vz] = match neb {
                    Some(neb) => {
                        // 如果有相邻原子
                        let vx = (atoms.get(neb).vec() - atom.vec()).normalize();
                        let vz = if vx.angle(&vz).abs() < 0.1 { Vec3::rand() } else { vz };
                        let vy = vz.cross(&vx).normalize();
                        let vz = vx.cross(&vy).normalize(); // 重新正交化
                        [vx, vy, vz]
                    }
                    None => {
                        let cent = atoms.cent(&nebs); // 相邻原子的中心

                        if nebs.len() == 1 {
                            let vx: Vec3 = atoms.vec(nebs[0], atm).normalize(); // 这个是确定的，键轴的方向
                            // 如果vx和vz在同一个方向，则重新整个vz
                            let vz = if vx.angle(&vz).abs() < 0.1 { Vec3::rand() } else { vz };
                            let vy = vz.cross(&vx).normalize();
                            let vz = vx.cross(&vy).normalize(); // 重新正交化
                            // log::info!("原子:{},vx:{},vy:{},vz:{}", atm + 1, vx, vy, vz);
                            [vx, vy, vz]
                        } else if nebs.len() == 2 {
                            let idx1 = nebs[0];
                            let idx2 = nebs[1];
                            let vx = atoms.vec(idx1, idx2).normalize();
                            let vy = vz.cross(&vx).normalize();
                            let vx = vy.cross(&vz).normalize(); // 重新正交化
                            [vx, vy, vz]
                        } else {
                            let vect = cent - atom.vec(); // 中心到原子的向量
                            let dist = vect.length();
                            let vr = if dist > 0.01 {
                                vect.normalize() // 还不一定垂直于z轴
                            } else {
                                Vec3::rand()
                            };
                            let vy = vz.cross(&vr).normalize();
                            let vx = vy.cross(&vz).normalize();
                            [vx, vy, vz]
                        }
                    }
                };
                let ex: [f64; 3] = vx.into_arr();
                let ey: [f64; 3] = vy.into_arr();
                let ez: [f64; 3] = vz.into_arr();
                let stm = Stm::new(ex, ey, ez);
                log::info!("atm:{}", atm);
                if !stm.check() {
                    log::error!("局部坐标系有误");
                }
                Some(stm)
            }
        };
        stm
    }

    /// 随机局部坐标系
    pub fn RLCS(&self, atm: usize) -> Stm {
        let vz = match self.normal_vector(atm) {
            Some(normal) => normal,
            None => {
                log::warn!("原子{}没有法向量", atm);
                Vec3::rand()
            }
        };
        let rand_vec = Vec3::rand();
        let vx = vz.cross(&rand_vec).normalize();
        let vy = vz.cross(&vx).normalize();
        let ex = vx.into_arr();
        let ey = vy.into_arr();
        let ez = vz.into_arr();
        Stm::new(ex, ey, ez)
    }

    /// 多重局部坐标系
    pub fn MLCS(&mut self, atm: usize) -> Vec<Stm> {
        let nebs = self.mole.atoms().nebs(atm);
        let mut mlcs = vec![];
        for idx in nebs.iter() {
            if let Some(stm) = self.LCS(*idx, Some(atm)) {
                mlcs.push(stm);
            }
        }
        mlcs
    }
    pub fn reactions(&self, atm: usize) -> Vec<Vec3> {
        let atoms = self.mole.atoms();
        let nebs = atoms.nebs(atm);
        let mut vectors = Vec::new();
        if nebs.len() == 1 {
            let vr: Vec3 = Vec3::rand(); // 随机向量
            let bond_vect: Vec3 = atoms.vec(atm, nebs[0]).normalize(); // 键向量
            let vert_vect: Vec3 = bond_vect.cross(&vr); // 垂直于键的向量
            let atom_cent = atoms.get(nebs[0]).vec(); // 原子的中心

            let point = &atom_cent + &vert_vect;

            let axis_i = bond_vect.cross(&vert_vect);
            let axis_j = bond_vect.clone();

            for i in 0..4 {
                let angle = i as f64 * PI / 8.0;
                let pi = point.rotate_point(&axis_i, angle);
                for j in 0..36 {
                    let angle = j as f64 * PI / 18.0;
                    let pj = pi.rotate_point(&axis_j, angle);
                    let vec = pj - &atom_cent;
                    vectors.push(vec);
                }
            }
            vectors.push(bond_vect.clone() * -1.0);
            return vectors;
            // return self.gt_bond_angle(atm);
        } else if nebs.len() == 2 {
            let ai = nebs[0];
            let aj = nebs[1];
            let vi = atoms.vec(atm, ai).normalize();
            let vj = atoms.vec(atm, aj).normalize();
            let angle = vi.angle(&vj);
            let cent = Vec3::zero();
            if (angle - PI).abs() < 1e-1 {
                // 线性分子
                let vr = Vec3::rand();
                let bond_vec = atoms.vec(ai, aj).normalize();
                let vert_vec = bond_vec.cross(&vr);

                for i in 0..36 {
                    // let vector = maths::rotate_point(&vert_vec, &cent, &bond_vec, i as f64 * PI / 18.0);
                    let vector = vert_vec.rotate_arrow(&cent, &bond_vec, i as f64 * PI / 18.0);
                    vectors.push(vector.normalize());
                }
            } else {
                let angle = PI - angle; // 有效角度
                let vert_vec = vi.cross(&vj).normalize();
                let mut count = (angle / (PI / 18.0)) as usize;
                if count < 3 {
                    count = 3;
                }
                let dang = angle / (count as f64);
                for i in 0..count + 1 {
                    // 水平旋转
                    let rang = dang * i as f64; // 旋转角度，跳过90度
                    let axis = vi.rotate_arrow(&cent, &vert_vec, -rang).normalize();
                    for i in 0..18 {
                        // 垂直旋转
                        // 垂直方向都是18个
                        // let vector = maths::rotate_point(&vert_vec, &cent, &axis, i as f64 * PI / 18.0);
                        let vector = vert_vec.rotate_arrow(&cent, &axis, i as f64 * PI / 18.0);
                        vectors.push(vector.normalize());
                    }
                }
            }
            return vectors;
        } else if nebs.len() == 3 {
            let atom = atoms.get(atm);
            let ai = nebs[0];
            let aj = nebs[1];
            let ak = nebs[2];
            let pi = atoms.get(ai).vec();
            let pj = atoms.get(aj).vec();
            let pk = atoms.get(ak).vec();

            let vi = atoms.vec(atm, ai);
            let vj = atoms.vec(atm, aj);
            let vk = atoms.vec(atm, ak);

            let plan_params = maths::get_plan_by_3points(&pi, &pj, &pk);
            let distance = atom.vec().dist_to_plan(&plan_params);
            let is_plane = distance.abs() < 1e-2;
            let [a, b, c, _d] = plan_params;
            let mut normal = Vec3::new(a, b, c).normalize();
            if is_plane {
                // 如果是平面
                vectors.push(normal.clone());
                vectors.push(normal * -1.0);
            } else {
                let mean_vec = (vi + vj + vk) / 3.0;
                let angle = mean_vec.angle(&normal);
                if angle < PI / 2.0 {
                    normal = normal * -1.0;
                }
                vectors.push(normal.clone());
                for neb in nebs.iter() {
                    let neb_vec = atoms.vec(atm, *neb).normalize(); // vx
                    let vy = normal.cross(&neb_vec).normalize(); // vy=vz@vx
                    let vz = neb_vec.cross(&vy).normalize();
                    vectors.push(vz);
                }
            }
            return vectors;
        } else {
            return vectors;
        }
    }

    pub fn gt_bond_angle(&self, atm: usize) -> Vec<Vec3> {
        // let dir_caler=atomprop::direction::Calculator::new(self.mole);
        let gridCaler = gridprop::dftgrid::Calculator::new(self.mole);
        let (grids, _weits) = gridCaler.atm_grids(atm, 1, 434);
        let atom = self.mole.atoms().get(atm);
        let nebs = self.mole.atoms().nebs(atm);
        let atm_xyz = atom.vec();
        let mut dirs = vec![];
        for g in 0..grids.len() {
            let grid = Vec3::from_arr(grids[g]);
            let mut meet = true; // 该点是否满足条件
            let vec_ag = (&grid - &atm_xyz).normalize();
            for neb in nebs.iter() {
                let neb_xyz = self.mole.atoms().get(*neb).vec();
                let vec_ab = (&neb_xyz - &atm_xyz).normalize();
                // let angle = maths::vec_angle(&vec_ab.into(), &vec_ag.into());
                let angle = vec_ab.angle(&vec_ag);
                if angle < 0.5 {
                    meet = false;
                }
            }
            if meet {
                dirs.push(vec_ag);
            }
        }
        dirs
    }

    pub fn linear_circle(&self, atm: usize) -> Vec<Vec3> {
        let dirs: Vec<Vec3> = Vec::new();
        let atoms = self.mole.atoms();
        let nebs = atoms.nebs(atm);
        if nebs.len() == 2 {
            let ai = nebs[0];
            let aj = nebs[1];
            let vi = atoms.vec(atm, ai);
            let vj = atoms.vec(atm, aj);
            let angle = vi.angle(&vj);
            if (angle - PI).abs() < 1e-3 {
                return self.reactions(atm);
            } else {
                return dirs;
            }
        } else {
            return dirs;
        }
    }
}
