use core::fmt;
use ndarray::Array2;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::{Mutex, OnceLock};

use crate::base::{Mole, Vec3};
use crate::datas::consts::BOHR;
use crate::datas::elems;
use crate::{datas, maths, utils};

// 原子局部坐标系，局部坐标系的每一列是一个基向量，从全局坐标系到局部坐标系的变换方式为 R_l = LT @ R_g
#[derive(Debug, Clone)]
pub struct Stm {
    pub ex: [f64; 3],
    pub ey: [f64; 3],
    pub ez: [f64; 3],
    pub mat: Array2<f64>, // 局部坐标系的矩阵，每一列是一个单位基向量
    pub non: bool,
}

impl Stm {
    pub fn check(&self) -> bool {
        let mat_t = self.mat.t().to_owned();
        let mat_i = maths::inv(&self.mat);
        let same = utils::same_matrix(&mat_i, &mat_t);
        same
    }

    pub fn new(ex: [f64; 3], ey: [f64; 3], ez: [f64; 3]) -> Self {
        let [xx, xy, xz] = ex;
        let [yx, yy, yz] = ey;
        let [zx, zy, zz] = ez;
        let mat = Array2::from_shape_vec((3, 3), vec![xx, yx, zx, xy, yy, zy, xz, yz, zz]).unwrap();
        let non = false;
        Self { ex, ey, ez, mat, non }
    }

    pub fn blank() -> Self {
        let mut stm = Stm::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        stm.non = true;
        stm
    }

    pub fn vals(&self) -> Vec<f64> {
        let [xx, xy, xz] = self.ex;
        let [yx, yy, yz] = self.ey;
        let [zx, zy, zz] = self.ez;
        vec![xx, xy, xz, yx, yy, yz, zx, zy, zz]
    }

    pub fn from_mat(mat: Array2<f64>) -> Self {
        let mut es = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                es[i][j] = mat[(j, i)];
            }
        }
        let ex = es[0];
        let ey = es[1];
        let ez = es[2];
        Self {
            ex,
            ey,
            ez,
            mat,
            non: false,
        }
    }

    /// 局部坐标系B在局部坐标系A中的朝向
    pub fn relative(&self, stmb: &Stm) -> Stm {
        let mata = self.mat.clone();
        let matb = stmb.mat.clone();
        let dmat = mata.t().dot(&matb);
        Stm::from_mat(dmat)
    }

    pub fn dirs(&self) -> [[f64; 3]; 3] {
        [self.ex, self.ey, self.ez]
    }

    pub fn copy(&mut self, stm: &Stm) {
        self.ex = stm.ex;
        self.ey = stm.ey;
        self.ez = stm.ex;
        self.mat = stm.mat.clone();
        self.non = stm.non;
    }
}

impl fmt::Display for Stm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let [xx, xy, xz] = self.ex;
        let [yx, yy, yz] = self.ey;
        let [zx, zy, zz] = self.ez;
        write!(f, "{}  ", if self.non { "F" } else { "T" })?;
        write!(f, "ex: ({:>8.4}{:>8.4}{:>8.4})  |  ", xx, xy, xz)?;
        write!(f, "ey: ({:>8.4}{:>8.4}{:>8.4})  |  ", yx, yy, yz)?;
        write!(f, "ez: ({:>8.4}{:>8.4}{:>8.4})", zx, zy, zz)?;
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub idx: usize,    // 原子在分子中的索引
    pub nuc: usize,    // 原子核电荷数
    pub crg: f64,      // 有效核电荷数，赝势基组需要
    pub xyz: [f64; 3], // 原子坐标
    pub rad: f64,      // 原子半径
    pub sym: String,   // 原子符号
    pub stm: Option<Stm>,
}

impl Atom {
    pub fn vec(&self) -> Vec3 {
        let [x, y, z] = self.xyz;
        Vec3::new(x, y, z)
    }
}

impl Atom {
    pub fn new(nuc: usize, xyz: [f64; 3], idx: usize, crg: f64) -> Self {
        let mut rad = 0.0;
        let mut sym = None;
        for elem in datas::elems::ELEMS.iter() {
            if nuc == elem.idx {
                rad = elems::atm2rad(nuc);
                sym = Some(elem.sym.clone());
                break;
            }
        }
        Self {
            nuc: nuc,
            crg: crg,
            xyz: xyz,
            idx: idx,
            rad: rad,
            sym: sym.unwrap(),
            stm: None,
        }
    }

    pub fn irow(&self) -> usize {
        elems::atm2row(self.nuc)
    }
}

#[derive(Debug)]
pub struct Atoms {
    pub atoms: Vec<Atom>,
    pub dists: Array2<f64>,                  // 距离矩阵
    pub bonds: Array2<usize>,                // 键矩阵，两个原子之间是否成键 0 表示不成键，1 表示成键
    pub neles: [usize; 2],                   // α电子数和β电子数
    nebs: Mutex<HashMap<usize, Vec<usize>>>, // 每个原子的相邻原子
    angs: Mutex<HashMap<[usize; 3], f64>>,   // 键角
    vecs: Mutex<HashMap<[usize; 2], Vec3>>,
    is_lines: Mutex<HashMap<usize, bool>>,
    is_bonds: Mutex<HashMap<[usize; 2], bool>>,
    nucs: OnceLock<Vec<usize>>,
    crgs: OnceLock<Vec<f64>>,
    idxs: OnceLock<Vec<usize>>,
    xyzs: OnceLock<Vec<[f64; 3]>>,
    stms: OnceLock<Vec<Option<Stm>>>,
}

impl Atoms {
    pub fn new(atms: &Vec<usize>, xyzs: &Vec<[f64; 3]>, crgs: &Vec<f64>) -> Self {
        let natm = atms.len();
        let mut atoms = Vec::with_capacity(natm);
        let mut nele = 0;
        for i in 0..natm {
            let nuc = atms[i];
            let xyz = xyzs[i];
            let idx = i;
            let crg = crgs[i];
            let atom = Atom::new(nuc, xyz, idx, crg);
            atoms.push(atom);
            nele += nuc;
        }
        let nela = nele / 2 + nele % 2;
        let nelb = nele / 2;
        let neles = [nela, nelb];
        let mut dists = Array2::zeros((natm, natm));
        let mut bonds = Array2::zeros((natm, natm));
        for i in 0..natm - 1 {
            for j in i + 1..natm {
                let dist = maths::calc_dist(&xyzs[i], &xyzs[j]);
                dists[(i, j)] = dist;
                dists[(j, i)] = dist;
                if dist < (atoms[i].rad + atoms[j].rad) * 1.1 {
                    bonds[(i, j)] = 1;
                    bonds[(j, i)] = 1;
                }
            }
        }
        Atoms {
            atoms: atoms,
            dists: dists,
            bonds: bonds,
            neles: neles,
            nebs: Mutex::new(HashMap::new()),
            angs: Mutex::new(HashMap::new()),
            vecs: Mutex::new(HashMap::new()),
            is_lines: Mutex::new(HashMap::new()),
            is_bonds: Mutex::new(HashMap::new()),
            nucs: OnceLock::new(),
            crgs: OnceLock::new(),
            idxs: OnceLock::new(),
            xyzs: OnceLock::new(),
            stms: OnceLock::new(),
        }
    }

    pub fn add(&mut self, atm: usize, xyz: [f64; 3], crg: f64) {
        let idx = self.atoms.len() as usize;
        self.atoms.push(Atom::new(atm, xyz, idx, crg));
    }

    pub fn get(&self, idx: usize) -> &Atom {
        &self.atoms.get(idx as usize).unwrap_or_else(|| {
            log::error!("原子索引超出范围:{}/{}", idx, self.atoms.len());
            panic!("原子索引超出范围")
        })
    }

    pub fn len(&self) -> usize {
        return self.atoms.len();
    }

    /// 计算两个原子之间的距离
    pub fn dist(&self, i: usize, j: usize) -> f64 {
        let [xi, yi, zi] = self.atoms[i as usize].xyz;
        let [xj, yj, zj] = self.atoms[j as usize].xyz;
        ((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2)).sqrt()
    }

    pub fn dmat(&self) -> Array2<f64> {
        let natm = self.len();
        let mut dmat = Array2::zeros((natm, natm));
        for i in 0..natm - 1 {
            for j in i + 1..natm {
                dmat[(i, j)] = self.dist(i, j);
                dmat[(j, i)] = dmat[(i, j)];
            }
        }
        dmat
    }

    pub fn xyzs(&self) -> &Vec<[f64; 3]> {
        self.xyzs
            .get_or_init(|| self.atoms.iter().map(|atom| atom.xyz).collect())
    }

    pub fn xyzs_ang(&self) -> Vec<[f64; 3]> {
        self.atoms
            .iter()
            .map(|atom| {
                let [x, y, z] = atom.xyz;
                [x / BOHR, y / BOHR, z / BOHR]
            })
            .collect()
    }

    pub fn rads(&self) -> Vec<f64> {
        let mut rads = Vec::new();
        for atom in self.atoms.iter() {
            rads.push(atom.rad);
        }
        rads
    }

    pub fn syms(&self) -> Vec<String> {
        let mut syms = Vec::new();
        for atom in self.atoms.iter() {
            syms.push(atom.sym.clone());
        }
        syms
    }

    /// 每个原子的核电荷数
    pub fn nucs(&self) -> &Vec<usize> {
        let nucs = self
            .nucs
            .get_or_init(|| self.atoms.iter().map(|atom| atom.nuc).collect::<Vec<usize>>());
        nucs
    }

    // 每个原子的有效核电荷数
    pub fn crgs(&self) -> &Vec<f64> {
        let crgs = self
            .crgs
            .get_or_init(|| self.atoms.iter().map(|atom| atom.crg).collect::<Vec<f64>>());
        crgs
    }

    pub fn stms(&self) -> &Vec<Option<Stm>> {
        self.stms
            .get_or_init(|| self.atoms.iter().map(|atom| atom.stm.clone()).collect())
    }

    /// 原子索引范围 (0..natm-1)
    pub fn idxs(&self) -> &Vec<usize> {
        let idxs = self
            .idxs
            .get_or_init(|| self.atoms.iter().map(|atom| atom.idx).collect::<Vec<usize>>());
        idxs
    }

    // 获取指定原子的邻接原子
    pub fn nebs(&self, atm: usize) -> Vec<usize> {
        let mut nebs = self.nebs.lock().expect("锁定失败");
        let nebs = nebs
            .entry(atm)
            .or_insert_with(|| {
                let mut nebs = Vec::new();
                for i in 0..self.len() {
                    if i == atm as usize {
                        continue;
                    }
                    if self.is_bond(atm, i) {
                        nebs.push(i);
                    }
                }
                nebs
            })
            .clone();
        // log::info!("获取邻接原子{}:{:?}", atm, nebs);
        nebs
    }

    /// ai->aj 的向量
    pub fn vec(&self, ai: usize, aj: usize) -> Vec3 {
        let mut vects = self.vecs.lock().expect("锁定失败");
        vects
            .entry([ai, aj])
            .or_insert_with(|| {
                let [xi, yi, zi] = self.get(ai).xyz;
                let [xj, yj, zj] = self.get(aj).xyz;
                Vec3::new(xj - xi, yj - yi, zj - zi)
            })
            .clone()
    }

    // ai-aj-ak 三个原子之间的键角
    pub fn ang(&self, ai: usize, aj: usize, ak: usize) -> f64 {
        let mut angles = self.angs.lock().expect("锁定失败");
        *angles.entry([ai, aj, ak]).or_insert_with(|| {
            let veci = self.vec(aj, ai);
            let vecj = self.vec(aj, ak);
            // maths::vec_angle(&veci, &vecj)
            veci.angle(&vecj)
        })
    }

    // 该原子是否为线性的
    pub fn is_line(&self, atm: usize) -> bool {
        let mut is_lines = self.is_lines.lock().expect("锁定失败");
        *is_lines.entry(atm).or_insert_with(|| {
            let nebs = self.nebs(atm);
            if nebs.len() == 1 {
                // // 单独的O原子直接算线性
                // if self.syms()[atm] == "O" {
                //     return true;
                // }
                let neb = nebs[0];
                if self.nebs(neb).len() == 1 {
                    return true;
                } else if self.nebs(neb).len() == 2 {
                    let nebs = self.nebs(neb);
                    let nebs = nebs.iter().filter(|idx| *idx != &atm).copied().collect::<Vec<usize>>();
                    let idx = nebs[0];
                    let angle = self.ang(atm, neb, idx);
                    return (angle - PI).abs() < PI / 18.0;
                } else {
                    return false;
                }
            }
            if nebs.len() == 2 {
                let ai = nebs[0];
                let aj = nebs[1];
                let vi = self.vec(atm, ai);
                let vj = self.vec(atm, aj);
                let angle = vi.angle(&vj);
                return (angle - PI).abs() < PI / 18.0;
            }
            return false;
        })
    }

    /// 判断两个原子之间是否成键
    pub fn is_bond(&self, ai: usize, aj: usize) -> bool {
        self.bonds[(ai, aj)] == 1
    }

    pub fn set_bonds(&self, bonds: HashMap<[usize; 2], bool>) {
        println!("设置成键:{:?}", bonds);
        let mut is_bonds = self.is_bonds.lock().expect("锁定失败");
        *is_bonds = bonds;
        let mut is_lines = self.is_lines.lock().expect("锁定失败");
        is_lines.clear();
        let mut nebs = self.nebs.lock().expect("锁定失败");
        nebs.clear();
    }

    // 获取指定原子的中心坐标
    pub fn cent(&self, atms: &Vec<usize>) -> Vec3 {
        let count = atms.len();
        let mut cent = [0.0, 0.0, 0.0];
        for atm in atms.iter() {
            let [px, py, pz] = self.get(*atm).xyz;
            cent[0] += px;
            cent[1] += py;
            cent[2] += pz;
        }
        cent[0] /= count as f64;
        cent[1] /= count as f64;
        cent[2] /= count as f64;
        Vec3::from_arr(cent)
    }

    // 判断两个结构是否相等
    pub fn equal(&self, other: &Atoms) -> bool {
        let atoms_0 = &self.atoms;
        let atoms_1 = &other.atoms;
        if atoms_0.len() != atoms_1.len() {
            log::info!("原子数量不相等");
            return false;
        }
        for i in 0..atoms_0.len() {
            let atom_0 = atoms_0.get(i as usize).unwrap();
            let atom_1 = atoms_1.get(i as usize).unwrap();
            if atom_0.sym != atom_1.sym {
                log::info!("元素符号不相等");
                return false;
            }
            let [x0, y0, z0] = atom_0.xyz;
            let [x1, y1, z1] = atom_1.xyz;
            let dist = ((x1 - x0).powi(2) + (y1 - y0).powi(2) + (z1 - z0).powi(2)).sqrt();
            if dist > 1e-2 {
                log::info!("原子坐标不相等");
                return false;
            }
        }
        return true;
    }
}

impl Clone for Atoms {
    fn clone(&self) -> Self {
        Atoms {
            atoms: self.atoms.clone(),
            dists: self.dists.clone(),
            bonds: self.bonds.clone(),
            neles: self.neles.clone(),
            angs: Mutex::new(HashMap::new()),
            nebs: Mutex::new(HashMap::new()),
            vecs: Mutex::new(HashMap::new()),
            is_lines: Mutex::new(HashMap::new()),
            is_bonds: Mutex::new(HashMap::new()),
            idxs: OnceLock::new(),
            crgs: OnceLock::new(),
            nucs: OnceLock::new(),
            xyzs: OnceLock::new(),
            stms: OnceLock::new(),
        }
    }
}

impl fmt::Display for Atoms {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let natm = self.atoms.len();
        writeln!(f, "{:>3}{:>14}{:>14}{:>14}{:>10}", "atm", "x", "y", "z", "crg")?;
        for i in 0..natm {
            let atom = self.get(i);
            let [x, y, z] = atom.xyz;
            writeln!(f, "{:>3}{:>14.8}{:>14.8}{:>14.8}{:>10.4}", atom.sym, x, y, z, atom.crg)?;
        }
        for i in 0..natm {
            let atom = self.get(i);
            if !atom.stm.is_none() {
                writeln!(f, "{}", atom.stm.as_ref().unwrap())?;
            }
        }
        Ok(())
    }
}

unsafe impl Send for Atoms {}
unsafe impl Sync for Atoms {}
