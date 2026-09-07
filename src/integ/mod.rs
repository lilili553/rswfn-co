#![allow(non_snake_case)]
// 此模块定义积分函数
// use ndarray::Array2;
use crate::base::{Ato, Atos, Mole, Stm};
use crate::datas::FACS;
use crate::datas::btrans::{ang2lmns, sym2lmn};
use crate::datas::integ::{RHM, WHM};
use crate::datas::rys_root;
use crate::datas::{btrans, mocv};
use crate::maths;
use ndarray::{Array2, Array4, Axis, s};
use std::collections::HashMap;
use std::f64::consts::PI;

mod cint;

/// 高斯型轨道(GTO)结构体
pub struct Gto {
    /// 角动量量子数 [l, m, n]
    lmn: [usize; 3],
    /// 指数参数
    alp: f64,
    /// 轨道中心坐标
    xyz: [f64; 3],
}

impl Gto {
    /// 创建新的GTO
    pub fn new(lmn: [usize; 3], alp: f64, xyz: [f64; 3]) -> Self {
        Self { lmn, alp, xyz }
    }

    /// 高斯乘积定理
    ///
    /// 两个高斯型轨道的乘积仍为高斯型轨道,返回乘积系数和新的轨道。
    pub fn prod(&self, gto_j: &Gto) -> (f64, Gto) {
        let gto_i = self;
        let [xi, yi, zi] = gto_i.xyz;
        let [xj, yj, zj] = gto_j.xyz;
        let ai = gto_i.alp;
        let aj = gto_j.alp;
        let am = ai + aj;
        let xm = (ai * xi + aj * xj) / am;
        let ym = (ai * yi + aj * yj) / am;
        let zm = (ai * zi + aj * zj) / am;
        let xyz = [xm, ym, zm];
        let [li, mi, ni] = gto_i.lmn;
        let [lj, mj, nj] = gto_j.lmn;
        let lmn = [li + lj, mi + mj, ni + nj];
        let gto_m = Gto::new(lmn, am, xyz);
        let r2 = (xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2);
        let K = f64::exp(-ai * aj / am * r2);
        (K, gto_m)
    }
}

struct IntVal1e {
    Ga: Gto,
    Gb: Gto,
    p: f64,
    u: f64,
    Vab: [f64; 3],
    Vpa: [f64; 3],
    Vpb: [f64; 3],
    sval_d: HashMap<(isize, isize, usize), f64>,
    tval_d: HashMap<(isize, isize, usize), f64>,
    eval_d: HashMap<(isize, isize, isize, usize), f64>,
}

impl IntVal1e {
    fn new(Ga: Gto, Gb: Gto) -> Self {
        let p = Ga.alp + Gb.alp;
        let u = Ga.alp * Gb.alp / p;
        let mut Vab = [0.0; 3];
        let mut Vpa = [0.0; 3];
        let mut Vpb = [0.0; 3];
        for d in 0..3 {
            Vab[d] = Ga.xyz[d] - Gb.xyz[d];
            Vpa[d] = -Gb.alp / p * Vab[d];
            Vpb[d] = Ga.alp / p * Vab[d];
        }

        Self {
            Ga,
            Gb,
            p,
            u,
            Vab,
            Vpa,
            Vpb,
            sval_d: HashMap::new(),
            tval_d: HashMap::new(),
            eval_d: HashMap::new(),
        }
    }

    pub fn sval_d(&mut self, i: isize, j: isize, d: usize) -> f64 {
        let key = (i, j, d);
        if self.sval_d.contains_key(&key) {
            return self.sval_d[&key];
        }
        if i < 0 || j < 0 {
            return 0.0;
        }
        let Xab = self.Vab[d];
        let Xpa = self.Vpa[d];
        let Xpb = self.Vpb[d];

        let fi = i as f64;
        let fj = j as f64;
        if i == 0 && j == 0 {
            let sval = f64::sqrt(PI / self.p) * f64::exp(-self.u * Xab.powi(2)); // 9.3.10
            self.sval_d.insert(key, sval);
            return sval;
        }
        if i > 0 {
            // 9.3.8
            let sval = Xpa * self.sval_d(i - 1, j, d)
                + ((fi - 1.0) * self.sval_d(i - 2, j, d) + fj * self.sval_d(i - 1, j - 1, d)) / (2.0 * self.p);
            self.sval_d.insert(key, sval);
            return sval;
        } else {
            // 9.3.9
            let sval = Xpb * self.sval_d(i, j - 1, d)
                + (fi * self.sval_d(i - 1, j - 1, d) + (fj - 1.0) * self.sval_d(i, j - 2, d)) / (2.0 * self.p);
            self.sval_d.insert(key, sval);
            return sval;
        }
    }

    pub fn tval_d(&mut self, i: isize, j: isize, d: usize) -> f64 {
        let key = (i, j, d);
        if self.tval_d.contains_key(&key) {
            return self.tval_d[&key];
        }
        if i < 0 || j < 0 {
            return 0.0;
        }
        let Xab = self.Vab[d];
        let Xpa = self.Vpa[d];
        let Xpb = self.Vpb[d];
        let a = self.Ga.alp;
        let b = self.Gb.alp;
        if i == 0 && j == 0 {
            // 9.3.42
            let S00 = self.sval_d(0, 0, d);
            // let term = a - 2.0 * a.powi(2) * (Xab.powi(2) + 2.0 / self.p);
            let term = a - 2.0 * a.powi(2) * (Xpa.powi(2) + 1.0 / (2.0 * self.p));
            let val = term * S00;
            self.tval_d.insert(key, val);
            return val;
        }
        if i > 0 {
            // 9.3.40
            let i = i - 1;
            let tzz = self.tval_d(i, j, d);
            let tnz = self.tval_d(i - 1, j, d);
            let tzn = self.tval_d(i, j - 1, d);
            let spz = self.sval_d(i + 1, j, d);
            let snz = self.sval_d(i - 1, j, d);
            let val = Xpa * tzz
                + (i as f64 * tnz + j as f64 * tzn) / (2.0 * self.p)
                + (b / self.p) * (2.0 * a * spz - i as f64 * snz);
            self.tval_d.insert(key, val);
            return val;
        } else {
            // 9.3.41
            let j = j - 1;
            let tzz = self.tval_d(i, j, d);
            let tnz = self.tval_d(i - 1, j, d);
            let tzn = self.tval_d(i, j - 1, d);
            let szp = self.sval_d(i, j + 1, d);
            let szn = self.sval_d(i, j - 1, d);
            let val = Xpb * tzz
                + (i as f64 * tnz + j as f64 * tzn) / (2.0 * self.p)
                + (a / self.p) * (2.0 * b * szp - j as f64 * szn);
            self.tval_d.insert(key, val);
            return val;
        }
    }
    pub fn sval(&mut self) -> f64 {
        let mut sval = 1.0;
        for d in 0..3 {
            let i = self.Ga.lmn[d] as isize;
            let j = self.Gb.lmn[d] as isize;
            sval *= self.sval_d(i, j, d);
        }
        sval
    }

    pub fn tval(&mut self) -> f64 {
        let i = self.Ga.lmn[0] as isize;
        let k = self.Ga.lmn[1] as isize;
        let m = self.Ga.lmn[2] as isize;
        let j = self.Gb.lmn[0] as isize;
        let l = self.Gb.lmn[1] as isize;
        let n = self.Gb.lmn[2] as isize;
        let V0 = self.tval_d(i, j, 0) * self.sval_d(k, l, 1) * self.sval_d(m, n, 2);
        let V1 = self.sval_d(i, j, 0) * self.tval_d(k, l, 1) * self.sval_d(m, n, 2);
        let V2 = self.sval_d(i, j, 0) * self.sval_d(k, l, 1) * self.tval_d(m, n, 2);

        let val = V0 + V1 + V2;
        val
    }

    pub fn eval_d(&mut self, i: isize, j: isize, t: isize, d: usize) -> f64 {
        let key = (i, j, t, d);
        if self.eval_d.contains_key(&key) {
            return self.eval_d[&key];
        }
        let Xab = self.Vab[d];
        let Xpa = self.Vpa[d];
        let Xpb = self.Vpb[d];
        let Kab = f64::exp(-self.u * Xab.powi(2)); // 9.2.15
        // 边界
        if i < 0 || j < 0 {
            return 0.0;
        }
        // 9.5.5
        if t > i + j || t < 0 {
            return 0.0;
        }
        // 9.5.8
        if i == 0 && j == 0 {
            return if t == 0 { Kab } else { 0.0 };
        }
        // McMurchie-Davidson递归关系
        if t == 0 {
            // s轨道
            if i > j {
                // 9.5.15
                let val = Xpa * self.eval_d(i - 1, j, 0, d) + self.eval_d(i - 1, j, 1, d);
                self.eval_d.insert(key, val);
                return val;
            } else {
                // 9.5.16
                let val = Xpb * self.eval_d(i, j - 1, 0, d) + self.eval_d(i, j - 1, 1, d);
                self.eval_d.insert(key, val);
                return val;
            }
        }
        if t > 0 {
            // 9.5.17
            let val = (i as f64 * self.eval_d(i - 1, j, t - 1, d) + j as f64 * self.eval_d(i, j - 1, t - 1, d))
                / (2.0 * self.p * t as f64);
            self.eval_d.insert(key, val);
            return val;
        }
        panic!("不符条件");
    }

    // 公式 9.11.36
    pub fn hval_d(&self, t: isize, d: usize, xc: f64, s2: f64) -> f64 {
        let xa = self.Ga.xyz[d];
        let xb = self.Gb.xyz[d];
        let a = self.Ga.alp;
        let b = self.Gb.alp;
        let xp = (a * xa + b * xb) / self.p; // 9.2.13 电荷中心坐标
        let xpc = xp - xc;
        if t < 0 {
            return 0.0;
        }
        if t == 0 {
            return 1.0;
        }
        if t == 1 {
            return -2.0 * self.p * s2 * xpc;
        }
        let term = xpc * self.hval_d(t - 1, d, xc, s2) + (t - 1) as f64 * self.hval_d(t - 2, d, xc, s2);
        -2.0 * self.p * s2 * term
    }

    // 公式 9.11.43
    pub fn ival_d(&mut self, i: isize, j: isize, d: usize, xc: f64, s2: f64) -> f64 {
        let mut sum = 0.0;
        for t in 0..=(i + j) {
            let val = self.eval_d(i, j, t, d) * self.hval_d(t, d, xc, s2);
            sum += val
        }
        sum
    }

    // 核引积分，单个原子核 公式9.11.42
    pub fn vval_c(&mut self, rc: &[f64; 3]) -> f64 {
        let [xa, ya, za] = self.Ga.xyz;
        let [xb, yb, zb] = self.Gb.xyz;
        let a = self.Ga.alp;
        let b = self.Gb.alp;
        let xp = (a * xa + b * xb) / (a + b);
        let yp = (a * ya + b * yb) / (a + b);
        let zp = (a * za + b * zb) / (a + b);
        let rp = [xp, yp, zp]; // 重叠中心坐标
        let Rpc = maths::calc_dist(&rp, rc); // 重叠中心到原子核距离
        let alpha = self.p * Rpc.powi(2); // 9.11.34
        let mut gamma = 0;
        for d in 0..3 {
            let i = self.Ga.lmn[d] as isize;
            let j = self.Gb.lmn[d] as isize;
            gamma += i + j;
        }
        let gamma = (gamma as f64 / 2.0).floor() as usize + 1;
        let mut sum = 0.0;
        let roots = rys_root::rys_roots(gamma, alpha);
        for kappa in 0..gamma {
            let [s, w] = roots[kappa];
            let mut val = w;
            assert!(w > 0.0, "权重需为正值");
            for d in 0..3 {
                let i = self.Ga.lmn[d] as isize;
                let j = self.Gb.lmn[d] as isize;
                let xc = rc[d];
                val *= self.ival_d(i, j, d, xc, s.powi(2));
            }
            sum += val;
        }
        let vval = 2.0 * PI / self.p * sum; // 9.11.42
        vval
    }

    pub fn vval(&mut self, xyzs: &Vec<[f64; 3]>, crgs: &Vec<f64>) -> f64 {
        let mut vval = 0.0;
        let natm = crgs.len();
        for i in 0..natm {
            let xyz = xyzs[i];
            let crg = crgs[i];
            let val = -crg * self.vval_c(&xyz);
            vval += val;
        }
        vval
    }
}

struct IntVal2e<'g> {
    Ga: &'g Gto,
    Gb: &'g Gto,
    Gc: &'g Gto,
    Gd: &'g Gto,
    Gp: Gto,
    Gq: Gto,
    Kp: f64,
    Kq: f64,
    hval_d: HashMap<(isize, usize, u64), f64>,
    eval_d: HashMap<(isize, isize, isize, usize, String), f64>,
}

impl<'g> IntVal2e<'g> {
    pub fn new(Ga: &'g Gto, Gb: &'g Gto, Gc: &'g Gto, Gd: &'g Gto) -> Self {
        let (Kp, Gp) = Ga.prod(&Gb);
        let (Kq, Gq) = Gc.prod(&Gd);
        Self {
            Ga,
            Gb,
            Gc,
            Gd,
            Gp,
            Gq,
            Kp,
            Kq,
            hval_d: HashMap::new(),
            eval_d: HashMap::new(),
        }
    }

    // 公式 9.11.39
    pub fn hval_d(&mut self, t: isize, d: usize, s2: f64) -> f64 {
        let bits = s2.to_bits();
        let key = (t, d, bits);
        if self.hval_d.contains_key(&key) {
            return self.hval_d[&key];
        }
        let p = self.Gp.alp; // 9.9.28
        let q = self.Gq.alp; // 9.9.30
        let alp = p * q / (p + q);
        // let alp = 1.0 / (p + q);

        let Xpq = self.Gp.xyz[d] - self.Gq.xyz[d]; // 9.2.13 两个新的基函数的中心距离
        if t < 0 {
            return 0.0;
        }
        if t == 0 {
            return 1.0;
        }
        if t == 1 {
            let val = -2.0 * alp * s2 * Xpq;
            self.hval_d.insert(key, val);
            return val;
        }
        // 9.11.39
        let term = Xpq * self.hval_d(t - 1, d, s2) + (t - 1) as f64 * self.hval_d(t - 2, d, s2);
        let val = -2.0 * alp * s2 * term;
        self.hval_d.insert(key, val);
        return val;
    }

    pub fn eval_d(&mut self, i: isize, j: isize, t: isize, d: usize, c: &str) -> f64 {
        let key = (i, j, t, d, c.to_string());
        if self.eval_d.contains_key(&key) {
            return self.eval_d[&key];
        }
        let (Ga, Gb, Gp) = match c {
            "ab" => (&self.Ga, &self.Gb, &self.Gp),
            "cd" => (&self.Gc, &self.Gd, &self.Gq),
            _ => panic!("c参数错误"),
        };
        let Xab = Ga.xyz[d] - Gb.xyz[d]; // 9.2.14
        let p = Gp.alp;
        let Xpa = Gp.xyz[d] - Ga.xyz[d];
        let Xpb = Gp.xyz[d] - Gb.xyz[d];
        let u = Ga.alp * Gb.alp / p;
        let Kab = f64::exp(-u * Xab.powi(2)); // 9.2.15
        // 边界
        if i < 0 || j < 0 {
            return 0.0;
        }
        // 9.5.5
        if t > i + j || t < 0 {
            return 0.0;
        }
        // 9.5.8
        if i == 0 && j == 0 && t == 0 {
            return if t == 0 { Kab } else { 0.0 };
        }
        // McMurchie-Davidson递归关系
        if t == 0 {
            // s轨道
            if i > j {
                // 9.5.15
                let val = Xpa * self.eval_d(i - 1, j, 0, d, c) + self.eval_d(i - 1, j, 1, d, c);
                self.eval_d.insert(key, val);
                return val;
            } else {
                // 9.5.16
                let val = Xpb * self.eval_d(i, j - 1, 0, d, c) + self.eval_d(i, j - 1, 1, d, c);
                self.eval_d.insert(key, val);
                return val;
            }
        }
        if t > 0 {
            // 9.5.17
            let val = (i as f64 * self.eval_d(i - 1, j, t - 1, d, c) + j as f64 * self.eval_d(i, j - 1, t - 1, d, c))
                / (2.0 * p * t as f64);
            self.eval_d.insert(key, val);
            return val;
        }
        panic!("不符条件");
    }

    // 公式 9.11.49 - 9.11.51
    pub fn ival_d(&mut self, i: isize, j: isize, k: isize, l: isize, d: usize, s2: f64) -> f64 {
        let mut sum = 0.0;
        for t in 0..=(i + j) {
            for u in 0..=(k + l) {
                let sign = if u % 2 == 0 { 1.0 } else { -1.0 };
                // 这里u是公式9.11.49中的tau
                sum += sign * self.eval_d(i, j, t, d, "ab") * self.eval_d(k, l, u, d, "cd") * self.hval_d(t + u, d, s2);
            }
        }
        sum
    }

    // 双电子积分 公式 9.11.48
    pub fn gval(&mut self) -> f64 {
        let p = self.Ga.alp + self.Gb.alp; // 9.9.28
        let q = self.Gc.alp + self.Gd.alp; // 9.9.30
        let mu = p * q / (p + q);
        let Rpq = maths::calc_dist(&self.Gp.xyz, &self.Gq.xyz);
        let alpha = mu * Rpq.powi(2);

        let mut L = 0;
        for d in 0..3 {
            let la = self.Ga.lmn[d] as isize;
            let lb = self.Gb.lmn[d] as isize;
            let lc = self.Gc.lmn[d] as isize;
            let ld = self.Gd.lmn[d] as isize;
            L += la + lb + lc + ld;
        }
        let gamma = (L as f64 / 2.0).floor() as usize + 1;
        let roots = rys_root::rys_roots(gamma, alpha); //rys多项式的根和权重
        assert_eq!(roots.len(), gamma);
        let mut sum = 0.0;
        for kappa in 0..gamma {
            let [s, w] = roots[kappa];
            let s2 = s.powi(2);
            let mut mul = w;
            for d in 0..3 {
                let i = self.Ga.lmn[d] as isize;
                let j = self.Gb.lmn[d] as isize;
                let k = self.Gc.lmn[d] as isize;
                let l = self.Gd.lmn[d] as isize;
                mul *= self.ival_d(i, j, k, l, d, s2);
            }
            sum += mul
        }
        let k = 2.0 * PI.powf(2.5) / (p * q * f64::sqrt(p + q));
        k * sum
    }
}

/// 双电子积分压缩存储结构体
///
/// 利用对称性将四维双电子积分矩阵压缩为一维存储。
pub struct Mat2e {
    /// 原子轨道数量
    pub nato: usize,
    /// 有效积分数量
    pub nval: usize,
    /// 标记位,记录积分是否已计算
    pub idxs: Vec<bool>,
    /// 积分值列表
    pub vals: Vec<f64>,
}

impl Mat2e {
    /// 创建新的Mat2e
    pub fn new(nato: usize) -> Self {
        let nval = nato * (nato + 1) * (nato * (nato + 1) + 2) / 8; // 有效电子积分数量
        let idxs = vec![false; nval];
        let vals = vec![0.0; nval];
        Self { nato, nval, idxs, vals }
    }

    /// 获取积分值
    pub fn get(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        let idx = get_2eid(i, j, k, l);
        self.vals[idx]
    }

    /// 设置积分值
    pub fn set(&mut self, i: usize, j: usize, k: usize, l: usize, val: f64) {
        let idx = get_2eid(i, j, k, l);
        self.vals[idx] = val;
    }

    /// 检查是否所有积分都已计算完成
    pub fn check(&self) {
        for idx in self.idxs.iter() {
            if *idx == false {
                panic!("没有遍历完成");
            }
        }
    }
    /// 从四维数组构建Mat2e
    pub fn from_mat4(mat4: Array4<f64>) -> Self {
        let nato = mat4.shape()[0];
        let mut mat2e = Mat2e::new(nato);
        for i in 0..nato {
            for j in i..nato {
                for k in 0..nato {
                    for l in k..nato {
                        mat2e.set(i, j, k, l, mat4[[i, j, k, l]]);
                    }
                }
            }
        }
        mat2e
    }

    /// 转为四维数组
    pub fn into_mat4(&self) -> Array4<f64> {
        let nato = self.nato;
        let mut mat4 = Array4::zeros((nato, nato, nato, nato));
        for i in 0..nato {
            for j in i..nato {
                for k in 0..nato {
                    for l in k..nato {
                        mat4[[i, j, k, l]] = self.get(i, j, k, l);
                    }
                }
            }
        }
        mat4
    }
}

/// Rust实现的分子积分计算器
///
/// 提供单电子积分(重叠、动能、核引)和双电子积分的计算功能,
/// 以及球谐与笛卡尔基函数的转换。
pub struct RsInt<'m> {
    /// 分子信息引用
    mole: &'m Mole,
}

impl<'m> RsInt<'m> {
    /// 创建新的RsInt实例
    ///
    /// # Parameters
    /// - `mole`: 分子信息，提供原子结构与基组用于积分计算
    pub fn new(mole: &'m Mole) -> Self {
        Self { mole }
    }

    /// 计算重叠矩阵
    pub fn smat(&self, form: &str) -> Array2<f64> {
        self.mat1e([1, 0, 0], form)[0].clone()
    }

    /// 计算动能矩阵
    pub fn tmat(&self, form: &str) -> Array2<f64> {
        self.mat1e([0, 1, 0], form)[1].clone()
    }

    /// 计算核吸引势矩阵
    pub fn vmat(&self, form: &str) -> Array2<f64> {
        self.mat1e([0, 0, 1], form)[2].clone()
    }

    /// 计算单电子积分(使用libcint)
    #[cfg(feature = "use-cint")]
    pub fn mat1e(&self, cals: [usize; 3], form: &str) -> [Array2<f64>; 3] {
        let cint = cint::Cint::new(self.mole);
        cint.calc_int1e(form, cals)
    }

    /// 计算单电子积分(Rust实现)
    #[cfg(not(feature = "use-cint"))]
    pub fn mat1e(&self, cals: [usize; 3], form: &str) -> [Array2<f64>; 3] {
        let atos = &self.mole.basis().atos().cars;
        let nato = atos.len(); // 原子轨道数量
        let xyzs = self.mole.atoms().xyzs();
        let crgs = self.mole.atoms().crgs();
        let mut smat = Array2::zeros((nato, nato)); // 核引矩阵
        let mut tmat = Array2::zeros((nato, nato)); // 核引矩阵
        let mut vmat = Array2::zeros((nato, nato)); // 核引矩阵

        for i in 0..nato {
            for j in 0..nato {
                let ni = atos[i].alps.len();
                let nj = atos[j].alps.len();

                let coes_i = atos[i].coes_norm();
                let coes_j = atos[j].coes_norm();

                let gi_lmn = atos[i].lmn();
                let gj_lmn = atos[j].lmn();
                let gi_xyz = xyzs[atos[i].atm];
                let gj_xyz = xyzs[atos[j].atm];
                for a in 0..ni {
                    for b in 0..nj {
                        let Ga = Gto::new(*gi_lmn, atos[i].alps[a], gi_xyz);
                        let Gb = Gto::new(*gj_lmn, atos[j].alps[b], gj_xyz);
                        let mut ints = IntVal1e::new(Ga, Gb);
                        let coef = coes_i[a] * coes_j[b];
                        if cals[0] == 1 {
                            let sval = ints.sval();
                            smat[(i, j)] += coef * sval;
                        }
                        if cals[1] == 1 {
                            let tval = ints.tval();
                            tmat[(i, j)] += coef * tval;
                        }
                        if cals[2] == 1 {
                            let vval = ints.vval(xyzs, crgs);
                            vmat[(i, j)] += coef * vval;
                        }
                    }
                }
            }
        }
        if form == "sph" {
            if cals[0] == 1 {
                smat = self.mat1e_sph(&smat);
            }
            if cals[1] == 1 {
                tmat = self.mat1e_sph(&tmat);
            }
            if cals[2] == 1 {
                vmat = self.mat1e_sph(&vmat);
            }
        }
        if form == "raw" {
            if cals[0] == 1 {
                smat = self.mat1e_raw(&smat);
            }
            if cals[1] == 1 {
                tmat = self.mat1e_raw(&tmat);
            }
            if cals[2] == 1 {
                vmat = self.mat1e_raw(&vmat);
            }
        }
        [smat, tmat, vmat]
    }

    /// 计算双电子积分(使用libcint)
    #[cfg(feature = "use-cint")]
    pub fn mat2e(&self, form: &str) -> Mat2e {
        println!("使用libcint");
        let cint = cint::Cint::new(self.mole);
        let mat4 = cint.calc_int2e(form);
        return Mat2e::from_mat4(mat4);
    }
    /// 计算双电子库伦积分(Rust实现)
    #[cfg(not(feature = "use-cint"))]
    pub fn mat2e(&self, form: &str) -> Mat2e {
        let atos = &self.mole.basis().atos().cars;
        let nato = atos.len(); // 原子轨道数量
        let xyzs = self.mole.atoms().xyzs();
        let mut emat = Mat2e::new(nato);
        // let count = nato * nato * nato * nato;
        // let mut index = 0;
        for i in 0..nato {
            for j in 0..nato {
                for k in 0..nato {
                    for l in 0..nato {
                        let idx = get_2eid(i, j, k, l);
                        if emat.idxs[idx] == false {
                            let ni = atos[i].alps.len();
                            let nj = atos[j].alps.len();
                            let nk = atos[k].alps.len();
                            let nl = atos[l].alps.len();

                            let gi_lmn = atos[i].lmn();
                            let gj_lmn = atos[j].lmn();
                            let gk_lmn = atos[k].lmn();
                            let gl_lmn = atos[l].lmn();

                            let gi_xyz = xyzs[atos[i].atm];
                            let gj_xyz = xyzs[atos[j].atm];
                            let gk_xyz = xyzs[atos[k].atm];
                            let gl_xyz = xyzs[atos[l].atm];

                            let coes_i = atos[i].coes_norm();
                            let coes_j = atos[j].coes_norm();
                            let coes_k = atos[k].coes_norm();
                            let coes_l = atos[l].coes_norm();

                            let mut sum = 0.0;
                            for a in 0..ni {
                                for b in 0..nj {
                                    for c in 0..nk {
                                        for d in 0..nl {
                                            let Ga = Gto::new(*gi_lmn, atos[i].alps[a], gi_xyz);
                                            let Gb = Gto::new(*gj_lmn, atos[j].alps[b], gj_xyz);
                                            let Gc = Gto::new(*gk_lmn, atos[k].alps[c], gk_xyz);
                                            let Gd = Gto::new(*gl_lmn, atos[l].alps[d], gl_xyz);
                                            let gval = IntVal2e::new(&Ga, &Gb, &Gc, &Gd).gval();
                                            sum += gval * coes_i[a] * coes_j[b] * coes_k[c] * coes_l[d];
                                        }
                                    }
                                }
                            }
                            emat.idxs[idx] = true;
                            emat.vals[idx] = sum;
                        }
                    }
                }
            }
        }
        emat
    }

    /// 将笛卡尔单电子积分转换为球谐形式
    pub fn mat1e_sph(&self, car_mat: &Array2<f64>) -> Array2<f64> {
        let shls = &self.mole.basis().shls;
        let nshl = shls.len();
        let uplo_sph = shls.uplo("sph");
        let uplo_car = shls.uplo("car");
        let nato_sph = shls.nato("sph");
        let mut sph_mat = Array2::zeros((nato_sph, nato_sph));
        for i in 0..nshl {
            let angi = shls.get(i).ang;
            let [cui, cli] = uplo_car[i];
            let [sui, sli] = uplo_sph[i];
            let tmati = btrans::car2sph_tmat(angi);
            for j in 0..nshl {
                let angj = shls.get(j).ang;
                let [cuj, clj] = uplo_car[j];
                let [suj, slj] = uplo_sph[j];
                let tmatj = btrans::sph2car_tmat(angj);
                let car_part = car_mat.slice(s![cui..cli, cuj..clj]).to_owned();
                let sph_part: Array2<f64> = tmati.dot(&car_part);
                let sph_part = sph_part.dot(&tmatj);

                sph_mat.slice_mut(s![sui..sli, suj..slj]).assign(&sph_part);
            }
        }
        sph_mat
    }

    /// 将笛卡尔单电子积分转换为原始(非收缩)形式
    pub fn mat1e_raw(&self, car_mat: &Array2<f64>) -> Array2<f64> {
        let shls = &self.mole.basis().shls;
        let nshl = shls.len();
        let uplo_raw = shls.uplo("raw");
        let uplo_car = shls.uplo("car");
        let nato_raw = shls.nato("raw");
        let mut smat_raw = Array2::zeros((nato_raw, nato_raw));
        for i in 0..nshl {
            let ang_i = shls.get(i).ang;
            let [up_car_i, lo_car_i] = uplo_car[i];
            let [up_raw_i, lo_raw_i] = uplo_raw[i];
            let n_car_i = lo_car_i - up_car_i;
            let n_raw_i = lo_raw_i - up_raw_i;
            let tmati = btrans::sph2car_tmat(ang_i);
            for j in 0..nshl {
                let angj = shls.get(j).ang;
                let [up_car_j, lo_car_j] = uplo_car[j];
                let [up_raw_j, lo_raw_j] = uplo_raw[j];
                let n_car_j = lo_car_j - up_car_j;
                let n_raw_j = lo_raw_j - up_raw_j;
                let tmatj = btrans::sph2car_tmat(angj);
                let mut raw_part = car_mat.slice(s![up_car_i..lo_car_i, up_car_j..lo_car_j]).to_owned();
                if n_car_i != n_raw_i {
                    raw_part = tmati.t().dot(&raw_part)
                }
                if n_car_j != n_raw_j {
                    raw_part = raw_part.dot(&tmatj)
                }
                smat_raw
                    .slice_mut(s![up_raw_i..lo_raw_i, up_raw_j..lo_raw_j])
                    .assign(&raw_part);
            }
        }
        smat_raw
    }

    /// 获取局部坐标系的单电子积分
    pub fn mat1e_lcs(&self, stms: &HashMap<usize, Stm>, cals: [usize; 3]) -> [Array2<f64>; 3] {
        let shls = &self.mole.basis().shls;
        let uplo = shls.uplo("car");
        let nshl = shls.len();
        // let atms = &map_shl.atms;

        let mut Tmats = HashMap::new();
        for i in 0..nshl {
            let num = uplo[i][1] - uplo[i][0];
            let atm = shls.get(i).atm;
            let tmat = mocv::get_tmat_car(num, &stms[&atm].mat);
            Tmats.insert(i, tmat); // 局部坐标系系数到全局坐标系系数的转换矩阵
        }

        // 全局坐标系的重叠积分
        let [gsmat, gtmat, gvmat] = self.mat1e(cals, "car");
        let mut lsmat = gsmat.clone();
        let mut ltmat = gtmat.clone();
        let mut lvmat = gvmat.clone();

        for i in 0..nshl {
            let Mt_i = Tmats.get(&i).unwrap();
            let [ui, li] = uplo[i];

            for j in 0..nshl {
                let Mt_j = Tmats.get(&j).unwrap();
                let [uj, lj] = uplo[j];

                let info = s![ui..li, uj..lj];
                let gs = gsmat.slice(info);
                let gt = gtmat.slice(info);
                let gv = gvmat.slice(info);
                let gs: Array2<f64> = Mt_i.dot(&gs).dot(&Mt_j.t());
                let gt: Array2<f64> = Mt_i.dot(&gt).dot(&Mt_j.t());
                let gv: Array2<f64> = Mt_i.dot(&gv).dot(&Mt_j.t());
                lsmat.slice_mut(info).assign(&gs);
                ltmat.slice_mut(info).assign(&gt);
                lvmat.slice_mut(info).assign(&gv);
            }
        }
        [lsmat, ltmat, lvmat]
    }

    /// 获取局部坐标系的双电子积分
    pub fn mat2e_lcs(&self, stms: &HashMap<usize, Stm>) -> Mat2e {
        let shls = &self.mole.basis().shls;
        let uplo = shls.uplo("car");
        let nshl = uplo.len();
        // let atms = &map_shl.atms;
        let gemat = self.mat2e("car"); // 全局坐标系的双电子积分
        let nato = gemat.nato;
        let mut lemat = Mat2e::new(nato);
        // 获取每个壳层的变换矩阵
        let mut Tmats = HashMap::new();
        for i in 0..nshl {
            let num = uplo[i][1] - uplo[i][0];
            let atm = shls.get(i).atm;
            let tmat = mocv::get_tmat_car(num, &stms[&atm].mat);
            Tmats.insert(i, tmat); // 局部坐标系系数到全局坐标系系数的转换矩阵
        }

        // 先对壳层循环
        for si in 0..nshl {
            for sj in 0..nshl {
                for sk in 0..nshl {
                    for sl in 0..nshl {
                        let [ui, li] = uplo[si];
                        let [uj, lj] = uplo[sj];
                        let [uk, lk] = uplo[sk];
                        let [ul, ll] = uplo[sl];
                        let Ti = Tmats.get(&si).unwrap();
                        let Tj = Tmats.get(&sj).unwrap();
                        let Tk = Tmats.get(&sk).unwrap();
                        let Tl = Tmats.get(&sl).unwrap();
                        // 再对壳层内的基函数循环
                        for (i, ai) in (ui..li).enumerate() {
                            for (j, aj) in (uj..lj).enumerate() {
                                for (k, ak) in (uk..lk).enumerate() {
                                    for (l, al) in (ul..ll).enumerate() {
                                        let mut sum = 0.0;
                                        // 每一个局部基函数是全局基函数的线性组合
                                        for (a, aa) in (ui..li).enumerate() {
                                            for (b, ab) in (uj..lj).enumerate() {
                                                for (c, ac) in (uk..lk).enumerate() {
                                                    for (d, ad) in (ul..ll).enumerate() {
                                                        sum += Ti[(i, a)]
                                                            * Tj[(j, b)]
                                                            * Tk[(k, c)]
                                                            * Tl[(l, d)]
                                                            * gemat.get(aa, ab, ac, ad);
                                                    }
                                                }
                                            }
                                        }
                                        lemat.set(ai, aj, ak, al, sum);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        lemat
    }
}

/// 将双电子积分的四维索引压缩为一维索引
///
/// 利用(ij|kl)积分的对称性,将(i,j,k,l)四维索引映射到一维数组。
pub fn get_2eid(i: usize, j: usize, k: usize, l: usize) -> usize {
    let p = if i >= j {
        i * (i + 1) / 2 + j
    } else {
        j * (j + 1) / 2 + i
    };
    let q = if k >= l {
        k * (k + 1) / 2 + l
    } else {
        l * (l + 1) / 2 + k
    };
    if p >= q {
        p * (p + 1) / 2 + q
    } else {
        q * (q + 1) / 2 + p
    }
}
