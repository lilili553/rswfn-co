#![allow(non_snake_case)]

use ndarray::Array2;
use std::collections::HashMap;

use crate::base::{Mole, Vec3};

use crate::atomprop;
use crate::bondprop;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    pub fn fukui(&self, mole_n: &Mole, mole_p: &Mole, ctype: &str) -> Vec<[f64; 7]> {
        let caler_0 = atomprop::charge::Calculator::new(self.mole);
        let caler_n = atomprop::charge::Calculator::new(&mole_n); // 多一个电子的分子
        let caler_p = atomprop::charge::Calculator::new(&mole_p); // 少一个电子的分子

        let eles_0 = caler_0.eles(ctype);
        let eles_n = caler_n.eles(ctype);
        let eles_p = caler_p.eles(ctype);

        let nele = eles_n.len();
        let mut vals = Vec::new();
        for i in 0..nele {
            let e0 = eles_0[i];
            let en = eles_n[i];
            let ep = eles_p[i];
            let _fp = en - e0;
            let _fn = e0 - ep;
            let f0 = (_fn + _fp) / 2.0;
            let df = _fp - _fn;
            vals.push([e0, en, ep, _fn, _fp, f0, df]);
        }
        vals
    }

    pub fn fukui_pi(&self, mole_n: &Mole, mole_p: &Mole, ctype: &str) -> (HashMap<usize, Vec3>, Vec<[f64; 7]>) {
        let caler_0 = atomprop::charge::Calculator::new(self.mole);
        let caler_n = atomprop::charge::Calculator::new(&mole_n);
        let caler_p = atomprop::charge::Calculator::new(&mole_p);

        let (_dirs, eles_n) = caler_n.pi_pocv(ctype);
        let (_dirs, eles_p) = caler_p.pi_pocv(ctype);
        let (dirs, eles_0) = caler_0.pi_pocv(ctype);

        let nele = eles_n.len();
        let mut vals = Vec::new();
        for i in 0..nele {
            let e0 = eles_0[i];
            let en = eles_n[i];
            let ep = eles_p[i];
            let _fp = en - e0;
            let _fn = e0 - ep;
            let f0 = (_fn + _fp) / 2.0;
            let df = _fp - _fn;
            vals.push([e0, en, ep, _fn, _fp, f0, df]);
        }
        (dirs, vals)
    }

    // 一个原子在某个方向的福井函数
    pub fn fukui_dir(&self, atm: usize, dir: &Vec3, mole_n: &Mole, mole_p: &Mole, ctype: &str) -> [f64; 7] {
        let caler_0 = atomprop::charge::Calculator::new(self.mole);
        let caler_n = atomprop::charge::Calculator::new(&mole_n);
        let caler_p = atomprop::charge::Calculator::new(&mole_p);

        let mut dirs_map = HashMap::new();
        dirs_map.insert(atm, dir.clone());

        let eles_0 = caler_0.pocv(&dirs_map, false, false, ctype);
        let eles_n = caler_n.pocv(&dirs_map, false, false, ctype);
        let eles_p = caler_p.pocv(&dirs_map, false, false, ctype);

        let e0 = eles_0[atm as usize];
        let en = eles_n[atm as usize];
        let ep = eles_p[atm as usize];
        let _fp = en - e0;
        let _fn = e0 - ep;
        let f0 = (_fn + _fp) / 2.0;
        let df = _fp - _fn;
        [e0, en, ep, _fn, _fp, f0, df]
    }

    // 自由价 freevalence
    pub fn freev(&self, otype: &str) -> Vec<[f64; 3]> {
        let caler = bondprop::order::Calculator::new(self.mole);
        let omat = caler.order(otype);
        let natm = self.mole.atoms.len();
        let mut crgs = self.mole.atoms.crgs();
        let mut vals = vec![[0.0; 3]; natm];
        for ai in 0..natm {
            vals[ai][0] = element_valence(crgs[ai] as usize)[0] as f64;
            for aj in self.mole.atoms.nebs(ai) {
                vals[ai][1] += omat[(ai, aj)];
            }
            vals[ai][2] = vals[ai][0] - vals[ai][1]
        }
        vals
    }

    pub fn freev_pi(&self, otype: &str) -> (HashMap<usize, Vec3>, Vec<[f64; 3]>) {
        let caler = bondprop::order::Calculator::new(&self.mole);
        let (dirs, omat) = caler.pi_pocv(otype);
        let natm = self.mole.atoms.len();
        let mut crgs = self.mole.atoms.crgs();
        let mut vals = vec![[0.0; 3]; natm];
        for ai in 0..natm {
            vals[ai][0] = element_valence(crgs[ai] as usize)[1] as f64;
            for aj in self.mole.atoms.nebs(ai) {
                vals[ai][1] += omat[(ai, aj)];
            }
            vals[ai][2] = vals[ai][0] - vals[ai][1]
        }
        (dirs, vals)
    }

    // 计算某个原子某个方向的活性 活性矢量/方向自由价
    pub fn freev_dir(&self, atm: usize, dir: &Vec3, nebs: &Vec<usize>, otype: &str) -> (f64, Vec<f64>, f64) {
        let crg = self.mole.atoms.crgs()[atm];
        let STAND = element_valence(crg as usize)[1];
        let caler = bondprop::order::Calculator::new(&self.mole);
        let vals = caler.bound(atm, dir, nebs, otype);
        let value = STAND - vals.iter().sum::<f64>();
        (STAND, vals, value)
    }
}

// 返回元素的价数和pi价数
pub fn element_valence(crg: usize) -> [f64; 2] {
    match crg {
        6 => [4.0, 3.0],
        _ => [4.0, 3.0],
    }
}
