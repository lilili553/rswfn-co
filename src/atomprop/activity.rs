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

    // 计算某个原子某个方向的活性
    pub fn vector(&self, atm: usize, dir: &Vec3, otype: &str) -> (Vec<usize>, Vec<f64>, f64) {
        let STAND = 3.0;
        let caler = bondprop::order::Calculator::new(&self.mole);
        let (nebs, vals) = caler.bound(atm, dir, otype);

        let value = STAND - vals.iter().sum::<f64>();
        (nebs, vals, value)
    }

    pub fn free_valence(&self) -> Vec<f64> {
        let caler = bondprop::order::Calculator::new(&self.mole);
        let (dirs, omat) = caler.pi_pocv();
        let atoms = self.mole.atoms();
        let natm = atoms.len();
        let mut vals = vec![0.0; natm];
        for i in 0..natm {
            for neb in atoms.nebs(i as usize) {
                vals[i] += omat[(i, neb as usize)]
            }
        }
        vals
    }
}
