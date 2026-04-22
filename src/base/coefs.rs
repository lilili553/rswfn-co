use crate::{
    base::mole::Mole,
    datas::{
        self,
        btrans::{self, len2syms},
    },
    maths,
};
use crate::{integ, utils};
use core::fmt;
use ndarray::{Array2, ArrayView2, Axis, concatenate, s};
use ndarray_linalg::{Eigh, Norm, UPLO};
use std::collections::HashMap;
use std::fmt::Debug;
use std::sync::OnceLock;

#[derive(Debug)]
pub struct Coef {
    pub atms: Vec<usize>,
    pub shls: Vec<usize>,
    pub syms: Vec<String>,
    cmat: Array2<f64>,
    pub pmat: OnceLock<Array2<f64>>,
    pub smat: OnceLock<Array2<f64>>,
}

impl Coef {
    pub fn new(atms: Vec<usize>, shls: Vec<usize>, syms: Vec<String>, cmat: Array2<f64>) -> Self {
        Self {
            atms: atms,
            shls: shls,
            syms: syms,
            cmat: cmat,
            pmat: OnceLock::new(),
            smat: OnceLock::new(),
        }
    }

    pub fn blank() -> Self {
        Self {
            atms: vec![],
            shls: vec![],
            syms: vec![],
            cmat: Array2::zeros((0, 0)),
            pmat: OnceLock::new(),
            smat: OnceLock::new(),
        }
    }

    pub fn build(&mut self, atms: Vec<usize>, shls: Vec<usize>, syms: Vec<String>, cmat: Array2<f64>) {
        self.atms = atms;
        self.shls = shls;
        self.syms = syms;
        self.cmat = cmat;
    }
}

impl Clone for Coef {
    fn clone(&self) -> Coef {
        Coef {
            atms: self.atms.clone(),
            shls: self.shls.clone(),
            syms: self.syms.clone(),
            cmat: self.cmat.clone(),
            pmat: OnceLock::new(),
            smat: OnceLock::new(),
        }
    }
}

#[derive(Debug)]
pub struct Coefs {
    pub coef_sph: Coef,
    pub coef_car: Coef,
    pub obt_engs: Vec<f64>,
    pub obt_occs: Vec<usize>,
    pub obt_idxs: OnceLock<(Vec<usize>, Vec<usize>)>, // 占据轨道和空轨道索引
}
#[rustfmt::skip]
impl Coefs {

    pub fn new(
        atms: Vec<usize>,
        shls: Vec<usize>,
        syms: Vec<String>,
        engs: Vec<f64>,
        occs: Vec<usize>,
        cmat: Array2<f64>,
    )->Self{
        // 将原始的统一转为球谐形式，要么使用纯sph，要么使用纯car，混合的没有意义
        
        let coef_car=raw2car(&atms, &shls, &syms, &cmat);
        let coef_sph=car2sph(&coef_car);
        Self{
            
            coef_sph:coef_sph,
            coef_car:coef_car,
            obt_engs:engs.clone(),
            obt_occs:occs.clone(),
            obt_idxs:OnceLock::new()
        }
    }
    pub fn blank() -> Self {
        Self {
            coef_sph: Coef::blank(),
            coef_car: Coef::blank(),
            obt_engs: vec![],
            obt_occs: vec![],
            obt_idxs: OnceLock::new(),
        }
    }

    // pub fn build(
    //     &mut self,
    //     atms: Vec<usize>,
    //     shls: Vec<usize>,
    //     syms: Vec<String>,
    //     engs: Vec<f64>,
    //     occs: Vec<usize>,
    //     cmat: Array2<f64>,
    // ){  
    //     self.coef_sph.build(ato_atms, ato_shls, ato_syms, cmat);

    //     self.obt_engs = obt_engs;
    //     self.obt_occs = obt_occs;
    // }

    pub fn is_open(&self) -> bool {
        let nato=self.coef_sph.cmat.nrows();
        let nobt=self.obt_engs.len();
        if (nobt as f64) < (nato as f64) * 1.5 {true} else {false}
    }

    pub fn is_none(&self) -> bool {
        let coef=&self.coef_sph;
        coef.atms.is_empty() || 
        coef.shls.is_empty() || 
        coef.syms.is_empty() || 
        coef.cmat.is_empty() || 
        self.obt_engs.is_empty() || 
        self.obt_occs.is_empty()
    }

    pub fn get_cmat(&self,form:&str)->&Array2<f64>{
        match form{
            "sph"=>&self.coef_sph.cmat,
            "car"=>&self.coef_car.cmat,
            _=>panic!("未知类型")
        }
    }

    pub fn set_cmat(&mut self,form:&str,cmat:&Array2<f64>){
        match form{
            "sph"=>{
                self.coef_sph.cmat=cmat.clone();
                self.coef_car=sph2car(&self.coef_sph);
            },
            "car"=>{
                self.coef_car.cmat=cmat.clone();
                self.coef_sph=car2sph(&self.coef_car);
            }
            _=>panic!("未知的form")
        };
    }

    pub fn ato_atms(&self,form:&str)->&Vec<usize>{
        match form{
            "sph"=>&self.coef_sph.atms,
            "car"=>&self.coef_car.atms,
            _=>panic!("未知类型")
        }
    }

    pub fn ato_shls(&self,form:&str)->&Vec<usize>{
        match form{
            "sph"=>&self.coef_sph.shls,
            "car"=>&self.coef_car.shls,
            _=>panic!("未知类型")
        }
    }

    pub fn ato_syms(&self,form:&str)->&Vec<String>{
        match form{
            "sph"=>&self.coef_sph.syms,
            "car"=>&self.coef_car.syms,
            _=>panic!("未知类型")
        }
    }

    pub fn ato_keys(&self,form:&str)->Vec<String>{
        let atms=self.ato_atms(form);
        let shls=self.ato_shls(form);
        let syms=self.ato_syms(form);
        let mut keys=Vec::new();
        for i in 0..atms.len(){ 
            keys.push(format!("{},{},{}",atms[i],shls[i],syms[i]));
        }
        return keys;
    }

    pub fn ato_xyzs(&self,form:&str,mole:&Mole)->Vec<[f64;3]>{
        let atms=self.ato_atms(form);
        let mut xyzs=Vec::new();
        for atm in atms{
            xyzs.push(mole.atoms.get(*atm).xyz);
        }
        xyzs
    }

    /// 获取占据轨道(otype=true)或空轨道(otype=false)的索引
    pub fn obt_idxs(&self,otype:bool)->&Vec<usize>{
        let obt_idxs=self.obt_idxs.get_or_init(||{
            let mut o_occs=Vec::new();
            let mut v_occs=Vec::new();
            for (i,occ) in self.obt_occs.iter().enumerate(){
                if *occ!=0{
                    o_occs.push(i as usize);
                }
                if *occ==0 && otype==false{
                    v_occs.push(i as usize);
                }
            }
            (o_occs,v_occs)
        });
        if otype==true{ 
            return &obt_idxs.0;
        }else {
            return &obt_idxs.1;
        }
    }

    // 获取α和β电子数
    pub fn obt_eles(&self)->[usize;2]{
        let cmat=self.get_cmat("sph");
        let nato=cmat.nrows(); // 原子轨道数
        let nobt=cmat.ncols(); // 分子轨道数
        let open=nobt as f64 > nato as f64 *1.5;
        let obt_occs=&self.obt_occs;
        // let obt_engs=&self.obt_engs;
        let mut eles=[0,0];
        if open{ // 开壳层
            for (i,ele) in obt_occs.iter().enumerate(){
                if i<nobt/2 && *ele!=0{
                    eles[0]+=1;
                }
                if i>=nobt/2 && *ele!=0{
                    eles[1]+=1;
                }
            }
        }else{
            for (_i,ele) in obt_occs.iter().enumerate(){
                if *ele!=0{
                    eles[0]+=1;
                    eles[1]+=1;
                }
            }
        }
        eles
    }

    pub fn get_dmat(&self,form:&str)->&Array2<f64>{
        let cmat=self.get_cmat(form);
        let idxs=self.obt_idxs(true);
        let open=if cmat.ncols() as f64 > cmat.nrows() as f64 *1.5{true}else{false};
        let dmat=match form {
            "sph"=>{
                let coef_raw=&self.coef_sph;
                coef_raw.pmat.get_or_init(||{
                    maths::calc_dmat(cmat, idxs,open)
                })
            },
            "car"=>{
                let coef_car=&self.coef_car;
                coef_car.pmat.get_or_init(||{
                    maths::calc_dmat(cmat, idxs,open)
                })
            },
            _=>{
                panic!("Invalid form for get_pmat");
            }
        };
        dmat
    }

    // pub fn cmat_norm(self,form:&str,mole:&Mole)->Array2<f64>{
    //     let smat=self.get_smat(form,mole);
    //     let cmat=self.get_cmat(form);
    //     let (
    //         eigval,
    //         eigvec
    //     )=smat.eigh(UPLO::Lower).expect("矩阵对角化失败");
    //     let diag=Array2::from_diag(&eigval);
    //     let xmat:Array2<f64>=eigvec.dot(&diag).dot(&eigvec.t());
    //     xmat.t().dot(cmat)
    // }

    pub fn equal(&self, other: &Coefs)->bool{
        if all_like(&self.obt_engs,&other.obt_engs)==false{
            log::info!("obt_engs不等");
            return false;
        }
        if utils::all_same(&self.obt_occs, &other.obt_occs)==false{
            log::info!("obt_occs不等");
            return false;
        }
        if utils::all_same(&self.coef_sph.atms, &other.coef_sph.atms)==false{
            log::info!("ato_atms不等");
            return false;
        }
        if utils::all_same(&self.coef_sph.shls, &other.coef_sph.shls)==false{
            log::info!("ato_shls不等");
            return false;
        }
        if utils::all_same(&self.coef_sph.syms, &other.coef_sph.syms)==false{
            log::info!("ato_syms不等");
            return false;
        }
        let cmat_0=&self.coef_sph.cmat;
        let cmat_1=&other.coef_sph.cmat;
        if cmat_0.shape()!=cmat_1.shape(){
            log::info!("cmat形状不等");
            return false;
        }
        let delta=cmat_0-cmat_1;
        let min_val = delta.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_val = delta.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        if (max_val-min_val)>1e-3{
            log::info!("cmat不等,min:{},max:{}",min_val,max_val);
            return false;
        }
        
        return true;
    }
}

fn all_like(vals0: &Vec<f64>, vals1: &Vec<f64>) -> bool {
    if vals0.len() != vals1.len() {
        return false;
    }
    for i in 0..vals0.len() {
        if (vals0[i] - vals1[i]).abs() > 1e-3 {
            return false;
        }
    }
    return true;
}

impl Clone for Coefs {
    fn clone(&self) -> Coefs {
        Coefs {
            coef_sph: self.coef_sph.clone(),
            coef_car: self.coef_car.clone(),
            obt_engs: self.obt_engs.clone(),
            obt_occs: self.obt_occs.clone(),
            obt_idxs: OnceLock::new(),
        }
    }
}

impl fmt::Display for Coefs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // let coef = &self.coef_sph;
        // let cmat = &coef.cmat;
        let cmat = self.get_cmat("car");
        let ato_atms = self.ato_atms("car");
        let ato_shls = self.ato_shls("car");
        let ato_syms = self.ato_syms("car");
        let obt_engs = &self.obt_engs;
        let obt_occs = &self.obt_occs;
        let shape = cmat.shape();
        let nato = shape[0];
        let nobt = shape[1];
        for b in (0..nobt).step_by(10) {
            write!(f, "               ")?;
            for i in b..(b + 10).min(nobt) {
                write!(f, "{:>7}{:>3}", i + 1, obt_occs[i])?;
            }
            writeln!(f, "")?;
            write!(f, "               ")?;
            for i in b..(b + 10).min(nobt) {
                write!(f, "{:>10.4}", obt_engs[i])?;
            }
            writeln!(f, "")?;
            for i in 0..nato {
                let atm = ato_atms[i];
                let shl = ato_shls[i];
                let sym = &ato_syms[i];
                write!(f, "{:>5}{:>5}{:<5}", atm + 1, shl + 1, sym)?;
                for j in b..(b + 10).min(nobt) {
                    write!(f, "{:>10.4}", cmat[(i, j)])?;
                }
                writeln!(f, "")?;
            }
        }

        Ok(())
    }
}

// 将球谐系数转为笛卡尔系数
pub fn sph2car(coef_sph: &Coef) -> Coef {
    let sph_atms = &coef_sph.atms;
    let sph_shls = &coef_sph.shls;
    let sph_syms = &coef_sph.syms;
    let sph_cmat = &coef_sph.cmat;
    let mut car_atms = Vec::new();
    let mut car_shls = Vec::new();
    let mut car_syms = Vec::new();
    let mut car_parts = Vec::new();
    let mut idx = 0;
    while idx < sph_syms.len() {
        let atm = sph_atms[idx];
        let shl = sph_shls[idx];
        let sym = &sph_syms[idx];
        // 先把球谐函数的系数取出
        let [ncar, nsph] = match sym.as_str() {
            "S" => [1, 1],
            "PX" => [3, 3],
            "D 0" => [6, 5],
            "F 0" => [10, 7],
            "G 0" => [16, 9],
            "H 0" => [21, 11],
            _ => panic!("不支持的球谐函数类型:{}", sym),
        };
        let part = sph_cmat.slice(s![idx..idx + nsph, ..]).to_owned();
        idx += nsph;
        let tmat = match nsph {
            1 => Array2::eye(1),
            3 => Array2::eye(3),
            5 => btrans::DMAT.clone(),
            7 => btrans::FMAT.clone(),
            9 => btrans::GMAT.clone(),
            11 => btrans::HMAT.clone(),
            _ => panic!("Unsupported spherical type"),
        };
        // println!("part1:{:?}", part);
        let part = tmat.dot(&part);
        // println!("part2:{:?}", part);
        car_atms.extend(vec![atm; ncar]);
        car_shls.extend(vec![shl; ncar]);
        car_syms.extend(len2syms(ncar));
        car_parts.push(part);
    }
    let car_views = car_parts.iter().map(|part| part.view()).collect::<Vec<_>>();
    let car_cmat = concatenate(Axis(0), car_views.as_slice()).unwrap();
    Coef::new(car_atms, car_shls, car_syms, car_cmat)
}

pub fn car2sph(coef_car: &Coef) -> Coef {
    let car_atms = &coef_car.atms;
    let car_shls = &coef_car.shls;
    let car_syms = &coef_car.syms;
    let car_cmat = &coef_car.cmat;
    let mut sph_atms = Vec::new();
    let mut sph_shls = Vec::new();
    let mut sph_syms = Vec::new();
    let mut sph_parts = Vec::new();
    let mut idx = 0;
    while idx < car_syms.len() {
        let atm = car_atms[idx];
        let shl = car_shls[idx];
        let sym = &car_syms[idx];
        // 先把球谐函数的系数取出
        let [num_car, num_sph] = match sym.as_str() {
            "S" => [1, 1],
            "PX" => [3, 3],
            "XX" => [6, 5],
            "XXX" => [10, 7],
            "ZZZZ" => [15, 9],
            "ZZZZZ" => [21, 7],
            _ => panic!("Unsupported spherical type"),
        };
        let part = car_cmat.slice(s![idx..idx + num_car, ..]).to_owned();
        idx += num_car;
        let tmat = match num_car {
            1 => Array2::eye(1),
            3 => Array2::eye(3),
            6 => btrans::DMAT.clone(),
            10 => btrans::FMAT.clone(),
            15 => btrans::GMAT.clone(),
            21 => btrans::HMAT.clone(),
            _ => panic!("Unsupported spherical type"),
        };
        let part = tmat.t().dot(&part);
        sph_atms.extend(vec![atm; num_sph]);
        sph_shls.extend(vec![shl; num_sph]);
        sph_syms.extend(len2syms(num_sph));
        sph_parts.push(part);
    }
    let sph_views = sph_parts.iter().map(|part| part.view()).collect::<Vec<_>>();
    let sph_cmat = concatenate(Axis(0), sph_views.as_slice()).unwrap();
    Coef::new(sph_atms, sph_shls, sph_syms, sph_cmat)
}

pub fn raw2car(atms: &Vec<usize>, shls: &Vec<usize>, syms: &Vec<String>, cmat: &Array2<f64>) -> Coef {
    assert_ne!(syms.len(), 0, "原始系数的符号列表不能为空");
    let mut car_atms = vec![];
    let mut car_shls = vec![];
    let mut car_syms = vec![];
    let mut car_parts = vec![];
    let mut idx = 0;
    while idx < syms.len() {
        let atm = atms[idx];
        let shl = shls[idx];
        let sym = &syms[idx];
        let [nraw, ncar] = match sym.as_str() {
            "S" => [1, 1],
            "PX" => [3, 3],
            "XX" => [6, 6],
            "XXX" => [10, 10],
            "ZZZZ" => [15, 15],
            "ZZZZZ" => [21, 21],
            "D 0" => [5, 6],
            "F 0" => [7, 10],
            "G 0" => [9, 15],
            "H 0" => [11, 21],
            _ => panic!("未知符号"),
        };
        let mut part = cmat.slice(s![idx..idx + nraw, ..]).to_owned();
        idx += nraw;
        // 不相等时说明需要转换
        if nraw != ncar {
            let tmat = btrans::get_mat(ncar);
            part = tmat.dot(&part);
        }
        car_parts.push(part);
        car_atms.extend(vec![atm; ncar]);
        car_shls.extend(vec![shl; ncar]);
        car_syms.extend(len2syms(ncar));
    }
    let car_views = car_parts.iter().map(|part| part.view()).collect::<Vec<_>>();
    let car_cmat = concatenate(Axis(0), car_views.as_slice()).expect("数组拼接失败");
    Coef::new(car_atms, car_shls, car_syms, car_cmat)
}
