use crate::{
    base::Basis,
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
use std::fmt::Debug;
use std::sync::OnceLock;
use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

const ANG_LEN_CAR: [usize; 5] = [1, 3, 6, 10, 15];
const ANG_LEN_SPH: [usize; 5] = [1, 3, 5, 7, 9];

#[derive(Debug)]
pub struct Coef {
    pub atms: Vec<usize>,
    pub shls: Vec<usize>,
    pub syms: Vec<String>,
    cmat: Array2<f64>,
    _dmat: OnceLock<Array2<f64>>,
    _uplo: OnceLock<Vec<[usize; 2]>>,
}

impl Coef {
    pub fn new(atms: &Vec<usize>, shls: &Vec<usize>, syms: &Vec<String>, cmat: &Array2<f64>) -> Self {
        Self {
            atms: atms.clone(),
            shls: shls.clone(),
            syms: syms.clone(),
            cmat: cmat.clone(),
            _dmat: OnceLock::new(),
            _uplo: OnceLock::new(),
        }
    }

    pub fn uplo(&self) -> &Vec<[usize; 2]> {
        self._uplo.get_or_init(|| {
            let mut uplo = vec![[0, 0]];
            let mut loc = self.atms[0];
            for (i, atm) in self.atms.iter().enumerate() {
                if *atm != loc {
                    uplo.push([i, 0]);
                    loc = *atm;
                }
            }
            let len = uplo.len();
            for i in 0..len - 1 {
                uplo[i][1] = uplo[i + 1][0];
            }
            uplo[len - 1][1] = self.atms.len();
            uplo
        })
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
            _dmat: OnceLock::new(),
            _uplo: OnceLock::new(),
        }
    }
}

#[derive(Debug)]
pub struct Coefs {
    pub _coef_raw: Coef,
    pub _coef_sph: OnceLock<Coef>,
    pub _coef_car: OnceLock<Coef>,
    pub obt_engs: Vec<f64>,
    pub obt_occs: Vec<f64>,
    pub obt_idxs: OnceLock<(Vec<usize>, Vec<usize>)>, // 占据轨道和空轨道索引
}

#[rustfmt::skip]
impl Coefs {
    pub fn new(
        ato_atms: Vec<usize>,
        ato_shls: Vec<usize>,
        ato_syms: Vec<String>,
        obt_engs: Vec<f64>,
        obt_occs: Vec<f64>,
        cmat: Array2<f64>,
    )->Self{
        // 将原始的统一转为球谐形式，要么使用纯sph，要么使用纯car，混合的没有意义
        let coef_raw=Coef::new(&ato_atms, &ato_shls, &ato_syms, &cmat);
        Self{
            _coef_raw:coef_raw,
            _coef_sph:OnceLock::new(),
            _coef_car:OnceLock::new(),
            obt_engs:obt_engs,
            obt_occs:obt_occs,
            obt_idxs:OnceLock::new(),
        }
    }

    pub fn new_ref(
        ato_atms: &Vec<usize>,
        ato_shls: &Vec<usize>,
        ato_syms: &Vec<String>,
        obt_engs: &Vec<f64>,
        obt_occs: &Vec<f64>,
        cmat: &Array2<f64>,
    ) -> Self {
        Self {
            _coef_raw: Coef::new(ato_atms, ato_shls, ato_syms, cmat),
            _coef_sph: OnceLock::new(),
            _coef_car: OnceLock::new(),
            obt_engs: obt_engs.clone(),
            obt_occs: obt_occs.clone(),
            obt_idxs: OnceLock::new(),
        }
    }

    pub fn coef_raw(&self) -> &Coef {
        &self._coef_raw
    }

    pub fn coef_sph(&self, basis: &Basis) -> &Coef {
        self._coef_sph.get_or_init(|| {
            let (cmat,atms,shls,syms)=from_raw(&self._coef_raw.cmat, basis, "sph");
            Coef::new(&atms, &shls, &syms, &cmat)
        })
    }

    pub fn coef_car(&self, basis: &Basis) -> &Coef {
        self._coef_car.get_or_init(|| {
            let (cmat,atms,shls,syms)=from_raw(&self._coef_raw.cmat, basis, "car");
            Coef::new(&atms, &shls, &syms, &cmat)
        })
    }

    pub fn is_none(&self) -> bool {
        let coef=&self._coef_raw;
        coef.atms.is_empty() ||
        coef.shls.is_empty() ||
        coef.syms.is_empty() ||
        coef.cmat.is_empty() ||
        self.obt_engs.is_empty() ||
        self.obt_occs.is_empty()
    }

    pub fn get_cmat(&self,form:&str,basis: &Basis)->&Array2<f64>{
        match form{
            "raw"=>&self.coef_raw().cmat,
            "sph"=>&self.coef_sph(basis).cmat,
            "car"=>&self.coef_car(basis).cmat,
            _=>panic!("未知类型")
        }
    }

    pub fn set_cmat(&mut self,form:&str,cmat:&Array2<f64>,basis: &Basis){
        self._coef_car.take();
        self._coef_sph.take();
        match form{
            "raw"=>{
                self._coef_raw.cmat=cmat.clone();
            }
            "sph"=>{
                self._coef_raw.cmat=into_raw(cmat, basis, "sph").0;
            },
            "car"=>{
                self._coef_raw.cmat=into_raw(cmat, basis, "car").0;
            }
            _=>panic!("未知的form")
        };
    }

    /// 获取占据轨道(otype=true)或空轨道(otype=false)的索引
    pub fn obt_idxs(&self,otype:bool)->&Vec<usize>{
        let obt_idxs=self.obt_idxs.get_or_init(||{
            let mut o_occs=Vec::new();
            let mut v_occs=Vec::new();
            for (i,occ) in self.obt_occs.iter().enumerate(){
                if *occ!=0.0{
                    o_occs.push(i);
                }
                if *occ==0.0{
                    v_occs.push(i);
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
    pub fn obt_eles(&self,basis: &Basis)->[usize;2]{
        let cmat=self.get_cmat("sph",basis);
        let nato=cmat.nrows(); // 原子轨道数
        let nobt=cmat.ncols(); // 分子轨道数
        let open=nobt as f64 > nato as f64 *1.5;
        let obt_occs=&self.obt_occs;
        // let obt_engs=&self.obt_engs;
        let mut eles=[0,0];
        if open{ // 开壳层
            for (i,ele) in obt_occs.iter().enumerate(){
                if i<nobt/2 && *ele!=0.0{
                    eles[0]+=1;
                }
                if i>=nobt/2 && *ele!=0.0{
                    eles[1]+=1;
                }
            }
        }else{
            for (_i,ele) in obt_occs.iter().enumerate(){
                if *ele!=0.0{
                    eles[0]+=1;
                    eles[1]+=1;
                }
            }
        }
        eles
    }

    pub fn get_dmat(&self,form:&str,basis: &Basis)->&Array2<f64>{
        // let mole=self.mole();
        let cmat=self.get_cmat(form,basis);
        let dmat=match form {
            "raw"=>{
                let coef_raw=self.coef_raw();
                coef_raw._dmat.get_or_init(||{
                    maths::calc_dmat(cmat, &self.obt_occs)
                })
            }
            "sph"=>{
                let coef_sph=self.coef_sph(basis);
                coef_sph._dmat.get_or_init(||{
                    maths::calc_dmat(cmat, &self.obt_occs)
                })
            },
            "car"=>{
                let coef_car=self.coef_car(basis);
                coef_car._dmat.get_or_init(||{
                    maths::calc_dmat(cmat, &self.obt_occs)
                })
            },
            _=>{
                panic!("Invalid form for get_dmat");
            }
        };
        dmat
    }

    pub fn atm_uplo(&self,form:&str,basis: &Basis)->&Vec<[usize;2]>{
        match form{
            "raw"=>self.coef_raw().uplo(),
            "sph"=>self.coef_sph(basis).uplo(),
            "car"=>self.coef_car(basis).uplo(),
            _=>{
                panic!("Invalid form for get_dmat");
            }
        }
    }

    pub fn equal(&self, other: &Coefs)->bool{
        if all_like(&self.obt_engs,&other.obt_engs)==false{
            log::info!("obt_engs不等");
            return false;
        }
        if utils::all_same(&self.obt_occs, &other.obt_occs)==false{
            log::info!("obt_occs不等");
            return false;
        }
        if utils::all_same(&self._coef_raw.atms, &other._coef_raw.atms)==false{
            log::info!("ato_atms不等");
            return false;
        }
        if utils::all_same(&self._coef_raw.shls, &other._coef_raw.shls)==false{
            log::info!("ato_shls不等");
            return false;
        }
        if utils::all_same(&self._coef_raw.syms, &other._coef_raw.syms)==false{
            log::info!("ato_syms不等");
            return false;
        }
        let cmat_0=&self._coef_raw.cmat;
        let cmat_1=&other._coef_raw.cmat;
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
            _coef_raw: self._coef_raw.clone(),
            _coef_sph: OnceLock::new(),
            _coef_car: OnceLock::new(),
            obt_engs: self.obt_engs.clone(),
            obt_occs: self.obt_occs.clone(),
            obt_idxs: OnceLock::new(),
        }
    }
}

impl fmt::Display for Coefs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let cmat = &self.coef_raw().cmat;
        // let ato_atms = self.ato_atms("raw");
        // let ato_shls = self.ato_shls("raw");
        // let ato_syms = self.ato_syms("raw");
        let obt_engs = &self.obt_engs;
        let obt_occs = &self.obt_occs;
        let shape = cmat.shape();
        let nato = shape[0];
        let nobt = shape[1];
        for b in (0..nobt).step_by(10) {
            // write!(f, "               ")?;
            for i in b..(b + 10).min(nobt) {
                write!(f, "{:>7}{:>3}", i + 1, obt_occs[i])?;
            }
            writeln!(f, "")?;
            // write!(f, "               ")?;
            for i in b..(b + 10).min(nobt) {
                write!(f, "{:>10.4}", obt_engs[i])?;
            }
            writeln!(f, "")?;
            for i in 0..nato {
                // let atm = ato_atms[i];
                // let shl = ato_shls[i];
                // let sym = &ato_syms[i];
                // write!(f, "{:>5}{:>5}{:<5}", atm + 1, shl + 1, sym)?;
                for j in b..(b + 10).min(nobt) {
                    write!(f, "{:>10.4}", cmat[(i, j)])?;
                }
                writeln!(f, "")?;
            }
        }

        Ok(())
    }
}

pub fn from_raw(
    cmat_raw: &Array2<f64>,
    basis: &Basis,
    form: &str,
) -> (Array2<f64>, Vec<usize>, Vec<usize>, Vec<String>) {
    let shls = &basis.shls;
    let mut atms_stp = vec![];
    let mut shls_stp = vec![];
    let mut syms_stp = vec![];
    let mut idx = 0;
    let uplo_raw = shls.uplo("raw");
    let uplo_stp = shls.uplo(form);
    let nato_raw = shls.nato("raw");
    let nato_stp = shls.nato(form);
    let mut cmat_stp = Array2::zeros((nato_stp, cmat_raw.ncols()));
    for i in 0..shls.len() {
        let shl = shls.get(i);
        let ang = shl.ang;
        let atm = shl.atm;
        let [up_raw, lo_raw] = uplo_raw[i];
        let [up_stp, lo_stp] = uplo_stp[i];
        let nraw = lo_raw - up_raw;
        let nstp = lo_stp - up_stp;
        let cmat_raw_part = if nraw != nstp {
            let tmat = match form {
                "car" => btrans::sph2car_mat(ang),
                "sph" => btrans::car2sph_mat(ang),
                _ => panic!("Unsupported form: {}", form),
            };
            let part = cmat_raw.slice(s![up_raw..lo_raw, ..]);
            tmat.dot(&part)
        } else {
            cmat_raw.slice(s![up_raw..lo_raw, ..]).to_owned()
        };
        cmat_stp.slice_mut(s![up_stp..lo_stp, ..]).assign(&cmat_raw_part);
        atms_stp.extend(vec![shl.atm; nstp]);
        shls_stp.extend(vec![shl.shl; nstp]);
        syms_stp.extend(len2syms(nstp));
    }
    (cmat_stp, atms_stp, shls_stp, syms_stp)
}

pub fn into_raw(
    cmat_stp: &Array2<f64>,
    basis: &Basis,
    form: &str,
) -> (Array2<f64>, Vec<usize>, Vec<usize>, Vec<String>) {
    let shls = &basis.shls;
    let mut atms_raw = vec![];
    let mut shls_raw = vec![];
    let mut syms_raw = vec![];
    let uplo_raw = shls.uplo("raw");
    let uplo_stp = shls.uplo(form);
    let nato_raw = shls.nato("raw");
    let nato_stp = shls.nato(form);
    let mut cmat_raw = Array2::zeros((nato_raw, cmat_stp.ncols()));
    for i in 0..shls.len() {
        let shl = shls.get(i);
        let ang = shl.ang;
        let atm = shl.atm;
        let [up_raw, lo_raw] = uplo_raw[i];
        let [up_stp, lo_stp] = uplo_stp[i];
        let nraw = lo_raw - up_raw;
        let nstp = lo_stp - up_stp;
        let cmat_stp_part = if nraw != nstp {
            let tmat = match form {
                "car" => btrans::car2sph_mat(ang),
                "sph" => btrans::sph2car_mat(ang),
                _ => panic!("Unsupported form: {}", form),
            };
            let part = cmat_stp.slice(s![up_stp..lo_stp, ..]);
            tmat.dot(&part)
        } else {
            cmat_stp.slice(s![up_stp..lo_stp, ..]).to_owned()
        };
        cmat_raw.slice_mut(s![up_raw..lo_raw, ..]).assign(&cmat_stp_part);
        atms_raw.extend(vec![shl.atm; nraw]);
        shls_raw.extend(vec![shl.shl; nraw]);
        syms_raw.extend(len2syms(nraw));
    }
    (cmat_raw, atms_raw, shls_raw, syms_raw)
}

unsafe impl Send for Coefs {}
unsafe impl Sync for Coefs {}
