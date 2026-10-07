#![allow(non_snake_case)]
use cached::proc_macro::cached;
use ndarray::{Array2, s};
use std::collections::HashMap;
use std::fmt;

use crate::atomprop;
use crate::base::{Mole, Stm, Vec3};
use crate::datas;
use crate::datas::mocv;
use crate::maths;

/// 分子轨道系数矩阵操作计算器，支持 POCV、MOCV 和轨道定域化等
pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub vobt: bool, // 是否包含空轨道
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator {
            mole: mole,
            vobt: false,
            atms: atms,
        }
    }
    /// 使用pocv法计算的投影分子轨道系数矩阵
    ///
    /// # Arguments
    ///
    /// * `dirs` - 每个原子对应的方向向量
    /// * `keep_other_atm` - 其它原子轨道系数是否保留
    /// * `keep_other_sym` - 是否保留其它类型的轨道
    pub fn pocv(
        &self,
        dirs: &HashMap<usize, Vec3>,
        keep_other_atm: bool,
        keep_other_sym: bool,
    ) -> Result<Array2<f64>, String> {
        let atoms = self.mole.atoms();
        let uls = self.mole.basis().atm_uplo("raw");
        let mut atms: Vec<usize> = dirs.keys().cloned().collect();
        atms.sort_by(|x, y| x.cmp(y));
        let cmat_raw = self.mole.coefs().get_cmat("raw", self.mole);
        let [nato, nobt] = cmat_raw.shape().try_into().expect("获取矩阵形状失败");
        // 初始化投影的分子轨道
        let mut cmat_pocv: Array2<f64> = Array2::zeros([nato, nobt]);

        let ato_atms = self.mole.basis().atos().atms("raw");
        let ato_syms = self.mole.basis().atos().syms("raw");

        let obts = if self.vobt {
            &(0..nobt).into_iter().collect()
        } else {
            self.mole.coefs().obt_idxs(true)
        }; // 占据轨道的索引
        for obt in obts.iter() {
            for ato in 0..nato {
                let atm = &ato_atms[ato];
                let sym = &ato_syms[ato];
                if atms.contains(atm) && sym.contains("P") {
                    cmat_pocv[(ato, *obt)] = cmat_raw[(ato, *obt)]
                }
                if !atms.contains(atm) && keep_other_atm {
                    cmat_pocv[(ato, *obt)] = cmat_raw[(ato, *obt)]
                }
                if !sym.contains("P") && keep_other_sym {
                    cmat_pocv[(ato, *obt)] = cmat_raw[(ato, *obt)]
                }
            }
        }

        let mut pidxs = Vec::new(); // p轨道的索引
        let ato_syms = self.mole.basis().atos().syms("raw");

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
        assert_eq!(pidxs.len() % 3, 0, "原子P轨道的数量必须是3的倍数:{}", pidxs.len());
        for i in (0..pidxs.len()).step_by(3) {
            let [atm, ato] = pidxs[i];
            let dir = dirs.get(&(atm as usize)).expect("dir获取atm失败");
            let [dx, dy, dz] = dir.into_arr();
            for obt in obts.iter() {
                let px = &cmat_pocv[(ato + 0, *obt as usize)];
                let py = &cmat_pocv[(ato + 1, *obt as usize)];
                let pz = &cmat_pocv[(ato + 2, *obt as usize)];
                let [tx, ty, tz] = maths::vec_project([*px, *py, *pz], [dx, dy, dz])?;
                cmat_pocv[(ato + 0, *obt as usize)] = tx;
                cmat_pocv[(ato + 1, *obt as usize)] = ty;
                cmat_pocv[(ato + 2, *obt as usize)] = tz;
            }
        }
        Ok(cmat_pocv)
    }

    /// 使用pocv方法计算得到pi分子轨道系数矩阵
    pub fn pi_pocv(&self) -> Result<(HashMap<usize, Vec3>, Array2<f64>), String> {
        let dir_caler = atomprop::direction::Calculator::new(self.mole);
        let mut dirs = HashMap::new();
        let atoms = self.mole.atoms();
        for atm in self.atms {
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
        let cmat = self.pocv(&dirs, false, false)?;
        Ok((dirs, cmat))
    }

    /// 轨道映射法，需要使用笛卡尔类型原子轨道
    ///
    /// # Arguments
    ///
    /// * `stms` - 原子局部坐标系
    /// * `atos` - 原子轨道索引
    pub fn mocv(&self, stms: &HashMap<usize, Stm>, atos: &Vec<usize>) -> Array2<f64> {
        let cmat_raw = self.mole.coefs().get_cmat("sph", self.mole);
        let nato = cmat_raw.nrows();
        let nobt = cmat_raw.ncols();
        if atos.len() == 0 {
            log::warn!("mocv:原子轨道长度为0");
            return Array2::zeros((nato, nobt));
        }
        let mut cmat_mocv = Array2::zeros([nato, nobt]);

        let mut ato_keeps = vec![0; nato];
        for ato in atos {
            ato_keeps[*ato] = 1;
        }

        let shls = &self.mole.basis().shls;
        let uplo = shls.uplo("sph");

        assert_eq!(shls.len(), uplo.len());
        for i in 0..shls.len() {
            let shl = shls.get(i);
            let atm = shl.atm;
            let ang = shl.ang;
            // 没有局部坐标系的原子跳过
            if !stms.contains_key(&atm) {
                continue;
            }

            let [up, lo] = uplo[i];
            let ato_keep = ato_keeps[up..lo].to_vec();
            let keeps: Vec<f64> = (0..lo - up).map(|i| if ato_keep[i] == 1 { 1.0 } else { 0.0 }).collect();
            let keeps = Array2::from_shape_vec((lo - up, 1), keeps).unwrap();
            let rcoefs = cmat_raw.slice(s![up..lo, ..]).to_owned(); // 原始的系数
            let L = &stms[&atm].mat;
            let tcoefs = mocv(L, &rcoefs, &keeps);
            cmat_mocv.slice_mut(s![up..lo, ..]).assign(&tcoefs);
        }
        cmat_mocv
    }

    /// 使用mocv方法计算得到pi分子轨道系数矩阵
    /// 使用 MOCV 方法计算 π 分子轨道系数矩阵
    pub fn pi_mocv(&self, dest: &MocvDest) -> (HashMap<usize, Stm>, Vec<usize>, Array2<f64>) {
        let mut dir_caler = atomprop::direction::Calculator::new(self.mole);
        let mut stms = HashMap::new();
        let atoms = self.mole.atoms();
        for i in 0..atoms.len() {
            if atoms.get(i as usize).nuc == 1 {
                continue;
            }
            let stm = match dir_caler.LCS(i as usize, None) {
                Some(stm) => stm,
                None => {
                    log::warn!("原子{}没有局部坐标系，使用随机坐标系", i);
                    dir_caler.RLCS(i as usize)
                } // 如果有局部坐标系的话使用局部坐标系否则使用随机座标系？
            };
            stms.insert(i as usize, stm);
        }
        let mut ato_idxs = mocv_atos(self.mole, dest);
        log::info!("选择原子轨道索引:{:?}", ato_idxs);
        let cmat = self.mocv(&stms, &ato_idxs);
        (stms, ato_idxs, cmat)
    }

    /// 原子局部坐标系的分子轨道系数矩阵
    /// 将分子轨道系数矩阵转换到原子局部坐标系
    pub fn ldao(&self, stms: &HashMap<usize, Stm>) -> Array2<f64> {
        let shls = &self.mole.basis().shls;
        let cmat_mole = self.mole.coefs().get_cmat("car", self.mole).clone();
        let mut cmat_ldao = cmat_mole.clone();
        let uplo = shls.uplo("car");
        let nshl = shls.len();
        for i in 0..nshl {
            let [up, lo] = uplo[i];
            let atm = shls.get(i).atm;
            if stms.contains_key(&atm) {
                let tmat = mocv::get_tmat_car(lo - up, &stms[&atm].mat).t().to_owned();
                let tmat = maths::inv(&tmat);
                let info = s![up..lo, ..];
                let coefs = cmat_mole.slice(info);
                cmat_ldao.slice_mut(info).assign(&tmat.dot(&coefs));
            } else {
                println!("没有计算:{}", atm);
            }
        }
        cmat_ldao
    }

    /// 轨道定域化分子轨道
    /// 使用 Pipek-Mezey 方法进行轨道定域化
    pub fn LMO(&self) {
        let mut cmat = self.mole.coefs().get_cmat("raw", self.mole).to_owned();
        let smat = self.mole.basis().get_smat("raw", self.mole).to_owned();
        let natm = self.mole.atoms().len();
        let uplo = self.mole.basis().atm_uplo("raw");
        // let [nato, nobt] = cmat.shape().try_into().expect("获取矩阵形状失败");
        let obts = self.mole.coefs().obt_idxs(true);
        let nobt = obts.len();
        let mut pval_old = 0.0;
        for _cyl in 0..30 {
            // 最多循环30次
            for i in 0..nobt - 1 {
                let oi = obts[i] as usize;
                let cmi = cmat.slice(s![.., oi]).insert_axis(ndarray::Axis(1)).to_owned(); // 第i个分子轨道系数
                for j in i + 1..nobt {
                    let oj = obts[j] as usize;
                    let cmj = cmat.slice(s![.., oj]).insert_axis(ndarray::Axis(1)).to_owned(); // 第j个分子轨道系数
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
                let cmi = cmat.slice(s![.., imo]).insert_axis(ndarray::Axis(1)).to_owned();
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

/// 使用 MOCV 方法将全局坐标系系数映射到局部坐标系再映射回来
pub fn mocv(L: &Array2<f64>, rcoefs: &Array2<f64>, keeps: &Array2<f64>) -> Array2<f64> {
    let nrow = rcoefs.nrows();
    let Mo = if [1, 3, 5, 7, 9, 11].contains(&nrow) {
        mocv::get_tmat_sph(nrow, L)
    } else {
        // mocv::get_tmat_car(nrow, L)
        panic!("只能使用球谐轨道系数")
    };
    // let Mi = mocv::get_tmat_sph(nrow, L);
    let Mt = Mo.t().to_owned();
    let Mi = maths::inv(&Mt);
    let mut tcoefs = Mi.dot(rcoefs); //全局系数映射到局部坐标系
    tcoefs = &tcoefs * keeps; // 保留某个维度
    let fcoefs = Mt.dot(&tcoefs); // 再映射回全局坐标系
    fcoefs
}

/// dest：计算目标，可以为 0:popul 或 1:order
/// 根据计算目标选择用于 MOCV 的原子轨道索引
pub fn mocv_atos(mole: &Mole, dest: &MocvDest) -> Vec<usize> {
    let mut atos = vec![];
    let ato_syms = mole.basis().atos().syms("sph");
    let ato_atms = mole.basis().atos().atms("sph");

    let mut atm_syms = HashMap::new(); // 每个原子保留的符号
    let natm = mole.atoms.len();
    for i in 0..natm {
        let geom = if mole.atoms.is_line(i) { 1 } else { 0 };
        let syms = mocv_syms(dest, geom);
        atm_syms.insert(i, syms);
    }
    let nato = ato_syms.len();
    let ato_locs = mole.basis().atos().locs("sph");
    for i in 0..nato {
        let atm = ato_atms[i];
        let sym = &ato_syms[i];
        let keeps = atm_syms.get(&atm).unwrap();
        if keeps.contains(sym) {
            atos.push(i);
        }
    }
    atos
}

pub enum OrderDType {
    SI, // σ键级
    PI, // π键级
    DE, // δ键级
    PH, // φ键级
}

pub enum MocvDest {
    PiProperty,             // pi电子性质
    OrderDecom(OrderDType), // 键级分解
    RoundOrder,             // 环轴键级
}

impl MocvDest {
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "PiProperty" => Some(MocvDest::PiProperty),
            "OrderDecom-SI" => Some(MocvDest::OrderDecom(OrderDType::SI)),
            "OrderDecom-PI" => Some(MocvDest::OrderDecom(OrderDType::PI)),
            "OrderDecom-DE" => Some(MocvDest::OrderDecom(OrderDType::DE)),
            "OrderDecom-PH" => Some(MocvDest::OrderDecom(OrderDType::PH)),
            "RoundOrder" => Some(MocvDest::RoundOrder),
            _ => None,
        }
    }
}

// 计算目标：广义pi性质-平面/直线、键级分解、环轴键级
/// dest：计算目标，可以为 MocvDest 枚举值
/// geom：几何类型，可以为 0:plan 或 1:line
// #[cached]
/// 根据计算目标和几何类型返回要保留的轨道对称性
fn mocv_syms(dest: &MocvDest, geom: usize) -> Vec<String> {
    let syms = match (dest, geom) {
        (MocvDest::PiProperty, 0) => {
            vec![
                "PX", // P轨道
                // "XZ", "XY", // 笛卡尔D轨道
                "D+1", "D-2", // 球谐D轨道
                // "XXY", "XXZ", "XYY", "XYZ", "XZZ", // 笛卡尔F轨道
                "F-2", "F+1", "F+3", // 球谐F轨道
            ] // 平面分子电子布局
        }
        (MocvDest::PiProperty, 1) => {
            vec![
                "PX", "PY", "D+1", "D-1", "D+2", "D-2", "F+1", "F-1", "F+2", "F-2", "F+3", "F-3",
            ] // 线性分子电子布居
        }
        (MocvDest::OrderDecom(OrderDType::SI), _) => {
            vec!["S", "PZ", "D 0", "F 0"]
        }
        (MocvDest::OrderDecom(OrderDType::PI), _) => {
            vec!["PX", "PY", "D+1", "D-1", "F+1", "F-1"]
        }
        (MocvDest::OrderDecom(OrderDType::DE), _) => {
            vec!["D+2", "D-2", "F+2", "F-2"]
        }
        (MocvDest::OrderDecom(OrderDType::PH), _) => {
            vec!["F+3", "F-3"]
        }
        (MocvDest::RoundOrder, _) => {
            vec![
                "PX", // P轨道
                // "XZ", "XY", // 笛卡尔D轨道
                "D+1", "D-2", // 球谐D轨道
                // "XXY", "XXZ", "XYY", "XYZ", "XZZ", // 笛卡尔F轨道
                "F-2", "F+1", "F+3", // 球谐F轨道
            ] // 环轴键级
        }
        _ => vec![],
    };
    syms.into_iter().map(|sym| sym.to_string()).collect()
}
