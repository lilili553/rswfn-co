// libcint的绑定
// 这个模块用来计算电子积分
use crate::{
    base::{Ato, Mole},
    datas::{FACS, btrans::sym2lmn},
    utils,
};
use log::info;
use ndarray::{Array, Array2, Array4, ArrayBase, Dim, Order, OwnedRepr, ShapeBuilder, s};
use std::{
    env,
    f64::consts::PI,
    fmt,
    os::raw::{c_double, c_int, c_void},
    sync::OnceLock,
};
#[cfg(feature = "use-cint")]
unsafe extern "C" {
    fn CINTcgto_cart(bas_id: c_int, bas: *const c_int) -> c_int;
    fn cint1e_ovlp_cart(
        buf: *mut c_double,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
    ) -> c_int;
    fn cint1e_kin_cart(
        buf: *mut c_double,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
    ) -> c_int;
    fn cint1e_nuc_cart(
        buf: *mut c_double,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
    ) -> c_int;

    fn CINTcgto_spheric(bas_id: c_int, bas: *const c_int) -> c_int;
    fn cint1e_ovlp_sph(
        buf: *mut c_double,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
    ) -> c_int;
    fn cint1e_kin_sph(
        buf: *mut c_double,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
    ) -> c_int;
    fn cint1e_nuc_sph(
        buf: *mut c_double,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
    ) -> c_int;

    fn cint2e_cart(
        buf: *mut c_double,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
        opt: *mut c_void,
    ) -> c_int;
    fn cint2e_cart_optimizer(
        opt: *mut *mut c_void,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
    ) -> c_int;

    fn cint2e_sph(
        buf: *mut c_double,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
        opt: *mut c_void,
    ) -> c_int;
    fn cint2e_sph_optimizer(
        opt: *mut *mut c_void,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const c_double,
    ) -> c_int;

    fn CINTdel_optimizer(opt: *mut *mut c_void);

    fn CINTgto_norm(ang: c_int, alp: c_double) -> c_double;
}

/// libcint绑定的积分计算器
///
/// 封装libcint库的C接口,提供高效的电子积分计算。
#[cfg(feature = "use-cint")]
pub struct Cint<'m> {
    /// 分子结构引用
    pub mole: &'m Mole,
    /// libcint的atm数组
    pub atm: Vec<i32>,
    /// libcint的bas数组
    pub bas: Vec<i32>,
    /// libcint的env数组
    pub env: Vec<f64>,
    /// 原子轨道数量 [球谐, 笛卡尔]
    pub nato: [usize; 2],
    /// 校正系数 [球谐, 笛卡尔]
    pub fixs: [Vec<f64>; 2],
}

#[cfg(feature = "use-cint")]
impl<'m> Cint<'m> {
    /// 创建Cint实例
    pub fn new(mole: &'m Mole) -> Self {
        let natm = mole.atoms.len();
        let map_shl = mole.basis.map_shl();
        let shls = &map_shl.shls;
        let nato = [map_shl.nato("sph"), map_shl.nato("car")];
        let nshl = map_shl.len();
        let mut atm = vec![0; natm * 6];
        let mut bas = vec![0; nshl * 8];
        let mut env = vec![0.0; 20];
        let mut idx = 20;
        for i in 0..natm {
            let atom = mole.atoms.get(i);
            let nuc = atom.nuc as c_int;
            let [x, y, z] = atom.xyz;
            atm[i * 6 + 0] = nuc;
            atm[i * 6 + 1] = idx;
            println!("xyz:{:>10.4},{:>10.4},{:>10.4}", x, y, z);
            env.push(x);
            env.push(y);
            env.push(z);
            idx += 3;
        }
        // let coes_norm = map_shl.coes_norm();

        for i in 0..nshl {
            bas[i * 8 + 0] = shls[i].atm as c_int;
            bas[i * 8 + 1] = shls[i].ang as c_int;
            let alps = &shls[i].alps;
            let coes = &shls[i].coes;
            // let coes = &coes_norm[i];
            let num = alps.len() as c_int;
            bas[i * 8 + 2] = num; // 原始gto的数量
            bas[i * 8 + 3] = 1; // 收缩gto的数量
            bas[i * 8 + 4] = 0;
            env.extend(alps.clone());
            bas[i * 8 + 5] = idx; // 该层基函数在env中的起始位置
            idx += num;
            env.extend(coes.clone());
            bas[i * 8 + 6] = idx;
            idx += num;
        }
        // let env_raw = env.clone();
        let env = norm_basis(&atm, &bas, &env); // 归一化

        // 计算校正系数

        let fixs = {
            let mut vals_sph = vec![];
            let mut vals_car = vec![];
            let smat_sph = &calc_int1e_sph(nato[0], &atm, &bas, &env, [1, 0, 0])[0];
            let smat_car = &calc_int1e_car(nato[1], &atm, &bas, &env, [1, 0, 0])[0];
            for i in 0..smat_sph.nrows() {
                let val = smat_sph[(i, i)];
                vals_sph.push(if val == 0.0 { 0.0 } else { 1.0 / val.sqrt() });
            }
            for i in 0..smat_car.nrows() {
                let val = smat_car[(i, i)];
                vals_car.push(if val == 0.0 { 0.0 } else { 1.0 / val.sqrt() });
            }
            [vals_sph, vals_car]
        };

        Self {
            mole,
            atm,
            bas,
            env,
            nato,
            fixs,
        }
    }

    /// 获取校正系数
    pub fn get_fixs(&self, form: &str) -> &Vec<f64> {
        match form {
            "sph" => &self.fixs[0],
            "car" => &self.fixs[1],
            _ => panic!("Invalid form"),
        }
    }

    /// 对单电子积分矩阵应用校正系数
    pub fn fix_int1e(&self, int1e: &mut Array2<f64>, form: &str) {
        let vals = self.get_fixs(form);
        let nato = vals.len();
        for i in 0..nato {
            for j in 0..nato {
                int1e[(i, j)] *= vals[i] * vals[j];
            }
        }
    }

    /// 对双电子积分张量应用校正系数
    pub fn fix_int2e(&self, int2e: &mut Array4<f64>, form: &str) {
        let vals = self.get_fixs(form);
        let nato = vals.len();
        for i in 0..nato {
            for j in 0..nato {
                for k in 0..nato {
                    for l in 0..nato {
                        int2e[(i, j, k, l)] *= vals[i] * vals[j] * vals[k] * vals[l];
                    }
                }
            }
        }
    }

    /// 计算单电子积分
    pub fn calc_int1e(&self, form: &str, cals: [usize; 3]) -> [Array2<f64>; 3] {
        let [mut smat, mut tmat, mut vmat] = match form {
            "sph" => calc_int1e_sph(self.nato[0], &self.atm, &self.bas, &self.env, cals),
            "car" => calc_int1e_car(self.nato[1], &self.atm, &self.bas, &self.env, cals),
            _ => panic!("Invalid form"),
        };
        // 系数矫正&顺序矫正
        if cals[0] == 1 {
            self.fix_int1e(&mut smat, form);
            smat = self.ord_int1e(&smat, form);
        }
        if cals[1] == 1 {
            self.fix_int1e(&mut tmat, form);
            tmat = self.ord_int1e(&tmat, form);
        }
        if cals[2] == 1 {
            self.fix_int1e(&mut vmat, form);
            vmat = self.ord_int1e(&vmat, form);
        }
        [smat, tmat, vmat]
    }
    /// 计算双电子积分
    pub fn calc_int2e(&self, form: &str) -> Array4<f64> {
        let mut mat2e = match form {
            "sph" => calc_int2e_sph(self.nato[0], &self.atm, &self.bas, &self.env),
            "car" => calc_int2e_car(self.nato[1], &self.atm, &self.bas, &self.env),
            _ => panic!("Invalid form"),
        };
        log::info!("双电子积分计算完成");

        self.fix_int2e(&mut mat2e, form);
        log::info!("双电子积分系数矫正完成");

        let mat2e = self.ord_int2e(mat2e, form);
        log::info!("双电子积分顺序矫正完成");

        mat2e
    }

    /// 将单电子积分的原子轨道顺序转为和高斯一致
    pub fn ord_int1e(&self, mat1e_raw: &Array2<f64>, form: &str) -> Array2<f64> {
        let mut mat1e_fix = mat1e_raw.clone();
        let uplo = self.mole.basis.map_shl().uplo(form);
        let nshl = uplo.len();
        for si in 0..nshl {
            let [ui, li] = uplo[si];
            let leni = li - ui;
            let idxi = len2idxs(leni);
            for sj in 0..nshl {
                let [uj, lj] = uplo[sj];
                let lenj = lj - uj;
                let idxj = len2idxs(lenj);
                if leni <= 3 && lenj <= 3 {
                    continue;
                }
                for i in 0..leni {
                    for j in 0..lenj {
                        mat1e_fix[(ui + i, uj + j)] = mat1e_raw[(ui + idxi[i], uj + idxj[j])]
                    }
                }
            }
        }
        mat1e_fix
    }

    /// 将双电子积分的原子轨道顺序转为和高斯一致
    pub fn ord_int2e(&self, mat2e_raw: Array4<f64>, form: &str) -> Array4<f64> {
        let mut mat2e_fix = mat2e_raw.clone();
        let uplo = self.mole.basis.map_shl().uplo(form);
        let nshl = uplo.len();
        for si in 0..nshl {
            let [ui, li] = uplo[si];
            let leni = li - ui;
            let idxi = len2idxs(leni);
            for sj in 0..nshl {
                let [uj, lj] = uplo[sj];
                let lenj = lj - uj;
                let idxj = len2idxs(lenj);
                for sk in 0..nshl {
                    let [uk, lk] = uplo[sk];
                    let lenk = lk - uk;
                    let idxk = len2idxs(lenk);
                    for sl in 0..nshl {
                        let [ul, ll] = uplo[sl];
                        let lenl = ll - ul;
                        let idxl = len2idxs(lenl);
                        for i in 0..leni {
                            for j in 0..lenj {
                                for k in 0..lenk {
                                    for l in 0..lenl {
                                        mat2e_fix[(ui + i, uj + j, uk + k, ul + l)] =
                                            mat2e_raw[(ui + idxi[i], uj + idxj[j], uk + idxk[k], ul + idxl[l])];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        mat2e_fix
    }
}

/// 归一化基函数
///
/// 对每个基函数应用CINTgto_norm归一化因子。
#[cfg(feature = "use-cint")]
pub fn norm_basis(atm: &Vec<i32>, bas: &Vec<i32>, env: &Vec<f64>) -> Vec<f64> {
    let mut env_new = env.clone();
    let nbas: usize = bas.len() / 8;
    for i in 0..nbas {
        let idx = i * 8;
        let ang = bas[idx + 1];
        let ialp = bas[idx + 5] as usize;
        let icoe = bas[idx + 6] as usize;
        let tot = icoe - ialp;
        for j in 0..tot {
            let alp = env[ialp + j];
            let coe = env[icoe + j];
            unsafe {
                env_new[icoe + j] = coe * CINTgto_norm(ang as c_int, alp as c_double);
            }
        }
    }
    env_new
}

/// 计算球谐单电子积分
///
/// 调用libcint的球谐积分函数,计算重叠、动能和核引积分。
#[cfg(feature = "use-cint")]
pub fn calc_int1e_sph(
    nato: usize,
    atm: &Vec<i32>,
    bas: &Vec<i32>,
    env: &Vec<f64>,
    cals: [usize; 3],
) -> [Array2<f64>; 3] {
    let nbas = bas.len() / 8;
    let natm = atm.len() / 6;
    let mut smat = Array2::<f64>::zeros((nato, nato).f()); // 重叠矩阵
    let mut tmat = Array2::<f64>::zeros((nato, nato).f()); // 动能矩阵
    let mut vmat = Array2::<f64>::zeros((nato, nato).f()); // 势能矩阵

    let mut locs = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
    let nums: [usize; 7] = [1, 3, 5, 7, 9, 11, 13]; // 每个角动量对应的数量
    let mut loc = 0usize;
    for i in 0..nbas {
        locs[i] = loc;
        let ang = bas[i * 8 + 1] as usize;
        loc += nums[ang];
    }
    // 获取第一个基函数的角动量对应的维度
    let mut dimi: usize;
    let mut dimj: usize;
    let mut loci: usize;
    let mut locj: usize;
    unsafe {
        let mut dims = vec![0 as usize; nbas]; // 每个基函数的原子轨道数量
        for i in 0..nbas {
            dims[i] = CINTcgto_spheric(i as c_int, bas.as_ptr()) as usize;
        }
        println!("dims_1e_sph:{:?}", dims);
        for i in 0..nbas {
            loci = locs[i];
            dimi = dims[i];
            assert!(loci + dimi <= nato, "数组越界,{i},{loci},{dimi},{nato}");
            for j in 0..nbas {
                locj = locs[j];
                dimj = dims[j];
                assert!(locj + dimj <= nato, "数组越界,{j},{locj},{dimj},{nato}");
                let mut buf = vec![0.0f64; dimi * dimj];
                let shls = [i as c_int, j as c_int];
                let shape = (dimi, dimj);
                if cals[0] == 1 {
                    // 计算重叠积分
                    cint1e_ovlp_sph(
                        buf.as_mut_ptr(),
                        shls.as_ptr(),
                        atm.as_ptr(),
                        natm as c_int,
                        bas.as_ptr(),
                        nbas as c_int,
                        env.as_ptr(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    smat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
                if cals[1] == 1 {
                    // 计算动能积分
                    cint1e_kin_sph(
                        buf.as_mut_ptr(),
                        shls.as_ptr(),
                        atm.as_ptr(),
                        natm as c_int,
                        bas.as_ptr(),
                        nbas as c_int,
                        env.as_ptr(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    tmat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
                if cals[2] == 1 {
                    // 计算核引积分
                    cint1e_nuc_sph(
                        buf.as_mut_ptr(),
                        shls.as_ptr(),
                        atm.as_ptr(),
                        natm as c_int,
                        bas.as_ptr(),
                        nbas as c_int,
                        env.as_ptr(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    vmat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
            }
        }
    }
    return [smat, tmat, vmat];
}

/// 计算笛卡尔单电子积分
///
/// 调用libcint的笛卡尔积分函数,计算重叠、动能和核引积分。
#[cfg(feature = "use-cint")]
pub fn calc_int1e_car(
    nato: usize,
    atm: &Vec<i32>,
    bas: &Vec<i32>,
    env: &Vec<f64>,
    cals: [usize; 3],
) -> [Array2<f64>; 3] {
    let nbas = bas.len() / 8;
    let natm = atm.len() / 6;
    let mut smat = Array2::<f64>::zeros((nato, nato).f()); // 重叠矩阵
    let mut tmat = Array2::<f64>::zeros((nato, nato).f()); // 动能矩阵
    let mut vmat = Array2::<f64>::zeros((nato, nato).f()); // 势能矩阵

    let mut locs = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
    let nums: [usize; 7] = [1, 3, 6, 10, 15, 21, 28]; // 每个角动量对应的数量
    let mut loc = 0usize;
    for i in 0..nbas {
        locs[i] = loc;
        let ang = bas[i * 8 + 1] as usize;
        loc += nums[ang];
    }
    // 获取第一个基函数的角动量对应的维度
    let mut dimi: usize;
    let mut dimj: usize;
    let mut loci: usize;
    let mut locj: usize;
    unsafe {
        let mut dims = vec![0 as usize; nbas]; // 每个基函数的原子轨道数量
        for i in 0..nbas {
            dims[i] = CINTcgto_cart(i as c_int, bas.as_ptr()) as usize;
        }
        println!("dims_1e_car:{:?}", dims);
        for i in 0..nbas {
            loci = locs[i];
            dimi = dims[i];
            assert!(loci + dimi <= nato, "数组越界,{i},{loci},{dimi},{nato}");
            for j in 0..nbas {
                locj = locs[j];
                dimj = dims[j];
                assert!(locj + dimj <= nato, "数组越界,{j},{locj},{dimj},{nato}");
                let mut buf = vec![0.0f64; dimi * dimj];
                let shls = [i as c_int, j as c_int];
                let shape = (dimi, dimj);
                if cals[0] == 1 {
                    // println!("计算重叠积分");
                    cint1e_ovlp_cart(
                        buf.as_mut_ptr(),
                        shls.as_ptr(),
                        atm.as_ptr(),
                        natm as c_int,
                        bas.as_ptr(),
                        nbas as c_int,
                        env.as_ptr(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    smat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
                if cals[1] == 1 {
                    // println!("计算动能积分");
                    cint1e_kin_cart(
                        buf.as_mut_ptr(),
                        shls.as_ptr(),
                        atm.as_ptr(),
                        natm as c_int,
                        bas.as_ptr(),
                        nbas as c_int,
                        env.as_ptr(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    tmat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
                if cals[2] == 1 {
                    // println!("计算核引积分");
                    cint1e_nuc_cart(
                        buf.as_mut_ptr(),
                        shls.as_ptr(),
                        atm.as_ptr(),
                        natm as c_int,
                        bas.as_ptr(),
                        nbas as c_int,
                        env.as_ptr(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    vmat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
            }
        }
    }
    return [smat, tmat, vmat];
}

/// 计算球谐双电子积分
///
/// 调用libcint的球谐双电子积分函数,使用优化器加速计算。
#[cfg(feature = "use-cint")]
pub fn calc_int2e_sph(nmat: usize, atm: &Vec<i32>, bas: &Vec<i32>, env: &Vec<f64>) -> Array4<f64> {
    let natm = atm.len() / 6; // 原子数量
    let nbas = bas.len() / 8; // 基函数数量

    let mut locs = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
    let nums: [usize; 7] = [1, 3, 5, 7, 9, 11, 13]; // 每个角动量对应的数量
    let mut loc = 0usize;
    for i in 0..nbas {
        locs[i] = loc;
        let ang = bas[i * 8 + 1] as usize;
        loc += nums[ang];
    }

    let mut dimi;
    let mut dimj;
    let mut dimk;
    let mut diml;
    let mut loci;
    let mut locj;
    let mut lock;
    let mut locl;

    let mut emat = Array4::<f64>::zeros((nmat, nmat, nmat, nmat).f()); // 双电子积分矩阵
    let mut opt: *mut c_void = std::ptr::null_mut();
    unsafe {
        // println!("开始双电子积分优化");
        cint2e_sph_optimizer(
            &mut opt,
            atm.as_ptr(),
            natm as c_int,
            bas.as_ptr(),
            nbas as c_int,
            env.as_ptr(),
        );

        let mut dims = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
        for i in 0..nbas {
            dims[i] = CINTcgto_spheric(i as c_int, bas.as_ptr()) as usize;
        }
        println!("dims_2e_sph:{:?}", dims);
        // println!("完成双电子积分优化");
        for i in 0..nbas {
            loci = locs[i];
            dimi = dims[i];
            assert!(loci + dimi <= nmat, "数组越界");
            for j in 0..nbas {
                locj = locs[j];
                dimj = dims[j];
                assert!(locj + dimj <= nmat, "数组越界");
                for k in 0..nbas {
                    lock = locs[k];
                    dimk = dims[k];
                    assert!(lock + dimk <= nmat, "数组越界");
                    for l in 0..nbas {
                        locl = locs[l];
                        diml = dims[l];
                        assert!(locl + diml <= nmat, "数组越界");
                        let mut buf = vec![0.0f64; dimi * dimj * dimk * diml];
                        let shls = [i as c_int, j as c_int, k as c_int, l as c_int];
                        // println!("双电子积分:{i}{j}{k}{l}");
                        cint2e_sph(
                            buf.as_mut_ptr(),
                            shls.as_ptr(),
                            atm.as_ptr(),
                            natm as c_int,
                            bas.as_ptr(),
                            nbas as c_int,
                            env.as_ptr(),
                            opt,
                        );
                        let shape = (dimi, dimj, dimk, diml);
                        let buf_array = Array::from_shape_vec(shape.f(), buf).unwrap();
                        emat.slice_mut(s![
                            loci..loci + dimi,
                            locj..locj + dimj,
                            lock..lock + dimk,
                            locl..locl + diml
                        ])
                        .assign(&buf_array);
                    }
                }
            }
        }
        CINTdel_optimizer(&mut opt);
    }
    println!("双电子积分计算完成");
    return emat;
}

/// 计算笛卡尔双电子积分
///
/// 调用libcint的笛卡尔双电子积分函数,使用优化器加速计算。
#[cfg(feature = "use-cint")]
pub fn calc_int2e_car(nmat: usize, atm: &Vec<i32>, bas: &Vec<i32>, env: &Vec<f64>) -> Array4<f64> {
    let natm = atm.len() / 6; // 原子数量
    let nbas = bas.len() / 8; // 基函数数量

    let mut locs = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
    let nums: [usize; 7] = [1, 3, 6, 10, 15, 21, 28]; // 每个角动量对应的数量
    let mut loc = 0usize;
    for i in 0..nbas {
        locs[i] = loc;
        let ang = bas[i * 8 + 1] as usize;
        loc += nums[ang];
    }

    let mut dimi;
    let mut dimj;
    let mut dimk;
    let mut diml;
    let mut loci;
    let mut locj;
    let mut lock;
    let mut locl;

    let mut emat = Array4::<f64>::zeros((nmat, nmat, nmat, nmat).f()); // 双电子积分矩阵
    let mut opt: *mut c_void = std::ptr::null_mut();
    unsafe {
        // println!("开始双电子积分优化");
        cint2e_cart_optimizer(
            &mut opt,
            atm.as_ptr(),
            natm as c_int,
            bas.as_ptr(),
            nbas as c_int,
            env.as_ptr(),
        );

        let mut dims = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
        for i in 0..nbas {
            dims[i] = CINTcgto_cart(i as c_int, bas.as_ptr()) as usize;
        }
        println!("dims_2e_car:{:?}", dims);
        // println!("完成双电子积分优化");
        for i in 0..nbas {
            loci = locs[i];
            dimi = dims[i];
            assert!(loci + dimi <= nmat, "数组越界");
            for j in 0..nbas {
                locj = locs[j];
                dimj = dims[j];
                assert!(locj + dimj <= nmat, "数组越界");
                for k in 0..nbas {
                    lock = locs[k];
                    dimk = dims[k];
                    assert!(lock + dimk <= nmat, "数组越界");
                    for l in 0..nbas {
                        locl = locs[l];
                        diml = dims[l];
                        assert!(locl + diml <= nmat, "数组越界");
                        let mut buf = vec![0.0f64; dimi * dimj * dimk * diml];
                        let shls = [i as c_int, j as c_int, k as c_int, l as c_int];
                        cint2e_cart(
                            buf.as_mut_ptr(),
                            shls.as_ptr(),
                            atm.as_ptr(),
                            natm as c_int,
                            bas.as_ptr(),
                            nbas as c_int,
                            env.as_ptr(),
                            opt,
                        );
                        let shape = (dimi, dimj, dimk, diml);
                        let buf_array = Array::from_shape_vec(shape.f(), buf).unwrap();
                        emat.slice_mut(s![
                            loci..loci + dimi,
                            locj..locj + dimj,
                            lock..lock + dimk,
                            locl..locl + diml
                        ])
                        .assign(&buf_array);
                    }
                }
            }
        }
        CINTdel_optimizer(&mut opt);
    }
    println!("双电子积分计算完成");
    return emat;
}

#[cfg(feature = "use-cint")]
impl<'m> fmt::Display for Cint<'m> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let env = &self.env;
        writeln!(f, "atm")?;
        let natm = self.atm.len() / 6;
        for i in 0..natm {
            let x = env[i * 3 + 0];
            let y = env[i * 3 + 1];
            let z = env[i * 3 + 2];
            write!(f, "{:>10.4},{:>10.4},{:>10.4} | ", x, y, z)?;
            writeln!(f, "{:?}", &self.atm[i * 6..i * 6 + 6])?;
        }
        let nbas = self.bas.len() / 8;
        writeln!(f, "bas")?;
        for i in 0..nbas {
            writeln!(f, "{:?}", &self.bas[i * 8..i * 8 + 8])?;
        }
        writeln!(f, "env")?;
        for i in 0..nbas {
            let idx = i * 8;
            let ialp = self.bas[idx + 5] as usize;
            let icoe = self.bas[idx + 6] as usize;
            let tot = icoe - ialp;
            write!(f, "alps:")?;
            for j in 0..tot {
                write!(f, "{:>16.10E}", env[ialp + j])?;
            }
            write!(f, "\n")?;
            write!(f, "coes:")?;
            for j in 0..tot {
                write!(f, "{:>16.10E}", env[icoe + j])?;
            }
            write!(f, "\n")?;
        }
        Ok(())
    }
}

fn len2idxs(len: usize) -> Vec<usize> {
    match len {
        1 => vec![0],
        3 => vec![0, 1, 2],
        5 => vec![2, 3, 1, 4, 0],
        // libcint d顺序[xx,xy,xz,yy,yz,zz] → 标准[xx,yy,zz,xy,xz,yz]
        6 => vec![0, 3, 5, 1, 2, 4],
        7 => vec![3, 4, 2, 5, 1, 6, 0],
        _ => panic!("未知长度:{}", len),
    }
}
