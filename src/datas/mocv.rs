#![allow(non_snake_case)]
use crate::datas::FACS;
use crate::datas::btrans;
/// 局部坐标系与全局坐标系的变换矩阵 Cg=Mt@Cl
use ndarray::Array2;

pub fn get_tmat_sph(num: usize, lmat: &Array2<f64>) -> Array2<f64> {
    match num {
        01 => tmat_s(lmat),
        03 => tmat_p(lmat),
        05 => btrans::DMAT.t().dot(&get_tmat_car(06, lmat)).dot(&*btrans::DMAT), // [5,6] [6,6] [6,5]
        07 => btrans::FMAT.t().dot(&get_tmat_car(10, lmat)).dot(&*btrans::FMAT),
        09 => btrans::GMAT.t().dot(&get_tmat_car(15, lmat)).dot(&*btrans::GMAT),
        _ => panic!("不支持的类型:{}", num),
    }
}

/// 获取变换矩阵,输入的矩阵是列为基向量，要转为行为基向量
pub fn get_tmat_car(num: usize, lmat: &Array2<f64>) -> Array2<f64> {
    match num {
        01 => tmat_s(lmat),
        03 => tmat_p(lmat),
        06 => tmat_d(lmat) * nmat_d(),
        10 => tmat_f(lmat) * nmat_f(),
        15 => tmat_g(lmat) * nmat_g(),
        _ => panic!("不支持的类型:{}", num),
    }
}

fn get_ts(L: &Array2<f64>) -> [f64; 9] {
    let L11 = L[(0, 0)];
    let L12 = L[(0, 1)];
    let L13 = L[(0, 2)];
    let L21 = L[(1, 0)];
    let L22 = L[(1, 1)];
    let L23 = L[(1, 2)];
    let L31 = L[(2, 0)];
    let L32 = L[(2, 1)];
    let L33 = L[(2, 2)];
    return [L11, L12, L13, L21, L22, L23, L31, L32, L33];
}

pub fn tmat_s(L: &Array2<f64>) -> Array2<f64> {
    Array2::ones((1, 1))
}

pub fn tmat_p(L: &Array2<f64>) -> Array2<f64> {
    let mut MP: Array2<f64> = Array2::zeros((3, 3));
    let [L11, L12, L13, L21, L22, L23, L31, L32, L33] = get_ts(L);
    MP[(0, 0)] = L11;
    MP[(0, 1)] = L21;
    MP[(0, 2)] = L31;
    MP[(1, 0)] = L12;
    MP[(1, 1)] = L22;
    MP[(1, 2)] = L32;
    MP[(2, 0)] = L13;
    MP[(2, 1)] = L23;
    MP[(2, 2)] = L33;
    MP
}

pub fn tmat_d(L: &Array2<f64>) -> Array2<f64> {
    let [L11, L12, L13, L21, L22, L23, L31, L32, L33] = get_ts(L);
    let mut MD: Array2<f64> = Array2::zeros((6, 6));
    MD[(0, 0)] = L11.powi(2);
    MD[(0, 1)] = L21.powi(2);
    MD[(0, 2)] = L31.powi(2);
    MD[(0, 3)] = 2.0 * L11 * L21;
    MD[(0, 4)] = 2.0 * L11 * L31;
    MD[(0, 5)] = 2.0 * L21 * L31;
    MD[(1, 0)] = L12.powi(2);
    MD[(1, 1)] = L22.powi(2);
    MD[(1, 2)] = L32.powi(2);
    MD[(1, 3)] = 2.0 * L12 * L22;
    MD[(1, 4)] = 2.0 * L12 * L32;
    MD[(1, 5)] = 2.0 * L22 * L32;
    MD[(2, 0)] = L13.powi(2);
    MD[(2, 1)] = L23.powi(2);
    MD[(2, 2)] = L33.powi(2);
    MD[(2, 3)] = 2.0 * L13 * L23;
    MD[(2, 4)] = 2.0 * L13 * L33;
    MD[(2, 5)] = 2.0 * L23 * L33;
    MD[(3, 0)] = L11 * L12;
    MD[(3, 1)] = L21 * L22;
    MD[(3, 2)] = L31 * L32;
    MD[(3, 3)] = L11 * L22 + L12 * L21;
    MD[(3, 4)] = L11 * L32 + L12 * L31;
    MD[(3, 5)] = L21 * L32 + L22 * L31;
    MD[(4, 0)] = L11 * L13;
    MD[(4, 1)] = L21 * L23;
    MD[(4, 2)] = L31 * L33;
    MD[(4, 3)] = L11 * L23 + L13 * L21;
    MD[(4, 4)] = L11 * L33 + L13 * L31;
    MD[(4, 5)] = L21 * L33 + L23 * L31;
    MD[(5, 0)] = L12 * L13;
    MD[(5, 1)] = L22 * L23;
    MD[(5, 2)] = L32 * L33;
    MD[(5, 3)] = L12 * L23 + L13 * L22;
    MD[(5, 4)] = L12 * L33 + L13 * L32;
    MD[(5, 5)] = L22 * L33 + L23 * L32;
    MD
}

// #[rustfmt::skip]
pub fn tmat_f(L: &Array2<f64>) -> Array2<f64> {
    let [L11, L12, L13, L21, L22, L23, L31, L32, L33] = get_ts(L);
    let mut MF: Array2<f64> = Array2::zeros((10, 10));
    MF[(0, 0)] = L11.powi(3);
    MF[(0, 1)] = L21.powi(3);
    MF[(0, 2)] = L31.powi(3);
    MF[(0, 3)] = 3.0 * L11 * L21.powi(2);
    MF[(0, 4)] = 3.0 * L11.powi(2) * L21;
    MF[(0, 5)] = 3.0 * L11.powi(2) * L31;
    MF[(0, 6)] = 3.0 * L11 * L31.powi(2);
    MF[(0, 7)] = 3.0 * L21 * L31.powi(2);
    MF[(0, 8)] = 3.0 * L21.powi(2) * L31;
    MF[(0, 9)] = 6.0 * L11 * L21 * L31;
    MF[(1, 0)] = L12.powi(3);
    MF[(1, 1)] = L22.powi(3);
    MF[(1, 2)] = L32.powi(3);
    MF[(1, 3)] = 3.0 * L12 * L22.powi(2);
    MF[(1, 4)] = 3.0 * L12.powi(2) * L22;
    MF[(1, 5)] = 3.0 * L12.powi(2) * L32;
    MF[(1, 6)] = 3.0 * L12 * L32.powi(2);
    MF[(1, 7)] = 3.0 * L22 * L32.powi(2);
    MF[(1, 8)] = 3.0 * L22.powi(2) * L32;
    MF[(1, 9)] = 6.0 * L12 * L22 * L32;
    MF[(2, 0)] = L13.powi(3);
    MF[(2, 1)] = L23.powi(3);
    MF[(2, 2)] = L33.powi(3);
    MF[(2, 3)] = 3.0 * L13 * L23.powi(2);
    MF[(2, 4)] = 3.0 * L13.powi(2) * L23;
    MF[(2, 5)] = 3.0 * L13.powi(2) * L33;
    MF[(2, 6)] = 3.0 * L13 * L33.powi(2);
    MF[(2, 7)] = 3.0 * L23 * L33.powi(2);
    MF[(2, 8)] = 3.0 * L23.powi(2) * L33;
    MF[(2, 9)] = 6.0 * L13 * L23 * L33;
    MF[(3, 0)] = L11 * L12.powi(2);
    MF[(3, 1)] = L21 * L22.powi(2);
    MF[(3, 2)] = L31 * L32.powi(2);
    MF[(3, 3)] = L11 * L22.powi(2) + 2.0 * L12 * L21 * L22;
    MF[(3, 4)] = 2.0 * L11 * L12 * L22 + L12.powi(2) * L21;
    MF[(3, 5)] = 2.0 * L11 * L12 * L32 + L12.powi(2) * L31;
    MF[(3, 6)] = L11 * L32.powi(2) + 2.0 * L12 * L31 * L32;
    MF[(3, 7)] = L21 * L32.powi(2) + 2.0 * L22 * L31 * L32;
    MF[(3, 8)] = 2.0 * L21 * L22 * L32 + L22.powi(2) * L31;
    MF[(3, 9)] = 2.0 * L11 * L22 * L32 + 2.0 * L12 * L21 * L32 + 2.0 * L12 * L22 * L31;
    MF[(4, 0)] = L11.powi(2) * L12;
    MF[(4, 1)] = L21.powi(2) * L22;
    MF[(4, 2)] = L31.powi(2) * L32;
    MF[(4, 3)] = 2.0 * L11 * L21 * L22 + L12 * L21.powi(2);
    MF[(4, 4)] = L11.powi(2) * L22 + 2.0 * L11 * L12 * L21;
    MF[(4, 5)] = L11.powi(2) * L32 + 2.0 * L11 * L12 * L31;
    MF[(4, 6)] = 2.0 * L11 * L31 * L32 + L12 * L31.powi(2);
    MF[(4, 7)] = 2.0 * L21 * L31 * L32 + L22 * L31.powi(2);
    MF[(4, 8)] = L21.powi(2) * L32 + 2.0 * L21 * L22 * L31;
    MF[(4, 9)] = 2.0 * L11 * L21 * L32 + 2.0 * L11 * L22 * L31 + 2.0 * L12 * L21 * L31;
    MF[(5, 0)] = L11.powi(2) * L13;
    MF[(5, 1)] = L21.powi(2) * L23;
    MF[(5, 2)] = L31.powi(2) * L33;
    MF[(5, 3)] = 2.0 * L11 * L21 * L23 + L13 * L21.powi(2);
    MF[(5, 4)] = L11.powi(2) * L23 + 2.0 * L11 * L13 * L21;
    MF[(5, 5)] = L11.powi(2) * L33 + 2.0 * L11 * L13 * L31;
    MF[(5, 6)] = 2.0 * L11 * L31 * L33 + L13 * L31.powi(2);
    MF[(5, 7)] = 2.0 * L21 * L31 * L33 + L23 * L31.powi(2);
    MF[(5, 8)] = L21.powi(2) * L33 + 2.0 * L21 * L23 * L31;
    MF[(5, 9)] = 2.0 * L11 * L21 * L33 + 2.0 * L11 * L23 * L31 + 2.0 * L13 * L21 * L31;
    MF[(6, 0)] = L11 * L13.powi(2);
    MF[(6, 1)] = L21 * L23.powi(2);
    MF[(6, 2)] = L31 * L33.powi(2);
    MF[(6, 3)] = L11 * L23.powi(2) + 2.0 * L13 * L21 * L23;
    MF[(6, 4)] = 2.0 * L11 * L13 * L23 + L13.powi(2) * L21;
    MF[(6, 5)] = 2.0 * L11 * L13 * L33 + L13.powi(2) * L31;
    MF[(6, 6)] = L11 * L33.powi(2) + 2.0 * L13 * L31 * L33;
    MF[(6, 7)] = L21 * L33.powi(2) + 2.0 * L23 * L31 * L33;
    MF[(6, 8)] = 2.0 * L21 * L23 * L33 + L23.powi(2) * L31;
    MF[(6, 9)] = 2.0 * L11 * L23 * L33 + 2.0 * L13 * L21 * L33 + 2.0 * L13 * L23 * L31;
    MF[(7, 0)] = L12 * L13.powi(2);
    MF[(7, 1)] = L22 * L23.powi(2);
    MF[(7, 2)] = L32 * L33.powi(2);
    MF[(7, 3)] = L12 * L23.powi(2) + 2.0 * L13 * L22 * L23;
    MF[(7, 4)] = 2.0 * L12 * L13 * L23 + L13.powi(2) * L22;
    MF[(7, 5)] = 2.0 * L12 * L13 * L33 + L13.powi(2) * L32;
    MF[(7, 6)] = L12 * L33.powi(2) + 2.0 * L13 * L32 * L33;
    MF[(7, 7)] = L22 * L33.powi(2) + 2.0 * L23 * L32 * L33;
    MF[(7, 8)] = 2.0 * L22 * L23 * L33 + L23.powi(2) * L32;
    MF[(7, 9)] = 2.0 * L12 * L23 * L33 + 2.0 * L13 * L22 * L33 + 2.0 * L13 * L23 * L32;
    MF[(8, 0)] = L12.powi(2) * L13;
    MF[(8, 1)] = L22.powi(2) * L23;
    MF[(8, 2)] = L32.powi(2) * L33;
    MF[(8, 3)] = 2.0 * L12 * L22 * L23 + L13 * L22.powi(2);
    MF[(8, 4)] = L12.powi(2) * L23 + 2.0 * L12 * L13 * L22;
    MF[(8, 5)] = L12.powi(2) * L33 + 2.0 * L12 * L13 * L32;
    MF[(8, 6)] = 2.0 * L12 * L32 * L33 + L13 * L32.powi(2);
    MF[(8, 7)] = 2.0 * L22 * L32 * L33 + L23 * L32.powi(2);
    MF[(8, 8)] = L22.powi(2) * L33 + 2.0 * L22 * L23 * L32;
    MF[(8, 9)] = 2.0 * L12 * L22 * L33 + 2.0 * L12 * L23 * L32 + 2.0 * L13 * L22 * L32;
    MF[(9, 0)] = L11 * L12 * L13;
    MF[(9, 1)] = L21 * L22 * L23;
    MF[(9, 2)] = L31 * L32 * L33;
    MF[(9, 3)] = L11 * L22 * L23 + L12 * L21 * L23 + L13 * L21 * L22;
    MF[(9, 4)] = L11 * L12 * L23 + L11 * L13 * L22 + L12 * L13 * L21;
    MF[(9, 5)] = L11 * L12 * L33 + L11 * L13 * L32 + L12 * L13 * L31;
    MF[(9, 6)] = L11 * L32 * L33 + L12 * L31 * L33 + L13 * L31 * L32;
    MF[(9, 7)] = L21 * L32 * L33 + L22 * L31 * L33 + L23 * L31 * L32;
    MF[(9, 8)] = L21 * L22 * L33 + L21 * L23 * L32 + L22 * L23 * L31;
    MF[(9, 9)] = L11 * L22 * L33 + L11 * L23 * L32 + L12 * L21 * L33 + L12 * L23 * L31 + L13 * L21 * L32 + L13 * L22 * L31;
    MF
}

#[rustfmt::skip]
pub fn tmat_g(L: &Array2<f64>) -> Array2<f64> {
    let [L11, L12, L13, L21, L22, L23, L31, L32, L33] = get_ts(L);
    let mut MG: Array2<f64> = Array2::zeros((10, 10));
    MG[(0,0)]=L33.powi(4);
    MG[(0,1)]=4.0*L23*L33.powi(3);
    MG[(0,2)]=6.0*L23.powi(2)*L33.powi(2);
    MG[(0,3)]=4.0*L23.powi(3)*L33;
    MG[(0,4)]=L23.powi(4);
    MG[(0,5)]=4.0*L13*L33.powi(3);
    MG[(0,6)]=12.0*L13*L23*L33.powi(2);
    MG[(0,7)]=12.0*L13*L23.powi(2)*L33;
    MG[(0,8)]=4.0*L13*L23.powi(3);
    MG[(0,9)]=6.0*L13.powi(2)*L33.powi(2);
    MG[(0,10)]=12.0*L13.powi(2)*L23*L33;
    MG[(0,11)]=6.0*L13.powi(2)*L23.powi(2);
    MG[(0,12)]=4.0*L13.powi(3)*L33;
    MG[(0,13)]=4.0*L13.powi(3)*L23;
    MG[(0,14)]=L13.powi(4);
    MG[(1,0)]=L32*L33.powi(3);
    MG[(1,1)]=L22*L33.powi(3) + 3.0*L23*L32*L33.powi(2);
    MG[(1,2)]=3.0*L22*L23*L33.powi(2) + 3.0*L23.powi(2)*L32*L33;
    MG[(1,3)]=3.0*L22*L23.powi(2)*L33 + L23.powi(3)*L32;
    MG[(1,4)]=L22*L23.powi(3);
    MG[(1,5)]=L12*L33.powi(3) + 3.0*L13*L32*L33.powi(2);
    MG[(1,6)]=3.0*L12*L23*L33.powi(2) + 3.0*L13*L22*L33.powi(2) + 6.0*L13*L23*L32*L33;
    MG[(1,7)]=3.0*L12*L23.powi(2)*L33 + 6.0*L13*L22*L23*L33 + 3.0*L13*L23.powi(2)*L32;
    MG[(1,8)]=L12*L23.powi(3) + 3.0*L13*L22*L23.powi(2);
    MG[(1,9)]=3.0*L12*L13*L33.powi(2) + 3.0*L13.powi(2)*L32*L33;
    MG[(1,10)]=6.0*L12*L13*L23*L33 + 3.0*L13.powi(2)*L22*L33 + 3.0*L13.powi(2)*L23*L32;
    MG[(1,11)]=3.0*L12*L13*L23.powi(2) + 3.0*L13.powi(2)*L22*L23;
    MG[(1,12)]=3.0*L12*L13.powi(2)*L33 + L13.powi(3)*L32;
    MG[(1,13)]=3.0*L12*L13.powi(2)*L23 + L13.powi(3)*L22;
    MG[(1,14)]=L12*L13.powi(3);
    MG[(2,0)]=L32.powi(2)*L33.powi(2);
    MG[(2,1)]=2.0*L22*L32*L33.powi(2) + 2.0*L23*L32.powi(2)*L33;
    MG[(2,2)]=L22.powi(2)*L33.powi(2) + 4.0*L22*L23*L32*L33 + L23.powi(2)*L32.powi(2);
    MG[(2,3)]=2.0*L22.powi(2)*L23*L33 + 2.0*L22*L23.powi(2)*L32;
    MG[(2,4)]=L22.powi(2)*L23.powi(2);
    MG[(2,5)]=2.0*L12*L32*L33.powi(2) + 2.0*L13*L32.powi(2)*L33;
    MG[(2,6)]=2.0*L12*L22*L33.powi(2) + 4.0*L12*L23*L32*L33 + 4.0*L13*L22*L32*L33 + 2.0*L13*L23*L32.powi(2);
    MG[(2,7)]=4.0*L12*L22*L23*L33 + 2.0*L12*L23.powi(2)*L32 + 2.0*L13*L22.powi(2)*L33 + 4.0*L13*L22*L23*L32;
    MG[(2,8)]=2.0*L12*L22*L23.powi(2) + 2.0*L13*L22.powi(2)*L23;
    MG[(2,9)]=L12.powi(2)*L33.powi(2) + 4.0*L12*L13*L32*L33 + L13.powi(2)*L32.powi(2);
    MG[(2,10)]=2.0*L12.powi(2)*L23*L33 + 4.0*L12*L13*L22*L33 + 4.0*L12*L13*L23*L32 + 2.0*L13.powi(2)*L22*L32;
    MG[(2,11)]=L12.powi(2)*L23.powi(2) + 4.0*L12*L13*L22*L23 + L13.powi(2)*L22.powi(2);
    MG[(2,12)]=2.0*L12.powi(2)*L13*L33 + 2.0*L12*L13.powi(2)*L32;
    MG[(2,13)]=2.0*L12.powi(2)*L13*L23 + 2.0*L12*L13.powi(2)*L22;
    MG[(2,14)]=L12.powi(2)*L13.powi(2);
    MG[(3,0)]=L32.powi(3)*L33;
    MG[(3,1)]=3.0*L22*L32.powi(2)*L33 + L23*L32.powi(3);
    MG[(3,2)]=3.0*L22.powi(2)*L32*L33 + 3.0*L22*L23*L32.powi(2);
    MG[(3,3)]=L22.powi(3)*L33 + 3.0*L22.powi(2)*L23*L32;
    MG[(3,4)]=L22.powi(3)*L23;
    MG[(3,5)]=3.0*L12*L32.powi(2)*L33 + L13*L32.powi(3);
    MG[(3,6)]=6.0*L12*L22*L32*L33 + 3.0*L12*L23*L32.powi(2) + 3.0*L13*L22*L32.powi(2);
    MG[(3,7)]=3.0*L12*L22.powi(2)*L33 + 6.0*L12*L22*L23*L32 + 3.0*L13*L22.powi(2)*L32;
    MG[(3,8)]=3.0*L12*L22.powi(2)*L23 + L13*L22.powi(3);
    MG[(3,9)]=3.0*L12.powi(2)*L32*L33 + 3.0*L12*L13*L32.powi(2);
    MG[(3,10)]=3.0*L12.powi(2)*L22*L33 + 3.0*L12.powi(2)*L23*L32 + 6.0*L12*L13*L22*L32;
    MG[(3,11)]=3.0*L12.powi(2)*L22*L23 + 3.0*L12*L13*L22.powi(2);
    MG[(3,12)]=L12.powi(3)*L33 + 3.0*L12.powi(2)*L13*L32;
    MG[(3,13)]=L12.powi(3)*L23 + 3.0*L12.powi(2)*L13*L22;
    MG[(3,14)]=L12.powi(3)*L13;
    MG[(4,0)]=L32.powi(4);
    MG[(4,1)]=4.0*L22*L32.powi(3);
    MG[(4,2)]=6.0*L22.powi(2)*L32.powi(2);
    MG[(4,3)]=4.0*L22.powi(3)*L32;
    MG[(4,4)]=L22.powi(4);
    MG[(4,5)]=4.0*L12*L32.powi(3);
    MG[(4,6)]=12.0*L12*L22*L32.powi(2);
    MG[(4,7)]=12.0*L12*L22.powi(2)*L32;
    MG[(4,8)]=4.0*L12*L22.powi(3);
    MG[(4,9)]=6.0*L12.powi(2)*L32.powi(2);
    MG[(4,10)]=12.0*L12.powi(2)*L22*L32;
    MG[(4,11)]=6.0*L12.powi(2)*L22.powi(2);
    MG[(4,12)]=4.0*L12.powi(3)*L32;
    MG[(4,13)]=4.0*L12.powi(3)*L22;
    MG[(4,14)]=L12.powi(4);
    MG[(5,0)]=L31*L33.powi(3);
    MG[(5,1)]=L21*L33.powi(3) + 3.0*L23*L31*L33.powi(2);
    MG[(5,2)]=3.0*L21*L23*L33.powi(2) + 3.0*L23.powi(2)*L31*L33;
    MG[(5,3)]=3.0*L21*L23.powi(2)*L33 + L23.powi(3)*L31;
    MG[(5,4)]=L21*L23.powi(3);
    MG[(5,5)]=L11*L33.powi(3) + 3.0*L13*L31*L33.powi(2);
    MG[(5,6)]=3.0*L11*L23*L33.powi(2) + 3.0*L13*L21*L33.powi(2) + 6.0*L13*L23*L31*L33;
    MG[(5,7)]=3.0*L11*L23.powi(2)*L33 + 6.0*L13*L21*L23*L33 + 3.0*L13*L23.powi(2)*L31;
    MG[(5,8)]=L11*L23.powi(3) + 3.0*L13*L21*L23.powi(2);
    MG[(5,9)]=3.0*L11*L13*L33.powi(2) + 3.0*L13.powi(2)*L31*L33;
    MG[(5,10)]=6.0*L11*L13*L23*L33 + 3.0*L13.powi(2)*L21*L33 + 3.0*L13.powi(2)*L23*L31;
    MG[(5,11)]=3.0*L11*L13*L23.powi(2) + 3.0*L13.powi(2)*L21*L23;
    MG[(5,12)]=3.0*L11*L13.powi(2)*L33 + L13.powi(3)*L31;
    MG[(5,13)]=3.0*L11*L13.powi(2)*L23 + L13.powi(3)*L21;
    MG[(5,14)]=L11*L13.powi(3);
    MG[(6,0)]=L31*L32*L33.powi(2);
    MG[(6,1)]=L21*L32*L33.powi(2) + L22*L31*L33.powi(2) + 2.0*L23*L31*L32*L33;
    MG[(6,2)]=L21*L22*L33.powi(2) + 2.0*L21*L23*L32*L33 + 2.0*L22*L23*L31*L33 + L23.powi(2)*L31*L32;
    MG[(6,3)]=2.0*L21*L22*L23*L33 + L21*L23.powi(2)*L32 + L22*L23.powi(2)*L31;
    MG[(6,4)]=L21*L22*L23.powi(2);
    MG[(6,5)]=L11*L32*L33.powi(2) + L12*L31*L33.powi(2) + 2.0*L13*L31*L32*L33;
    MG[(6,6)]=L11*L22*L33.powi(2) + 2.0*L11*L23*L32*L33 + L12*L21*L33.powi(2) + 2.0*L12*L23*L31*L33 + 2.0*L13*L21*L32*L33 + 2.0*L13*L22*L31*L33 + 2.0*L13*L23*L31*L32;
    MG[(6,7)]=2.0*L11*L22*L23*L33 + L11*L23.powi(2)*L32 + 2.0*L12*L21*L23*L33 + L12*L23.powi(2)*L31 + 2.0*L13*L21*L22*L33 + 2.0*L13*L21*L23*L32 + 2.0*L13*L22*L23*L31;
    MG[(6,8)]=L11*L22*L23.powi(2) + L12*L21*L23.powi(2) + 2.0*L13*L21*L22*L23;
    MG[(6,9)]=L11*L12*L33.powi(2) + 2.0*L11*L13*L32*L33 + 2.0*L12*L13*L31*L33 + L13.powi(2)*L31*L32;
    MG[(6,10)]=2.0*L11*L12*L23*L33 + 2.0*L11*L13*L22*L33 + 2.0*L11*L13*L23*L32 + 2.0*L12*L13*L21*L33 + 2.0*L12*L13*L23*L31 + L13.powi(2)*L21*L32 + L13.powi(2)*L22*L31;
    MG[(6,11)]=L11*L12*L23.powi(2) + 2.0*L11*L13*L22*L23 + 2.0*L12*L13*L21*L23 + L13.powi(2)*L21*L22;
    MG[(6,12)]=2.0*L11*L12*L13*L33 + L11*L13.powi(2)*L32 + L12*L13.powi(2)*L31;
    MG[(6,13)]=2.0*L11*L12*L13*L23 + L11*L13.powi(2)*L22 + L12*L13.powi(2)*L21;
    MG[(6,14)]=L11*L12*L13.powi(2);
    MG[(7,0)]=L31*L32.powi(2)*L33;
    MG[(7,1)]=L21*L32.powi(2)*L33 + 2.0*L22*L31*L32*L33 + L23*L31*L32.powi(2);
    MG[(7,2)]=2.0*L21*L22*L32*L33 + L21*L23*L32.powi(2) + L22.powi(2)*L31*L33 + 2.0*L22*L23*L31*L32;
    MG[(7,3)]=L21*L22.powi(2)*L33 + 2.0*L21*L22*L23*L32 + L22.powi(2)*L23*L31;
    MG[(7,4)]=L21*L22.powi(2)*L23;
    MG[(7,5)]=L11*L32.powi(2)*L33 + 2.0*L12*L31*L32*L33 + L13*L31*L32.powi(2);
    MG[(7,6)]=2.0*L11*L22*L32*L33 + L11*L23*L32.powi(2) + 2.0*L12*L21*L32*L33 + 2.0*L12*L22*L31*L33 + 2.0*L12*L23*L31*L32 + L13*L21*L32.powi(2) + 2.0*L13*L22*L31*L32;
    MG[(7,7)]=L11*L22.powi(2)*L33 + 2.0*L11*L22*L23*L32 + 2.0*L12*L21*L22*L33 + 2.0*L12*L21*L23*L32 + 2.0*L12*L22*L23*L31 + 2.0*L13*L21*L22*L32 + L13*L22.powi(2)*L31;
    MG[(7,8)]=L11*L22.powi(2)*L23 + 2.0*L12*L21*L22*L23 + L13*L21*L22.powi(2);
    MG[(7,9)]=2.0*L11*L12*L32*L33 + L11*L13*L32.powi(2) + L12.powi(2)*L31*L33 + 2.0*L12*L13*L31*L32;
    MG[(7,10)]=2.0*L11*L12*L22*L33 + 2.0*L11*L12*L23*L32 + 2.0*L11*L13*L22*L32 + L12.powi(2)*L21*L33 + L12.powi(2)*L23*L31 + 2.0*L12*L13*L21*L32 + 2.0*L12*L13*L22*L31;
    MG[(7,11)]=2.0*L11*L12*L22*L23 + L11*L13*L22.powi(2) + L12.powi(2)*L21*L23 + 2.0*L12*L13*L21*L22;
    MG[(7,12)]=L11*L12.powi(2)*L33 + 2.0*L11*L12*L13*L32 + L12.powi(2)*L13*L31;
    MG[(7,13)]=L11*L12.powi(2)*L23 + 2.0*L11*L12*L13*L22 + L12.powi(2)*L13*L21;
    MG[(7,14)]=L11*L12.powi(2)*L13;
    MG[(8,0)]=L31*L32.powi(3);
    MG[(8,1)]=L21*L32.powi(3) + 3.0*L22*L31*L32.powi(2);
    MG[(8,2)]=3.0*L21*L22*L32.powi(2) + 3.0*L22.powi(2)*L31*L32;
    MG[(8,3)]=3.0*L21*L22.powi(2)*L32 + L22.powi(3)*L31;
    MG[(8,4)]=L21*L22.powi(3);
    MG[(8,5)]=L11*L32.powi(3) + 3.0*L12*L31*L32.powi(2);
    MG[(8,6)]=3.0*L11*L22*L32.powi(2) + 3.0*L12*L21*L32.powi(2) + 6.0*L12*L22*L31*L32;
    MG[(8,7)]=3.0*L11*L22.powi(2)*L32 + 6.0*L12*L21*L22*L32 + 3.0*L12*L22.powi(2)*L31;
    MG[(8,8)]=L11*L22.powi(3) + 3.0*L12*L21*L22.powi(2);
    MG[(8,9)]=3.0*L11*L12*L32.powi(2) + 3.0*L12.powi(2)*L31*L32;
    MG[(8,10)]=6.0*L11*L12*L22*L32 + 3.0*L12.powi(2)*L21*L32 + 3.0*L12.powi(2)*L22*L31;
    MG[(8,11)]=3.0*L11*L12*L22.powi(2) + 3.0*L12.powi(2)*L21*L22;
    MG[(8,12)]=3.0*L11*L12.powi(2)*L32 + L12.powi(3)*L31;
    MG[(8,13)]=3.0*L11*L12.powi(2)*L22 + L12.powi(3)*L21;
    MG[(8,14)]=L11*L12.powi(3);
    MG[(9,0)]=L31.powi(2)*L33.powi(2);
    MG[(9,1)]=2.0*L21*L31*L33.powi(2) + 2.0*L23*L31.powi(2)*L33;
    MG[(9,2)]=L21.powi(2)*L33.powi(2) + 4.0*L21*L23*L31*L33 + L23.powi(2)*L31.powi(2);
    MG[(9,3)]=2.0*L21.powi(2)*L23*L33 + 2.0*L21*L23.powi(2)*L31;
    MG[(9,4)]=L21.powi(2)*L23.powi(2);
    MG[(9,5)]=2.0*L11*L31*L33.powi(2) + 2.0*L13*L31.powi(2)*L33;
    MG[(9,6)]=2.0*L11*L21*L33.powi(2) + 4.0*L11*L23*L31*L33 + 4.0*L13*L21*L31*L33 + 2.0*L13*L23*L31.powi(2);
    MG[(9,7)]=4.0*L11*L21*L23*L33 + 2.0*L11*L23.powi(2)*L31 + 2.0*L13*L21.powi(2)*L33 + 4.0*L13*L21*L23*L31;
    MG[(9,8)]=2.0*L11*L21*L23.powi(2) + 2.0*L13*L21.powi(2)*L23;
    MG[(9,9)]=L11.powi(2)*L33.powi(2) + 4.0*L11*L13*L31*L33 + L13.powi(2)*L31.powi(2);
    MG[(9,10)]=2.0*L11.powi(2)*L23*L33 + 4.0*L11*L13*L21*L33 + 4.0*L11*L13*L23*L31 + 2.0*L13.powi(2)*L21*L31;
    MG[(9,11)]=L11.powi(2)*L23.powi(2) + 4.0*L11*L13*L21*L23 + L13.powi(2)*L21.powi(2);
    MG[(9,12)]=2.0*L11.powi(2)*L13*L33 + 2.0*L11*L13.powi(2)*L31;
    MG[(9,13)]=2.0*L11.powi(2)*L13*L23 + 2.0*L11*L13.powi(2)*L21;
    MG[(9,14)]=L11.powi(2)*L13.powi(2);
    MG[(10,0)]=L31.powi(2)*L32*L33;
    MG[(10,1)]=2.0*L21*L31*L32*L33 + L22*L31.powi(2)*L33 + L23*L31.powi(2)*L32;
    MG[(10,2)]=L21.powi(2)*L32*L33 + 2.0*L21*L22*L31*L33 + 2.0*L21*L23*L31*L32 + L22*L23*L31.powi(2);
    MG[(10,3)]=L21.powi(2)*L22*L33 + L21.powi(2)*L23*L32 + 2.0*L21*L22*L23*L31;
    MG[(10,4)]=L21.powi(2)*L22*L23;
    MG[(10,5)]=2.0*L11*L31*L32*L33 + L12*L31.powi(2)*L33 + L13*L31.powi(2)*L32;
    MG[(10,6)]=2.0*L11*L21*L32*L33 + 2.0*L11*L22*L31*L33 + 2.0*L11*L23*L31*L32 + 2.0*L12*L21*L31*L33 + L12*L23*L31.powi(2) + 2.0*L13*L21*L31*L32 + L13*L22*L31.powi(2);
    MG[(10,7)]=2.0*L11*L21*L22*L33 + 2.0*L11*L21*L23*L32 + 2.0*L11*L22*L23*L31 + L12*L21.powi(2)*L33 + 2.0*L12*L21*L23*L31 + L13*L21.powi(2)*L32 + 2.0*L13*L21*L22*L31;
    MG[(10,8)]=2.0*L11*L21*L22*L23 + L12*L21.powi(2)*L23 + L13*L21.powi(2)*L22;
    MG[(10,9)]=L11.powi(2)*L32*L33 + 2.0*L11*L12*L31*L33 + 2.0*L11*L13*L31*L32 + L12*L13*L31.powi(2);
    MG[(10,10)]=L11.powi(2)*L22*L33 + L11.powi(2)*L23*L32 + 2.0*L11*L12*L21*L33 + 2.0*L11*L12*L23*L31 + 2.0*L11*L13*L21*L32 + 2.0*L11*L13*L22*L31 + 2.0*L12*L13*L21*L31;
    MG[(10,11)]=L11.powi(2)*L22*L23 + 2.0*L11*L12*L21*L23 + 2.0*L11*L13*L21*L22 + L12*L13*L21.powi(2);
    MG[(10,12)]=L11.powi(2)*L12*L33 + L11.powi(2)*L13*L32 + 2.0*L11*L12*L13*L31;
    MG[(10,13)]=L11.powi(2)*L12*L23 + L11.powi(2)*L13*L22 + 2.0*L11*L12*L13*L21;
    MG[(10,14)]=L11.powi(2)*L12*L13;
    MG[(11,0)]=L31.powi(2)*L32.powi(2);
    MG[(11,1)]=2.0*L21*L31*L32.powi(2) + 2.0*L22*L31.powi(2)*L32;
    MG[(11,2)]=L21.powi(2)*L32.powi(2) + 4.0*L21*L22*L31*L32 + L22.powi(2)*L31.powi(2);
    MG[(11,3)]=2.0*L21.powi(2)*L22*L32 + 2.0*L21*L22.powi(2)*L31;
    MG[(11,4)]=L21.powi(2)*L22.powi(2);
    MG[(11,5)]=2.0*L11*L31*L32.powi(2) + 2.0*L12*L31.powi(2)*L32;
    MG[(11,6)]=2.0*L11*L21*L32.powi(2) + 4.0*L11*L22*L31*L32 + 4.0*L12*L21*L31*L32 + 2.0*L12*L22*L31.powi(2);
    MG[(11,7)]=4.0*L11*L21*L22*L32 + 2.0*L11*L22.powi(2)*L31 + 2.0*L12*L21.powi(2)*L32 + 4.0*L12*L21*L22*L31;
    MG[(11,8)]=2.0*L11*L21*L22.powi(2) + 2.0*L12*L21.powi(2)*L22;
    MG[(11,9)]=L11.powi(2)*L32.powi(2) + 4.0*L11*L12*L31*L32 + L12.powi(2)*L31.powi(2);
    MG[(11,10)]=2.0*L11.powi(2)*L22*L32 + 4.0*L11*L12*L21*L32 + 4.0*L11*L12*L22*L31 + 2.0*L12.powi(2)*L21*L31;
    MG[(11,11)]=L11.powi(2)*L22.powi(2) + 4.0*L11*L12*L21*L22 + L12.powi(2)*L21.powi(2);
    MG[(11,12)]=2.0*L11.powi(2)*L12*L32 + 2.0*L11*L12.powi(2)*L31;
    MG[(11,13)]=2.0*L11.powi(2)*L12*L22 + 2.0*L11*L12.powi(2)*L21;
    MG[(11,14)]=L11.powi(2)*L12.powi(2);
    MG[(12,0)]=L31.powi(3)*L33;
    MG[(12,1)]=3.0*L21*L31.powi(2)*L33 + L23*L31.powi(3);
    MG[(12,2)]=3.0*L21.powi(2)*L31*L33 + 3.0*L21*L23*L31.powi(2);
    MG[(12,3)]=L21.powi(3)*L33 + 3.0*L21.powi(2)*L23*L31;
    MG[(12,4)]=L21.powi(3)*L23;
    MG[(12,5)]=3.0*L11*L31.powi(2)*L33 + L13*L31.powi(3);
    MG[(12,6)]=6.0*L11*L21*L31*L33 + 3.0*L11*L23*L31.powi(2) + 3.0*L13*L21*L31.powi(2);
    MG[(12,7)]=3.0*L11*L21.powi(2)*L33 + 6.0*L11*L21*L23*L31 + 3.0*L13*L21.powi(2)*L31;
    MG[(12,8)]=3.0*L11*L21.powi(2)*L23 + L13*L21.powi(3);
    MG[(12,9)]=3.0*L11.powi(2)*L31*L33 + 3.0*L11*L13*L31.powi(2);
    MG[(12,10)]=3.0*L11.powi(2)*L21*L33 + 3.0*L11.powi(2)*L23*L31 + 6.0*L11*L13*L21*L31;
    MG[(12,11)]=3.0*L11.powi(2)*L21*L23 + 3.0*L11*L13*L21.powi(2);
    MG[(12,12)]=L11.powi(3)*L33 + 3.0*L11.powi(2)*L13*L31;
    MG[(12,13)]=L11.powi(3)*L23 + 3.0*L11.powi(2)*L13*L21;
    MG[(12,14)]=L11.powi(3)*L13;
    MG[(13,0)]=L31.powi(3)*L32;
    MG[(13,1)]=3.0*L21*L31.powi(2)*L32 + L22*L31.powi(3);
    MG[(13,2)]=3.0*L21.powi(2)*L31*L32 + 3.0*L21*L22*L31.powi(2);
    MG[(13,3)]=L21.powi(3)*L32 + 3.0*L21.powi(2)*L22*L31;
    MG[(13,4)]=L21.powi(3)*L22;
    MG[(13,5)]=3.0*L11*L31.powi(2)*L32 + L12*L31.powi(3);
    MG[(13,6)]=6.0*L11*L21*L31*L32 + 3.0*L11*L22*L31.powi(2) + 3.0*L12*L21*L31.powi(2);
    MG[(13,7)]=3.0*L11*L21.powi(2)*L32 + 6.0*L11*L21*L22*L31 + 3.0*L12*L21.powi(2)*L31;
    MG[(13,8)]=3.0*L11*L21.powi(2)*L22 + L12*L21.powi(3);
    MG[(13,9)]=3.0*L11.powi(2)*L31*L32 + 3.0*L11*L12*L31.powi(2);
    MG[(13,10)]=3.0*L11.powi(2)*L21*L32 + 3.0*L11.powi(2)*L22*L31 + 6.0*L11*L12*L21*L31;
    MG[(13,11)]=3.0*L11.powi(2)*L21*L22 + 3.0*L11*L12*L21.powi(2);
    MG[(13,12)]=L11.powi(3)*L32 + 3.0*L11.powi(2)*L12*L31;
    MG[(13,13)]=L11.powi(3)*L22 + 3.0*L11.powi(2)*L12*L21;
    MG[(13,14)]=L11.powi(3)*L12;
    MG[(14,0)]=L31.powi(4);
    MG[(14,1)]=4.0*L21*L31.powi(3);
    MG[(14,2)]=6.0*L21.powi(2)*L31.powi(2);
    MG[(14,3)]=4.0*L21.powi(3)*L31;
    MG[(14,4)]=L21.powi(4);
    MG[(14,5)]=4.0*L11*L31.powi(3);
    MG[(14,6)]=12.0*L11*L21*L31.powi(2);
    MG[(14,7)]=12.0*L11*L21.powi(2)*L31;
    MG[(14,8)]=4.0*L11*L21.powi(3);
    MG[(14,9)]=6.0*L11.powi(2)*L31.powi(2);
    MG[(14,10)]=12.0*L11.powi(2)*L21*L31;
    MG[(14,11)]=6.0*L11.powi(2)*L21.powi(2);
    MG[(14,12)]=4.0*L11.powi(3)*L31;
    MG[(14,13)]=4.0*L11.powi(3)*L21;
    MG[(14,14)]=L11.powi(4);
    MG
}

fn get_nmat(lmns: &Vec<[usize; 3]>) -> Array2<f64> {
    let nval = |lmn: &[usize; 3]| {
        let [l, m, n] = lmn;
        let fac = FACS[*l] * FACS[*m] * FACS[*n];
        1.0 / fac.sqrt()
    };
    let size = lmns.len();
    let mut nmat = Array2::zeros([size, size]);
    for (i, lmn_i) in lmns.iter().enumerate() {
        let val_i = nval(lmn_i);
        for (j, lmn_j) in lmns.iter().enumerate() {
            let val_j = nval(lmn_j);
            nmat[(i, j)] = val_i / val_j;
        }
    }
    nmat
}

#[rustfmt::skip]
pub fn nmat_d()->Array2<f64>{
    let lmns=btrans::CAR_D_LMNS.to_vec();
    get_nmat(&lmns)
}

#[rustfmt::skip]
pub fn nmat_f()->Array2<f64>{
    let lmns=btrans::CAR_F_LMNS.to_vec();
    get_nmat(&lmns)
}

#[rustfmt::skip]
pub fn nmat_g()->Array2<f64>{
    let lmns=btrans::CAR_G_LMNS.to_vec();
    get_nmat(&lmns)
}
