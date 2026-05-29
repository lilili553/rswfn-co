use crate::maths;
use ndarray::{Array1, Array2};
use std::path::Path;

pub fn has_nan(matrix: &Array2<f64>) -> bool {
    matrix.iter().any(|&x| x.is_nan())
}

pub fn show_matrix(mat: &Array2<f64>) {
    let nrow = mat.nrows();
    let ncol = mat.ncols();
    println!("shape: ({}, {})", nrow, ncol);
    for i in 0..nrow {
        for j in 0..ncol {
            print!("{:>10.5}", mat[(i, j)]);
        }
        println!("");
    }
}

// 比较两个矩阵
pub fn comp_matrix(mat1: &Array2<f64>, mat2: &Array2<f64>, show: isize) {
    let mut del_max = f64::MIN;
    let mut del_min = f64::MAX;
    assert_eq!(mat1.shape(), mat2.shape());
    println!("矩阵对比:");
    let nrow = mat1.nrows();
    let ncol = mat1.ncols();
    let mut errs = vec![];
    for i in 0..nrow {
        for j in 0..ncol {
            let val1 = mat1[(i, j)];
            let val2 = mat2[(i, j)];
            let del = val1 - val2;
            if show >= 1 {
                print!("{:>8.4}-{:>8.4}={:>8.4}    |", val1, val2, del);
            }

            if del.abs() > 1e-3 {
                errs.push((i, j, val1, val2));
            }
            if del > del_max {
                del_max = del;
            }
            if del < del_min {
                del_min = del;
            }
        }
        if show >= 1 {
            println!("");
        }
    }
    if show >= 2 {
        for (i, j, v1, v2) in errs {
            println!("({:>2},{:>2}):{:>8.4}-{:>8.4}={:>8.4}", i + 1, j + 1, v1, v2, v1 - v2);
        }
    }

    println!("最大差:{},最小差:{}", del_max, del_min);
}

// 判断两个两个矩阵是否相同
pub fn same_matrix(mat1: &Array2<f64>, mat2: &Array2<f64>) -> bool {
    let nrow = mat1.nrows();
    let ncol = mat1.ncols();
    for i in 0..nrow {
        for j in 0..ncol {
            let val1 = mat1[(i, j)];
            let val2 = mat2[(i, j)];
            let del = val1 - val2;
            if del.abs() > 1e-4 {
                return false;
            }
        }
    }
    true
}

pub fn symm_matrix(mat: &Array2<f64>) -> bool {
    let nrow = mat.nrows();
    let ncol = mat.nrows();
    if nrow != ncol {
        return false;
    }
    for i in 0..nrow {
        for j in 0..i {
            let del = mat[(i, j)] - mat[(j, i)];
            if del.abs() > 1e-4 {
                return false;
            }
        }
    }
    true
}

// 判断矩阵是否为正定矩阵
pub fn posi_matrix(mat: &Array2<f64>) -> bool {
    let nrow = mat.nrows();
    let ncol = mat.nrows();
    if nrow != ncol {
        return false;
    }
    let mat_t = mat.t();
    let mat_i = maths::inv(mat);
    for i in 0..nrow {
        for j in 0..ncol {
            let del = mat_t[(i, j)] - mat_i[(i, j)];
            if del.abs() > 1e-4 {
                return false;
            }
        }
    }
    true
}

pub fn fort_form(val: f64, pec: usize) -> String {
    let raw_text = format!("{:>.*E}", pec, val);
    let idx = raw_text.find("E").expect("解析E失败");
    let head = &raw_text[..idx];
    let tail: isize = raw_text[idx + 1..].to_string().parse().expect("解析失败");
    let tail = tail.abs();
    let unit = if val.abs() >= 1.0 { "+" } else { "-" };
    let fmt_text = format!("{head}E{unit}{tail:0>2}");
    fmt_text
}

pub fn file_type(path: &str) -> String {
    let path = Path::new(path);
    path.extension().map(|ext| ext.to_string_lossy().into_owned()).expect("文件拓展名获取失败")
}

pub fn file_name(path: &str) -> String {
    let path = Path::new(path);
    path.file_stem().map(|ext| ext.to_string_lossy().into_owned()).expect("文件名获取失败")
}

// 判断两个数组中的元素是否完全一样
pub fn all_same<T: PartialEq>(vals0: &Vec<T>, vals1: &Vec<T>) -> bool {
    if vals0.len() != vals1.len() {
        return false;
    }
    for i in 0..vals0.len() {
        if vals0[i] != vals1[i] {
            return false;
        }
    }
    return true;
}
