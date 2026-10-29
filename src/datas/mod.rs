pub mod btrans;
pub mod consts;
pub mod dftfc;
pub mod dhmo2;
pub mod elems;
pub mod integ;
pub mod lebedev;
pub mod rdens;
pub mod roots;
pub mod rys_root;
pub mod temps;

// use std::collections::HashMap;

// let mut cache: HashMap<String, CacheTypes> = HashMap::new();

// enum CacheTypes {
//     Int(usize),
//     Float(f64),
// }

// pub static RAND_VECT: LazyLock<Vector3<f64>> = LazyLock::new(|| {
//     //随机向量
//     let x = rand::random();
//     let y = rand::random();
//     let z = rand::random();
//     Vector3::new(x, y, z).normalize()
// });

// 双阶乘 (2i-1)!!
pub const FACS: [f64; 6] = [1., 1., 3., 15., 105., 945.];
