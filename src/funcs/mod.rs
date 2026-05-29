#[cfg(feature = "use-func")]
mod bind_funcs;
#[cfg(feature = "use-func")]
mod bind_xc;

pub mod cxcs;
pub mod rxcs;
use std::mem::MaybeUninit;

#[cfg(feature = "use-func")]
pub struct RsFun {
    pub xc_func: MaybeUninit<bind_xc::xc_func_type>,
    pub xc_family: u32,
}

#[cfg(not(feature = "use-func"))]
pub struct RsFun {}

impl RsFun {
    #[cfg(feature = "use-func")]
    pub fn new() -> Self {
        let xc_func = cxcs::init_func(false);
        let xc_family = cxcs::get_family(xc_func);
        println!("libxc版本{},xc_family:{}", cxcs::xc_version(), xc_family);
        Self { xc_func, xc_family }
    }

    #[cfg(not(feature = "use-func"))]
    pub fn new() -> Self {
        Self {}
    }

    #[cfg(feature = "use-func")]
    pub fn run(&self, rho: &Vec<f64>, sigma: &Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let p = self.xc_func.as_ptr();
        let np = rho.len();
        let mut zk = vec![0.0; np];
        let mut vrho = vec![0.0; np];
        let mut vsigma = vec![0.0; np];
        unsafe {
            match self.xc_family {
                bind_xc::XC_FAMILY_LDA => {
                    bind_xc::xc_lda_exc_vxc(p, np, rho.as_ptr(), zk.as_mut_ptr(), vrho.as_mut_ptr());
                }
                bind_xc::XC_FAMILY_GGA => {
                    bind_xc::xc_gga_exc_vxc(p, np, rho.as_ptr(), sigma.as_ptr(), zk.as_mut_ptr(), vrho.as_mut_ptr(), vsigma.as_mut_ptr());
                }
                bind_xc::XC_FAMILY_HYB_GGA => {
                    bind_xc::xc_gga_exc_vxc(p, np, rho.as_ptr(), sigma.as_ptr(), zk.as_mut_ptr(), vrho.as_mut_ptr(), vsigma.as_mut_ptr());
                }
                _ => panic!("未知泛函类型:{}", self.xc_family),
            }
        }
        (zk, vrho, vsigma)
    }

    #[cfg(not(feature = "use-func"))]
    pub fn run(&self, rho: &Vec<f64>, sigma: &Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let np = rho.len();
        let zk = vec![0.0; np];
        let vrho = vec![0.0; np];
        let vsigma = vec![0.0; np];
        (zk, vrho, vsigma)
    }
}
