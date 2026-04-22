use std::mem::MaybeUninit;
use std::os::raw::c_int;

use crate::funcs::bind_funcs;
use crate::funcs::bind_xc;

// 调用libxc

#[cfg(feature = "use-c")]
pub fn xc_version() -> String {
    let mut major = 0;
    let mut minor = 0;
    let mut micro = 0;
    unsafe {
        bind_xc::xc_version(&mut major, &mut minor, &mut micro);
    }
    format!("{}.{}.{}", major, minor, micro)
}

#[cfg(feature = "use-c")]
pub fn init_func(open: bool) -> MaybeUninit<bind_xc::xc_func_type> {
    let mut xc_func = MaybeUninit::<bind_xc::xc_func_type>::uninit();
    let nspin = if open { bind_xc::XC_POLARIZED } else { bind_xc::XC_UNPOLARIZED };
    let func_id = bind_funcs::XC_HYB_GGA_XC_B3LYP;
    // let func_id = bind_funcs::XC_LDA_X;
    unsafe {
        let res = bind_xc::xc_func_init(xc_func.as_mut_ptr(), func_id as c_int, nspin as c_int);
        println!("xc_init:{}", res);
        assert!(res == 0, "未发现的泛函ID");
    }
    xc_func
}

#[cfg(feature = "use-c")]
// 获取泛函的类型
pub fn get_family(mut xc_func: MaybeUninit<bind_xc::xc_func_type>) -> u32 {
    let xc_family = unsafe {
        let xc_info = bind_xc::xc_func_get_info(xc_func.as_mut_ptr());
        let xc_family = bind_xc::xc_func_info_get_family(xc_info);
        println!("获取泛函类型:{}", xc_family);
        assert!(xc_family != -1, "未知的泛函类型");
        xc_family
    };
    return xc_family as u32;
}
