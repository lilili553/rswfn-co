fn main() {
    #[cfg(all(feature = "use-cint", feature = "use-func"))]
    {
        let root_dir = std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR 未设置");
        let libs_dir = format!("{}/libs", root_dir);
        let cint_lib = format!("{}/libcint.a", libs_dir);
        let cxcs_lib = format!("{}/libxc.a", libs_dir);

        println!("cargo:warning=Link search path: {}", libs_dir);
        println!("cargo:rustc-link-search=native={}", libs_dir);

        println!("cargo:rustc-link-lib=static=cint");
        println!("cargo:warning=Linking static library: {}", cint_lib);
        println!("cargo:rustc-link-arg={}", cint_lib);
        println!("cargo:rerun-if-changed={}", cint_lib);

        println!("cargo:rustc-link-lib=static=xc");
        println!("cargo:warning=Linking static library: {}", cxcs_lib);
        println!("cargo:rustc-link-arg={}", cxcs_lib);
        println!("cargo:rerun-if-changed={}", cxcs_lib);
    }
}
