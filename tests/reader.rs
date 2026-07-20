use rswfn::base::Mole;

#[test]
/// 测试读取fch文件
fn test1() {
    let path = String::from(r"E:\code\rswfn\moles\C6H6.fch");
    let mole = Mole::from_file(&path);
}

#[test]
/// 测试带赝势基组的log文件的读取
fn test2() {
    let mole1 = Mole::from_file(r"E:\code\rswfn\moles\dui-Nb.log_scan_14.fchk");
    let mole2 = Mole::from_file(r"E:\code\rswfn\moles\dui-Nb.log_scan_14.log");
    let same = mole1.equal(&mole2);
    assert!(same);
}
