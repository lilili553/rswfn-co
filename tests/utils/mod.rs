pub fn allclose(res: Vec<f64>, ans: Vec<f64>, thread: f64) -> bool {
    if res.len() != ans.len() {
        return false;
    }
    for i in 0..res.len() {
        if (res[i] - ans[i]).abs() > thread {
            return false;
        }
    }
    return true;
}
