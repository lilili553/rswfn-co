// 编辑分子结构，改变结构后的波函数信息肯定就没用了，因此只需要返回一个只包含结构的分子
use crate::base::{Atoms, Mole};
use ndarray::Array2;

/// 分子结构编辑器，提供键的旋转和伸缩等几何操作
pub struct Editor<'m> {
    pub mole: &'m Mole,
}

impl<'m> Editor<'m> {
    /// 创建一个新的编辑器实例
    pub fn new(mole: &'m Mole) -> Self {
        Editor { mole }
    }

    /// 获取与指定原子相连的原子分组（排除指定原子）
    pub fn get_group(&self, atm: usize, excludes: &Vec<usize>) -> Vec<usize> {
        let mut searched = vec![];
        let mut searches = vec![atm];
        while !searches.is_empty() {
            let atm = searches.pop().unwrap();

            if excludes.contains(&atm) {
                continue;
            }
            searched.push(atm);
            let nebs = self.mole.atoms().nebs(atm as usize);
            for neb in nebs {
                let neb = neb as usize;
                if searched.contains(&neb) {
                    continue;
                } else {
                    if !searches.contains(&neb) {
                        searches.push(neb as usize);
                    }
                }
            }
        }
        searched
    }

    /// 旋转一根化学键，返回旋转后的新分子
    /// - a1: 第一个原子索引
    /// - a2: 第二个原子索引
    /// - angle: 旋转角度
    /// - rtype: 旋转类型（0:旋转基团1, 1:旋转基团2, 2:两个基团都旋转）
    pub fn rotate_bond(&self, a1: usize, a2: usize, angle: f64, rtype: usize) -> Mole {
        let mut xyzs = self.mole.atoms().xyzs().clone();
        let g1 = self.get_group(a1, &vec![a2]);
        let g2 = self.get_group(a2, &vec![a1]);
        let axis = self.mole.atoms().vec(a1 as usize, a2 as usize);
        let [dx, dy, dz] = axis.into_arr();
        let axis = [dx, dy, dz];
        let [ang1, ang2] = match rtype {
            0 => [-angle, 0.0],
            1 => [0.0, angle],
            2 => [-angle / 2.0, angle / 2.0],
            _ => panic!("不支持"),
        };
        let angs = [ang1, ang2];
        let idxs = [a1, a2];
        let gs = [g1, g2];
        for i in 0..2 {
            let cent = xyzs[idxs[i]];
            let angle = angs[i];
            let points: Vec<[f64; 3]> = gs[i].iter().map(|atm| xyzs[*atm]).collect();
            let points = rotate_points(&points, &cent, &axis, angle);
            for (j, atm) in gs[i].iter().enumerate() {
                xyzs[*atm as usize] = points[j].clone();
            }
        }
        let mut mole = Mole::new(self.mole.atoms.clone(), None, None);
        let atms = self.mole.atoms().nucs();
        let crgs = self.mole.atoms().crgs();
        let stms = self.mole.atoms().stms().clone();
        let atoms = Atoms::new(atms, &xyzs, crgs);
        mole.atoms = atoms;
        mole
    }

    /// 伸缩一根化学键，返回新的分子
    /// - a1: 第一个原子索引
    /// - a2: 第二个原子索引
    /// - length: 伸缩长度
    /// - stype: 伸缩类型（0:移动基团1, 1:移动基团2, 2:两个基团都移动）
    pub fn stretch_bond(&self, a1: usize, a2: usize, length: f64, stype: usize) -> Mole {
        let mut xyzs = self.mole.atoms().xyzs().clone();
        let g1 = self.get_group(a1, &vec![a2]);
        let g2 = self.get_group(a2, &vec![a1]);
        let [len1, len2] = match stype {
            0 => [-length, 0.0],
            1 => [0.0, length],
            2 => [-length / 2.0, length / 2.0],
            _ => panic!("不支持"),
        };
        let axis = self.mole.atoms().vec(a1 as usize, a2 as usize).normalize();
        let [dx, dy, dz] = axis.into_arr();
        let lens = [len1, len2];
        let gs = [g1, g2];
        for i in 0..2 {
            for atm in gs[i].iter() {
                xyzs[*atm][0] += dx * lens[i];
                xyzs[*atm][1] += dy * lens[i];
                xyzs[*atm][2] += dz * lens[i];
            }
        }
        let mut mole = Mole::new(self.mole.atoms.clone(), None, None);
        let atms = self.mole.atoms().nucs();
        let crgs = self.mole.atoms().crgs();
        let stms = self.mole.atoms().stms().clone();
        mole.atoms = Atoms::new(atms, &xyzs, crgs);
        mole
    }
}

/// 绕指定轴旋转一组点
/// - points: 待旋转的点集
/// - cent: 旋转中心
/// - axis: 旋转轴
/// - angle: 旋转角度
pub fn rotate_points(points: &Vec<[f64; 3]>, cent: &[f64; 3], axis: &[f64; 3], angle: f64) -> Vec<[f64; 3]> {
    let [dx, dy, dz] = axis;
    let length = f64::sqrt(dx * dx + dy * dy + dz * dz);
    let [dx, dy, dz] = [dx / length, dy / length, dz / length];
    let c = f64::cos(angle);
    let s = f64::sin(angle);
    #[rustfmt::skip]
    let tmat=Array2::<f64>::from_shape_vec((3,3), vec![
        dx*dx*(1.0-c)+c,dx*dy*(1.0-c)+dz*s,dx*dz*(1.0-c)-dy*s,
        dx*dy*(1.0-c)-dz*s,dy*dy*(1.0-c)+c,dy*dz*(1.0-c)+dx*s,
        dx*dz*(1.0-c)+dy*s,dy*dz*(1.0-c)-dx*s,dz*dz*(1.0-c)+c
    ]).unwrap();
    let mut new_points = vec![];
    let cent = Array2::<f64>::from_shape_vec((3, 1), cent.to_vec()).unwrap();
    for point in points.iter() {
        let mut point = Array2::<f64>::from_shape_vec((3, 1), point.to_vec()).unwrap();
        point -= &cent;
        point = tmat.dot(&point);
        point += &cent;
        let px = point[(0, 0)];
        let py = point[(1, 0)];
        let pz = point[(2, 0)];
        new_points.push([px, py, pz]);
    }
    new_points
}
