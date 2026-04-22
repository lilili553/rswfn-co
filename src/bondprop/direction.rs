use crate::base::{Mole, Vec3};
use crate::maths;
use std::f64::consts::PI;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    pub fn verts(&self, ai: usize, aj: usize) -> Vec<Vec3> {
        let mut vectors = vec![];
        let vr = Vec3::rand();
        let bond_vec = self.mole.atoms().vec(ai, aj).normalize(); // 键轴单位向量
        let norm_vec = bond_vec.cross(&vr);
        let bond_cent = self.mole.atoms.cent(&vec![ai, aj]);
        for i in 0..36 {
            let vector = norm_vec.rotate_arrow(&bond_cent, &bond_vec, i as f64 * PI / 18.0);
            vectors.push(vector.normalize());
        }
        vectors
    }
}
