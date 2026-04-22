```rust
use crate::base::atoms::{Atoms, Stm};
use core::fmt;

// use crate::base::Mole;

#[derive(Debug, Clone)]
pub struct Geome {
    pub atoms: Atoms,
    pub neles: [usize; 2], // α电子数和β电子数
}

impl Geome {
    pub fn new(atms: &Vec<usize>, xyzs: &Vec<[f64; 3]>, crgs: &Vec<f64>, stms: Option<Vec<Stm>>) -> Self {
        assert_eq!(atms.len(), xyzs.len());
        assert_eq!(atms.len(), crgs.len());
        let mut atoms = Atoms::new();
        let natm = atms.len();
        let stms = match stms {
            Some(stms) => stms,
            None => (0..natm).into_iter().map(|_| Stm::blank()).collect(),
        };
        let mut nele = 0;
        for i in 0..natm {
            let stm = stms[i].clone();
            atoms.add(atms[i], xyzs[i], crgs[i], stm);
            nele += atms[i];
        }

        let nela = nele / 2 + nele % 2;
        let nelb = nele / 2;
        let neles = [nela, nelb];
        Geome { atoms, neles }
    }

    pub fn blank() -> Self {
        Geome { atoms: Atoms::new(), neles: [0, 0] }
    }

    pub fn build(&mut self, atms: &Vec<usize>, xyzs: &Vec<[f64; 3]>, crgs: &Vec<f64>, stms: &Vec<Stm>) {
        let natm = atms.len();
        for i in 0..natm {
            let stm = stms[i].clone();
            self.atoms.add(atms[i], xyzs[i], crgs[i], stm);
        }
    }

    // 判断两个结构是否相等
    pub fn equal(&self, other: &Geome) -> bool {
        let atoms_0 = &self.atoms;
        let atoms_1 = &other.atoms;
        if atoms_0.len() != atoms_1.len() {
            log::info!("原子数量不相等");
            return false;
        }
        for i in 0..atoms_0.len() {
            let atom_0 = atoms_0.get(i as usize);
            let atom_1 = atoms_1.get(i as usize);
            if atom_0.sym != atom_1.sym {
                log::info!("元素符号不相等");
                return false;
            }
            let [x0, y0, z0] = atom_0.xyz;
            let [x1, y1, z1] = atom_1.xyz;
            let dist = ((x1 - x0).powi(2) + (y1 - y0).powi(2) + (z1 - z0).powi(2)).sqrt();
            if dist > 1e-2 {
                log::info!("原子坐标不相等");
                return false;
            }
        }
        return true;
    }
}

impl fmt::Display for Geome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let natm = self.atoms.len();
        writeln!(f, "{:>3}{:>14}{:>14}{:>14}{:>10}", "atm", "x", "y", "z", "crg")?;
        for i in 0..natm {
            let atom = self.atoms.get(i);
            let [x, y, z] = atom.xyz;
            writeln!(f, "{:>3}{:>14.8}{:>14.8}{:>14.8}{:>10.4}", atom.sym, x, y, z, atom.crg)?;
        }
        for i in 0..natm {
            let atom = self.atoms.get(i);
            if !atom.stm.non {
                writeln!(f, "{}", atom.stm)?;
            }
        }
        Ok(())
    }
}
```
