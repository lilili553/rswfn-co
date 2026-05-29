use crate::base::{Atoms, Basis, Coefs, Stm};
use crate::reader::Reader;

pub struct NonReader {
    path: String,
    lines: Vec<String>,
    atoms: Atoms,
    basis: Option<Basis>,
    coefs: Option<Coefs>,
    neles: [usize; 2],
}

impl NonReader {
    pub fn new(path: &str) -> Self {
        Self {
            path: path.to_string(),
            lines: vec![],
            atoms: Atoms::new(&vec![], &vec![], &vec![]),
            basis: None,
            coefs: None,
            neles: [0, 0],
        }
    }

    pub fn set_atoms(&mut self, atms: &Vec<usize>, xyzs: &Vec<[f64; 3]>) {
        let crgs = atms.iter().map(|atm| *atm as f64).collect();
        let atoms = Atoms::new(atms, xyzs, &crgs);
        self.atoms = atoms;
    }
}

impl Reader for NonReader {
    fn get_atoms(&self) -> Atoms {
        self.atoms.clone()
    }

    fn get_basis(&self) -> Option<Basis> {
        self.basis.clone()
    }

    fn get_coefs(&self) -> Option<Coefs> {
        self.coefs.clone()
    }
}
