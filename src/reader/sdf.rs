use crate::base::{Atoms, Bas, Basis, Coefs};
use crate::datas::elems::sym2atm;
use crate::reader::{Reader, read_lines};

pub struct SdfReader {
    path: String,
    lines: Vec<String>,
}

impl SdfReader {
    pub fn new(path: &str) -> Self {
        let lines = read_lines(path);
        Self { path: path.to_string(), lines }
    }

    pub fn read_geome(&self) -> (Vec<usize>, Vec<[f64; 3]>) {
        let mut atms = vec![];
        let mut xyzs = vec![];
        let natm = &self.lines[3][..3].trim().parse::<usize>().expect("usize parse error");
        let lines = &self.lines[4..4 + natm];
        for line in lines {
            let x = line[00..10].trim().parse::<f64>().expect("f64 parse error");
            let y = line[10..20].trim().parse::<f64>().expect("f64 parse error");
            let z = line[20..30].trim().parse::<f64>().expect("f64 parse error");
            let sym = line[30..32].trim().parse::<String>().expect("String parse error");
            // log::info!("{},{}", line, s);
            let atm = sym2atm(&sym);
            atms.push(atm);
            xyzs.push([x, y, z]);
        }
        (atms, xyzs)
    }
}

impl Reader for SdfReader {
    fn get_atoms(&self) -> Atoms {
        let (atms, xyzs) = self.read_geome();
        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
        Atoms::new(&atms, &xyzs, &crgs)
    }

    fn get_basis(&self) -> Option<Basis> {
        None
    }

    fn get_coefs(&self) -> Option<Coefs> {
        None
    }
}
