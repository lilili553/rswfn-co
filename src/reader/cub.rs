use std::path;

use crate::{
    base::{Atoms, Basis, Coefs},
    datas::consts::BOHR,
    gridprop::CubeGrid,
    reader::{Reader, read_lines},
};

pub struct CubReader {
    pub path: String,
    pub lines: Vec<String>,
    pub mobt: bool,
    pub natm: usize,
}

impl CubReader {
    pub fn new(path: &str) -> Self {
        let lines = read_lines(path);
        let path = path.to_string();
        let natm = lines[2][..5].trim().parse::<isize>().unwrap();
        let mobt = if natm < 0 { true } else { false };
        let natm = natm.abs() as usize;
        Self { path, lines, mobt, natm }
    }

    pub fn read_cube(&self) -> CubeGrid {
        let mut start = [0.0; 3];
        start[0] = self.lines[2][5..17].trim().parse::<f64>().unwrap();
        start[1] = self.lines[2][17..29].trim().parse::<f64>().unwrap();
        start[2] = self.lines[2][29..41].trim().parse::<f64>().unwrap();
        let mut shape = [0; 3];
        for i in 0..3 {
            let line = &self.lines[i + 3];
            shape[i] = line[..5].trim().parse::<usize>().unwrap();
        }
        let sx = self.lines[3][5..17].trim().parse::<f64>().unwrap();
        let sy = self.lines[4][17..29].trim().parse::<f64>().unwrap();
        let sz = self.lines[5][29..41].trim().parse::<f64>().unwrap();
        let steps = [sx, sy, sz];
        CubeGrid { shape, start, steps }
    }

    pub fn read_atms(&self) -> (Vec<usize>, Vec<f64>, Vec<[f64; 3]>) {
        let natm = self.natm;
        let mut atms = Vec::with_capacity(natm);
        let mut crgs = Vec::with_capacity(natm);
        let mut xyzs = Vec::with_capacity(natm);
        for l in 6..6 + natm {
            let line = &self.lines[l];
            let atm = line[..5].trim().parse::<usize>().unwrap();
            let crg = line[5..17].trim().parse::<f64>().unwrap();
            let x = line[17..29].trim().parse::<f64>().unwrap() / BOHR;
            let y = line[29..41].trim().parse::<f64>().unwrap() / BOHR;
            let z = line[41..53].trim().parse::<f64>().unwrap() / BOHR;
            atms.push(atm);
            crgs.push(crg);
            xyzs.push([x, y, z]);
        }
        (atms, crgs, xyzs)
    }

    pub fn read_vals(&self) -> (Vec<usize>, Vec<Vec<f64>>) {
        let nobt = self.lines[6 + self.natm][..5].trim().parse::<usize>().unwrap();
        let obts: Vec<usize> = self.lines[6 + self.natm][5..].trim().split_whitespace().map(|e| e.parse::<usize>().unwrap()).collect();
        let [nx, ny, nz] = self.read_cube().shape;
        let count = nx * ny * nz;
        let mut vals = vec![vec![0.0; count]; nobt];
        let mut idx = 0;
        for line in self.lines.iter().skip(6 + self.natm + 1) {
            for i in (0..line.len()).step_by(13) {
                let val = line[i..i + 13].trim().parse::<f64>().unwrap();
                // println!("{}->{}", line[i..i + 13].to_string(), val);
                vals[idx % nobt][idx / nobt] = val;
                idx += 1;
            }
        }
        (obts, vals)
    }
}

impl Reader for CubReader {
    fn get_atoms(&self) -> Atoms {
        let (atms, crgs, xyzs) = self.read_atms();
        Atoms::new(&atms, &xyzs, &crgs, None)
    }

    fn get_basis(&self) -> Basis {
        Basis::blank()
    }

    fn get_coefs(&self) -> Coefs {
        Coefs::blank()
    }
}
