use ndarray::Array2;

use crate::base::atoms::Atoms;
use crate::base::basis::Basis;
use crate::base::coefs::Coefs;
use crate::datas::elems;
use crate::reader::{FileType, Reader, get_fileType, read_lines};
use crate::{maths, reader};

use core::fmt;
use std::collections::{HashMap, HashSet};
use std::path::Path;

// use mem_dbg::MemDbg;

#[derive(Debug, Clone)]
pub struct Mole {
    pub atoms: Atoms,
    pub basis: Basis,
    pub coefs: Coefs,
}

impl Mole {
    pub fn new(atoms: Atoms, basis: Basis, coefs: Coefs) -> Mole {
        let mole = Self { atoms, basis, coefs };
        mole
    }

    pub fn blank() -> Mole {
        Mole { atoms: Atoms::blank(), basis: Basis::blank(), coefs: Coefs::blank() }
    }

    pub fn from_read(reader: &dyn Reader) -> Mole {
        let atoms = reader.get_atoms();
        let basis = reader.get_basis();
        let coefs = reader.get_coefs();
        let mole = Mole { atoms, basis, coefs };
        mole
    }

    pub fn from_file(path: &str) -> Mole {
        let file_type = get_fileType(path);
        match file_type {
            FileType::Gjf => return Mole::from_read(&reader::GjfReader::new(path)),
            FileType::Log => return Mole::from_read(&reader::LogReader::new(path)),
            FileType::Fch => return Mole::from_read(&reader::FchReader::new(path)),
            FileType::Mol => return Mole::from_read(&reader::MolReader::new(path)),
            FileType::Sdf => return Mole::from_read(&reader::SdfReader::new(path)),
            FileType::Xyz => return Mole::from_read(&reader::XyzReader::new(path)),
            FileType::Mde => return Mole::from_read(&reader::MdeReader::new(path)),
            FileType::Dmo => return Mole::from_read(&reader::DmoReader::new(path)),
            FileType::Inp => return Mole::from_read(&reader::InpReader::new(path)),
            FileType::Cub => return Mole::from_read(&reader::CubReader::new(path)),
        }
    }

    pub fn atoms(&self) -> &Atoms {
        &self.atoms
    }

    pub fn basis(&self) -> &Basis {
        return &self.basis;
    }

    pub fn coefs(&self) -> &Coefs {
        return &self.coefs;
    }

    pub fn is_open(&self) -> bool {
        self.multip() != 1
    }

    pub fn formula(&self) -> String {
        let mut formula = String::new();
        let syms = self.atoms().syms();
        for elem in elems::ELEMS.iter() {
            let count = syms.iter().filter(|sym| **sym == elem.sym).count();
            if count > 0 {
                formula += &format!("{}{}", elem.sym, count);
            }
        }
        formula
    }

    // 分子总核电荷数
    pub fn ncrg(&self) -> usize {
        let mut ncrg = 0.0;
        let crgs = self.atoms().crgs();
        for crg in crgs.iter() {
            ncrg += *crg;
        }
        ncrg as usize
    }

    pub fn charge(&self) -> isize {
        let [nela, nelb] = self.atoms.neles;
        let ncrg = self.ncrg();
        nela as isize + nelb as isize - ncrg as isize
    }

    pub fn multip(&self) -> usize {
        let [nela, nelb] = self.atoms.neles;
        nela - nelb + 1
    }

    pub fn dist(&self, grid: &[f64; 3]) -> f64 {
        // 格点到分子的距离
        let [gx, gy, gz] = *grid;
        let atoms = self.atoms();
        let mut min_dist = f64::MAX;
        for i in 0..atoms.len() {
            let [ax, ay, az] = atoms.get(i as usize).xyz;
            let dist = ((gx - ax).powi(2) + (gy - ay).powi(2) + (gz - az).powi(2)).sqrt();
            if dist < min_dist {
                min_dist = dist;
            }
        }
        min_dist
    }

    pub fn rings(&self) -> Vec<Vec<usize>> {
        let atoms = self.atoms();
        let natm = atoms.len();
        let mut adjacency_list = HashMap::new();
        for i in 0..natm {
            for j in 0..natm {
                if atoms.is_bond(i as usize, j as usize) {
                    adjacency_list.entry(i as usize).or_insert(Vec::new()).push(j as usize);
                }
            }
        }
        let rings = find_sssr(&adjacency_list);
        // 如果环内有相邻的原子，则排除
        let rings = rings
            .iter()
            .filter(|ring| {
                let natm = ring.len();
                let mut nbond = 0;
                for i in 0..natm - 1 {
                    for j in i + 1..natm {
                        if atoms.is_bond(ring[i], ring[j]) {
                            nbond += 1;
                        }
                    }
                }
                nbond == natm // 原子的数量等于键的数量，说明是小环
            })
            .cloned()
            .collect();
        rings
    }

    pub fn equal(&self, other: &Mole) -> bool {
        if self.atoms().equal(&other.atoms()) == false {
            log::info!("分子结构不相等");
            return false;
        }
        if self.basis().equal(&other.basis()) == false {
            log::info!("分子基组不相等");
            return false;
        }
        if self.coefs().equal(&other.coefs()) == false {
            log::info!("分子系数不相等");
            return false;
        }
        return true;
    }

    pub fn border(&self) -> [[f64; 3]; 2] {
        let xyzs = self.atoms().xyzs();
        let xs: Vec<f64> = xyzs.iter().map(|xyz| xyz[0]).collect();
        let ys: Vec<f64> = xyzs.iter().map(|xyz| xyz[1]).collect();
        let zs: Vec<f64> = xyzs.iter().map(|xyz| xyz[2]).collect();
        let minx = maths::arg_min(&xs).1;
        let miny = maths::arg_min(&ys).1;
        let minz = maths::arg_min(&zs).1;
        let maxx = maths::arg_max(&xs).1;
        let maxy = maths::arg_max(&ys).1;
        let maxz = maths::arg_max(&zs).1;
        let p0 = [minx, miny, minz];
        let p1 = [maxx, maxy, maxz];
        [p0, p1]
    }
}

impl fmt::Display for Mole {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let charge = self.charge();
        let multip = self.multip();
        writeln!(f, "{:>3}{:>3}", charge, multip)?;
        write!(f, "{}", self.atoms())?;
        write!(f, "{}", self.basis())?;
        write!(f, "{}", self.coefs())?;
        Ok(())
    }
}

pub fn find_sssr(adjacency_list: &HashMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
    let n = adjacency_list.len();
    let mut rings: Vec<Vec<usize>> = Vec::new();

    fn dfs_find_rings(current: usize, start: usize, path: &mut Vec<usize>, visited: &mut HashSet<usize>, adjacency_list: &HashMap<usize, Vec<usize>>, rings: &mut Vec<Vec<usize>>) {
        path.push(current);
        visited.insert(current);

        if let Some(neighbors) = adjacency_list.get(&current) {
            for &neighbor in neighbors {
                if neighbor == start && path.len() >= 3 {
                    // 找到环
                    let ring = path.clone();
                    if is_new_ring(&ring, rings) {
                        rings.push(ring);
                    }
                } else if !visited.contains(&neighbor) {
                    dfs_find_rings(neighbor, start, path, visited, adjacency_list, rings);
                }
            }
        }

        path.pop();
        visited.remove(&current);
    }

    fn is_new_ring(candidate_ring: &[usize], existing_rings: &[Vec<usize>]) -> bool {
        let candidate_set: HashSet<usize> = candidate_ring.iter().cloned().collect();
        for ring in existing_rings {
            let ring_set: HashSet<usize> = ring.iter().cloned().collect();
            if candidate_set == ring_set {
                return false;
            }
        }
        true
    }

    for start_atom in 0..n as usize {
        let mut path = Vec::new();
        let mut visited = HashSet::new();
        dfs_find_rings(start_atom, start_atom, &mut path, &mut visited, adjacency_list, &mut rings);
    }

    rings
}
