/// Fch格式文件(.fch)写入器
use ndarray::{Array2, s};
use std::collections::HashMap;

use crate::base::{Atoms, Basis, Coefs, Mole};
use crate::datas::elems;
use crate::utils;
use crate::writer::Writer;

pub struct FchWriter {
    /// 分子电荷
    pub charge: isize,
    /// 自旋多重度
    pub multip: usize,
    /// 原子信息
    pub atoms: Atoms,
    /// 基组信息
    pub basis: Basis,
    /// 波函数信息
    pub coefs: Coefs,
}

impl FchWriter {
    /// 创建FchWriter实例
    pub fn new() -> Self {
        Self {
            charge: 0,
            multip: 1,
            atoms: Atoms::new(&vec![], &vec![], &vec![]),
            basis: Basis::new(vec![]),
            coefs: Coefs::new(vec![], vec![], vec![], vec![], vec![], Array2::zeros((0, 0))),
        }
    }

    /// 从Mole构造FchWriter
    pub fn from_mole(mole: &Mole) -> Self {
        Self {
            charge: mole.charge(),
            multip: mole.multip(),
            atoms: mole.atoms().clone(),
            basis: mole.basis().clone(),
            coefs: mole.coefs().clone(),
        }
    }

    #[rustfmt::skip]
    /// 生成fch文件内容字符串
    pub fn build(&self)->String {
        let mut text = String::new();
        let atoms=&self.atoms;
        let natm=atoms.len();
        let nele=self.coefs.obt_eles(&self.basis);
        let syms=atoms.syms();
        let xyzs=atoms.xyzs();
        let cmat=self.coefs.get_cmat("sph", &self.basis);
        let [nato,nobt]=cmat.shape().try_into().expect("获取分子轨道形状失败");
        // let nbas=self.basis.map_ato().len();
        let nbas=nato;

        text.push_str("generate by rswfn");
        text.push_str(&format!("Number of atoms                            I{:>17}\n",natm));
        text.push_str(&format!("Charge                                     I{:>17}\n",self.charge));
        text.push_str(&format!("Multiplicity                               I{:>17}\n",self.multip));
        text.push_str(&format!("Number of electrons                        I{:>17}\n",nele[0] + nele[1]));
        text.push_str(&format!("Number of alpha electrons                  I{:>17}\n",nele[0]));
        text.push_str(&format!("Number of beta electrons                   I{:>17}\n",nele[1]));
        text.push_str(&format!("Number of basis functions                  I{:>17}\n",nbas));
        text.push_str(&format!("Number of independent functions            I{:>17}\n",nobt));
        text.push_str(&format!("Atomic numbers                             I   N={:>12}\n",natm));
        let atms:Vec<isize>=syms.iter().map(|sym|elems::sym2atm(sym) as isize).collect();
        text.push_str(&i2str(&atms));
        text.push_str(&format!("Nuclear charges                            R   N={:>12}\n",natm));
        let crgs:Vec<f64>=atms.iter().map(|atm|*atm as f64).collect();
        text.push_str(&f2str(&crgs));
        text.push_str(&format!("Current cartesian coordinates              R   N={:>12}\n",natm*3));
        let xyz_vec:Vec<f64>=xyzs.iter().flat_map(|xyz|xyz.iter().copied()).collect();
        text.push_str(&f2str(&xyz_vec));

        text.push_str("Translation vectors                        R   N=           0\n");
        // 基组信息
        let (nshl,shl_types,shl_nbass,shl_atoms,pm_alps,pm_coes,sp_coes)=basis_info(&self.basis);
        let shl_nbass=shl_nbass.iter().map(|v|*v as isize).collect();
        let shl_atoms=shl_atoms.iter().map(|v|*v as isize).collect();
        text.push_str(&format!("Shell types                                I   N={:>12}\n",nshl));
        text.push_str(&i2str(&shl_types));
        text.push_str(&format!("Number of primitives per shell             I   N={:>12}\n",nshl));
        text.push_str(&i2str(&shl_nbass));
        text.push_str(&format!("Shell to atom map                          I   N={:>12}\n",nshl));
        text.push_str(&i2str(&shl_atoms));
        text.push_str(&format!("Primitive exponents                        R   N={:>12}\n",pm_alps.len()));
        text.push_str(&f2str(&pm_alps));
        text.push_str(&format!("Contraction coefficients                   R   N={:>12}\n",pm_coes.len()));
        text.push_str(&f2str(&pm_coes));
        text.push_str(&format!("P(S=P) Contraction coefficients            R   N={:>12}\n",sp_coes.len()));
        text.push_str(&f2str(&sp_coes));
        text.push_str(&format!("Coordinates of each shell                  R   N={:>12}\n",nshl*3));
        let xyz_vec:Vec<f64>=shl_atoms.iter().flat_map(|atm|xyzs[*atm as usize-1].iter().copied()).collect();
        text.push_str(&f2str(&xyz_vec));
        text.push_str(&format!("Virial Ratio                               R{:>27.15E}\n",0.0));
        text.push_str(&format!("Total Energy                               R{:>27.15E}\n",0.0));

        // 分子轨道信息
        let obt_engs=&self.coefs.obt_engs;
        text.push_str(&format!("Alpha Orbital Energies                     R   N={:>12}\n",obt_engs.len()));
        text.push_str(&f2str(obt_engs));
        let nato=cmat.nrows();
        let nobt=cmat.ncols();
        if nobt<=nato{
            let coefs=cmat.t().flatten().to_vec();
            text.push_str(&format!("Alpha MO coefficients                      R   N={:>12}\n",coefs.len()));
            text.push_str(&f2str(&coefs));
        }else{
            let cmat_a=cmat.slice(s![..,..nobt/2]).to_owned();
            let cmat_b=cmat.slice(s![..,nobt/2..]).to_owned();
            let coefs_a=cmat_a.t().flatten().to_vec();
            let coefs_b=cmat_b.t().flatten().to_vec();
            text.push_str(&format!("Alpha MO coefficients                      R   N={:>12}\n",coefs_a.len()));
            text.push_str(&f2str(&coefs_a));
            text.push_str(&format!("Beta MO coefficients                       R   N={:>12}\n",coefs_b.len()));
            text.push_str(&f2str(&coefs_b));
        }

        let dmat=self.coefs.get_dmat("sph", &self.basis);
        let mut dvec=Vec::new();
        for i in 0..nato{
            for j in 0..i{
                dvec.push(dmat[(i,j)]);
            }
        }
        text.push_str(&format!("Total SCF Density                          R   N={:>12}\n",dvec.len()));
        text.push_str(&f2str(&dvec));
        text
    }

    /// 保存到文件
    pub fn save(&self, path: &str) {
        let text = self.build();
        std::fs::write(path, text).expect("Unable to write file");
    }
}

fn i2str(nums: &Vec<isize>) -> String {
    let mut text = String::new();
    for (i, num) in nums.iter().enumerate() {
        text.push_str(&format!("{:>12}", num));
        if (i + 1) % 6 == 0 {
            text.push_str("\n");
        }
    }
    if nums.len() % 6 != 0 {
        text.push_str("\n");
    }
    text
}

fn f2str(nums: &Vec<f64>) -> String {
    let mut text = String::new();
    for (i, num) in nums.iter().enumerate() {
        let fort = utils::fort_form(*num, 8);
        text.push_str(&format!("{:>16}", fort));

        if (i + 1) % 5 == 0 {
            text.push_str("\n");
        }
    }
    if nums.len() % 5 != 0 {
        text.push_str("\n");
    }
    text
}

struct ShlInfo {
    atom: usize,
    alps: Vec<f64>,
    coes: Vec<f64>,
    pm_alps: Vec<f64>,
    pm_coes: Vec<f64>,
    sp_coes: Vec<f64>,
    angs: Vec<usize>,
}

impl ShlInfo {
    pub fn blank() -> Self {
        Self {
            atom: 0,
            alps: Vec::new(),
            coes: Vec::new(),
            pm_alps: Vec::new(),
            pm_coes: Vec::new(),
            sp_coes: Vec::new(),
            angs: Vec::new(),
        }
    }
}

fn basis_info(basis: &Basis) -> (usize, Vec<isize>, Vec<usize>, Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut shl_types = Vec::new();
    let mut shl_nbass = Vec::new();
    let mut shl_atoms = Vec::new();
    // let mut keys=Vec::new();
    let mut nbas: usize = 1;
    let mut pm_alps = Vec::new();
    let mut pm_coes = Vec::new();
    let mut sp_coes = Vec::new();
    let nshl = basis.shls.len();
    for i in 0..nshl {
        let shl = basis.shls.get(i);
        shl_types.push(shl.ang as isize);
        shl_nbass.push(shl.alps.len());
        pm_alps.extend(shl.alps.clone());
        pm_coes.extend(shl.coes.clone());
    }
    shl_nbass.push(nbas);
    (nshl, shl_types, shl_nbass, shl_atoms, pm_alps, pm_coes, sp_coes)
}

impl Writer for FchWriter {
    fn read_mole(&mut self, mole: &crate::base::Mole) {
        self.atoms = mole.atoms().clone();
        self.basis = mole.basis().clone();
        self.coefs = mole.coefs().clone();
    }
}

fn u2i(nums: &mut Vec<usize>) {
    for num in nums {
        *num as isize;
    }
}
