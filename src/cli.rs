use std::io;

use crate::base::Mole;

use crate::atomprop;
use crate::bondprop;
use crate::gridprop;
use crate::writer::CubWriter;

fn input(tip: &str) -> String {
    print!("{}", tip);
    let mut txt = String::new();
    io::stdin().read_line(&mut txt).expect("无法读取输入");
    txt
}

pub struct Main {
    path: String,
}

impl Main {
    pub fn new() -> Self {
        println!("欢迎使用rswfn: 基于rust的波函数分析工具");
        let path = input("请输入文件路径: ");
        Self { path }
    }

    pub fn atom_prop() {
        AtomProp {}.main();
    }

    pub fn bond_prop() {}

    pub fn grid_prop() {}
}

struct AtomProp {}

impl AtomProp {
    pub fn main(&self) {
        loop {
            println!("1.电子布局");
            println!("2.原子活性");
            let opt = input("请输入选项: ");
            match opt.trim() {
                "1" => self.atom_charge(),
                "2" => {}
                _ => break,
            }
        }
    }

    pub fn atom_charge(&self) {}
}

struct BondProp {}

struct GridProp {}
impl GridProp {
    pub fn main(&self, mole: &Mole) {
        let mut grid_caler = gridprop::CubeGrid::new();
        let [p0, p1] = mole.border();
        grid_caler.set_v1(&p0, &p1, 0.2, 4.0);
        let (shape, grids) = grid_caler.get();
        println!("1.分子轨道");
        println!("2.电子密度");
        let opt = input("请输入选项: ");
        let vals = match opt.trim() {
            "1" => {
                let obt: usize = input("输入轨道序号: ").parse().expect("轨道序号错误");
                let caler = gridprop::wfnfunc::Calculator::new(mole);
                let vals = grids.iter().map(|grid| caler.obt_wfn(obt, grid, 0).0).collect();
                vals
            }
            "2" => {
                let caler = gridprop::density::Calculator::new(mole);
                let vals: Vec<f64> = grids.iter().map(|grid| caler.mol_rho_cm(grid, 0).0).collect();
                vals
            }
            _ => return,
        };
        let mut cub_writer = CubWriter::new("new.cub");
        cub_writer.pos0 = p0;
        cub_writer.step = [0.2; 3];
        cub_writer.size = shape;
        cub_writer.xyzs = mole.atoms().xyzs().clone();
        cub_writer.syms = mole.atoms().syms().clone();
        cub_writer.obts = vec![-1];
        cub_writer.vals = Vec::from([vals]);
        cub_writer.save();
    }
}
