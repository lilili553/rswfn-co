use std::io;

use crate::base::Mole;

use crate::atomprop;
use crate::bondprop;
use crate::gridprop;
use crate::writer::CubWriter;

fn input(tip: &str) -> String {
    let mut txt = String::new();
    io::stdin().read_line(&mut txt).expect("无法读取输入");
    txt
}

/// TUI 主界面，提供交互式功能选择
pub struct Main {
    path: String,
}

impl Main {
    /// 创建新的 TUI 主界面并启动用户交互
    pub fn new() -> Self {
        println!("欢迎使用rswfn: 基于rust的波函数分析工具");
        let path = input("请输入文件路径: ");
        Self { path }
    }

    /// 进入原子属性子菜单
    pub fn atom_prop() {
        AtomProp {}.main();
    }

    /// 进入化学键属性子菜单
    pub fn bond_prop() {}

    /// 进入格点属性子菜单
    pub fn grid_prop() {}
}

/// 原子属性计算器
struct AtomProp {}

impl AtomProp {
    /// 进入原子属性交互菜单，循环处理用户输入
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

    /// 计算原子电荷/电子布局
    pub fn atom_charge(&self) {}
}

struct BondProp {}

/// 格点属性计算器
struct GridProp {}
impl GridProp {
    /// 进入格点属性交互菜单，支持输出分子轨道和电子密度到cube文件
    pub fn main(&self, mole: &Mole) {
        let mut grid_caler = gridprop::CubeGrid::new();
        let [p0, p1] = mole.border(4.0);
        grid_caler.set_v1(&p0, &p1, 0.2, 0.0);
        let (shape, grids) = grid_caler.get();
        println!("1.分子轨道");
        println!("2.电子密度");
        let opt = input("请输入选项: ");
        let vals = match opt.trim() {
            "1" => {
                let obt: usize = input("输入轨道序号: ").parse().expect("轨道序号错误");
                let caler = gridprop::wfnfunc::Calculator::new(mole);
                let vals = grids.iter().map(|grid| caler.obt_wfn(obt, grid, 0).val0).collect();
                vals
            }
            "2" => {
                let caler = gridprop::density::Calculator::new(mole);
                let vals: Vec<f64> = grids.iter().map(|grid| caler.mol_rho_cm(grid, 0).val0).collect();
                vals
            }
            _ => return,
        };
        let mut cub_writer = CubWriter::new();
        cub_writer.pos0 = p0;
        cub_writer.step = [0.2; 3];
        cub_writer.size = shape;
        cub_writer.xyzs = mole.atoms().xyzs().clone();
        cub_writer.syms = mole.atoms().syms().clone();
        cub_writer.obts = vec![-1];
        cub_writer.vals = Vec::from([vals]);
        cub_writer.save("new.cub");
    }
}
