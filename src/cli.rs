use clap::{Parser, Subcommand, command};
use rayon::prelude::*;
use std::path::{PathBuf, absolute};

use crate::base::Mole;
use crate::dev;
use crate::gridprop::{CubeGrid, density, potential, topolog};
use crate::reader::{FchReader, GjfReader, LogReader};
use crate::writer::{CubWriter, DmoWriter, GjfWriter, Writer};
use crate::{atomprop, qmcalc, utils};

#[derive(Parser, Debug)]
#[command(author = "Xiaofei Shi", about = "Weave Function Analysis Tool")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// 原子属性分析
    #[command(subcommand)]
    Atomprop(AtompropCmds),
    /// 键属性分析
    #[command(subcommand)]
    Bondprop(BondpropCmds),
    /// 分子属性分析
    #[command(subcommand)]
    Moleprop(MolepropCmds),
    /// 轨道属性分析
    #[command(subcommand)]
    Orbtprop(OrbtpropCmds),
    /// 格点属性分析，导出 cub 文件
    #[command(subcommand)]
    Gridprop(GridpropCmds),

    /// 文件格式转换
    Trans { input: PathBuf, output: PathBuf },
    /// DHMO 计算
    Dmo { input: PathBuf, output: Option<PathBuf> },
    /// 开发时测试
    Dev,
}

// ===== atomprop 子命令 =====

#[derive(Subcommand, Debug)]
enum AtompropCmds {
    /// 原子电荷/布居分析
    #[command(subcommand)]
    Charge(ChargeCmds),
}

#[derive(Subcommand, Debug)]
enum ChargeCmds {
    /// Mulliken 布居分析
    Mulliken {
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk/.gjf）
        filepath: PathBuf,
    },
    /// Löwdin 布居分析
    Lowdin {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// Hirshfeld 电荷分析
    Hirshfeld {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// 自旋密度布居分析
    Spin {
        /// 电荷类型: mulliken | lowdin | hirshfeld
        #[arg(default_value = "mulliken")]
        ctype: String,
        /// 波函数文件路径
        filepath: PathBuf,
    },
}

// ===== moleprop 子命令 =====

#[derive(Subcommand, Debug)]
enum MolepropCmds {
    /// 芳香性分析
    Aromacity {
        /// 波函数文件路径
        filepath: PathBuf,
    },
}

// ===== orbtprop 子命令 =====

#[derive(Subcommand, Debug)]
enum OrbtpropCmds {
    /// 分子轨道 π 电子数（POCV 方法）
    PiElesPocv {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// 分子轨道 π 电子数（MOCV 方法）
    PiElesMocv {
        /// 波函数文件路径
        filepath: PathBuf,
    },
}

// ===== bondprop 子命令 =====

#[derive(Subcommand, Debug)]
enum BondpropCmds {
    /// 键级分析
    #[command(subcommand)]
    Order(OrderCmds),
}

#[derive(Subcommand, Debug)]
enum OrderCmds {
    /// Mayer 键级
    Mayer {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// Wiberg 键级
    Wiberg {
        /// 波函数文件路径
        filepath: PathBuf,
    },
}

// ===== gridprop 子命令 =====

#[derive(Subcommand, Debug)]
enum GridpropCmds {
    /// 电子密度，输出 cub 文件
    Density {
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk）
        filepath: PathBuf,
        /// 输出 cub 文件路径（默认: 输入文件名_density.cub）
        output: Option<PathBuf>,
        /// 格点步长，单位 Bohr（默认 0.2）
        #[arg(default_value = "0.2", long = "step")]
        step: f64,
    },
    /// 静电势，输出 cub 文件
    Potential {
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk）
        filepath: PathBuf,
        /// 输出 cub 文件路径（默认: 输入文件名_potential.cub）
        output: Option<PathBuf>,
        /// 格点步长，单位 Bohr（默认 0.2）
        #[arg(default_value = "0.2", long = "step")]
        step: f64,
    },
    /// 拓扑分析，寻找临界点
    Topology {
        /// 函数名（如 mol_rho_dm, pi_rho）
        #[arg(short = 'f')]
        func: String,
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk）
        filepath: PathBuf,
    },
}

// ===== 读取分子文件 =====

/// 根据文件路径自动选择合适的读取器加载分子
fn load_mole(filepath: &PathBuf) -> Mole {
    let path = absolute(filepath).expect("无法获取文件绝对路径");
    let path_str = path.to_str().expect("路径转字符串失败");
    let ftype = utils::file_type(path_str);
    match ftype.as_str() {
        "log" | "out" => Mole::from_read(&LogReader::new(path_str)),
        "gjf" => Mole::from_read(&GjfReader::new(path_str)),
        "fch" | "fchk" => Mole::from_read(&FchReader::new(path_str)),
        _ => panic!("不支持的波函数文件类型: {ftype}，支持: log/out/gjf/fch/fchk"),
    }
}

/// 打印原子布居结果表
fn print_atom_population(mole: &Mole, pops: &[f64], title: &str) {
    let atoms = mole.atoms();
    println!("\n{}", title);
    println!("{}", "=".repeat(title.len()));
    println!("{:>4}  {:>6}  {:>14}", "原子", "符号", "电子布居");
    println!("{:-<28}", "");
    for (i, &pop) in pops.iter().enumerate() {
        println!("{:>4}  {:>6}  {:>14.8}", i + 1, atoms.atoms[i].sym, pop);
    }
    let total: f64 = pops.iter().sum();
    println!("{:-<28}", "");
    println!("{:>12}  {:>14.8}", "总电子数", total);
}

/// 打印键级矩阵（仅上三角，筛选显著键级）
fn print_bond_order(mole: &Mole, omat: &ndarray::Array2<f64>, title: &str) {
    let natm = mole.atoms().len();
    let atoms = mole.atoms();
    println!("\n{}", title);
    println!("{}", "=".repeat(title.len()));
    println!("{:>6}  {:>6}  {:>12}", "原子对", "序号", "键级");
    println!("{:-<30}", "");
    for i in 0..natm {
        for j in (i + 1)..natm {
            let order = omat[(i, j)];
            if order.abs() > 0.01 {
                println!(
                    "{:>6}  {:>3}-{:<3}  {:>12.6}",
                    format!("{}-{}", atoms.atoms[i].sym, atoms.atoms[j].sym),
                    i + 1,
                    j + 1,
                    order
                );
            }
        }
    }
}

/// 根据分子坐标计算含边框的立方体空间范围（单位: Bohr）
fn molecule_bounding_box(mole: &Mole, border: f64) -> ([f64; 3], [f64; 3]) {
    let xyzs = mole.atoms().xyzs();
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for xyz in xyzs.iter() {
        for i in 0..3 {
            if xyz[i] < min[i] {
                min[i] = xyz[i];
            }
            if xyz[i] > max[i] {
                max[i] = xyz[i];
            }
        }
    }
    let p0 = min.map(|v| v - border);
    let p1 = max.map(|v| v + border);
    (p0, p1)
}

/// 在立方体格点上计算标量场值，并写入 cub 文件
fn compute_and_write_cube<F>(mole: &Mole, output_path: &str, step: f64, title: &str, field_calc: F)
where
    F: Fn(&[f64; 3]) -> f64 + Send + Sync,
{
    // 计算分子边界 + 4 Bohr 边框
    let border = 4.0;
    let (p0, p1) = molecule_bounding_box(mole, border);

    // 创建立方体格点
    let mut cube_grid = CubeGrid::new();
    cube_grid.set_v1(&p0, &p1, step, 0.0);
    let (shape, grids) = cube_grid.get();
    let [nx, ny, nz] = shape;
    let start = cube_grid.start;
    let steps = cube_grid.steps;

    let npts = grids.len();
    println!("格点数: {} ({} x {} x {}), 步长: {:.3} Bohr", npts, nx, ny, nz, step);
    println!(
        "空间范围: [{:.3}, {:.3}] x [{:.3}, {:.3}] x [{:.3}, {:.3}] (Bohr)",
        p0[0], p1[0], p0[1], p1[1], p0[2], p1[2]
    );
    println!("正在计算 {} ...", title);

    // 使用 rayon 并行计算每个格点的值
    let vals: Vec<f64> = grids.par_iter().map(|grid| field_calc(grid)).collect();

    // 写入 cub 文件
    let mut writer = CubWriter::new();
    writer.read_mole(mole);
    writer.title = title.to_string();
    writer.pos0 = start;
    writer.size = [nx, ny, nz];
    writer.step = steps;
    writer.obts = vec![0];
    writer.vals = vec![vals];
    writer.save(output_path);

    println!("cub 文件已保存: {}", output_path);
}

// ===== 主入口 =====

/// CLI 模式入口函数，解析命令行参数执行对应操作
pub fn main() {
    let cli = Cli::parse();
    match &cli.command {
        // ==================== Atomprop ====================
        Commands::Atomprop(cmd) => match cmd {
            AtompropCmds::Charge(cmd) => match cmd {
                ChargeCmds::Mulliken { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = atomprop::charge::Calculator::new(&mole);
                    let pops = caler.mulliken();
                    print_atom_population(&mole, &pops, "Mulliken 电子布居分析");
                }
                ChargeCmds::Lowdin { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = atomprop::charge::Calculator::new(&mole);
                    let pops = caler.lowdin();
                    print_atom_population(&mole, &pops, "Löwdin 电子布居分析");
                }
                ChargeCmds::Hirshfeld { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = atomprop::charge::Calculator::new(&mole);
                    let pops = caler.hirshfeld();
                    print_atom_population(&mole, &pops, "Hirshfeld 电子布居分析");
                }
                ChargeCmds::Spin { ctype, filepath } => {
                    let mole = load_mole(filepath);
                    if !mole.is_open() {
                        panic!("非开壳层分子（多重度=1）无法计算自旋布居");
                    }
                    let caler = atomprop::charge::Calculator::new(&mole);
                    let spin = caler.spin(ctype);
                    println!("\n自旋密度布居分析（ctype: {ctype}）");
                    println!("{}", "=".repeat(58));
                    println!(
                        "{:>4}  {:>6}  {:>14}  {:>14}  {:>14}",
                        "原子", "符号", "α电子", "β电子", "α - β"
                    );
                    println!("{:-<58}", "");
                    for (i, &[a, b, d]) in spin.iter().enumerate() {
                        println!(
                            "{:>4}  {:>6}  {:>14.8}  {:>14.8}  {:>14.8}",
                            i + 1,
                            mole.atoms().atoms[i].sym,
                            a,
                            b,
                            d
                        );
                    }
                }
            },
        },

        // ==================== Moleprop ====================
        Commands::Moleprop(cmd) => match cmd {
            MolepropCmds::Aromacity { filepath } => {
                let mole = load_mole(filepath);
                let caler = crate::moleprop::aromacity::Calculator::new(&mole);
                // 芳香性分析的具体方法待补全
                println!("\n芳香性分析（待完善）");
            }
        },

        // ==================== Orbtprop ====================
        Commands::Orbtprop(cmd) => match cmd {
            OrbtpropCmds::PiElesPocv { filepath } => {
                let mole = load_mole(filepath);
                let caler = crate::orbtprop::obtval::Calculator::new(&mole);
                let pi_eles = caler.pi_eles_pocv();
                println!("\n分子轨道 π 电子数（POCV 方法）");
                println!("{}", "=".repeat(40));
                for (i, &val) in pi_eles.iter().enumerate() {
                    println!("轨道 {:>4}: {:>12.8}", i + 1, val);
                }
            }
            OrbtpropCmds::PiElesMocv { filepath } => {
                let mole = load_mole(filepath);
                let caler = crate::orbtprop::obtval::Calculator::new(&mole);
                let pi_eles = caler.pi_eles_mocv();
                println!("\n分子轨道 π 电子数（MOCV 方法）");
                println!("{}", "=".repeat(40));
                for (i, &val) in pi_eles.iter().enumerate() {
                    println!("轨道 {:>4}: {:>12.8}", i + 1, val);
                }
            }
        },

        // ==================== Bondprop ====================
        Commands::Bondprop(cmd) => match cmd {
            BondpropCmds::Order(cmd) => match cmd {
                OrderCmds::Mayer { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = crate::bondprop::order::Calculator::new(&mole);
                    let omat = caler.mayer();
                    print_bond_order(&mole, &omat, "Mayer 键级分析");
                }
                OrderCmds::Wiberg { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = crate::bondprop::order::Calculator::new(&mole);
                    let omat = caler.wiberg();
                    print_bond_order(&mole, &omat, "Wiberg 键级分析");
                }
            },
        },

        // ==================== Gridprop ====================
        Commands::Gridprop(cmd) => match cmd {
            GridpropCmds::Density { filepath, output, step } => {
                let mole = load_mole(filepath);
                let caler = density::Calculator::new(&mole);
                // 确定输出路径
                let output = match output {
                    Some(path) => path.to_str().expect("输出路径转换失败").to_owned(),
                    None => {
                        let path = absolute(filepath).expect("无法获取文件绝对路径");
                        let path_str = path.to_str().expect("路径转字符串失败");
                        let name = utils::file_name(path_str);
                        format!("{}_density.cub", name)
                    }
                };
                compute_and_write_cube(&mole, &output, *step, "Electron Density", |grid| {
                    caler.mol_rho_dm(grid, 0).val0
                });
            }
            GridpropCmds::Topology { func, filepath } => {
                let mole = load_mole(filepath);
                let caler = topolog::Calculator::new(&mole);
                let cpts = caler.analysis(func);
                println!("\n拓扑分析结果（函数: {func}）");
                println!("{:=<90}", "");
                println!(
                    "{:>4}  {:>20}  {:>20}  {:>20}  {:>12}  {:>8}",
                    "序号", "X (Bohr)", "Y (Bohr)", "Z (Bohr)", "函数值", "临界点类型"
                );
                println!("{:-<90}", "");
                for (i, (pos, val, cpt)) in cpts.iter().enumerate() {
                    let cpt_name = match cpt {
                        -3 => "(3, -3) 核临界点",
                        -1 => "(3, -1) 键临界点",
                        1 => "(3, +1) 环临界点",
                        3 => "(3, +3) 笼临界点",
                        _ => "其他",
                    };
                    println!(
                        "{:>4}  {:>20.12}  {:>20.12}  {:>20.12}  {:>12.8}  {:>8}",
                        i + 1,
                        pos[0],
                        pos[1],
                        pos[2],
                        val,
                        cpt_name,
                    );
                }
                println!("{:=<90}", "");
            }
            GridpropCmds::Potential { filepath, output, step } => {
                let mole = load_mole(filepath);
                let caler = potential::Calculator::new(&mole);
                // 确定输出路径
                let output = match output {
                    Some(path) => path.to_str().expect("输出路径转换失败").to_owned(),
                    None => {
                        let path = absolute(filepath).expect("无法获取文件绝对路径");
                        let path_str = path.to_str().expect("路径转字符串失败");
                        let name = utils::file_name(path_str);
                        format!("{}_potential.cub", name)
                    }
                };
                compute_and_write_cube(&mole, &output, *step, "Electrostatic Potential", |grid| {
                    let nuc = caler.nuc_potential(grid);
                    let ele = caler.ele_potential(grid);
                    ele - nuc
                });
            }
        },

        // ==================== 原有命令 ====================
        Commands::Trans { input, output } => {
            let from = absolute(input).expect("输入文件转换失败");
            let dest = absolute(output).expect("输出文件转换失败");
            let from_path = from.to_str().expect("准换字符串失败");
            let dest_path = dest.to_str().expect("准换字符串失败");
            let from_type = utils::file_type(from_path);
            let dest_type = utils::file_type(dest_path);
            let dest_name = utils::file_name(dest_path);
            println!("{} -> {}", from_type, dest_type);
            match (from_type.as_str(), dest_type.as_str()) {
                ("log" | "out", "gjf") => {
                    let reader = LogReader::new(&from_path);
                    let mole = Mole::from_read(&reader);
                    let job = reader.read_job();
                    let mut writer = GjfWriter::new();
                    writer.read_mole(&mole);
                    writer.job = job.to_owned();
                    writer.chk = dest_name;
                    writer.charge = reader.charge as isize;
                    writer.multip = reader.multip as usize;
                    writer.save(&dest_path);
                    println!("文件转换成功 >_< ");
                }
                ("gjf", "dmi") => {
                    let reader = GjfReader::new(&from_path);
                    let mole = Mole::from_read(&reader);
                    let mut writer = DmoWriter::blank();
                    writer.read_mole(&mole);
                    writer.save(&dest_path);
                    println!("文件转换成功 >_< ");
                }
                _ => panic!("不支持的转换类型!!"),
            }
        }
        Commands::Dmo { input, output } => {
            let from = absolute(input)
                .expect("输入文件转换失败")
                .to_str()
                .expect("转换字符串失败")
                .to_owned();
            let dest = match output {
                Some(output) => {
                    let dest = absolute(output).expect("输出文件转换失败");
                    let dest = dest.to_str().expect("转换字符串失败").to_owned();
                    Some(dest)
                }
                None => None,
            };
            qmcalc::dhmo::output(&from, dest);
            println!("DHMO计算完成 >_< ");
        }
        Commands::Dev => {
            dev::test();
        }
    }
}
