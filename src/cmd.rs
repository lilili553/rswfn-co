use clap::{Parser, Subcommand, command};
use std::path::{PathBuf, absolute};

use crate::base::Mole;
use crate::reader::{GjfReader, LogReader};
use crate::tests;
use crate::writer::{DmoWriter, GjfWriter, Writer};
use crate::{qmcalc, utils};

#[derive(Parser, Debug)]
#[command(author = "石小飞", about = "波函数分析程序")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Trans { input: PathBuf, output: PathBuf },
    Dmo { input: PathBuf, output: Option<PathBuf> },
    Test,
}

pub fn main() {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Trans { input, output } => {
            let from = absolute(input).expect("输入文件转换失败");
            let dest = absolute(output).expect("输出文件转换失败");
            let from_path = from.to_str().expect("准换字符串失败");
            let dest_path = dest.to_str().expect("准换字符串失败");
            let from_type = utils::file_type(from_path);
            let dest_type = utils::file_type(dest_path);
            // let from_name = utils::file_name(from_path);
            let dest_name = utils::file_name(dest_path);
            println!("{} -> {}", from_type, dest_type);
            match (from_type.as_str(), dest_type.as_str()) {
                ("log" | "out", "gjf") => {
                    let reader = LogReader::new(&from_path);
                    let mole = Mole::from_read(&reader);
                    let job = reader.read_job();
                    let mut writer = GjfWriter::new();
                    writer.from_mole(&mole);
                    writer.job = job;
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
                    writer.from_mole(&mole);
                    writer.save(&dest_path);
                    println!("文件转换成功 >_< ");
                }
                _ => panic!("不支持的转换类型!!"),
            }
        }
        Commands::Dmo { input, output } => {
            let from = absolute(input).expect("输入文件转换失败").to_str().expect("转换字符串失败").to_owned();
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
        Commands::Test => {
            tests::test();
        }
    }
}
