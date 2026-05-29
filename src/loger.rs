use simplelog;

#[rustfmt::skip]
pub fn init_log() {
    let file = std::fs::File::create("rswfn日志.txt").expect("创建日志文件失败");
    simplelog::CombinedLogger::init(vec![
        simplelog::TermLogger::new(
            simplelog::LevelFilter::Info,
            simplelog::Config::default(),
            simplelog::TerminalMode::Mixed,
            simplelog::ColorChoice::Auto,
        ),
        simplelog::WriteLogger::new(
            simplelog::LevelFilter::Info,
            simplelog::Config::default(),
            file,
        ),
    ])
    .expect("无法创建日志系统");
    log::info!("日志系统已启动");
}
