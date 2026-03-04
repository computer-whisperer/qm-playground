mod app;
mod colormap;

fn main() -> eframe::Result {
    env_logger::init();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_title("Emergence Sim — Quantum Mechanics"),
        ..Default::default()
    };

    eframe::run_native(
        "emergence-sim",
        options,
        Box::new(|cc| Ok(Box::new(app::App::new(cc)))),
    )
}
