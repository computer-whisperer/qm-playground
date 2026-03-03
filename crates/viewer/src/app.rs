use egui::*;
use sim_core::potential::Potential;
use sim_core::Simulation;

/// Preset scenarios for exploring QM phenomena.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Scenario {
    FreeParticle,
    InfiniteWell,
    TunnelingBarrier,
    HarmonicOscillator,
    DoubleWell,
}

impl Scenario {
    const ALL: &[Scenario] = &[
        Scenario::FreeParticle,
        Scenario::InfiniteWell,
        Scenario::TunnelingBarrier,
        Scenario::HarmonicOscillator,
        Scenario::DoubleWell,
    ];

    fn label(&self) -> &'static str {
        match self {
            Scenario::FreeParticle => "Free Particle",
            Scenario::InfiniteWell => "Infinite Well",
            Scenario::TunnelingBarrier => "Tunneling (Barrier)",
            Scenario::HarmonicOscillator => "Harmonic Oscillator",
            Scenario::DoubleWell => "Double Well",
        }
    }

    fn description(&self) -> &'static str {
        match self {
            Scenario::FreeParticle => "Gaussian wave packet in free space. Watch it spread — this IS the uncertainty principle.",
            Scenario::InfiniteWell => "Particle confined between walls. Energy is quantized: only certain wavelengths fit.",
            Scenario::TunnelingBarrier => "Wave packet hits a barrier taller than its energy. Part reflects, part tunnels through.",
            Scenario::HarmonicOscillator => "Quantum spring. The ground state Gaussian is an eigenstate — it doesn't spread.",
            Scenario::DoubleWell => "Two wells separated by a barrier. The particle tunnels back and forth between them.",
        }
    }
}

pub struct App {
    sim: Simulation,
    scenario: Scenario,
    running: bool,
    steps_per_frame: usize,

    // Initial condition parameters
    x0: f64,
    sigma: f64,
    k0: f64,

    // Display toggles
    show_real: bool,
    show_imag: bool,
    show_probability: bool,
    show_potential: bool,

    // Grid/sim parameters
    grid_points: usize,
    dt: f64,
    x_range: f64,
}

impl App {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        let scenario = Scenario::TunnelingBarrier;
        let grid_points = 1024;
        let x_range = 30.0;
        let dt = 0.005;

        let x0 = -8.0;
        let sigma = 1.5;
        let k0 = 5.0;

        let mut app = Self {
            sim: Self::build_sim(scenario, grid_points, x_range, dt),
            scenario,
            running: false,
            steps_per_frame: 10,
            x0,
            sigma,
            k0,
            show_real: false,
            show_imag: false,
            show_probability: true,
            show_potential: true,
            grid_points,
            dt,
            x_range,
        };
        app.reset_wavefunction();
        app
    }

    fn build_sim(scenario: Scenario, n: usize, x_range: f64, dt: f64) -> Simulation {
        let potential = Self::potential_for(scenario);
        Simulation::new(n, -x_range, x_range, dt, potential)
    }

    fn potential_for(scenario: Scenario) -> Potential {
        match scenario {
            Scenario::FreeParticle => Potential::free(),
            Scenario::InfiniteWell => Potential::well(-10.0, 10.0, 1000.0),
            Scenario::TunnelingBarrier => Potential::barrier(-0.5, 0.5, 15.0),
            Scenario::HarmonicOscillator => Potential::harmonic(0.0, 1.0),
            Scenario::DoubleWell => Potential::double_well(-4.0, 4.0, 1.0, 3.0),
        }
    }

    fn reset_wavefunction(&mut self) {
        self.sim = Self::build_sim(self.scenario, self.grid_points, self.x_range, self.dt);
        self.sim.wf.set_gaussian(self.x0, self.sigma, self.k0);
        self.sim.time = 0.0;
    }

    fn defaults_for_scenario(&mut self) {
        match self.scenario {
            Scenario::FreeParticle => {
                self.x0 = 0.0;
                self.sigma = 2.0;
                self.k0 = 3.0;
            }
            Scenario::InfiniteWell => {
                self.x0 = 0.0;
                self.sigma = 2.0;
                self.k0 = 0.0;
            }
            Scenario::TunnelingBarrier => {
                self.x0 = -8.0;
                self.sigma = 1.5;
                self.k0 = 5.0;
            }
            Scenario::HarmonicOscillator => {
                self.x0 = -3.0;
                self.sigma = 1.0;
                self.k0 = 0.0;
            }
            Scenario::DoubleWell => {
                self.x0 = -4.0;
                self.sigma = 1.0;
                self.k0 = 0.0;
            }
        }
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        // Step simulation if running
        if self.running {
            self.sim.step_n(self.steps_per_frame);
            ctx.request_repaint();
        }

        // Left panel: controls
        SidePanel::left("controls").min_width(250.0).show(ctx, |ui| {
            ui.heading("Scenario");
            for &s in Scenario::ALL {
                if ui
                    .selectable_label(self.scenario == s, s.label())
                    .clicked()
                {
                    self.scenario = s;
                    self.defaults_for_scenario();
                    self.reset_wavefunction();
                    self.running = false;
                }
            }

            ui.separator();
            ui.label(self.scenario.description());

            ui.separator();
            ui.heading("Initial State");

            ui.horizontal(|ui| {
                ui.label("x₀:");
                ui.add(Slider::new(&mut self.x0, -self.x_range..=self.x_range).step_by(0.1));
            });
            ui.horizontal(|ui| {
                ui.label("σ:");
                ui.add(Slider::new(&mut self.sigma, 0.1..=5.0).step_by(0.1));
            });
            ui.horizontal(|ui| {
                ui.label("k₀:");
                ui.add(Slider::new(&mut self.k0, -10.0..=10.0).step_by(0.1));
            });

            if ui.button("Reset").clicked() {
                self.reset_wavefunction();
                self.running = false;
            }

            ui.separator();
            ui.heading("Simulation");

            ui.horizontal(|ui| {
                if ui
                    .button(if self.running { "⏸ Pause" } else { "▶ Run" })
                    .clicked()
                {
                    self.running = !self.running;
                }
                if ui.button("⏭ Step").clicked() {
                    self.sim.step_n(self.steps_per_frame);
                }
            });

            ui.horizontal(|ui| {
                ui.label("Steps/frame:");
                ui.add(Slider::new(&mut self.steps_per_frame, 1..=100).logarithmic(true));
            });
            ui.horizontal(|ui| {
                ui.label("dt:");
                let old_dt = self.dt;
                ui.add(Slider::new(&mut self.dt, 0.0001..=0.05).logarithmic(true));
                if (self.dt - old_dt).abs() > f64::EPSILON {
                    self.sim.dt = self.dt;
                    self.sim.rebuild_integrator();
                }
            });

            ui.separator();
            ui.heading("Display");
            ui.checkbox(&mut self.show_probability, "|ψ|² probability density");
            ui.checkbox(&mut self.show_real, "Re(ψ)");
            ui.checkbox(&mut self.show_imag, "Im(ψ)");
            ui.checkbox(&mut self.show_potential, "V(x) potential");

            ui.separator();
            ui.heading("Observables");
            ui.label(format!("t = {:.3}", self.sim.time));
            ui.label(format!("‖ψ‖² = {:.6}", self.sim.norm()));
            ui.label(format!("⟨x⟩ = {:.3}", self.sim.expected_x()));
            ui.label(format!("⟨p⟩ = {:.3}", self.sim.expected_p()));
            ui.label(format!("⟨E⟩ = {:.3}", self.sim.expected_energy()));
        });

        // Central panel: plot
        CentralPanel::default().show(ctx, |ui| {
            self.draw_plot(ui);
        });
    }
}

impl App {
    fn draw_plot(&self, ui: &mut Ui) {
        let rect = ui.available_rect_before_wrap();
        let painter = ui.painter_at(rect);

        // Dark background
        painter.rect_filled(rect, 0.0, Color32::from_rgb(20, 22, 28));

        let xs = self.sim.wf.xs();
        let density = self.sim.wf.probability_density();
        let potential_vals = self.sim.potential.values_on_grid(&xs);

        // Find Y scale for wavefunction
        let max_density = density.iter().cloned().fold(0.0_f64, f64::max).max(0.01);
        let max_psi = self
            .sim
            .wf
            .psi
            .iter()
            .map(|c| c.re.abs().max(c.im.abs()))
            .fold(0.0_f64, f64::max)
            .max(0.01);
        let y_scale = max_density.max(max_psi) * 1.2;

        // Potential scale: normalize to fit in the top portion of the plot
        let max_potential = potential_vals
            .iter()
            .cloned()
            .filter(|v| *v < 1e6) // ignore "infinite" walls for scaling
            .fold(0.0_f64, f64::max)
            .max(0.01);
        let potential_scale = y_scale / max_potential * 0.8;

        let margin = 40.0;
        let plot_rect = Rect::from_min_max(
            pos2(rect.min.x + margin, rect.min.y + margin),
            pos2(rect.max.x - margin, rect.max.y - margin),
        );

        // Map from simulation coordinates to screen coordinates
        let x_to_screen = |x: f64| -> f32 {
            let t = (x - self.sim.wf.x_min) / (self.sim.wf.x_max - self.sim.wf.x_min);
            plot_rect.min.x + t as f32 * plot_rect.width()
        };
        let y_to_screen = |y: f64| -> f32 {
            let t = y / y_scale;
            plot_rect.max.y - t as f32 * plot_rect.height() * 0.5 - plot_rect.height() * 0.15
        };

        // Draw grid lines
        let zero_y = y_to_screen(0.0);
        painter.line_segment(
            [pos2(plot_rect.min.x, zero_y), pos2(plot_rect.max.x, zero_y)],
            Stroke::new(1.0, Color32::from_gray(60)),
        );

        // Draw potential
        if self.show_potential {
            let points: Vec<Pos2> = xs
                .iter()
                .zip(potential_vals.iter())
                .map(|(&x, &v)| {
                    let v_clamped = (v * potential_scale).min(y_scale);
                    pos2(x_to_screen(x), y_to_screen(v_clamped))
                })
                .collect();

            // Fill under potential curve
            let mut fill_points = vec![pos2(points[0].x, zero_y)];
            fill_points.extend_from_slice(&points);
            fill_points.push(pos2(points.last().unwrap().x, zero_y));

            painter.add(Shape::convex_polygon(
                fill_points,
                Color32::from_rgba_premultiplied(80, 80, 120, 30),
                Stroke::NONE,
            ));

            if points.len() >= 2 {
                painter.add(Shape::line(
                    points,
                    Stroke::new(1.5, Color32::from_rgb(100, 100, 160)),
                ));
            }
        }

        // Draw |ψ|²
        if self.show_probability {
            let points: Vec<Pos2> = xs
                .iter()
                .zip(density.iter())
                .map(|(&x, &d)| pos2(x_to_screen(x), y_to_screen(d)))
                .collect();

            // Fill under curve
            let mut fill_points = vec![pos2(points[0].x, zero_y)];
            fill_points.extend_from_slice(&points);
            fill_points.push(pos2(points.last().unwrap().x, zero_y));

            painter.add(Shape::convex_polygon(
                fill_points,
                Color32::from_rgba_premultiplied(0, 120, 255, 40),
                Stroke::NONE,
            ));

            if points.len() >= 2 {
                painter.add(Shape::line(
                    points,
                    Stroke::new(2.0, Color32::from_rgb(40, 160, 255)),
                ));
            }
        }

        // Draw Re(ψ)
        if self.show_real {
            let points: Vec<Pos2> = xs
                .iter()
                .enumerate()
                .map(|(i, &x)| pos2(x_to_screen(x), y_to_screen(self.sim.wf.psi[i].re)))
                .collect();
            if points.len() >= 2 {
                painter.add(Shape::line(
                    points,
                    Stroke::new(1.0, Color32::from_rgb(255, 100, 100)),
                ));
            }
        }

        // Draw Im(ψ)
        if self.show_imag {
            let points: Vec<Pos2> = xs
                .iter()
                .enumerate()
                .map(|(i, &x)| pos2(x_to_screen(x), y_to_screen(self.sim.wf.psi[i].im)))
                .collect();
            if points.len() >= 2 {
                painter.add(Shape::line(
                    points,
                    Stroke::new(1.0, Color32::from_rgb(100, 255, 100)),
                ));
            }
        }

        // Axis labels
        painter.text(
            pos2(plot_rect.center().x, plot_rect.max.y + 20.0),
            Align2::CENTER_TOP,
            "x (Bohr radii)",
            FontId::proportional(14.0),
            Color32::from_gray(160),
        );

        // Legend
        let mut legend_y = plot_rect.min.y + 5.0;
        let legend_x = plot_rect.max.x - 160.0;
        let legend_spacing = 18.0;

        if self.show_probability {
            painter.text(
                pos2(legend_x, legend_y),
                Align2::LEFT_TOP,
                "— |ψ|² probability",
                FontId::proportional(12.0),
                Color32::from_rgb(40, 160, 255),
            );
            legend_y += legend_spacing;
        }
        if self.show_real {
            painter.text(
                pos2(legend_x, legend_y),
                Align2::LEFT_TOP,
                "— Re(ψ)",
                FontId::proportional(12.0),
                Color32::from_rgb(255, 100, 100),
            );
            legend_y += legend_spacing;
        }
        if self.show_imag {
            painter.text(
                pos2(legend_x, legend_y),
                Align2::LEFT_TOP,
                "— Im(ψ)",
                FontId::proportional(12.0),
                Color32::from_rgb(100, 255, 100),
            );
            legend_y += legend_spacing;
        }
        if self.show_potential {
            painter.text(
                pos2(legend_x, legend_y),
                Align2::LEFT_TOP,
                "— V(x) potential",
                FontId::proportional(12.0),
                Color32::from_rgb(100, 100, 160),
            );
        }

        // Allocate the rect so egui knows it's used
        ui.allocate_rect(rect, Sense::hover());
    }
}
