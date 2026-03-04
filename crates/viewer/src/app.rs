use egui::*;
use sim_core::potential::Potential;
use sim_core::two_particle::potential::{Interaction, Potential2D};
use sim_core::two_particle::wavefunction::ParticleSymmetry;
use sim_core::two_particle::TwoParticleSimulation;
use sim_core::Simulation;

use crate::colormap;

/// Preset scenarios for exploring QM phenomena.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Scenario {
    // Single-particle (1D)
    FreeParticle,
    InfiniteWell,
    TunnelingBarrier,
    HarmonicOscillator,
    DoubleWell,
    // Two-particle (2D)
    TwoParticleFree,
    TwoParticleScattering,
    TwoParticleTrapped,
    // Identical particles (Chapter 3)
    BosonBunching,
    FermionAntibunching,
    PauliExclusion,
}

impl Scenario {
    const ONE_PARTICLE: &[Scenario] = &[
        Scenario::FreeParticle,
        Scenario::InfiniteWell,
        Scenario::TunnelingBarrier,
        Scenario::HarmonicOscillator,
        Scenario::DoubleWell,
    ];

    const TWO_PARTICLE: &[Scenario] = &[
        Scenario::TwoParticleFree,
        Scenario::TwoParticleScattering,
        Scenario::TwoParticleTrapped,
    ];

    const IDENTICAL: &[Scenario] = &[
        Scenario::BosonBunching,
        Scenario::FermionAntibunching,
        Scenario::PauliExclusion,
    ];

    fn label(&self) -> &'static str {
        match self {
            Scenario::FreeParticle => "Free Particle",
            Scenario::InfiniteWell => "Infinite Well",
            Scenario::TunnelingBarrier => "Tunneling (Barrier)",
            Scenario::HarmonicOscillator => "Harmonic Oscillator",
            Scenario::DoubleWell => "Double Well",
            Scenario::TwoParticleFree => "Free (no interaction)",
            Scenario::TwoParticleScattering => "Scattering",
            Scenario::TwoParticleTrapped => "Trapped + interaction",
            Scenario::BosonBunching => "Boson bunching",
            Scenario::FermionAntibunching => "Fermion antibunching",
            Scenario::PauliExclusion => "Pauli exclusion",
        }
    }

    fn description(&self) -> &'static str {
        match self {
            Scenario::FreeParticle => "Gaussian wave packet in free space. Watch it spread — this IS the uncertainty principle.",
            Scenario::InfiniteWell => "Particle confined between walls. Energy is quantized: only certain wavelengths fit.",
            Scenario::TunnelingBarrier => "Wave packet hits a barrier taller than its energy. Part reflects, part tunnels through.",
            Scenario::HarmonicOscillator => "Quantum spring. The ground state Gaussian is an eigenstate — it doesn't spread.",
            Scenario::DoubleWell => "Two wells separated by a barrier. The particle tunnels back and forth between them.",
            Scenario::TwoParticleFree => "Two particles with opposite momenta, no interaction. Stays a product state — no entanglement.",
            Scenario::TwoParticleScattering => "Approaching particles repel via soft-Coulomb. Watch entanglement develop as they scatter.",
            Scenario::TwoParticleTrapped => "Both particles in a harmonic well with interaction. Bound-state entanglement builds up.",
            Scenario::BosonBunching => "Two bosons in a harmonic well. Symmetric wavefunction enhances probability on the x₁=x₂ diagonal — bosons cluster.",
            Scenario::FermionAntibunching => "Two fermions in a harmonic well. Antisymmetric wavefunction has a node on x₁=x₂ — fermions repel without any force.",
            Scenario::PauliExclusion => "Two fermions initialized in the SAME state. Antisymmetry forces ψ=0 everywhere — the state cannot exist.",
        }
    }

    fn is_two_particle(&self) -> bool {
        matches!(
            self,
            Scenario::TwoParticleFree
                | Scenario::TwoParticleScattering
                | Scenario::TwoParticleTrapped
                | Scenario::BosonBunching
                | Scenario::FermionAntibunching
                | Scenario::PauliExclusion
        )
    }
}

enum SimMode {
    OneParticle(Simulation),
    TwoParticle(TwoParticleSimulation),
}

pub struct App {
    mode: SimMode,
    scenario: Scenario,
    running: bool,
    steps_per_frame: usize,

    // Particle 1 / single-particle initial condition
    x0: f64,
    sigma: f64,
    k0: f64,

    // Particle 2 initial condition (two-particle mode)
    x0_2: f64,
    sigma_2: f64,
    k0_2: f64,

    // Interaction parameters
    interaction_g: f64,
    interaction_epsilon: f64,

    // Particle symmetry (for identical particle scenarios)
    symmetry: ParticleSymmetry,

    // Display toggles
    show_real: bool,
    show_imag: bool,
    show_probability: bool,
    show_potential: bool,

    // Grid/sim parameters
    grid_points: usize,
    grid_points_2p: usize,
    dt: f64,
    x_range: f64,

    // Two-particle visualization cache
    purity_cache: f64,
    purity_frame_counter: usize,
}

impl App {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        let scenario = Scenario::TunnelingBarrier;
        let grid_points = 1024;
        let grid_points_2p = 256;
        let x_range = 30.0;
        let dt = 0.005;

        let mut app = Self {
            mode: SimMode::OneParticle(Self::build_1p_sim(scenario, grid_points, x_range, dt)),
            scenario,
            running: false,
            steps_per_frame: 10,
            x0: -8.0,
            sigma: 1.5,
            k0: 5.0,
            x0_2: 0.0,
            sigma_2: 1.5,
            k0_2: 0.0,
            interaction_g: 1.0,
            interaction_epsilon: 0.5,
            symmetry: ParticleSymmetry::Distinguishable,
            show_real: false,
            show_imag: false,
            show_probability: true,
            show_potential: true,
            grid_points,
            grid_points_2p,
            dt,
            x_range,
            purity_cache: 1.0,
            purity_frame_counter: 0,
        };
        app.defaults_for_scenario();
        app.reset_wavefunction();
        app
    }

    fn build_1p_sim(scenario: Scenario, n: usize, x_range: f64, dt: f64) -> Simulation {
        let potential = Self::potential_1p_for(scenario);
        Simulation::new(n, -x_range, x_range, dt, potential)
    }

    fn potential_1p_for(scenario: Scenario) -> Potential {
        match scenario {
            Scenario::FreeParticle => Potential::free(),
            Scenario::InfiniteWell => Potential::well(-10.0, 10.0, 1000.0),
            Scenario::TunnelingBarrier => Potential::barrier(-0.5, 0.5, 15.0),
            Scenario::HarmonicOscillator => Potential::harmonic(0.0, 1.0),
            Scenario::DoubleWell => Potential::double_well(-4.0, 4.0, 1.0, 3.0),
            _ => Potential::free(),
        }
    }

    fn build_2p_sim(&self, scenario: Scenario) -> TwoParticleSimulation {
        let (v1, v2, interaction) = match scenario {
            Scenario::TwoParticleFree => (
                Potential::free(),
                Potential::free(),
                Interaction::None,
            ),
            Scenario::TwoParticleScattering => (
                Potential::free(),
                Potential::free(),
                Interaction::SoftCoulomb {
                    g: self.interaction_g,
                    epsilon: self.interaction_epsilon,
                },
            ),
            Scenario::TwoParticleTrapped => (
                Potential::harmonic(0.0, 1.0),
                Potential::harmonic(0.0, 1.0),
                Interaction::SoftCoulomb {
                    g: self.interaction_g,
                    epsilon: self.interaction_epsilon,
                },
            ),
            Scenario::BosonBunching | Scenario::FermionAntibunching | Scenario::PauliExclusion => (
                Potential::harmonic(0.0, 0.5),
                Potential::harmonic(0.0, 0.5),
                Interaction::SoftCoulomb {
                    g: self.interaction_g,
                    epsilon: self.interaction_epsilon,
                },
            ),
            _ => (Potential::free(), Potential::free(), Interaction::None),
        };
        let pot = Potential2D::new(v1, v2, interaction);
        TwoParticleSimulation::new(self.grid_points_2p, -self.x_range, self.x_range, self.dt, pot)
    }

    fn reset_wavefunction(&mut self) {
        if self.scenario.is_two_particle() {
            let mut sim = self.build_2p_sim(self.scenario);
            sim.wf.set_symmetrized_gaussian(
                self.x0, self.sigma, self.k0,
                self.x0_2, self.sigma_2, self.k0_2,
                self.symmetry,
            );
            self.purity_cache = if self.symmetry == ParticleSymmetry::Distinguishable {
                1.0
            } else if sim.norm() > 1e-10 {
                sim.purity()
            } else {
                0.0
            };
            self.purity_frame_counter = 0;
            self.mode = SimMode::TwoParticle(sim);
        } else {
            let mut sim =
                Self::build_1p_sim(self.scenario, self.grid_points, self.x_range, self.dt);
            sim.wf.set_gaussian(self.x0, self.sigma, self.k0);
            self.mode = SimMode::OneParticle(sim);
        }
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
            Scenario::TwoParticleFree => {
                self.x0 = -5.0;
                self.sigma = 1.5;
                self.k0 = 3.0;
                self.x0_2 = 5.0;
                self.sigma_2 = 1.5;
                self.k0_2 = -3.0;
                self.interaction_g = 0.0;
                self.interaction_epsilon = 0.5;
                self.symmetry = ParticleSymmetry::Distinguishable;
            }
            Scenario::TwoParticleScattering => {
                self.x0 = -8.0;
                self.sigma = 1.5;
                self.k0 = 3.0;
                self.x0_2 = 8.0;
                self.sigma_2 = 1.5;
                self.k0_2 = -3.0;
                self.interaction_g = 1.0;
                self.interaction_epsilon = 0.5;
                self.symmetry = ParticleSymmetry::Distinguishable;
            }
            Scenario::TwoParticleTrapped => {
                self.x0 = -3.0;
                self.sigma = 1.0;
                self.k0 = 0.0;
                self.x0_2 = 3.0;
                self.sigma_2 = 1.0;
                self.k0_2 = 0.0;
                self.interaction_g = 1.0;
                self.interaction_epsilon = 0.5;
                self.symmetry = ParticleSymmetry::Distinguishable;
            }
            Scenario::BosonBunching => {
                self.x0 = -3.0;
                self.sigma = 1.0;
                self.k0 = 0.0;
                self.x0_2 = 3.0;
                self.sigma_2 = 1.0;
                self.k0_2 = 0.0;
                self.interaction_g = 0.5;
                self.interaction_epsilon = 0.5;
                self.symmetry = ParticleSymmetry::Boson;
            }
            Scenario::FermionAntibunching => {
                self.x0 = -3.0;
                self.sigma = 1.0;
                self.k0 = 0.0;
                self.x0_2 = 3.0;
                self.sigma_2 = 1.0;
                self.k0_2 = 0.0;
                self.interaction_g = 0.5;
                self.interaction_epsilon = 0.5;
                self.symmetry = ParticleSymmetry::Fermion;
            }
            Scenario::PauliExclusion => {
                self.x0 = 0.0;
                self.sigma = 1.5;
                self.k0 = 0.0;
                self.x0_2 = 0.0;
                self.sigma_2 = 1.5;
                self.k0_2 = 0.0;
                self.interaction_g = 0.0;
                self.interaction_epsilon = 0.5;
                self.symmetry = ParticleSymmetry::Fermion;
            }
        }
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        // Step simulation if running
        if self.running {
            match &mut self.mode {
                SimMode::OneParticle(sim) => sim.step_n(self.steps_per_frame),
                SimMode::TwoParticle(sim) => sim.step_n(self.steps_per_frame),
            }
            // Update purity cache periodically
            if self.scenario.is_two_particle() {
                self.purity_frame_counter += 1;
                if self.purity_frame_counter >= 10 {
                    if let SimMode::TwoParticle(sim) = &self.mode {
                        self.purity_cache = sim.purity();
                    }
                    self.purity_frame_counter = 0;
                }
            }
            ctx.request_repaint();
        }

        // Left panel: controls
        SidePanel::left("controls").min_width(250.0).show(ctx, |ui| {
            self.draw_controls(ui);
        });

        // Central panel: plot
        CentralPanel::default().show(ctx, |ui| {
            match &self.mode {
                SimMode::OneParticle(_) => self.draw_1p_plot(ui),
                SimMode::TwoParticle(_) => self.draw_2p_plot(ui, ctx),
            }
        });
    }
}

// ---- Controls ----

impl App {
    fn draw_controls(&mut self, ui: &mut Ui) {
        ui.heading("One Particle");
        for &s in Scenario::ONE_PARTICLE {
            if ui.selectable_label(self.scenario == s, s.label()).clicked() {
                self.scenario = s;
                self.defaults_for_scenario();
                self.reset_wavefunction();
                self.running = false;
            }
        }

        ui.separator();
        ui.heading("Two Particles");
        for &s in Scenario::TWO_PARTICLE {
            if ui.selectable_label(self.scenario == s, s.label()).clicked() {
                self.scenario = s;
                self.defaults_for_scenario();
                self.reset_wavefunction();
                self.running = false;
            }
        }

        ui.separator();
        ui.heading("Identical Particles");
        for &s in Scenario::IDENTICAL {
            if ui.selectable_label(self.scenario == s, s.label()).clicked() {
                self.scenario = s;
                self.defaults_for_scenario();
                self.reset_wavefunction();
                self.running = false;
            }
        }

        ui.separator();
        ui.label(self.scenario.description());

        ui.separator();
        if self.scenario.is_two_particle() {
            self.draw_2p_controls(ui);
        } else {
            self.draw_1p_controls(ui);
        }

        ui.separator();
        ui.heading("Simulation");
        self.draw_sim_controls(ui);

        if !self.scenario.is_two_particle() {
            ui.separator();
            ui.heading("Display");
            ui.checkbox(&mut self.show_probability, "|ψ|² probability density");
            ui.checkbox(&mut self.show_real, "Re(ψ)");
            ui.checkbox(&mut self.show_imag, "Im(ψ)");
            ui.checkbox(&mut self.show_potential, "V(x) potential");
        }

        ui.separator();
        ui.heading("Observables");
        self.draw_observables(ui);
    }

    fn draw_1p_controls(&mut self, ui: &mut Ui) {
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
    }

    fn draw_2p_controls(&mut self, ui: &mut Ui) {
        ui.heading("Particle 1");
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

        ui.separator();
        ui.heading("Particle 2");
        ui.horizontal(|ui| {
            ui.label("x₀:");
            ui.add(Slider::new(&mut self.x0_2, -self.x_range..=self.x_range).step_by(0.1));
        });
        ui.horizontal(|ui| {
            ui.label("σ:");
            ui.add(Slider::new(&mut self.sigma_2, 0.1..=5.0).step_by(0.1));
        });
        ui.horizontal(|ui| {
            ui.label("k₀:");
            ui.add(Slider::new(&mut self.k0_2, -10.0..=10.0).step_by(0.1));
        });

        if !matches!(self.scenario, Scenario::TwoParticleFree) {
            ui.separator();
            ui.heading("Interaction");
            ui.horizontal(|ui| {
                ui.label("g:");
                ui.add(Slider::new(&mut self.interaction_g, 0.0..=5.0).step_by(0.1));
            });
            ui.horizontal(|ui| {
                ui.label("ε:");
                ui.add(Slider::new(&mut self.interaction_epsilon, 0.1..=3.0).step_by(0.1));
            });
        }

        if ui.button("Reset").clicked() {
            self.reset_wavefunction();
            self.running = false;
        }
    }

    fn draw_sim_controls(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            if ui
                .button(if self.running { "⏸ Pause" } else { "▶ Run" })
                .clicked()
            {
                self.running = !self.running;
            }
            if ui.button("⏭ Step").clicked() {
                match &mut self.mode {
                    SimMode::OneParticle(sim) => sim.step_n(self.steps_per_frame),
                    SimMode::TwoParticle(sim) => sim.step_n(self.steps_per_frame),
                }
                if self.scenario.is_two_particle() {
                    if let SimMode::TwoParticle(sim) = &self.mode {
                        self.purity_cache = sim.purity();
                    }
                }
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
                match &mut self.mode {
                    SimMode::OneParticle(sim) => {
                        sim.dt = self.dt;
                        sim.rebuild_integrator();
                    }
                    SimMode::TwoParticle(sim) => {
                        sim.dt = self.dt;
                        sim.rebuild_integrator();
                    }
                }
            }
        });
    }

    fn draw_observables(&self, ui: &mut Ui) {
        match &self.mode {
            SimMode::OneParticle(sim) => {
                ui.label(format!("t = {:.3}", sim.time));
                ui.label(format!("‖ψ‖² = {:.6}", sim.norm()));
                ui.label(format!("⟨x⟩ = {:.3}", sim.expected_x()));
                ui.label(format!("⟨p⟩ = {:.3}", sim.expected_p()));
                ui.label(format!("⟨E⟩ = {:.3}", sim.expected_energy()));
            }
            SimMode::TwoParticle(sim) => {
                if self.symmetry != ParticleSymmetry::Distinguishable {
                    let sym_label = match self.symmetry {
                        ParticleSymmetry::Boson => "Bosons (symmetric)",
                        ParticleSymmetry::Fermion => "Fermions (antisymmetric)",
                        _ => "",
                    };
                    ui.label(sym_label);
                }
                ui.label(format!("t = {:.3}", sim.time));
                ui.label(format!("‖ψ‖² = {:.6}", sim.norm()));
                ui.label(format!("⟨x₁⟩ = {:.3}", sim.expected_x1()));
                ui.label(format!("⟨x₂⟩ = {:.3}", sim.expected_x2()));
                ui.label(format!("Purity = {:.4}", self.purity_cache));
                ui.label(format!(
                    "Schmidt K = {:.2}",
                    1.0 / self.purity_cache.max(1e-10)
                ));
            }
        }
    }
}

// ---- 1D Plot ----

impl App {
    fn draw_1p_plot(&self, ui: &mut Ui) {
        let sim = match &self.mode {
            SimMode::OneParticle(s) => s,
            _ => return,
        };

        let rect = ui.available_rect_before_wrap();
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, Color32::from_rgb(20, 22, 28));

        let xs = sim.wf.xs();
        let density = sim.wf.probability_density();
        let potential_vals = sim.potential.values_on_grid(&xs);

        let max_density = density.iter().cloned().fold(0.0_f64, f64::max).max(0.01);
        let max_psi = sim
            .wf
            .psi
            .iter()
            .map(|c| c.re.abs().max(c.im.abs()))
            .fold(0.0_f64, f64::max)
            .max(0.01);
        let y_scale = max_density.max(max_psi) * 1.2;

        let max_potential = potential_vals
            .iter()
            .cloned()
            .filter(|v| *v < 1e6)
            .fold(0.0_f64, f64::max)
            .max(0.01);
        let potential_scale = y_scale / max_potential * 0.8;

        let margin = 40.0;
        let plot_rect = Rect::from_min_max(
            pos2(rect.min.x + margin, rect.min.y + margin),
            pos2(rect.max.x - margin, rect.max.y - margin),
        );

        let x_to_screen = |x: f64| -> f32 {
            let t = (x - sim.wf.x_min) / (sim.wf.x_max - sim.wf.x_min);
            plot_rect.min.x + t as f32 * plot_rect.width()
        };
        let y_to_screen = |y: f64| -> f32 {
            let t = y / y_scale;
            plot_rect.max.y - t as f32 * plot_rect.height() * 0.5 - plot_rect.height() * 0.15
        };

        let zero_y = y_to_screen(0.0);
        painter.line_segment(
            [pos2(plot_rect.min.x, zero_y), pos2(plot_rect.max.x, zero_y)],
            Stroke::new(1.0, Color32::from_gray(60)),
        );

        if self.show_potential {
            let points: Vec<Pos2> = xs
                .iter()
                .zip(potential_vals.iter())
                .map(|(&x, &v)| {
                    let v_clamped = (v * potential_scale).min(y_scale);
                    pos2(x_to_screen(x), y_to_screen(v_clamped))
                })
                .collect();
            fill_under_curve(
                &painter,
                &points,
                zero_y,
                Color32::from_rgba_premultiplied(80, 80, 120, 30),
            );
            if points.len() >= 2 {
                painter.add(Shape::line(
                    points,
                    Stroke::new(1.5, Color32::from_rgb(100, 100, 160)),
                ));
            }
        }

        if self.show_probability {
            let points: Vec<Pos2> = xs
                .iter()
                .zip(density.iter())
                .map(|(&x, &d)| pos2(x_to_screen(x), y_to_screen(d)))
                .collect();
            fill_under_curve(
                &painter,
                &points,
                zero_y,
                Color32::from_rgba_premultiplied(0, 120, 255, 40),
            );
            if points.len() >= 2 {
                painter.add(Shape::line(
                    points,
                    Stroke::new(2.0, Color32::from_rgb(40, 160, 255)),
                ));
            }
        }

        if self.show_real {
            let points: Vec<Pos2> = xs
                .iter()
                .enumerate()
                .map(|(i, &x)| pos2(x_to_screen(x), y_to_screen(sim.wf.psi[i].re)))
                .collect();
            if points.len() >= 2 {
                painter.add(Shape::line(
                    points,
                    Stroke::new(1.0, Color32::from_rgb(255, 100, 100)),
                ));
            }
        }

        if self.show_imag {
            let points: Vec<Pos2> = xs
                .iter()
                .enumerate()
                .map(|(i, &x)| pos2(x_to_screen(x), y_to_screen(sim.wf.psi[i].im)))
                .collect();
            if points.len() >= 2 {
                painter.add(Shape::line(
                    points,
                    Stroke::new(1.0, Color32::from_rgb(100, 255, 100)),
                ));
            }
        }

        // Axis label
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

        let legend_entry = |legend_y: &mut f32, text: &str, color: Color32| {
            painter.text(
                pos2(legend_x, *legend_y),
                Align2::LEFT_TOP,
                text,
                FontId::proportional(12.0),
                color,
            );
            *legend_y += legend_spacing;
        };

        if self.show_probability {
            legend_entry(&mut legend_y, "— |ψ|² probability", Color32::from_rgb(40, 160, 255));
        }
        if self.show_real {
            legend_entry(&mut legend_y, "— Re(ψ)", Color32::from_rgb(255, 100, 100));
        }
        if self.show_imag {
            legend_entry(&mut legend_y, "— Im(ψ)", Color32::from_rgb(100, 255, 100));
        }
        if self.show_potential {
            legend_entry(&mut legend_y, "— V(x) potential", Color32::from_rgb(100, 100, 160));
        }

        ui.allocate_rect(rect, Sense::hover());
    }
}

// ---- 2D Heatmap + Marginals ----

impl App {
    fn draw_2p_plot(&self, ui: &mut Ui, ctx: &Context) {
        let sim = match &self.mode {
            SimMode::TwoParticle(s) => s,
            _ => return,
        };

        let rect = ui.available_rect_before_wrap();
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, Color32::from_rgb(20, 22, 28));

        let n = sim.wf.n;
        let margin = 40.0;
        let marginal_h = 80.0; // height for top marginal plot
        let marginal_w = 80.0; // width for right marginal plot

        // Compute heatmap region (square)
        let avail_w = rect.width() - 2.0 * margin - marginal_w;
        let avail_h = rect.height() - 2.0 * margin - marginal_h;
        let heatmap_size = avail_w.min(avail_h).max(100.0);

        let heatmap_left = rect.min.x + margin;
        let heatmap_top = rect.min.y + margin + marginal_h;
        let heatmap_rect = Rect::from_min_size(
            pos2(heatmap_left, heatmap_top),
            vec2(heatmap_size, heatmap_size),
        );

        // Compute density for heatmap
        let density = sim.wf.probability_density();
        let max_density = density.iter().cloned().fold(0.0_f64, f64::max).max(1e-10);

        // Build RGBA pixel data for the heatmap texture
        // Image row r, col c: x-axis = x₁ (col), y-axis = x₂ (row, flipped)
        let mut pixels = vec![0u8; n * n * 4];
        for r in 0..n {
            for c in 0..n {
                let i1 = c; // x₁ index = column
                let i2 = n - 1 - r; // x₂ index = flipped row (y up)
                let val = density[i1 * n + i2] / max_density;
                let [cr, cg, cb] = colormap::inferno(val.sqrt()); // sqrt for better dynamic range
                let idx = (r * n + c) * 4;
                pixels[idx] = cr;
                pixels[idx + 1] = cg;
                pixels[idx + 2] = cb;
                pixels[idx + 3] = 255;
            }
        }

        let image = ColorImage::from_rgba_unmultiplied([n, n], &pixels);
        let texture = ctx.load_texture("heatmap_2p", image, TextureOptions::LINEAR);

        // Draw heatmap
        painter.image(
            texture.id(),
            heatmap_rect,
            Rect::from_min_max(pos2(0.0, 0.0), pos2(1.0, 1.0)),
            Color32::WHITE,
        );

        // Axis labels
        painter.text(
            pos2(heatmap_rect.center().x, heatmap_rect.max.y + 20.0),
            Align2::CENTER_TOP,
            "x₁ (Bohr radii)",
            FontId::proportional(14.0),
            Color32::from_gray(160),
        );

        // Rotated y-axis label (draw vertically with individual characters)
        let label = "x₂";
        painter.text(
            pos2(heatmap_rect.min.x - 25.0, heatmap_rect.center().y),
            Align2::CENTER_CENTER,
            label,
            FontId::proportional(14.0),
            Color32::from_gray(160),
        );

        // x₁ = x₂ diagonal — where particles occupy the same position
        painter.line_segment(
            [
                pos2(heatmap_rect.min.x, heatmap_rect.max.y),  // bottom-left: x₁=x_min, x₂=x_min
                pos2(heatmap_rect.max.x, heatmap_rect.min.y),  // top-right: x₁=x_max, x₂=x_max
            ],
            Stroke::new(1.0, Color32::from_rgba_premultiplied(255, 255, 255, 80)),
        );

        // Draw marginal plots
        let m1 = sim.wf.marginal_1();
        let m2 = sim.wf.marginal_2();

        // Top marginal: ρ₁(x₁) — same x range as heatmap
        let top_rect = Rect::from_min_max(
            pos2(heatmap_rect.min.x, rect.min.y + margin),
            pos2(heatmap_rect.max.x, heatmap_rect.min.y - 4.0),
        );
        self.draw_marginal_top(&painter, &m1, top_rect);

        // Right marginal: ρ₂(x₂) — same y range as heatmap
        let right_rect = Rect::from_min_max(
            pos2(heatmap_rect.max.x + 4.0, heatmap_rect.min.y),
            pos2(heatmap_rect.max.x + 4.0 + marginal_w, heatmap_rect.max.y),
        );
        self.draw_marginal_right(&painter, &m2, right_rect);

        ui.allocate_rect(rect, Sense::hover());
    }

    /// Draw marginal ρ₁(x₁) as a filled curve in the given rect (x-axis matches heatmap).
    fn draw_marginal_top(&self, painter: &Painter, marginal: &[f64], rect: Rect) {
        let n = marginal.len();
        if n < 2 {
            return;
        }
        let max_val = marginal.iter().cloned().fold(0.0_f64, f64::max).max(1e-10);

        let points: Vec<Pos2> = (0..n)
            .map(|i| {
                let t = (i as f32 + 0.5) / n as f32;
                let x = rect.min.x + t * rect.width();
                let y = rect.max.y - (marginal[i] / max_val) as f32 * rect.height() * 0.9;
                pos2(x, y)
            })
            .collect();

        fill_under_curve(
            painter,
            &points,
            rect.max.y,
            Color32::from_rgba_premultiplied(0, 120, 255, 40),
        );
        painter.add(Shape::line(
            points,
            Stroke::new(1.5, Color32::from_rgb(40, 160, 255)),
        ));
    }

    /// Draw marginal ρ₂(x₂) as a filled curve in the given rect (y-axis matches heatmap, density extends right).
    fn draw_marginal_right(&self, painter: &Painter, marginal: &[f64], rect: Rect) {
        let n = marginal.len();
        if n < 2 {
            return;
        }
        let max_val = marginal.iter().cloned().fold(0.0_f64, f64::max).max(1e-10);

        // y-axis: x₂ increases upward (row 0 = top = x₂_max)
        let points: Vec<Pos2> = (0..n)
            .map(|i| {
                let i2 = n - 1 - i; // flip so top = x₂_max
                let t = (i as f32 + 0.5) / n as f32;
                let y = rect.min.y + t * rect.height();
                let x = rect.min.x + (marginal[i2] / max_val) as f32 * rect.width() * 0.9;
                pos2(x, y)
            })
            .collect();

        // Fill between curve and left edge
        fill_right_curve(painter, &points, rect.min.x, Color32::from_rgba_premultiplied(0, 120, 255, 40));
        painter.add(Shape::line(
            points,
            Stroke::new(1.5, Color32::from_rgb(40, 160, 255)),
        ));
    }
}

/// Fill the area between a curve and a horizontal baseline using a triangle mesh.
fn fill_under_curve(painter: &Painter, points: &[Pos2], baseline_y: f32, color: Color32) {
    if points.len() < 2 {
        return;
    }

    let mut mesh = epaint::Mesh::default();

    for i in 0..points.len() {
        mesh.vertices.push(epaint::Vertex {
            pos: points[i],
            uv: epaint::WHITE_UV,
            color,
        });
        mesh.vertices.push(epaint::Vertex {
            pos: pos2(points[i].x, baseline_y),
            uv: epaint::WHITE_UV,
            color,
        });
    }

    for i in 0..(points.len() - 1) {
        let top_left = (i * 2) as u32;
        let bot_left = (i * 2 + 1) as u32;
        let top_right = (i * 2 + 2) as u32;
        let bot_right = (i * 2 + 3) as u32;

        mesh.indices
            .extend_from_slice(&[top_left, bot_left, top_right]);
        mesh.indices
            .extend_from_slice(&[top_right, bot_left, bot_right]);
    }

    painter.add(Shape::mesh(mesh));
}

/// Fill the area between a vertical curve and a vertical baseline (for right marginal).
fn fill_right_curve(painter: &Painter, points: &[Pos2], baseline_x: f32, color: Color32) {
    if points.len() < 2 {
        return;
    }

    let mut mesh = epaint::Mesh::default();

    for i in 0..points.len() {
        // Curve point
        mesh.vertices.push(epaint::Vertex {
            pos: points[i],
            uv: epaint::WHITE_UV,
            color,
        });
        // Baseline point (same y, x = baseline)
        mesh.vertices.push(epaint::Vertex {
            pos: pos2(baseline_x, points[i].y),
            uv: epaint::WHITE_UV,
            color,
        });
    }

    for i in 0..(points.len() - 1) {
        let curve_top = (i * 2) as u32;
        let base_top = (i * 2 + 1) as u32;
        let curve_bot = (i * 2 + 2) as u32;
        let base_bot = (i * 2 + 3) as u32;

        mesh.indices
            .extend_from_slice(&[curve_top, base_top, curve_bot]);
        mesh.indices
            .extend_from_slice(&[curve_bot, base_top, base_bot]);
    }

    painter.add(Shape::mesh(mesh));
}
