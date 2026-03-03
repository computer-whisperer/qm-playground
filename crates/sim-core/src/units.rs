/// Physical constants in SI units, and conversion factors to/from atomic units.
///
/// In atomic units: ℏ = mₑ = e = 4πε₀ = 1
///
/// We work exclusively in atomic units inside the simulation.
/// These constants are provided for reference and for UI display conversion.

/// Reduced Planck constant (J·s)
pub const HBAR_SI: f64 = 1.054_571_817e-34;

/// Electron mass (kg)
pub const M_ELECTRON_SI: f64 = 9.109_383_702e-31;

/// Elementary charge (C)
pub const E_CHARGE_SI: f64 = 1.602_176_634e-19;

/// Bohr radius (m) — atomic unit of length
pub const BOHR_RADIUS_SI: f64 = 5.291_772_109e-11;

/// Hartree energy (J) — atomic unit of energy
pub const HARTREE_SI: f64 = 4.359_744_722e-18;

/// Hartree energy (eV)
pub const HARTREE_EV: f64 = 27.211_386_246;

/// Atomic unit of time (s) = ℏ / Eₕ
pub const TIME_AU_SI: f64 = 2.418_884_326e-17;

/// Convert atomic units of length to Ångströms.
pub fn bohr_to_angstrom(x: f64) -> f64 {
    x * 0.529_177_211
}

/// Convert atomic units of energy to eV.
pub fn hartree_to_ev(e: f64) -> f64 {
    e * HARTREE_EV
}

/// Convert atomic units of time to femtoseconds.
pub fn au_time_to_fs(t: f64) -> f64 {
    t * TIME_AU_SI * 1e15
}
