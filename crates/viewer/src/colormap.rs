/// Inferno colormap: black → purple → orange → yellow.
///
/// 8-point lookup table with linear interpolation.
const INFERNO: [[u8; 3]; 8] = [
    [0, 0, 4],        // 0.000
    [31, 12, 72],     // 0.143
    [85, 15, 109],    // 0.286
    [136, 34, 106],   // 0.429
    [186, 54, 85],    // 0.571
    [227, 89, 51],    // 0.714
    [249, 149, 21],   // 0.857
    [252, 255, 164],  // 1.000
];

/// Map a value in [0, 1] to an (R, G, B) color via the inferno colormap.
pub fn inferno(t: f64) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    let n = INFERNO.len() - 1;
    let scaled = t * n as f64;
    let i = (scaled as usize).min(n - 1);
    let frac = scaled - i as f64;

    let c0 = INFERNO[i];
    let c1 = INFERNO[i + 1];

    [
        (c0[0] as f64 + (c1[0] as f64 - c0[0] as f64) * frac) as u8,
        (c0[1] as f64 + (c1[1] as f64 - c0[1] as f64) * frac) as u8,
        (c0[2] as f64 + (c1[2] as f64 - c0[2] as f64) * frac) as u8,
    ]
}
