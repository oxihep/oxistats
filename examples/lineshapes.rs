//! The HEP lineshapes, tabulated: peaked shapes normalised to height 1 at the
//! peak, densities normalised to unit area.
//!
//! ```sh
//! cargo run --example lineshapes
//! ```

use oxistats::{
    argus, breit_wigner, crystal_ball, double_crystal_ball, gaussian, landau,
    relativistic_breit_wigner, voigtian,
};

fn main() {
    // A Z peak: mass 91.19 GeV, width 2.50 GeV, detector resolution 1.5 GeV.
    let (m, gamma, res) = (91.19, 2.50, 1.5);
    println!("    m |  gauss | CB(α=1.5,n=3) | double CB |    BW (pdf) | rel. BW | Voigt (pdf)");
    for i in 0..=12 {
        let x = 80.0 + 2.0 * f64::from(i);
        println!(
            "{x:5.1} | {:6.4} | {:13.4} | {:9.4} | {:11.5} | {:7.5} | {:11.5}",
            gaussian(x, m, res),
            crystal_ball(x, m, res, 1.5, 3.0),
            double_crystal_ball(x, m, res, 1.5, 3.0, 2.0, 5.0),
            breit_wigner(x, m, gamma),
            relativistic_breit_wigner(x, m, gamma),
            voigtian(x, m, res, gamma),
        );
    }

    // The densities integrate to one: check the Voigt profile numerically.
    let dx = 0.01;
    let area: f64 = (0..10_000)
        .map(|i| voigtian(41.19 + dx * (f64::from(i) + 0.5), m, res, gamma) * dx)
        .sum();
    println!("∫ Voigt over m ± 50 GeV = {area:.4} (the rest is in the Lorentzian tails)");

    // An endpoint background (B-meson beam-constrained mass) and an energy-loss
    // distribution.
    println!(
        "ARGUS at 5.20, 5.27, 5.289 GeV: {:.4}, {:.4}, {:.4}",
        argus(5.20, 5.29, -20.0, 0.5),
        argus(5.27, 5.29, -20.0, 0.5),
        argus(5.289, 5.29, -20.0, 0.5)
    );
    let peak = (0..2000)
        .map(|i| -5.0 + 0.01 * f64::from(i))
        .max_by(|a, b| landau(*a, 0.0, 1.0).total_cmp(&landau(*b, 0.0, 1.0)))
        .unwrap();
    println!("Landau(0, 1) peaks near x = {peak:.2}");
}
