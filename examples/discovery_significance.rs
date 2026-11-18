//! How significant is an excess? The counting-experiment formulas behind a
//! "3σ evidence" or "5σ discovery".
//!
//! For `s` expected signal events over `b` background:
//!
//! - the naive `s/√b` holds only for `s ≪ b`;
//! - the Asimov significance `Z_A = √(2[(s+b)ln(1+s/b) − s])` (Cowan et al.,
//!   Eur. Phys. J. C 71 (2011) 1554) is the median significance of the
//!   profile-likelihood test, valid at any `s/b`;
//! - an uncertainty `σ_b` on the background lowers it (the same paper's
//!   formula with a nuisance parameter).
//!
//! Then the look-elsewhere effect: a local excess found by scanning many mass
//! points is less significant globally.
//!
//! ```sh
//! cargo run --example discovery_significance
//! ```

use oxistats::{poisson_sf, pvalue_from_significance, significance_from_pvalue};

/// The median discovery significance for `s` over `b` (no background
/// uncertainty).
fn asimov(s: f64, b: f64) -> f64 {
    (2.0 * ((s + b) * (1.0 + s / b).ln() - s)).sqrt()
}

/// The same with an absolute background uncertainty `sigma_b`.
fn asimov_with_uncertainty(s: f64, b: f64, sigma_b: f64) -> f64 {
    let v = sigma_b * sigma_b;
    let t1 = (s + b) * (((s + b) * (b + v)) / (b * b + (s + b) * v)).ln();
    let t2 = b * b / v * (1.0 + v * s / (b * (b + v))).ln();
    (2.0 * (t1 - t2)).sqrt()
}

fn main() {
    println!("Expected significance of a counting experiment");
    println!(
        "{:>8} {:>8} {:>8} {:>9} {:>14} {:>14}",
        "s", "b", "s/√b", "Z_Asimov", "Z (σ_b = 10%)", "Z (σ_b = 30%)"
    );
    for (s, b) in [
        (10.0_f64, 1000.0_f64),
        (50.0, 1000.0),
        (10.0, 10.0),
        (20.0, 5.0),
        (5.0, 0.5),
    ] {
        println!(
            "{s:>8} {b:>8} {:>8.2} {:>9.2} {:>14.2} {:>14.2}",
            s / b.sqrt(),
            asimov(s, b),
            asimov_with_uncertainty(s, b, 0.1 * b),
            asimov_with_uncertainty(s, b, 0.3 * b),
        );
    }
    println!("At large s/b, s/√b overstates the significance; a background uncertainty");
    println!("can matter more than the statistics.\n");

    // An observed count: n events where b were expected. The p-value is the
    // chance of b alone fluctuating to n or more.
    let (b, n) = (12.5, 27.0);
    let p = poisson_sf(n - 1.0, b); // P(N ≥ n | b)
    let z = significance_from_pvalue(p);
    println!("Observed {n} events on an expected background of {b}:");
    println!("  local p-value {p:.3e}, i.e. {z:.2}σ");

    // Look-elsewhere effect: searching N independent mass windows, the chance
    // that *some* window fluctuates this far is larger.
    for trials in [10.0_f64, 100.0, 1000.0] {
        let global = 1.0 - (1.0 - p).powf(trials);
        println!(
            "  searched {trials:>5} windows: global p-value {global:.3e}, i.e. {:.2}σ",
            significance_from_pvalue(global)
        );
    }

    // How much more data reaches 5σ? With s and b both scaling with the
    // integrated luminosity, Z grows roughly as √L.
    let (s, b) = (8.0, 20.0);
    let z_now = asimov(s, b);
    let mut scale = 1.0_f64;
    while asimov(s * scale, b * scale) < 5.0 {
        scale *= 1.01;
    }
    println!(
        "\nWith s = {s}, b = {b} today (Z = {z_now:.2}), 5σ needs {scale:.1}× the luminosity."
    );
    println!(
        "For reference: 3σ is p = {:.2e}, 5σ is p = {:.2e}.",
        pvalue_from_significance(3.0),
        pvalue_from_significance(5.0)
    );
}
