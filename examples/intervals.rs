//! Significances, confidence intervals and combinations: the everyday numbers
//! of a counting analysis.
//!
//! ```sh
//! cargo run --example intervals
//! ```

use oxistats::{
    agresti_coull_interval, bootstrap_ci, clopper_pearson, combine_measurements, feldman_cousins,
    median, poisson_conf_interval, pvalue_from_significance, significance_from_pvalue,
    wilson_interval,
};

fn main() {
    // Significance ↔ p-value (one-sided): the 3σ and 5σ conventions.
    for z in [3.0, 5.0] {
        println!("{z}σ ↔ p = {:.3e}", pvalue_from_significance(z));
    }
    println!("p = 1e-4 ↔ {:.2}σ", significance_from_pvalue(1e-4));

    // An efficiency: 45 of 50 events pass, 68.27 % (1σ) intervals.
    let (k, n, cl) = (45.0, 50.0, 0.6827);
    let show = |name: &str, (lo, hi): (f64, f64)| println!("{name:>15}: [{lo:.4}, {hi:.4}]");
    println!("efficiency {k}/{n} = {:.3}", k / n);
    show("Clopper–Pearson", clopper_pearson(k, n, cl));
    show("Wilson", wilson_interval(k, n, cl));
    show("Agresti–Coull", agresti_coull_interval(k, n, cl));

    // A count: 4 events observed.
    show("Garwood (k=4)", poisson_conf_interval(4.0, cl));

    // A signal rate over an expected background of 1.2, at 90 % CL: the
    // Feldman–Cousins unified interval goes from an upper limit (μ from 0) to
    // a two-sided interval as the observed count rises.
    for n_obs in [0, 3, 8] {
        let (lo, hi) = feldman_cousins(n_obs, 1.2, 0.90);
        println!("Feldman–Cousins, n = {n_obs}, b = 1.2: μ ∈ [{lo:.2}, {hi:.2}]");
    }

    // Combining three measurements of the same quantity (inverse-variance).
    let (value, error) = combine_measurements(&[91.187, 91.192, 91.180], &[0.007, 0.006, 0.010])
        .expect("same length");
    println!("combined: {value:.4} ± {error:.4}");

    // A bootstrap interval for the median of a small sample.
    let sample = [2.1, 2.4, 1.9, 3.3, 2.8, 2.2, 5.0, 2.6, 2.0, 2.9];
    let (lo, hi) = bootstrap_ci(&sample, median, 5000, 0.95, 7);
    println!(
        "median {:.2}, 95 % bootstrap interval [{lo:.2}, {hi:.2}]",
        median(&sample)
    );
}
