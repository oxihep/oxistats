//! Distributions and hypothesis tests, the `scipy.stats` way.
//!
//! ```sh
//! cargo run --example tests_and_distributions
//! ```

use oxistats::{ChiSquared, Normal, StudentT, chisquare, ks_1samp, ks_2samp, ttest_ind};

/// A reproducible stream of standard-normal draws (Box–Muller over a small
/// LCG), so the example needs no RNG crate.
fn normal_draws(n: usize, seed: u64, mean: f64, sigma: f64) -> Vec<f64> {
    let mut s = seed;
    let mut uniform = move || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            let (u1, u2) = (uniform(), uniform());
            mean + sigma * (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
        })
        .collect()
}

fn main() {
    // Distributions: pdf, cdf, sf (1 − cdf, accurate in the tail) and ppf.
    let n = Normal::new(0.0, 1.0);
    println!(
        "Normal: Φ(1.96) = {:.4}, 97.5 % quantile = {:.4}",
        n.cdf(1.96),
        n.ppf(0.975)
    );
    let t = StudentT::new(10.0);
    println!(
        "Student t(10): two-sided p for t = 2.5 is {:.4}",
        2.0 * t.sf(2.5)
    );

    // A χ² goodness of fit: χ² = 35.2 on 25 degrees of freedom.
    let chi2 = ChiSquared::new(25.0);
    println!("χ²(25): p-value of 35.2 is {:.4}", chi2.sf(35.2));

    // Observed counts against expected ones (scipy.stats.chisquare).
    let observed = [18.0, 22.0, 31.0, 29.0];
    let expected = [25.0, 25.0, 25.0, 25.0];
    let (stat, p) = chisquare(&observed, &expected).expect("same length");
    println!("chisquare: χ² = {stat:.3}, p = {p:.4}");

    // Two samples: do they share a mean? a distribution?
    let a = normal_draws(400, 1, 0.0, 1.0);
    let b = normal_draws(400, 2, 0.15, 1.0);
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    println!("sample means: {:.3} and {:.3}", mean(&a), mean(&b));
    let (t_stat, p_t) = ttest_ind(&a, &b);
    println!("ttest_ind: t = {t_stat:.3}, p = {p_t:.4}");
    let (d, p_ks) = ks_2samp(&a, &b);
    println!("ks_2samp: D = {d:.3}, p = {p_ks:.4}");

    // One sample against a model's CDF.
    let (d1, p1) = ks_1samp(&a, |x| n.cdf(x));
    println!("ks_1samp against N(0, 1): D = {d1:.3}, p = {p1:.4}");
}
