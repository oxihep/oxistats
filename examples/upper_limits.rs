//! Upper limits on a signal in a counting experiment, three ways.
//!
//! An analysis sees `n` events where `b` background events were expected.
//! How large a signal `s` does that exclude at 95 % confidence?
//!
//! - **Classical (Neyman)**: the largest `s` with `P(N ≤ n | s + b) ≥ 5 %`.
//!   With a downward fluctuation it can exclude even `s = 0`.
//! - **CLs** (the LHC convention): divide by `P(N ≤ n | b)`, so a deficit
//!   does not exclude signals the experiment has no sensitivity to.
//! - **Feldman–Cousins**: the unified interval, which turns from an upper
//!   limit into a two-sided interval as the excess grows.
//!
//! ```sh
//! cargo run --example upper_limits
//! ```

use oxistats::{feldman_cousins, poisson_cdf};

/// The `s` where `f(s)` crosses `target`, for `f` decreasing in `s`.
fn solve(f: impl Fn(f64) -> f64, target: f64) -> f64 {
    let (mut lo, mut hi) = (0.0_f64, 100.0_f64);
    if f(lo) < target {
        return 0.0;
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > target { lo = mid } else { hi = mid }
    }
    0.5 * (lo + hi)
}

fn main() {
    let cl = 0.95;
    let alpha = 1.0 - cl;
    let b = 3.0;
    println!("95 % CL upper limits on the signal s, with b = {b} expected background");
    println!(
        "{:>4} {:>12} {:>10} {:>18}",
        "n", "classical", "CLs", "Feldman–Cousins"
    );
    for n in [0u32, 1, 2, 3, 5, 8, 12] {
        let nf = f64::from(n);
        let classical = solve(|s| poisson_cdf(nf, s + b), alpha);
        let cls = solve(|s| poisson_cdf(nf, s + b) / poisson_cdf(nf, b), alpha);
        let (fc_lo, fc_hi) = feldman_cousins(n as usize, b, cl);
        let fc = if fc_lo > 0.0 {
            format!("[{fc_lo:.2}, {fc_hi:.2}]")
        } else {
            format!("< {fc_hi:.2}")
        };
        println!("{n:>4} {classical:>12.2} {cls:>10.2} {fc:>18}");
    }
    println!();
    println!("With n = 0 on b = 3, the classical limit is about 0: it \"excludes\"");
    println!("almost any signal because the background itself fluctuated low. CLs");
    println!("stays at the familiar 3.0, the limit with no background at all.");
    println!("Feldman–Cousins quotes a two-sided interval once the excess is large.");
}
