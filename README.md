# oxistats

[![CI](https://github.com/oxihep/oxistats/actions/workflows/ci.yml/badge.svg)](https://github.com/oxihep/oxistats/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Rust 1.95+](https://img.shields.io/badge/rust-1.95%2B-orange.svg)](https://www.rust-lang.org)

Pure-Rust statistics for physics: a dependency-free `f64` core whose numerics
track [`scipy.stats`](https://docs.scipy.org/doc/scipy/reference/stats.html)
(and `scipy.special`) **function-for-function**, verified against it, plus the
HEP pieces scipy does not have: fit lineshapes, significance ↔ p-value,
efficiency intervals and the Feldman–Cousins unified interval.

Part of the [oxihep](https://github.com/oxihep) family of Rust crates for
high-energy physics, alongside [oxiroot](https://github.com/oxihep/oxiroot)
(ROOT file IO) and [oxifit](https://github.com/oxihep/oxifit) (fitting).

```toml
[dependencies]
oxistats = { git = "https://github.com/oxihep/oxistats" }
```

```rust
use oxistats::*;

let z = Normal::standard().ppf(0.975); // 1.959963984540054  (scipy.stats.norm.ppf)
let (t, p) = ttest_1samp(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0], 4.0);
```

| Topic | What it holds |
| --- | --- |
| Special functions | `erf`/`erfc`, `gammaln`, incomplete gamma & beta, inverse normal CDF |
| Distributions | `Normal`, `StudentT`, `ChiSquared`, `FisherF`; Poisson/Binomial |
| Descriptive statistics | `describe`, moments, `skew`/`kurtosis`, quantiles, `entropy`, … |
| Correlation | `pearsonr`, `spearmanr` |
| Hypothesis tests | t-tests, `normaltest`, `chisquare`, KS, Mann–Whitney, Wilcoxon |
| HEP lineshapes | Crystal Ball, Breit–Wigner, Voigt, Novosibirsk, ARGUS, Landau, … |
| Physics helpers | significance ↔ p-value, weighted means, efficiency & Feldman–Cousins intervals |
| Resampling | seeded percentile `bootstrap_ci` |

The test suite pins every reference value against scipy. The
[statistics guide](https://oxihep.github.io/oxiroot/guide/statistics/) walks
through each area with examples (it uses the `oxiroot::stat` re-export, which is
this crate).

## History

This crate was `oxiroot-stat` inside the
[oxiroot](https://github.com/oxihep/oxiroot) workspace; its history was carried
over when it moved here.

## License

MIT — see [LICENSE](LICENSE).
