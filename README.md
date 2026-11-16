# oxistats

[![CI](https://github.com/oxihep/oxistats/actions/workflows/ci.yml/badge.svg)](https://github.com/oxihep/oxistats/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Rust 1.100+](https://img.shields.io/badge/rust-1.100%2B-orange.svg)](https://www.rust-lang.org)

Pure-Rust statistics for physics: a dependency-free `f64` core whose numerics
track [`scipy.stats`](https://docs.scipy.org/doc/scipy/reference/stats.html)
(and `scipy.special`) **function-for-function**, verified against it, plus the
HEP pieces scipy does not have: fit lineshapes, significance ↔ p-value,
efficiency intervals and the Feldman–Cousins unified interval.

```toml
[dependencies]
oxistats = { git = "https://github.com/oxihep/oxistats", tag = "v0.1.0" }
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

## The oxihep family

Rust crates for high-energy physics, split the way Python's
[Scikit-HEP](https://scikit-hep.org) is. Each works on its own; together they
cover reading data, histogramming, fitting and plotting.

| Crate | What it does | Python counterpart |
| --- | --- | --- |
| [oxiroot](https://github.com/oxihep/oxiroot) | Read and write ROOT files: `TTree`, RNTuple, ROOT histograms and graphs | uproot |
| [oxihist](https://github.com/oxihep/oxihist) | Histograms: typed axes, weighted and profile storages | boost-histogram, hist |
| [oxiplot](https://github.com/oxihep/oxiplot) | Plotting to SVG/PNG/PDF with an mplhep look | matplotlib, mplhep |
| [oxifit](https://github.com/oxihep/oxifit) | Curve fitting with Minuit2, and a `TFormula` engine | iminuit |
| **oxistats** (this crate) | Statistics, HEP lineshapes and intervals | scipy.stats, hepstats |
| [oxiparticle](https://github.com/oxihep/oxiparticle) | PDG particle data and ID decoding | particle |

The crates are not on crates.io yet: depend on them from GitHub, pinned to a
release tag, as shown above.

## History

This crate was `oxiroot-stat` inside the
[oxiroot](https://github.com/oxihep/oxiroot) workspace; its history was carried
over when it moved here.

## License

MIT — see [LICENSE](LICENSE).
