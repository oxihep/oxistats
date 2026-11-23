# Changelog

Notable changes to oxistats, by release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.2.1]

### Changed

- **`erf` and `erfc` are about 16× faster** (7 ns, from 111) and more
  accurate: fdlibm's rational approximations, as the C library and ROOT's
  `TMath::Erf` use, in place of the incomplete gamma. Against the C library
  they are now within 2 ulps (`erf`) and 6 ulps (`erfc`), where they were
  within 39 ulps, and `erfc` follows its tail to x = 28 where it flushed to
  0 from x ≈ 26.7.

## [0.2.0]

API review toward 1.0 (breaking).

### Added

- Runnable examples: `tests_and_distributions`, `lineshapes` and
  `intervals` (`cargo run --example …`), listed in the README.

### Changed

- **`#[non_exhaustive]`** on `Describe`, `Normal`, `StudentT`, `ChiSquared`
  and `FisherF`: build them with their constructors (`Normal::new`, …);
  their fields stay readable.
- **`Normal::std` is now `Normal::std_dev`**, the name the oxihep crates use
  for a standard deviation everywhere.
- **Rust 1.100 and edition 2024.** The minimum supported Rust version
  (`rust-version`) moves from 1.95 to 1.100, the current stable release that
  CI already builds with, and the crates move to the 2024 edition (migrated
  with `cargo fix --edition`; formatting follows the 2024 style).

## [0.1.0]

First release under the `oxistats` name. The code is `oxiroot-stat` 0.3.0 from
the [oxiroot](https://github.com/oxihep/oxiroot) workspace, moved to its own
repository with its history.

### Changed

- Renamed from `oxiroot-stat` to `oxistats`: replace `oxiroot_stat::` with
  `oxistats::` (oxiroot users keep `oxiroot::stat`, which re-exports this crate).
- The Euler–Mascheroni constant now comes from `std::f64::consts::EULER_GAMMA`.
