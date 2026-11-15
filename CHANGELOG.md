# Changelog

Notable changes to oxistats, by release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0]

First release under the `oxistats` name. The code is `oxiroot-stat` 0.3.0 from
the [oxiroot](https://github.com/oxihep/oxiroot) workspace, moved to its own
repository with its history.

### Changed

- Renamed from `oxiroot-stat` to `oxistats`: replace `oxiroot_stat::` with
  `oxistats::` (oxiroot users keep `oxiroot::stat`, which re-exports this crate).
- The Euler–Mascheroni constant now comes from `std::f64::consts::EULER_GAMMA`.
