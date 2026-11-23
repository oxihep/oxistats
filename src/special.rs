//! Special functions (pure `f64`, no dependencies): the log-gamma and
//! regularized incomplete gamma/beta functions, the error function, and the
//! inverse normal CDF. These underpin the distributions and tests in the rest of
//! the crate; names follow `scipy.special` where there is a direct counterpart.

mod temme;

const MACHEP: f64 = 1.1102230246251565e-16;
const BIG: f64 = 4.503599627370496e15;
const BIG_INV: f64 = 2.220446049250313e-16;
/// `ln(f64::MAX)`; an exponent below `-MAX_LOG` underflows to 0.
const MAX_LOG: f64 = 709.782712893384;
/// Above 2⁵³ a `f64` cannot represent `a + 1` exactly, so the incomplete-gamma
/// series and continued fraction, which step a counter from `a` by 1, stop
/// making progress. They are never reached there (see [`uses_temme`]); if they
/// were, they would return `NaN` rather than spin.
const MAX_EXACT_STEP: f64 = 9_007_199_254_740_992.0;
/// The Euler–Mascheroni constant γ.
const EULER: f64 = std::f64::consts::EULER_GAMMA;

/// Iteration budget for the incomplete-gamma series and continued fraction.
///
/// Near `x ≈ a` the series terms only fall off like `exp(−k²/2a)`, but that
/// region is left to Temme's expansion for `a > 20`, so the series and the
/// continued fraction converge within a few hundred steps wherever they are
/// used. The budget exists only to guarantee termination; running out yields
/// `NaN`, never a truncated sum.
fn iteration_budget(a: f64) -> usize {
    // Float-to-int casts saturate, so a huge `a` cannot overflow here.
    (2000.0 + 20.0 * a.sqrt()) as usize
}

/// Natural log of the absolute gamma function, `ln|Γ(x)|` — `scipy.special.gammaln`.
/// Lanczos approximation (g = 7), with the reflection formula for `x < 0.5`.
#[must_use]
pub fn gammaln(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.9999999999998099,
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        let pi = std::f64::consts::PI;
        pi.ln() - (pi * x).sin().abs().ln() - gammaln(1.0 - x)
    } else {
        let x = x - 1.0;
        let t = x + 7.5;
        let mut a = C[0];
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Regularized lower incomplete gamma `P(a, x)` — `scipy.special.gammainc`.
/// The CDF of a Gamma(`a`) at `x`.
///
/// As in scipy: a `NaN` argument, a negative one, and `(0, 0)` yield `NaN`;
/// `P(0, x) == 1` for `x > 0`, `P(a, 0) == 0`, `P(a, +∞) == 1`, `P(+∞, x) == 0`,
/// and `(+∞, +∞)` is `NaN`. Near `x ≈ a` with `a > 20` it uses Temme's uniform
/// asymptotic expansion, which stays accurate however large `a` is.
#[must_use]
pub fn gammainc(a: f64, x: f64) -> f64 {
    if let Some(edge) = gamma_edge(a, x, true) {
        return edge;
    }
    if uses_temme(a, x) {
        return temme(a, x, true);
    }
    if x > 1.0 && x > a {
        return 1.0 - gammaincc(a, x);
    }
    igam_series(a, x)
}

/// Regularized upper incomplete gamma `Q(a, x) = 1 - P(a, x)` —
/// `scipy.special.gammaincc`.
///
/// As in scipy: a `NaN` argument, a negative one, and `(0, 0)` yield `NaN`;
/// `Q(0, x) == 0` for `x > 0`, `Q(a, 0) == 1`, `Q(a, +∞) == 0`, `Q(+∞, x) == 1`,
/// and `(+∞, +∞)` is `NaN`. Near `x ≈ a` with `a > 20` it uses Temme's uniform
/// asymptotic expansion, which stays accurate however large `a` is.
#[must_use]
pub fn gammaincc(a: f64, x: f64) -> f64 {
    if let Some(edge) = gamma_edge(a, x, false) {
        return edge;
    }
    if uses_temme(a, x) {
        return temme(a, x, false);
    }
    // scipy's choice among the series for P, its complement, and the
    // continued fraction, by where each converges fast without cancelling.
    if x > 1.1 {
        if x < a {
            1.0 - igam_series(a, x)
        } else {
            igamc_continued_fraction(a, x)
        }
    } else if x <= 0.5 {
        if -0.4 / x.ln() < a {
            1.0 - igam_series(a, x)
        } else {
            igamc_series(a, x)
        }
    } else if x * 1.1 < a {
        1.0 - igam_series(a, x)
    } else {
        igamc_series(a, x)
    }
}

/// `P(a, x)` (`lower`) or `Q(a, x)` where one of the arguments is `NaN`,
/// negative, zero or infinite; `None` for finite, positive arguments.
fn gamma_edge(a: f64, x: f64, lower: bool) -> Option<f64> {
    // `P` for the lower function, `1 - P` for the upper one.
    let p = |value: f64| Some(if lower { value } else { 1.0 - value });
    if a.is_nan() || x.is_nan() || a < 0.0 || x < 0.0 {
        return Some(f64::NAN);
    }
    if a == 0.0 {
        return if x > 0.0 { p(1.0) } else { Some(f64::NAN) };
    }
    if x == 0.0 {
        return p(0.0);
    }
    match (a.is_infinite(), x.is_infinite()) {
        (true, true) => Some(f64::NAN),
        (true, false) => p(0.0),
        (false, true) => p(1.0),
        (false, false) => None,
    }
}

/// Whether `(a, x)` is in the region left to Temme's expansion: `a > 20` and `x`
/// within 30% of `a`. scipy uses it there for `a < 200`, and only within
/// `4.5/√a` of the peak above that; using it throughout keeps the series (which
/// cannot step past `a ≥ 2⁵³`) away from huge `a`, and the expansion is at least
/// as accurate there, since its terms fall off like `a⁻ᵏ`.
fn uses_temme(a: f64, x: f64) -> bool {
    a > 20.0 && (x - a).abs() / a < 0.3
}

/// `ln(1 + x) - x`, accurate for small `x` (scipy's `log1pmx`).
fn log1pmx(x: f64) -> f64 {
    if x.abs() >= 0.5 {
        return x.ln_1p() - x;
    }
    // -x²/2 + x³/3 - …: at |x| < 1/2 each term is under half the previous one.
    let mut power = x;
    let mut sum = 0.0;
    for n in 2..100 {
        power *= -x;
        let term = power / f64::from(n);
        sum += term;
        if term.abs() <= MACHEP * sum.abs() {
            break;
        }
    }
    sum
}

/// `ln Γ*(a)`, where `Γ*(a) = Γ(a) / (√(2π/a) (a/e)^a)` tends to 1 as `a` grows.
/// For `a ≥ 10` it is the sum of Stirling's series, `Σ B₂ⱼ / (2j(2j−1) a^(2j−1))`,
/// whose eighth term is below 3·10⁻¹⁷ there.
fn ln_gammastar(a: f64) -> f64 {
    if a < 10.0 {
        return gammaln(a) + a - (a - 0.5) * a.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln();
    }
    const STIRLING: [f64; 8] = [
        1.0 / 12.0,
        -1.0 / 360.0,
        1.0 / 1260.0,
        -1.0 / 1680.0,
        1.0 / 1188.0,
        -691.0 / 360_360.0,
        1.0 / 156.0,
        -3617.0 / 122_400.0,
    ];
    let inv = 1.0 / a;
    let inv2 = inv * inv;
    let mut power = inv;
    let mut sum = 0.0;
    for c in STIRLING {
        sum += c * power;
        power *= inv2;
    }
    sum
}

/// `x^a e^(-x) / Γ(a)`, the factor in front of both series. For `a ≥ 10` it is
/// `exp(a·log1pmx((x−a)/a)) · √(a/2π) / Γ*(a)`: the direct form subtracts terms
/// of size `a ln a` and loses about `a·10⁻¹⁶` in the exponent, which spoils the
/// result for large `a`.
fn igam_fac(a: f64, x: f64) -> f64 {
    let log = if a < 10.0 {
        a * x.ln() - x - gammaln(a)
    } else {
        let sigma = (x - a) / a;
        // ln((x/a)^a e^(a−x)); below x = a/2, 1 + sigma loses the digits of x/a.
        let main = if sigma > -0.5 {
            a * log1pmx(sigma)
        } else {
            a * (x / a).ln() - (x - a)
        };
        main + 0.5 * (a / (2.0 * std::f64::consts::PI)).ln() - ln_gammastar(a)
    };
    if log < -MAX_LOG {
        return 0.0;
    }
    log.exp()
}

/// `P(a, x)` by its power series (DLMF 8.11.4).
fn igam_series(a: f64, x: f64) -> f64 {
    let fac = igam_fac(a, x);
    if fac == 0.0 {
        return 0.0;
    }
    if !fac.is_finite() || a >= MAX_EXACT_STEP {
        return f64::NAN;
    }
    let mut r = a;
    let mut c = 1.0;
    let mut sum = 1.0;
    for _ in 0..iteration_budget(a) {
        r += 1.0;
        c *= x / r;
        sum += c;
        if c <= MACHEP * sum {
            return sum * fac / a;
        }
    }
    f64::NAN
}

/// `Q(a, x)` by its continued fraction (DLMF 8.9.2), for `x > a`.
fn igamc_continued_fraction(a: f64, x: f64) -> f64 {
    let fac = igam_fac(a, x);
    if fac == 0.0 {
        return 0.0;
    }
    if !fac.is_finite() || a >= MAX_EXACT_STEP {
        return f64::NAN;
    }
    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut c = 0.0;
    let mut pkm2 = 1.0;
    let mut qkm2 = x;
    let mut pkm1 = x + 1.0;
    let mut qkm1 = z * x;
    let mut ans = pkm1 / qkm1;
    for _ in 0..iteration_budget(a) {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        let yc = y * c;
        let pk = pkm1 * z - pkm2 * yc;
        let qk = qkm1 * z - qkm2 * yc;
        if qk != 0.0 {
            let r = pk / qk;
            let t = ((ans - r) / r).abs();
            ans = r;
            if t <= MACHEP {
                return ans * fac;
            }
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if pk.abs() > BIG {
            pkm2 *= BIG_INV;
            pkm1 *= BIG_INV;
            qkm2 *= BIG_INV;
            qkm1 *= BIG_INV;
        }
    }
    f64::NAN
}

/// `Q(a, x)` for small `x` by DLMF 8.7.3, which avoids the cancellation in
/// `1 − P(a, x)` when `P` is close to 1.
fn igamc_series(a: f64, x: f64) -> f64 {
    let mut fac = 1.0;
    let mut sum = 0.0;
    for n in 1..2000 {
        let n = f64::from(n);
        fac *= -x / n;
        let term = fac / (a + n);
        sum += term;
        if term.abs() <= MACHEP * sum.abs() {
            break;
        }
    }
    let logx = x.ln();
    -(a * logx - lgam1p(a)).exp_m1() - (a * logx - gammaln(a)).exp() * sum
}

/// `ln Γ(1 + x)`, accurate near `x = 0` and `x = 1` (scipy's `lgam1p`).
fn lgam1p(x: f64) -> f64 {
    if x.abs() <= 0.5 {
        lgam1p_taylor(x)
    } else if (x - 1.0).abs() < 0.5 {
        x.ln() + lgam1p_taylor(x - 1.0)
    } else {
        gammaln(x + 1.0)
    }
}

/// The Taylor series of `ln Γ(1 + x)` about 0: `−γx + Σ_{n≥2} ζ(n) (−x)ⁿ / n`.
fn lgam1p_taylor(x: f64) -> f64 {
    /// ζ(2) … ζ(41).
    const ZETA: [f64; 40] = [
        1.644_934_066_848_226_4,
        1.202_056_903_159_594_2,
        1.082_323_233_711_138_1,
        1.036_927_755_143_37,
        1.017_343_061_984_449_2,
        1.008_349_277_381_923,
        1.004_077_356_197_944_4,
        1.002_008_392_826_082_1,
        1.000_994_575_127_818,
        1.000_494_188_604_119_4,
        1.000_246_086_553_308,
        1.000_122_713_347_578_5,
        1.000_061_248_135_058_8,
        1.000_030_588_236_307,
        1.000_015_282_259_408_6,
        1.000_007_637_197_637_9,
        1.000_003_817_293_265,
        1.000_001_908_212_716_5,
        1.000_000_953_962_033_8,
        1.000_000_476_932_986_9,
        1.000_000_238_450_502_7,
        1.000_000_119_219_926,
        1.000_000_059_608_189,
        1.000_000_029_803_503_4,
        1.000_000_014_901_554_9,
        1.000_000_007_450_711_8,
        1.000_000_003_725_334,
        1.000_000_001_862_659_8,
        1.000_000_000_931_327_5,
        1.000_000_000_465_662_8,
        1.000_000_000_232_831,
        1.000_000_000_116_415_5,
        1.000_000_000_058_207_7,
        1.000_000_000_029_103_8,
        1.000_000_000_014_552,
        1.000_000_000_007_276,
        1.000_000_000_003_638,
        1.000_000_000_001_819,
        1.000_000_000_000_909_5,
        1.000_000_000_000_454_7,
    ];
    if x == 0.0 {
        return 0.0;
    }
    let mut sum = -EULER * x;
    let mut power = -x;
    for (i, zeta) in ZETA.iter().enumerate() {
        power *= -x;
        let term = zeta * power / (i + 2) as f64;
        sum += term;
        if term.abs() < MACHEP * sum.abs() {
            break;
        }
    }
    sum
}

/// `P(a, x)` (`lower`) or `Q(a, x)` by Temme's uniform asymptotic expansion
/// (DLMF 8.12.3/8.12.4), as scipy evaluates it.
fn temme(a: f64, x: f64, lower: bool) -> f64 {
    let sigma = (x - a) / a;
    // η² / 2 = λ − 1 − ln λ with λ = x/a, and η has the sign of λ − 1.
    let eta = (-2.0 * log1pmx(sigma)).sqrt().copysign(sigma);
    let eta = if sigma == 0.0 { 0.0 } else { eta };
    let sign = if lower { -1.0 } else { 1.0 };
    let head = 0.5 * erfc(sign * eta * (a / 2.0).sqrt());

    let mut eta_powers = [0.0; 25];
    eta_powers[0] = 1.0;
    let mut known_powers = 0;
    let mut sum = 0.0;
    let mut a_power = 1.0;
    let mut previous = f64::INFINITY;
    for row in &temme::D {
        let mut ck = row[0];
        for (n, &d) in row.iter().enumerate().skip(1) {
            if n > known_powers {
                eta_powers[n] = eta * eta_powers[n - 1];
                known_powers = n;
            }
            let part = d * eta_powers[n];
            ck += part;
            if part.abs() < MACHEP * ck.abs() {
                break;
            }
        }
        let term = ck * a_power;
        if term.abs() > previous {
            break; // the asymptotic series has started to diverge
        }
        sum += term;
        if term.abs() < MACHEP * sum.abs() {
            break;
        }
        previous = term.abs();
        a_power /= a;
    }
    head + sign * (-0.5 * a * eta * eta).exp() * sum / (2.0 * std::f64::consts::PI * a).sqrt()
}

/// The error function `erf(x)` — `scipy.special.erf`, to about 1 ulp.
///
/// fdlibm's rational approximations (as in musl and most C libraries, and so
/// in ROOT's `TMath::Erf`): a few dozen floating-point operations, where the
/// incomplete-gamma route took several hundred.
#[must_use]
pub fn erf(x: f64) -> f64 {
    let bits = x.to_bits();
    let ix = high_word(bits) & 0x7fff_ffff;
    let negative = bits >> 63 != 0;
    if ix >= 0x7ff0_0000 {
        // erf(NaN) = NaN, erf(±inf) = ±1.
        return if x.is_nan() {
            x
        } else if negative {
            -1.0
        } else {
            1.0
        };
    }
    if ix < 0x3feb_0000 {
        // |x| < 0.84375
        if ix < 0x3e30_0000 {
            // |x| < 2⁻²⁸: avoid underflow
            return 0.125 * (8.0 * x + EFX8 * x);
        }
        let z = x * x;
        return x + x * (erf_small_p(z) / erf_small_q(z));
    }
    let y = if ix < 0x4018_0000 {
        // |x| < 6
        1.0 - erfc_tail(ix, x)
    } else {
        1.0 - f64::MIN_POSITIVE
    };
    if negative { -y } else { y }
}

/// The complementary error function `erfc(x) = 1 - erf(x)` —
/// `scipy.special.erfc`, to about 1 ulp, accurate in the tail (fdlibm's, as
/// [`erf`]).
#[must_use]
pub fn erfc(x: f64) -> f64 {
    let bits = x.to_bits();
    let ix = high_word(bits) & 0x7fff_ffff;
    let negative = bits >> 63 != 0;
    if ix >= 0x7ff0_0000 {
        // erfc(NaN) = NaN, erfc(+inf) = 0, erfc(-inf) = 2.
        return if x.is_nan() {
            x
        } else if negative {
            2.0
        } else {
            0.0
        };
    }
    if ix < 0x3feb_0000 {
        // |x| < 0.84375
        if ix < 0x3c70_0000 {
            // |x| < 2⁻⁵⁶
            return 1.0 - x;
        }
        let z = x * x;
        let y = erf_small_p(z) / erf_small_q(z);
        return if negative || ix < 0x3fd0_0000 {
            1.0 - (x + x * y)
        } else {
            0.5 - (x - 0.5 + x * y)
        };
    }
    if ix < 0x403c_0000 {
        // |x| < 28
        let t = erfc_tail(ix, x);
        return if negative { 2.0 - t } else { t };
    }
    if negative { 2.0 } else { 0.0 }
}

/// The high 32 bits of an `f64`'s representation (fdlibm's `GET_HIGH_WORD`).
fn high_word(bits: u64) -> u32 {
    u32::try_from(bits >> 32).unwrap_or(u32::MAX)
}

// fdlibm's coefficients (Sun Microsystems, 1993; freely redistributable).
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const ERX: f64 = 8.450_629_115_104_675_292_97e-1;
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const EFX8: f64 = 1.027_033_336_764_100_690_53;
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const PP: [f64; 5] = [
    1.283_791_670_955_125_585_61e-1,
    -3.250_421_072_470_014_993_70e-1,
    -2.848_174_957_559_851_047_66e-2,
    -5.770_270_296_489_441_591_57e-3,
    -2.376_301_665_665_016_260_84e-5,
];
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const QQ: [f64; 5] = [
    3.979_172_239_591_553_528_19e-1,
    6.502_224_998_876_729_444_85e-2,
    5.081_306_281_875_765_627_76e-3,
    1.324_947_380_043_216_445_26e-4,
    -3.960_228_278_775_368_123_20e-6,
];
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const PA: [f64; 7] = [
    -2.362_118_560_752_659_440_77e-3,
    4.148_561_186_837_483_316_66e-1,
    -3.722_078_760_357_013_238_47e-1,
    3.183_466_199_011_617_536_74e-1,
    -1.108_946_942_823_966_774_76e-1,
    3.547_830_432_561_823_593_71e-2,
    -2.166_375_594_868_790_843_00e-3,
];
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const QA: [f64; 6] = [
    1.064_208_804_008_442_282_86e-1,
    5.403_979_177_021_710_489_37e-1,
    7.182_865_441_419_626_628_68e-2,
    1.261_712_198_087_616_421_12e-1,
    1.363_708_391_202_905_073_62e-2,
    1.198_449_984_679_910_741_70e-2,
];
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const RA: [f64; 8] = [
    -9.864_944_034_847_148_227_05e-3,
    -6.938_585_727_071_817_643_72e-1,
    -1.055_862_622_532_329_098_14e1,
    -6.237_533_245_032_600_603_96e1,
    -1.623_966_694_625_734_703_55e2,
    -1.846_050_929_067_110_359_94e2,
    -8.128_743_550_630_659_342_46e1,
    -9.814_329_344_169_145_485_92,
];
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const SA: [f64; 8] = [
    1.965_127_166_743_925_712_92e1,
    1.376_577_541_435_190_426_00e2,
    4.345_658_774_752_292_288_21e2,
    6.453_872_717_332_678_803_36e2,
    4.290_081_400_275_678_333_86e2,
    1.086_350_055_417_794_351_34e2,
    6.570_249_770_319_281_701_35,
    -6.042_441_521_485_809_874_38e-2,
];
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const RB: [f64; 7] = [
    -9.864_942_924_700_099_285_97e-3,
    -7.992_832_376_805_230_065_74e-1,
    -1.775_795_491_775_475_198_89e1,
    -1.606_363_848_558_219_160_62e2,
    -6.375_664_433_683_896_277_22e2,
    -1.025_095_131_611_077_249_54e3,
    -4.835_191_916_086_513_970_19e2,
];
#[allow(
    clippy::excessive_precision,
    reason = "fdlibm's coefficients, as published"
)]
const SB: [f64; 7] = [
    3.033_806_074_348_245_829_24e1,
    3.257_925_129_965_739_188_26e2,
    1.536_729_586_084_436_959_94e3,
    3.199_858_219_508_595_539_08e3,
    2.553_050_406_433_164_425_83e3,
    4.745_285_412_069_553_672_15e2,
    -2.244_095_244_658_581_833_62e1,
];

/// `c[0] + z·(c[1] + z·(…))`.
fn horner(c: &[f64], z: f64) -> f64 {
    c.iter().rev().fold(0.0, |acc, &k| acc * z + k)
}

/// `1 + z·(c[0] + z·(c[1] + …))`.
fn horner1(c: &[f64], z: f64) -> f64 {
    1.0 + z * horner(c, z)
}

fn erf_small_p(z: f64) -> f64 {
    horner(&PP, z)
}

fn erf_small_q(z: f64) -> f64 {
    horner1(&QQ, z)
}

/// `erfc(|x|)` for `0.84375 ≤ |x| < 28` (fdlibm's `erfc1`/`erfc2`).
fn erfc_tail(ix: u32, x: f64) -> f64 {
    if ix < 0x3ff4_0000 {
        // |x| < 1.25
        let s = x.abs() - 1.0;
        return 1.0 - ERX - horner(&PA, s) / horner1(&QA, s);
    }
    let x = x.abs();
    let s = 1.0 / (x * x);
    let (r, big_s) = if ix < 0x4006_db6d {
        // |x| < 1/0.35
        (horner(&RA, s), horner1(&SA, s))
    } else {
        (horner(&RB, s), horner1(&SB, s))
    };
    // x with its low 32 bits cleared, so exp(-z²) is exact enough.
    let z = f64::from_bits(x.to_bits() & 0xffff_ffff_0000_0000);
    (-z * z - 0.5625).exp() * ((z - x) * (z + x) + r / big_s).exp() / x
}

/// Natural log of the beta function, `ln B(a, b)` — `scipy.special.betaln`.
#[must_use]
pub fn betaln(a: f64, b: f64) -> f64 {
    gammaln(a) + gammaln(b) - gammaln(a + b)
}

/// The beta function `B(a, b) = Γ(a)Γ(b)/Γ(a+b)` — `scipy.special.beta`
/// (for `a, b > 0`).
#[must_use]
pub fn beta(a: f64, b: f64) -> f64 {
    betaln(a, b).exp()
}

/// The regularized incomplete beta `I_x(a, b)` — `scipy.special.betainc`. The CDF
/// of a Beta(`a`, `b`) at `x ∈ [0, 1]`. Numerical-Recipes continued fraction with
/// the `I_x(a,b) = 1 - I_{1-x}(b,a)` symmetry for fast convergence.
#[must_use]
pub fn betainc(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front = (gammaln(a + b) - gammaln(a) - gammaln(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * betacf(a, b, x) / a
    } else {
        1.0 - front * betacf(b, a, 1.0 - x) / b
    }
}

/// Lentz's continued fraction for the incomplete beta (Numerical Recipes `betacf`).
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const FPMIN: f64 = 1e-300;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=200 {
        let m = f64::from(m);
        let m2 = 2.0 * m;
        // Even step.
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        // Odd step.
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 3.0e-15 {
            break;
        }
    }
    h
}

/// Inverse of the regularized incomplete beta in `x`: the `y`-quantile of a
/// Beta(`a`, `b`) — `scipy.special.betaincinv`. Bisection on the monotone
/// [`betainc`].
#[must_use]
pub fn betaincinv(a: f64, b: f64, y: f64) -> f64 {
    if y <= 0.0 {
        return 0.0;
    }
    if y >= 1.0 {
        return 1.0;
    }
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if betainc(a, b, mid) < y {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo <= 1e-15 {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// Inverse of the standard normal CDF (the quantile / probit function) —
/// `scipy.special.ndtri`. Acklam's rational approximation refined by one Halley
/// step, accurate to full `f64` precision on `(0, 1)`.
// Acklam's published coefficients are kept verbatim.
#[allow(clippy::excessive_precision)]
#[must_use]
pub fn ndtri(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    const P_LOW: f64 = 0.02425;
    let mut x = if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= 1.0 - P_LOW {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    };
    // One Halley refinement to reach machine precision.
    let e = 0.5 * erfc(-x / std::f64::consts::SQRT_2) - p;
    let u = e * (2.0 * std::f64::consts::PI).sqrt() * (x * x / 2.0).exp();
    x -= u / (1.0 + x * u / 2.0);
    x
}
