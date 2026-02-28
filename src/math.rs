/// Computes `atan2(y, x)` in radians via a two-stage argument reduction and
/// a degree-7 Taylor polynomial, with no calls to libm or any platform math
/// library. The implementation is pure `core` with zero external dependencies.
///
/// # Algorithm
///
/// 1. Reduce to `t ∈ [0, 1]` via `min(|x|,|y|)/max(|x|,|y|)` — single division.
/// 2. Further reduce to `s ∈ [0, tan(π/8)]` via:
///    `atan(t) = π/4 + atan((t − 1) / (t + 1))`
/// 3. Evaluate the degree-7 Taylor polynomial (exact rational coefficients):
///    `atan(s) = s·(1 − s²/3 + s⁴/5 − s⁶/7)`
/// 4. Undo both reductions; restore sign and quadrant via `copysign`.
///
/// # Precondition
///
/// `(x, y) ≠ (0, 0)`.
#[inline(always)]
pub fn atan2f(y: f32, x: f32) -> f32 {
    use core::f32::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    let ax = x.abs();
    let ay = y.abs();

    let swapped = ay > ax;
    let t = ax.min(ay) / ax.max(ay);

    const TAN_PI_8: f32 = core::f32::consts::SQRT_2 - 1.0;
    let (s, offset) = if t > TAN_PI_8 {
        ((t - 1.0) / (t + 1.0), FRAC_PI_4)
    } else {
        (t, 0.0_f32)
    };

    let s2 = s * s;
    let p = s * (1.0 + s2 * (-1.0 / 3.0 + s2 * (1.0 / 5.0 + s2 * (-1.0 / 7.0))));

    let r = p + offset;
    let r = if swapped { FRAC_PI_2 - r } else { r };

    f32::copysign(if x < 0.0 { PI - r } else { r }, y)
}
