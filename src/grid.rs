/// Minimum angle on the φ/ψ grid, in degrees (−180°).
///
/// The grid spans from `GRID_MIN` (−180°) to `GRID_MIN + (GRID_COUNT - 1) × GRID_STEP`
/// (180°), inclusive on both endpoints.
pub const GRID_MIN: f32 = -180.0;

/// Step size between consecutive grid points, in degrees (10°).
pub const GRID_STEP: f32 = 10.0;

/// Number of grid points along each axis, φ and ψ (37).
///
/// Indices 0 through 36 cover −180° to +180°; both endpoints carry identical
/// data due to angular periodicity, eliminating the need for modular
/// arithmetic during interpolation.
pub const GRID_COUNT: usize = 37;

/// Map a backbone dihedral angle (in degrees) to a grid index and fractional offset.
///
/// The input is clamped to \[−180.0, 180.0\] before conversion. Returns
/// `(lo, frac)` where `lo ∈ [0, 35]` and `frac ∈ [0.0, 1.0]`, suitable
/// for bilinear interpolation between `table[lo]` and `table[lo + 1]`.
///
/// Both `clamp` and `.min(35)` compile to branchless instructions
/// (`minss` / `CMOV`).
#[inline(always)]
pub fn angle_to_grid(deg: f32) -> (usize, f32) {
    // Map [-180, 180] → [0.0, 36.0].
    let shifted = (deg.clamp(GRID_MIN, -GRID_MIN) + (-GRID_MIN)) * (1.0 / GRID_STEP);
    // Clamp low index so that lo+1 ≤ 36 (last valid table index).
    let lo = (shifted as usize).min(GRID_COUNT - 2);
    let frac = shifted - lo as f32;
    (lo, frac)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_angle_to_grid_at_minus_180() {
        let (lo, frac) = angle_to_grid(-180.0);
        assert_eq!(lo, 0);
        assert_relative_eq!(frac, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_angle_to_grid_at_plus_180() {
        let (lo, frac) = angle_to_grid(180.0);
        assert_eq!(lo, 35);
        assert_relative_eq!(frac, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_angle_to_grid_at_zero() {
        let (lo, frac) = angle_to_grid(0.0);
        assert_eq!(lo, 18);
        assert_relative_eq!(frac, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_angle_to_grid_midpoint() {
        let (lo, frac) = angle_to_grid(-175.0);
        assert_eq!(lo, 0);
        assert_relative_eq!(frac, 0.5, epsilon = 1e-6);
    }

    #[test]
    fn test_angle_to_grid_clamp_below() {
        let (lo, frac) = angle_to_grid(-200.0);
        assert_eq!(lo, 0);
        assert_relative_eq!(frac, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_angle_to_grid_clamp_above() {
        let (lo, frac) = angle_to_grid(200.0);
        assert_eq!(lo, 35);
        assert_relative_eq!(frac, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_angle_to_grid_exact_grid_point() {
        let (lo, frac) = angle_to_grid(-60.0);
        assert_eq!(lo, 12);
        assert_relative_eq!(frac, 0.0, epsilon = 1e-6);
    }
}
