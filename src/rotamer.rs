/// A backbone-dependent rotamer entry.
///
/// Each rotamer describes a discrete side-chain conformation at a specific
/// backbone (φ, ψ) grid point, or the bilinearly interpolated result at an
/// arbitrary (φ, ψ) query. The const generic `N` encodes the number of χ
/// dihedral angles for the residue type.
///
/// # Layout
///
/// `#[repr(C)]` guarantees a stable, deterministic field layout across
/// compiler versions. This is beneficial for reproducible binary output
/// and cache-line alignment predictability.
///
/// # Examples
///
/// ```
/// use dunbrack::{Residue, Rotamer, Val};
///
/// let rots: Vec<Rotamer<1>> = Val::rotamers(-60.0, -40.0).collect();
/// assert_eq!(rots.len(), 3);
/// assert!(rots[0].prob > 0.0);
/// assert!(rots[0].chi_sigma[0] > 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(C)]
pub struct Rotamer<const N: usize> {
    /// Rotamer bin indices (1-based).
    ///
    /// For a residue with `N` χ angles, only `r[0..N]` carry meaningful
    /// values. Each element identifies, for its corresponding χ angle,
    /// which discrete bin (gauche+, trans, gauche−, etc.) the rotamer
    /// belongs to.
    pub r: [u8; N],

    /// Prior probability P(r | φ, ψ).
    ///
    /// At exact grid points this is the raw value from the Dunbrack 2010
    /// library. After bilinear interpolation the probabilities are
    /// re-normalized so that Σ prob = 1.0 across all rotamers.
    pub prob: f32,

    /// Mean χ dihedral angles in degrees.
    ///
    /// After interpolation, these are computed via circular weighted mean
    /// (using sin/cos decomposition) to correctly handle the ±180°
    /// wraparound.
    pub chi_mean: [f32; N],

    /// Standard deviations of the χ angles in degrees.
    ///
    /// Always positive. Interpolated via standard (linear) bilinear
    /// weighting — no circular treatment is needed for positive scalars.
    pub chi_sigma: [f32; N],
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_rotamer_copy_clone() {
        let rot = Rotamer {
            r: [1, 2],
            prob: 0.5,
            chi_mean: [60.0, -60.0],
            chi_sigma: [10.0, 12.0],
        };
        let copy = rot;
        let clone = rot.clone();

        assert_eq!(rot, copy);
        assert_eq!(rot, clone);
    }

    #[test]
    fn test_rotamer_fields() {
        let rot = Rotamer {
            r: [3],
            prob: 0.42,
            chi_mean: [175.0],
            chi_sigma: [8.5],
        };

        assert_eq!(rot.r, [3]);
        assert_relative_eq!(rot.prob, 0.42, epsilon = 1e-6);
        assert_relative_eq!(rot.chi_mean[0], 175.0, epsilon = 1e-6);
        assert_relative_eq!(rot.chi_sigma[0], 8.5, epsilon = 1e-6);
    }
}
