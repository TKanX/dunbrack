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
/// use dunbrack::{Rotamer, Val};
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
