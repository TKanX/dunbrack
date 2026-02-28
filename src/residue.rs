use crate::sealed::Sealed;

/// Backbone-dependent rotamer library interface.
///
/// Each of the 22 amino acid residue types implements this trait, providing
/// type-safe, zero-allocation access to the Dunbrack 2010 rotamer data.
/// The trait is `sealed` — downstream crates cannot add new implementations.
///
/// # Examples
///
/// ```
/// use dunbrack::{Residue, Val};
///
/// // Query rotamers generically.
/// fn count<R: Residue>(phi: f32, psi: f32) -> usize {
///     R::rotamers(phi, psi).count()
/// }
///
/// assert_eq!(count::<Val>(-60.0, -40.0), 3);
/// ```
pub trait Residue: Sealed + Copy + 'static {
    /// Number of χ dihedral angles for this residue type.
    const N_CHI: usize;

    /// Number of distinct rotamers per (φ, ψ) grid cell.
    const N_ROTAMERS: usize;

    /// Concrete rotamer type, always `Rotamer<{N_CHI}>`.
    type Rot: Copy + 'static;

    /// Concrete iterator type, always `RotamerIter<{N_CHI}, {N_ROTAMERS}>`.
    type Iter: Iterator<Item = Self::Rot> + ExactSizeIterator;

    /// Return an iterator of bilinearly interpolated rotamers at the given
    /// backbone dihedral angles.
    ///
    /// Both `phi` and `psi` are clamped to \[−180.0, 180.0\] before use.
    /// The iterator yields exactly [`N_ROTAMERS`](Self::N_ROTAMERS) items
    /// whose probabilities sum to 1.0.
    ///
    /// # Examples
    ///
    /// ```
    /// use dunbrack::{Residue, Val};
    ///
    /// let rots: Vec<_> = Val::rotamers(-60.0, -40.0).collect();
    /// assert_eq!(rots.len(), 3);
    /// assert!(rots[0].prob > 0.0);
    /// ```
    fn rotamers(phi: f32, psi: f32) -> Self::Iter;
}
