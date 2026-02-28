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

    /// Three-letter residue name (uppercase ASCII, e.g. `"ARG"`, `"VAL"`).
    const NAME: &'static str;

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

/// Arginine (4 χ angles, 75 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Arg;

/// Asparagine (2 χ angles, 36 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Asn;

/// Aspartate (2 χ angles, 18 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Asp;

/// Cis-proline (3 χ angles, 2 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Cpr;

/// Disulfide-bonded cysteine (1 χ angle, 3 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Cyd;

/// Free (non-disulfide) cysteine (1 χ angle, 3 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Cyh;

/// Combined cysteine pool (1 χ angle, 3 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Cys;

/// Glutamine (3 χ angles, 108 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Gln;

/// Glutamate (3 χ angles, 54 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Glu;

/// Histidine (2 χ angles, 36 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct His;

/// Isoleucine (2 χ angles, 9 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Ile;

/// Leucine (2 χ angles, 9 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Leu;

/// Lysine (4 χ angles, 73 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Lys;

/// Methionine (3 χ angles, 27 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Met;

/// Phenylalanine (2 χ angles, 18 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Phe;

/// Combined proline pool (3 χ angles, 2 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Pro;

/// Serine (1 χ angle, 3 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Ser;

/// Threonine (1 χ angle, 3 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Thr;

/// Trans-proline (3 χ angles, 2 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Tpr;

/// Tryptophan (2 χ angles, 36 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Trp;

/// Tyrosine (2 χ angles, 18 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Tyr;

/// Valine (1 χ angle, 3 rotamers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Val;
