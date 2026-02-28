use crate::rotamer::Rotamer;

/// Eagerly constructed iterator over bilinearly interpolated rotamers.
///
/// All `R` interpolated values (including probability re-normalization)
/// are computed eagerly in the constructor. Calling [`Iterator::next`]
/// simply returns the next pre-computed entry — no floating-point work,
/// no branches.
///
/// This type is produced by [`Residue::rotamers`](crate::Residue::rotamers)
/// and should not be constructed directly.
///
/// # Examples
///
/// ```
/// use dunbrack::{Residue, Val};
///
/// let mut iter = Val::rotamers(-60.0, -40.0);
/// assert_eq!(iter.len(), 3);
/// let first = iter.next().unwrap();
/// assert!(first.prob > 0.0);
/// ```
pub struct RotamerIter<const N: usize, const R: usize> {
    items: [Rotamer<N>; R],
    idx: usize,
}

impl<const N: usize, const R: usize> Iterator for RotamerIter<N, R> {
    type Item = Rotamer<N>;

    #[inline]
    fn next(&mut self) -> Option<Rotamer<N>> {
        if self.idx < R {
            let item = self.items[self.idx];
            self.idx += 1;
            Some(item)
        } else {
            None
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = R - self.idx;
        (remaining, Some(remaining))
    }
}

impl<const N: usize, const R: usize> ExactSizeIterator for RotamerIter<N, R> {
    #[inline]
    fn len(&self) -> usize {
        R - self.idx
    }
}

impl<const N: usize, const R: usize> core::iter::FusedIterator for RotamerIter<N, R> {}
