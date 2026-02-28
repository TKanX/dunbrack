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

/// Degrees-to-radians conversion factor.
const DEG_TO_RAD: f32 = core::f32::consts::PI / 180.0;

/// Radians-to-degrees conversion factor.
const RAD_TO_DEG: f32 = 180.0 / core::f32::consts::PI;

/// Computes the circular weighted mean of four angles (in degrees).
///
/// Uses sin/cos decomposition to correctly handle the ±180° discontinuity.
/// Returns a value in (−180°, 180°].
fn circular_mean(weights: [f32; 4], angles: [f32; 4]) -> f32 {
    let mut sin_sum = 0.0_f32;
    let mut cos_sum = 0.0_f32;
    for i in 0..4 {
        let rad = angles[i] * DEG_TO_RAD;
        sin_sum += weights[i] * libm::sinf(rad);
        cos_sum += weights[i] * libm::cosf(rad);
    }
    libm::atan2f(sin_sum, cos_sum) * RAD_TO_DEG
}

/// Computes a bilinear combination of four scalar values.
#[inline]
fn bilinear(weights: [f32; 4], values: [f32; 4]) -> f32 {
    weights[0] * values[0]
        + weights[1] * values[1]
        + weights[2] * values[2]
        + weights[3] * values[3]
}
