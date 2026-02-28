#[cfg(feature = "std")]
mod backend {
    #[inline(always)]
    pub fn sinf(x: f32) -> f32 {
        x.sin()
    }

    #[inline(always)]
    pub fn cosf(x: f32) -> f32 {
        x.cos()
    }

    #[inline(always)]
    pub fn atan2f(y: f32, x: f32) -> f32 {
        y.atan2(x)
    }
}

#[cfg(not(feature = "std"))]
mod backend {
    #[inline(always)]
    pub fn sinf(x: f32) -> f32 {
        ::libm::sinf(x)
    }

    #[inline(always)]
    pub fn cosf(x: f32) -> f32 {
        ::libm::cosf(x)
    }

    #[inline(always)]
    pub fn atan2f(y: f32, x: f32) -> f32 {
        ::libm::atan2f(y, x)
    }
}

pub use backend::{atan2f, cosf, sinf};
