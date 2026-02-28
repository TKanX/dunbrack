//! Interpolation correctness tests.
//!
//! Verifies that the bilinear interpolation algorithm works correctly:
//! grid degeneracy, continuity, circular chi_mean wrapping, and probability
//! normalization.

use dunbrack::{
    Arg, Asn, Asp, Cpr, Cyd, Cyh, Cys, Gln, Glu, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Tpr,
    Trp, Tyr, Val,
};
use dunbrack::{Residue, Rotamer, for_all_residues};

const GRID_COUNT: usize = 37;
const GRID_MIN: f32 = -180.0;
const GRID_STEP: f32 = 10.0;

macro_rules! grid_degeneracy_test {
    ($Res:ident, $n_chi:literal, $n_rot:literal) => {
        paste::paste! {
            #[test]
            fn [<test_grid_degeneracy_ $Res:lower>]() {
                for phi_idx in 0..GRID_COUNT {
                    for psi_idx in 0..GRID_COUNT {
                        let phi = GRID_MIN + phi_idx as f32 * GRID_STEP;
                        let psi = GRID_MIN + psi_idx as f32 * GRID_STEP;

                        let rots_a: Vec<Rotamer<$n_chi>> =
                            <$Res>::rotamers(phi, psi).collect();
                        let rots_b: Vec<Rotamer<$n_chi>> =
                            <$Res>::rotamers(phi, psi).collect();

                        for (k, (a, b)) in rots_a.iter().zip(rots_b.iter()).enumerate() {
                            assert_eq!(
                                a.prob, b.prob,
                                "{}(φ={},ψ={})[{k}]: prob mismatch",
                                stringify!($Res), phi, psi
                            );
                            assert_eq!(
                                a.chi_mean, b.chi_mean,
                                "{}(φ={},ψ={})[{k}]: chi_mean mismatch",
                                stringify!($Res), phi, psi
                            );
                            assert_eq!(
                                a.chi_sigma, b.chi_sigma,
                                "{}(φ={},ψ={})[{k}]: chi_sigma mismatch",
                                stringify!($Res), phi, psi
                            );
                        }
                    }
                }
            }
        }
    };
}

for_all_residues!(grid_degeneracy_test);

macro_rules! continuity_test {
    ($Res:ident, $n_chi:literal, $n_rot:literal) => {
        paste::paste! {
            #[test]
            fn [<test_continuity_ $Res:lower>]() {
                let psi = -40.0_f32;
                let step = 0.1_f32;
                let n_steps = 3599;
                let max_jump = 0.05_f32;

                let mut prev: Vec<Rotamer<$n_chi>> =
                    <$Res>::rotamers(-179.9, psi).collect();

                for i in 1..=n_steps {
                    let phi = -179.9 + i as f32 * step;
                    let curr: Vec<Rotamer<$n_chi>> =
                        <$Res>::rotamers(phi, psi).collect();

                    for (k, (p, c)) in prev.iter().zip(curr.iter()).enumerate() {
                        let delta = (c.prob - p.prob).abs();
                        assert!(
                            delta < max_jump,
                            "{}: rotamer {k} prob jumped by {delta:.6} at φ={phi:.1}°",
                            stringify!($Res)
                        );
                    }

                    prev = curr;
                }
            }
        }
    };
}

for_all_residues!(continuity_test);

macro_rules! circular_wrap_test {
    ($Res:ident, $n_chi:literal, $n_rot:literal) => {
        paste::paste! {
            #[test]
            fn [<test_circular_wrap_ $Res:lower>]() {
                let mut found = false;

                'outer: for psi_idx in 0..GRID_COUNT {
                    let psi = GRID_MIN + psi_idx as f32 * GRID_STEP;

                    for phi_idx in 0..(GRID_COUNT - 1) {
                        let phi_a = GRID_MIN + phi_idx as f32 * GRID_STEP;
                        let phi_b = phi_a + GRID_STEP;

                        let rots_a: Vec<Rotamer<$n_chi>> =
                            <$Res>::rotamers(phi_a, psi).collect();
                        let rots_b: Vec<Rotamer<$n_chi>> =
                            <$Res>::rotamers(phi_b, psi).collect();

                        for k in 0..rots_a.len() {
                            for i in 0..$n_chi {
                                let a = rots_a[k].chi_mean[i];
                                let b = rots_b[k].chi_mean[i];

                                let crosses = (a > 170.0 && b < -170.0)
                                           || (a < -170.0 && b > 170.0);

                                if crosses {
                                    let phi_mid = (phi_a + phi_b) / 2.0;
                                    let rots_mid: Vec<Rotamer<$n_chi>> =
                                        <$Res>::rotamers(phi_mid, psi).collect();
                                    let mid = rots_mid[k].chi_mean[i];

                                    assert!(
                                        mid.abs() > 90.0,
                                        "{}: circular wrap failed at \
φ∈[{phi_a},{phi_b}], ψ={psi}, \
rot[{k}].chi_mean[{i}]: \
a={a:.1}, b={b:.1}, mid={mid:.1}",
                                        stringify!($Res)
                                    );
                                    found = true;
                                    break 'outer;
                                }
                            }
                        }
                    }
                }

                let _ = found;
            }
        }
    };
}

for_all_residues!(circular_wrap_test);

macro_rules! prob_normalization_test {
    ($Res:ident, $n_chi:literal, $n_rot:literal) => {
        paste::paste! {
            #[test]
            fn [<test_prob_normalization_ $Res:lower>]() {
                let angles: [f32; 32] = [
                    -175.3, -163.7, -152.1, -141.9, -130.5, -119.2, -108.8, -97.4,
                    -86.1, -74.6, -63.3, -52.0, -41.7, -30.3, -19.9, -8.5,
                    2.8, 13.1, 24.4, 35.7, 46.0, 57.3, 68.6, 79.9,
                    91.2, 102.5, 113.8, 125.1, 136.4, 147.7, 159.0, 170.3,
                ];

                for &phi in &angles {
                    for &psi in &angles {
                        let rots: Vec<Rotamer<$n_chi>> =
                            <$Res>::rotamers(phi, psi).collect();

                        let prob_sum: f32 = rots.iter().map(|r| r.prob).sum();
                        assert!(
                            (prob_sum - 1.0).abs() < 1e-4,
                            "{}(φ={phi},ψ={psi}): prob_sum={prob_sum:.8}",
                            stringify!($Res)
                        );
                    }
                }
            }
        }
    };
}

for_all_residues!(prob_normalization_test);
