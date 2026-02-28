//! Domain completeness tests.
//!
//! Verifies that the library never panics and always returns well-formed
//! data across the entire defined domain.

use dunbrack::{
    Arg, Asn, Asp, Cpr, Cyd, Cyh, Cys, Gln, Glu, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Tpr,
    Trp, Tyr, Val,
};
use dunbrack::{Residue, Rotamer, for_all_residues};
use proptest::prelude::*;

const GRID_COUNT: usize = 37;
const GRID_MIN: f32 = -180.0;
const GRID_STEP: f32 = 10.0;

#[inline]
fn grid_angle(idx: usize) -> f32 {
    GRID_MIN + idx as f32 * GRID_STEP
}

fn assert_well_formed<const N: usize>(rots: &[Rotamer<N>], label: &str, n_rotamers: usize) {
    assert_eq!(
        rots.len(),
        n_rotamers,
        "{label}: expected {n_rotamers} rotamers, got {}",
        rots.len()
    );

    let prob_sum: f32 = rots.iter().map(|r| r.prob).sum();
    assert!(
        (prob_sum - 1.0).abs() < 1e-4,
        "{label}: prob_sum={prob_sum:.6}, expected ≈ 1.0"
    );

    for (k, rot) in rots.iter().enumerate() {
        assert!(
            (0.0..=1.0).contains(&rot.prob),
            "{label}[{k}]: prob={} out of [0,1]",
            rot.prob
        );

        for i in 0..N {
            assert!(
                rot.chi_sigma[i] > 0.0,
                "{label}[{k}]: chi_sigma[{i}]={} must be > 0",
                rot.chi_sigma[i]
            );

            assert!(
                rot.chi_mean[i].abs() < 360.0,
                "{label}[{k}]: chi_mean[{i}]={} out of bounds",
                rot.chi_mean[i]
            );

            assert!(rot.r[i] >= 1, "{label}[{k}]: r[{i}]=0 (must be 1-based)");
        }
    }
}

macro_rules! grid_coverage_test {
    ($Res:ident, $n_chi:literal, $n_rot:literal) => {
        paste::paste! {
            #[test]
            fn [<test_grid_coverage_ $Res:lower>]() {
                for phi_idx in 0..GRID_COUNT {
                    for psi_idx in 0..GRID_COUNT {
                        let phi = grid_angle(phi_idx);
                        let psi = grid_angle(psi_idx);
                        let rots: Vec<Rotamer<$n_chi>> =
                            <$Res>::rotamers(phi, psi).collect();
                        let label = format!(
                            "{}(φ={},ψ={})",
                            stringify!($Res), phi, psi
                        );
                        assert_well_formed(&rots, &label, $n_rot);
                    }
                }
            }
        }
    };
}

for_all_residues!(grid_coverage_test);

macro_rules! fuzz_no_panic_test {
    ($Res:ident, $n_chi:literal, $n_rot:literal) => {
        paste::paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(10_000))]

                #[test]
                fn [<test_fuzz_no_panic_ $Res:lower>](
                    phi in -200.0_f32..200.0,
                    psi in -200.0_f32..200.0
                ) {
                    let rots: Vec<Rotamer<$n_chi>> =
                        <$Res>::rotamers(phi, psi).collect();

                    prop_assert_eq!(
                        rots.len(), $n_rot,
                        "{}: expected {} rotamers, got {}",
                        stringify!($Res), $n_rot, rots.len()
                    );

                    let prob_sum: f32 = rots.iter().map(|r| r.prob).sum();
                    prop_assert!(
                        (prob_sum - 1.0).abs() < 1e-3,
                        "{}: prob_sum={:.6} (expected ≈ 1.0)",
                        stringify!($Res), prob_sum
                    );
                }
            }
        }
    };
}

for_all_residues!(fuzz_no_panic_test);
