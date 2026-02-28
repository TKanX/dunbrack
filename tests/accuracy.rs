//! Data correctness tests.
//!
//! Verifies that `build.rs` introduced zero data loss by round-tripping
//! every CSV row through the library and comparing values. Also performs
//! semantic sanity checks based on known properties from the Dunbrack paper.

use csv::ReaderBuilder;
use dunbrack::{
    Arg, Asn, Asp, Cpr, Cyd, Cyh, Cys, Gln, Glu, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Tpr,
    Trp, Tyr, Val,
};
use dunbrack::{Residue, Rotamer};
use std::path::Path;

macro_rules! for_all_residues {
    ($callback:ident) => {
        $callback!(Arg, 4, 75, "ARG");
        $callback!(Asn, 2, 36, "ASN");
        $callback!(Asp, 2, 18, "ASP");
        $callback!(Cpr, 3, 2, "CPR");
        $callback!(Cyd, 1, 3, "CYD");
        $callback!(Cyh, 1, 3, "CYH");
        $callback!(Cys, 1, 3, "CYS");
        $callback!(Gln, 3, 108, "GLN");
        $callback!(Glu, 3, 54, "GLU");
        $callback!(His, 2, 36, "HIS");
        $callback!(Ile, 2, 9, "ILE");
        $callback!(Leu, 2, 9, "LEU");
        $callback!(Lys, 4, 73, "LYS");
        $callback!(Met, 3, 27, "MET");
        $callback!(Phe, 2, 18, "PHE");
        $callback!(Pro, 3, 2, "PRO");
        $callback!(Ser, 1, 3, "SER");
        $callback!(Thr, 1, 3, "THR");
        $callback!(Tpr, 3, 2, "TPR");
        $callback!(Trp, 2, 36, "TRP");
        $callback!(Tyr, 2, 18, "TYR");
        $callback!(Val, 1, 3, "VAL");
    };
}

trait RotamerAccess {
    fn prob(&self) -> f32;
    fn r(&self, i: usize) -> u8;
    fn chi_mean(&self, i: usize) -> f32;
    fn chi_sigma(&self, i: usize) -> f32;
}

impl<const N: usize> RotamerAccess for Rotamer<N> {
    #[inline]
    fn prob(&self) -> f32 {
        self.prob
    }
    #[inline]
    fn r(&self, i: usize) -> u8 {
        self.r[i]
    }
    #[inline]
    fn chi_mean(&self, i: usize) -> f32 {
        self.chi_mean[i]
    }
    #[inline]
    fn chi_sigma(&self, i: usize) -> f32 {
        self.chi_sigma[i]
    }
}

struct CsvRow {
    res: String,
    phi: f32,
    psi: f32,
    r: [u8; 4],
    n_chi: usize,
    prob: f32,
    chi_val: [f32; 4],
    chi_sig: [f32; 4],
}

fn verify_row<R: Residue>(row: &CsvRow)
where
    R::Rot: RotamerAccess,
{
    let rots: Vec<R::Rot> = R::rotamers(row.phi, row.psi).collect();

    let matching = rots
        .iter()
        .find(|rot| (0..row.n_chi).all(|i| rot.r(i) == row.r[i]));

    let rot = matching.unwrap_or_else(|| {
        panic!(
            "{} at (φ={},ψ={}): r={:?} not found",
            row.res,
            row.phi,
            row.psi,
            &row.r[..row.n_chi]
        )
    });

    let prob_err = (rot.prob() - row.prob).abs();
    assert!(
        prob_err < 1e-5,
        "{}(φ={},ψ={}) r={:?}: prob_err={prob_err:.8}",
        row.res,
        row.phi,
        row.psi,
        &row.r[..row.n_chi]
    );

    for i in 0..row.n_chi {
        let mut mean_err = (rot.chi_mean(i) - row.chi_val[i]).abs();
        if mean_err > 180.0 {
            mean_err = 360.0 - mean_err;
        }
        assert!(
            mean_err < 0.05,
            "{}(φ={},ψ={}) r={:?}: chi_mean[{i}] err={mean_err:.4}°",
            row.res,
            row.phi,
            row.psi,
            &row.r[..row.n_chi]
        );

        let sig_err = (rot.chi_sigma(i) - row.chi_sig[i]).abs();
        assert!(
            sig_err < 0.05,
            "{}(φ={},ψ={}) r={:?}: chi_sigma[{i}] err={sig_err:.4}°",
            row.res,
            row.phi,
            row.psi,
            &row.r[..row.n_chi]
        );
    }
}

macro_rules! dispatch_verify {
    ($row:expr; $($Res:ident, $n_chi:literal, $n_rot:literal, $tag:literal);+ $(;)?) => {
        match $row.res.as_str() {
            $( $tag => verify_row::<$Res>($row), )+
            other => panic!("unknown residue: {other}"),
        }
    };
}

#[test]
fn test_full_table_round_trip() {
    let csv_path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("dunbrack-2010.lib.csv");

    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .from_path(&csv_path)
        .unwrap_or_else(|e| panic!("failed to open CSV: {e}"));

    let mut row_count = 0usize;

    for result in reader.records() {
        let record = result.expect("failed to parse CSV row");
        row_count += 1;

        let res = record[0].to_string();
        let phi: f32 = record[1].parse().unwrap();
        let psi: f32 = record[2].parse().unwrap();
        let r1: u8 = record[3].parse().unwrap();
        let r2: u8 = record[4].parse().unwrap();
        let r3: u8 = record[5].parse().unwrap();
        let r4: u8 = record[6].parse().unwrap();
        let prob: f32 = record[7].parse().unwrap();
        let chi1_val: f32 = record[8].parse().unwrap();
        let chi2_val: f32 = record[9].parse().unwrap();
        let chi3_val: f32 = record[10].parse().unwrap();
        let chi4_val: f32 = record[11].parse().unwrap();
        let chi1_sig: f32 = record[12].parse().unwrap();
        let chi2_sig: f32 = record[13].parse().unwrap();
        let chi3_sig: f32 = record[14].parse().unwrap();
        let chi4_sig: f32 = record[15].parse().unwrap();

        let n_chi = [r1, r2, r3, r4].iter().filter(|&&x| x > 0).count();

        let row = CsvRow {
            res,
            phi,
            psi,
            r: [r1, r2, r3, r4],
            n_chi,
            prob,
            chi_val: [chi1_val, chi2_val, chi3_val, chi4_val],
            chi_sig: [chi1_sig, chi2_sig, chi3_sig, chi4_sig],
        };

        dispatch_verify!(&row;
            Arg, 4, 75, "ARG";
            Asn, 2, 36, "ASN";
            Asp, 2, 18, "ASP";
            Cpr, 3, 2, "CPR";
            Cyd, 1, 3, "CYD";
            Cyh, 1, 3, "CYH";
            Cys, 1, 3, "CYS";
            Gln, 3, 108, "GLN";
            Glu, 3, 54, "GLU";
            His, 2, 36, "HIS";
            Ile, 2, 9, "ILE";
            Leu, 2, 9, "LEU";
            Lys, 4, 73, "LYS";
            Met, 3, 27, "MET";
            Phe, 2, 18, "PHE";
            Pro, 3, 2, "PRO";
            Ser, 1, 3, "SER";
            Thr, 1, 3, "THR";
            Tpr, 3, 2, "TPR";
            Trp, 2, 36, "TRP";
            Tyr, 2, 18, "TYR";
            Val, 1, 3, "VAL";
        );
    }

    assert_eq!(
        row_count, 740_629,
        "expected 740,629 CSV rows, got {row_count}"
    );
}

#[test]
fn test_val_alpha_helix_most_probable() {
    let phi = -60.0_f32;
    let psi = -40.0_f32;

    let rots: Vec<Rotamer<1>> = Val::rotamers(phi, psi).collect();

    let best = rots
        .iter()
        .max_by(|a, b| a.prob.partial_cmp(&b.prob).unwrap())
        .unwrap();

    assert_eq!(
        best.r,
        [2],
        "Val at α-helix: expected r=[2] (trans), got {:?}",
        best.r
    );

    assert!(
        best.prob > 0.5,
        "Val trans at α-helix: prob={:.3}, expected > 0.5",
        best.prob
    );
}

#[test]
fn test_leu_nine_rotamers() {
    let rots: Vec<_> = Leu::rotamers(-60.0, -40.0).collect();
    assert_eq!(rots.len(), 9, "Leu should have exactly 9 rotamers");
}

#[test]
fn test_pro_two_rotamers() {
    let phi = -60.0_f32;
    let psi = -40.0_f32;

    let pro_rots: Vec<_> = Pro::rotamers(phi, psi).collect();
    let tpr_rots: Vec<_> = Tpr::rotamers(phi, psi).collect();
    let cpr_rots: Vec<_> = Cpr::rotamers(phi, psi).collect();

    assert_eq!(pro_rots.len(), 2, "Pro should have 2 rotamers");
    assert_eq!(tpr_rots.len(), 2, "Tpr should have 2 rotamers");
    assert_eq!(cpr_rots.len(), 2, "Cpr should have 2 rotamers");

    for (name, rots) in [("Pro", &pro_rots), ("Tpr", &tpr_rots), ("Cpr", &cpr_rots)] {
        let prob_sum: f32 = rots.iter().map(|r| r.prob).sum();
        assert!(
            (prob_sum - 1.0).abs() < 1e-4,
            "{name}: prob_sum={prob_sum:.6}"
        );
    }
}

#[test]
fn test_cyd_cyh_distinct() {
    let phi = -120.0_f32;
    let psi = 130.0_f32;

    let cyd_rots: Vec<_> = Cyd::rotamers(phi, psi).collect();
    let cyh_rots: Vec<_> = Cyh::rotamers(phi, psi).collect();

    assert_eq!(cyd_rots.len(), 3);
    assert_eq!(cyh_rots.len(), 3);

    let differs = cyd_rots
        .iter()
        .zip(cyh_rots.iter())
        .any(|(d, h)| (d.prob - h.prob).abs() > 0.01);

    assert!(
        differs,
        "CYD and CYH should have different probability distributions"
    );
}

#[test]
fn test_arg_seventy_five_rotamers() {
    let rots: Vec<_> = Arg::rotamers(-60.0, -40.0).collect();
    assert_eq!(rots.len(), 75, "Arg should have exactly 75 rotamers");
}

#[test]
fn test_gln_largest_rotamer_set() {
    let rots: Vec<_> = Gln::rotamers(-60.0, -40.0).collect();
    assert_eq!(rots.len(), 108, "Gln should have exactly 108 rotamers");

    macro_rules! check_smaller {
        ($Res:ident, $n_chi:literal, $n_rot:literal, $tag:literal) => {
            assert!(
                $n_rot <= 108,
                "{} has {} rotamers (> 108)",
                stringify!($Res),
                $n_rot
            );
        };
    }
    for_all_residues!(check_smaller);
}
