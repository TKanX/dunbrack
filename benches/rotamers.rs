use criterion::{Criterion, black_box, criterion_group, criterion_main};
use dunbrack::Residue;
use dunbrack::for_all_residues;
use dunbrack::{
    Arg, Asn, Asp, Cpr, Cyd, Cyh, Cys, Gln, Glu, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Tpr,
    Trp, Tyr, Val,
};

const PHI: f32 = -65.0;
const PSI: f32 = -43.0;

const GRID_COUNT: usize = 37;
const GRID_MIN: f32 = -180.0;
const GRID_STEP: f32 = 10.0;

fn bench_query(c: &mut Criterion) {
    let mut group = c.benchmark_group("query");

    macro_rules! add_bench {
        ($Res:ident, $n_chi:literal, $n_rot:literal) => {
            paste::paste! {
                group.bench_function(stringify!([<$Res:lower>]), |b| {
                    b.iter(|| {
                        <$Res>::rotamers(black_box(PHI), black_box(PSI))
                            .for_each(|r| { black_box(r); })
                    });
                });
            }
        };
    }
    for_all_residues!(add_bench);

    group.finish();
}

fn bench_sweep(c: &mut Criterion) {
    let mut group = c.benchmark_group("sweep");

    macro_rules! add_bench {
        ($Res:ident, $n_chi:literal, $n_rot:literal) => {
            paste::paste! {
                group.bench_function(stringify!([<$Res:lower>]), |b| {
                    b.iter(|| {
                        for phi_idx in 0..GRID_COUNT {
                            for psi_idx in 0..GRID_COUNT {
                                let phi = GRID_MIN + phi_idx as f32 * GRID_STEP;
                                let psi = GRID_MIN + psi_idx as f32 * GRID_STEP;
                                <$Res>::rotamers(phi, psi)
                                    .for_each(|r| { black_box(r); });
                            }
                        }
                    });
                });
            }
        };
    }
    for_all_residues!(add_bench);

    group.finish();
}

criterion_group!(benches, bench_query, bench_sweep);
criterion_main!(benches);
