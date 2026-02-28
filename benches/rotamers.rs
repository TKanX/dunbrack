use criterion::{Criterion, black_box, criterion_group, criterion_main};
use dunbrack::Residue;
use dunbrack::{
    Arg, Asn, Asp, Cpr, Cyd, Cyh, Cys, Gln, Glu, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Tpr,
    Trp, Tyr, Val,
};

const PHI: f32 = -65.0;
const PSI: f32 = -43.0;

const GRID_COUNT: usize = 37;
const GRID_MIN: f32 = -180.0;
const GRID_STEP: f32 = 10.0;

macro_rules! for_all_residues {
    ($callback:ident) => {
        $callback!(Arg, 4, 75);
        $callback!(Asn, 2, 36);
        $callback!(Asp, 2, 18);
        $callback!(Cpr, 3, 2);
        $callback!(Cyd, 1, 3);
        $callback!(Cyh, 1, 3);
        $callback!(Cys, 1, 3);
        $callback!(Gln, 3, 108);
        $callback!(Glu, 3, 54);
        $callback!(His, 2, 36);
        $callback!(Ile, 2, 9);
        $callback!(Leu, 2, 9);
        $callback!(Lys, 4, 73);
        $callback!(Met, 3, 27);
        $callback!(Phe, 2, 18);
        $callback!(Pro, 3, 2);
        $callback!(Ser, 1, 3);
        $callback!(Thr, 1, 3);
        $callback!(Tpr, 3, 2);
        $callback!(Trp, 2, 36);
        $callback!(Tyr, 2, 18);
        $callback!(Val, 1, 3);
    };
}

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
