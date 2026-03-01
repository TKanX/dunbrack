# Dunbrack

**A zero-cost Rust interface to the Dunbrack 2010 backbone-dependent rotamer library.**

Provides bilinearly interpolated side-chain rotamer probabilities, mean χ angles, and standard deviations for 22 amino acid types at any (φ, ψ) backbone conformation. All 740,629 source rows are baked into `.rodata` at compile time; queries touch zero heap memory and link zero runtime dependencies.

[Features](#features) • [Installation](#installation) • [Usage](#usage) • [Residue Types](#residue-types) • [Performance](#performance) • [Verification](#verification) • [License](#license)

---

## Features

- **Zero startup latency.** The entire ~28 MB rotamer database is embedded in `.rodata` at compile time via `build.rs`. No file I/O, no deserialization, no lazy initialization.
- **Zero heap allocation.** Every query returns a `RotamerIter<N, R>` — a stack-allocated array of exactly `R` `Rotamer<N>` values. No `Vec`, no `Box`, no allocator required.
- **`#![no_std]` compatible.** No standard library, no libm linkage. Usable in embedded firmware, OS kernels, and WASM environments.
- **Type-safe χ dimensionality.** The number of χ angles per residue is a compile-time constant `N` encoded in `Rotamer<N>` and `RotamerIter<N, R>`. There are no padding zeros, no runtime bounds checks, no wrong-length arrays.
- **Bilinear interpolation with circular χ means.** `Residue::rotamers(phi, psi)` bilinearly interpolates across the four surrounding grid cells. χ means are computed via circular weighted mean (sin/cos decomposition), correctly handling the ±180° wraparound. Probabilities are re-normalized to Σ = 1.0 after interpolation.
- **Precomputed (sin χ, cos χ) in the static table.** `build.rs` stores sin/cos pairs rather than raw angles, eliminating 8N trigonometric calls per query (4 sin + 4 cos per χ angle, per corner cell).
- **Custom branchless `atan2f`.** A two-stage argument-reduction + degree-7 Taylor polynomial implementation with zero conditional branches and ±0.002° maximum error — 25× more accurate than the 0.05° precision requirement, with no libm dependency.
- **Compile-time data integrity.** `build.rs` asserts five invariants before emitting any code: rotamer count, probability sums, φ/ψ = ±180° periodicity, and bin index consistency across all 1,369 grid cells. Compilation fails loudly on data corruption.
- **`for_all_residues!` macro.** A generated declarative macro for writing generic code over all 22 residue types without runtime dispatch.

---

## Installation

```toml
[dependencies]
dunbrack = "0.1.0"
```

**Note:** `build.rs` reads `data/dunbrack-2010.lib.csv` (740,629 rows) and generates ~28 MB of static Rust source. Initial compilation takes 15–30 seconds depending on hardware.

---

## Usage

### Basic Query

```rust
use dunbrack::{Residue, Val};

// Bilinearly interpolated rotamers for Val at α-helical backbone.
for rot in Val::rotamers(-60.0, -40.0) {
    // rot.r:         [u8; 1]  — rotamer bin index (1-based)
    // rot.prob:      f32      — probability (Σ = 1.0 across all rotamers)
    // rot.chi_mean:  [f32; 1] — mean χ angle in degrees, ±180° range
    // rot.chi_sigma: [f32; 1] — standard deviation in degrees
    println!("r={:?}  p={:.4}  χ₁={:.1}°±{:.1}°",
        rot.r, rot.prob, rot.chi_mean[0], rot.chi_sigma[0]);
}
```

Output (Val at φ=−60°, ψ=−40°):

```text
r=[1]  p=0.0414  χ₁=68.0°±7.0°
r=[2]  p=0.9391  χ₁=171.5°±5.0°
r=[3]  p=0.0194  χ₁=-61.0°±9.6°
```

### Generic Usage

The `Residue` trait exposes compile-time constants usable in fully generic code:

```rust
use dunbrack::Residue;

fn rotamer_count<R: Residue>() -> usize {
    R::N_ROTAMERS
}

fn residue_name<R: Residue>() -> &'static str {
    R::NAME
}
```

Accessing rotamer fields (`.prob`, `.chi_mean`, `.chi_sigma`, `.r`) requires a concrete type or a monomorphized context, since `Residue::Rot` carries no field bounds:

```rust
use dunbrack::{Residue, Val};

// Collect and find the most probable rotamer for Val.
let best = Val::rotamers(-60.0, -40.0)
    .max_by(|a, b| a.prob.partial_cmp(&b.prob).unwrap())
    .unwrap();
```

### `for_all_residues!` Macro

This macro invokes `$callback!(Type, N_CHI, N_ROTAMERS)` for all 22 residue types. It drives generic infrastructure like benchmarks, coverage tests, and per-type dispatch with zero boilerplate.

```rust
use dunbrack::*;

macro_rules! print_info {
    ($Res:ident, $n_chi:literal, $n_rot:literal) => {
        println!("{}: {} χ angles, {} rotamers",
            <$Res as Residue>::NAME, $n_chi, $n_rot);
    };
}

for_all_residues!(print_info);
```

---

## Residue Types

All 22 residue types from the Dunbrack 2010 library, including separated cysteine and proline variants:

| Type  | `N_CHI` | `N_ROTAMERS` | Notes                              |
| :---- | :-----: | :----------: | :--------------------------------- |
| `Arg` |    4    |      75      |                                    |
| `Asn` |    2    |      36      |                                    |
| `Asp` |    2    |      18      |                                    |
| `Gln` |    3    |     108      | Largest table                      |
| `Glu` |    3    |      54      |                                    |
| `His` |    2    |      36      |                                    |
| `Ile` |    2    |      9       |                                    |
| `Leu` |    2    |      9       |                                    |
| `Lys` |    4    |      73      |                                    |
| `Met` |    3    |      27      |                                    |
| `Phe` |    2    |      18      |                                    |
| `Ser` |    1    |      3       |                                    |
| `Thr` |    1    |      3       |                                    |
| `Trp` |    2    |      36      |                                    |
| `Tyr` |    2    |      18      |                                    |
| `Val` |    1    |      3       |                                    |
| `Cyh` |    1    |      3       | Free (non-disulfide) cysteine      |
| `Cyd` |    1    |      3       | Disulfide-bonded cysteine          |
| `Cys` |    1    |      3       | Combined cysteine pool (CYH + CYD) |
| `Tpr` |    3    |      2       | Trans-proline                      |
| `Cpr` |    3    |      2       | Cis-proline                        |
| `Pro` |    3    |      2       | Combined proline pool (TPR + CPR)  |

Each type implements `Residue + Copy + PartialEq + Eq + Hash + Debug`.

---

## Performance

Benchmarked with `Criterion.rs` on an Intel® Core™ i7-13620H (Raptor Lake, 4.90 GHz turbo, AVX2), Linux, `opt-level=3, lto=true, codegen-units=1`.

**Single-point query** — time to call `Residue::rotamers(phi, psi)` and consume the full iterator:

| Residue | N_CHI | N_ROTAMERS |           Time |      Throughput |
| :------ | :---: | :--------: | -------------: | --------------: |
| `Val`   |   1   |     3      |    **31.1 ns** | **32.1 MOps/s** |
| `Ser`   |   1   |     3      |        32.0 ns |     31.2 MOps/s |
| `Pro`   |   3   |     2      |        43.7 ns |     22.9 MOps/s |
| `Leu`   |   2   |     9      |       143.5 ns |     6.97 MOps/s |
| `Phe`   |   2   |     18     |       268.0 ns |     3.73 MOps/s |
| `Met`   |   3   |     27     |       358.0 ns |     2.79 MOps/s |
| `Asn`   |   2   |     36     |       551.7 ns |     1.81 MOps/s |
| `Glu`   |   3   |     54     |       673.4 ns |     1.49 MOps/s |
| `Arg`   |   4   |     75     |     1,142.6 ns |     0.88 MOps/s |
| `Lys`   |   4   |     73     |     1,158.9 ns |     0.86 MOps/s |
| `Gln`   |   3   |    108     | **1,316.7 ns** | **0.76 MOps/s** |

Query time scales linearly with `N_ROTAMERS` at ~12–16 ns per rotamer, dominated by `atan2f` calls (one per χ angle per rotamer).

**Full grid sweep** (all 37×37 = 1,369 cells, sustained throughput):

| Residue |       Time | Per-query | Table size |
| :------ | ---------: | --------: | ---------: |
| `Val`   |    38.2 µs |   27.9 ns |     64 KiB |
| `Gln`   | 1,834.8 µs |  1,340 ns |  5,776 KiB |

Per-query time drops ~10% in sweep mode due to cache warmth across adjacent cells.

For full data including all 22 residues and methodology, see [BENCHMARKS.md](BENCHMARKS.md).

### Why it's fast

| Optimization                            | Savings                                                                                       |
| :-------------------------------------- | :-------------------------------------------------------------------------------------------- |
| Precomputed `(sin χ, cos χ)` in table   | Eliminates 8N trig calls per query (e.g. 32 calls → 4 for Arg, N=4)                           |
| Custom branchless `atan2f`              | Eliminates libm overhead; zero branch-prediction penalties                                    |
| Compile-time static tables (`build.rs`) | Zero startup cost; OS can share read-only pages across processes                              |
| `KEYS` deduplication                    | _bin indices_ stored once per residue (in `KEYS`), not per cell; saves ~401 KiB for Arg alone |
| Stack-only `RotamerIter<N, R>`          | No allocator, no pointer indirection; `next()` is a single array read + increment             |

---

## Verification

The library is verified at three levels:

**Compile time (`build.rs` assertions)** — compilation aborts if any of the following fail:

- Rotamer count per cell matches the registered `N_ROTAMERS`
- Probability sum per cell ∈ [0.99, 1.01]
- φ = −180° and φ = +180° cells are bitwise identical (periodic boundary)
- ψ = −180° and ψ = +180° cells are bitwise identical
- bin index key sets are identical across all 1,369 cells for each residue

**Unit tests** (21 tests in `src/`):

- `atan2f` accuracy: maximum error 3.5×10⁻⁵ rad (±0.002°) over a dense grid
- Circular mean with ±180° wraparound
- `angle_to_grid` at boundaries, midpoints, and out-of-range inputs

**Integration tests** (140 tests in `tests/`):

| File               | Tests | What is verified                                                                                                                                  |
| :----------------- | ----: | :------------------------------------------------------------------------------------------------------------------------------------------------ |
| `accuracy.rs`      |     8 | Full 740,629-row CSV round-trip; `\|prob_err\| < 1e-5`, `\|chi_mean_err\| < 0.05°` at every grid point                                            |
| `coverage.rs`      |    44 | Every (residue, φ, ψ) combination on the 37×37 grid: correct count, Σprob ≈ 1.0, valid ranges; 220,000 random (φ, ψ) fuzz inputs per residue type |
| `interpolation.rs` |    88 | Determinacy, continuity (Δprob < 0.05 per 0.1° step), circular χ wrap correctness, normalization at 32×32 off-grid angles                         |

The `atan2f` error of ±0.002° is 25× below the 0.05° accuracy threshold, meaning the precision ceiling is the source data (CSV values are stored to one decimal place), not the implementation.

Run the full suite:

```bash
cargo test
```

---

## License

MIT — see [LICENSE](LICENSE).
