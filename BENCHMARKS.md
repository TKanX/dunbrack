# Performance Benchmarks

This document details the performance metrics for the `dunbrack` crate. All benchmarks were conducted using the `Criterion.rs` framework with 100 samples per benchmark, ensuring statistical significance.

## Executive Summary

The Dunbrack 2010 rotamer library achieves extreme throughput by leveraging compile-time data embedding, zero-allocation iterators, precomputed sin/cos tables, and a custom branchless `atan2f` implementation.

- **Fastest Query**: Val (1 χ angle, 3 rotamers) — **31.1 ns** (~32 million queries/s)
- **Largest Table**: Gln (3 χ angles, 108 rotamers) — **1.32 µs** (~759k queries/s)
- **Most Complex**: Arg (4 χ angles, 75 rotamers) — **1.14 µs** (~875k queries/s)
- **Zero Allocations**: All results returned via stack-allocated iterators
- **Zero Runtime Dependencies**: Pure `#![no_std]` with no libm linkage

## Detailed Results

### Single-Point Query (`query/*`)

Time to call `Residue::rotamers(phi, psi)` and iterate through all rotamers at a single (φ, ψ) coordinate. This measures the complete hot path: grid indexing, bilinear interpolation, circular mean computation, and probability normalization.

| Residue | N_CHI | N_ROTAMERS |   Time (ns) | Throughput (MOps/s) |
| :------ | :---: | :--------: | ----------: | ------------------: |
| Val     |   1   |     3      |    **31.1** |            **32.1** |
| Ser     |   1   |     3      |        32.0 |                31.2 |
| Cys     |   1   |     3      |        32.0 |                31.2 |
| Cyd     |   1   |     3      |        32.0 |                31.2 |
| Cyh     |   1   |     3      |        32.2 |                31.1 |
| Thr     |   1   |     3      |        32.1 |                31.1 |
| Pro     |   3   |     2      |        43.7 |                22.9 |
| Tpr     |   3   |     2      |        43.8 |                22.9 |
| Cpr     |   3   |     2      |        44.3 |                22.6 |
| Leu     |   2   |     9      |       143.5 |                6.97 |
| Ile     |   2   |     9      |       149.9 |                6.67 |
| Phe     |   2   |     18     |       268.0 |                3.73 |
| Tyr     |   2   |     18     |       273.5 |                3.66 |
| Asp     |   2   |     18     |       281.4 |                3.55 |
| Met     |   3   |     27     |       358.0 |                2.79 |
| Trp     |   2   |     36     |       535.0 |                1.87 |
| His     |   2   |     36     |       547.4 |                1.83 |
| Asn     |   2   |     36     |       551.7 |                1.81 |
| Glu     |   3   |     54     |       673.4 |                1.49 |
| Arg     |   4   |     75     |     1,142.6 |                0.88 |
| Lys     |   4   |     73     |     1,158.9 |                0.86 |
| Gln     |   3   |    108     | **1,316.7** |            **0.76** |

### Full Grid Sweep (`sweep/*`)

Time to query all 37 × 37 = 1,369 grid points for a single residue type. This stress-tests cache behavior and measures sustained throughput.

| Residue | N_CHI | N_ROTAMERS | Time (µs) | Per-Query (ns) | Table Size |
| :------ | :---: | :--------: | --------: | -------------: | ---------: |
| Thr     |   1   |     3      |      37.6 |           27.5 |     64 KiB |
| Cyh     |   1   |     3      |      37.4 |           27.3 |     64 KiB |
| Val     |   1   |     3      |      38.2 |           27.9 |     64 KiB |
| Ser     |   1   |     3      |      38.4 |           28.1 |     64 KiB |
| Cyd     |   1   |     3      |      38.2 |           27.9 |     64 KiB |
| Cys     |   1   |     3      |      38.1 |           27.9 |     64 KiB |
| Cpr     |   3   |     2      |      61.3 |           44.8 |    107 KiB |
| Pro     |   3   |     2      |      61.7 |           45.1 |    107 KiB |
| Tpr     |   3   |     2      |      62.0 |           45.3 |    107 KiB |
| Leu     |   2   |     9      |     204.7 |          149.5 |    337 KiB |
| Ile     |   2   |     9      |     208.4 |          152.3 |    337 KiB |
| Asp     |   2   |     18     |     391.8 |          286.3 |    674 KiB |
| Phe     |   2   |     18     |     394.2 |          288.0 |    674 KiB |
| Tyr     |   2   |     18     |     389.5 |          284.6 |    674 KiB |
| Met     |   3   |     27     |     498.5 |          364.2 |  1,443 KiB |
| Trp     |   2   |     36     |     752.8 |          550.0 |  1,348 KiB |
| His     |   2   |     36     |     772.8 |          564.6 |  1,348 KiB |
| Asn     |   2   |     36     |     788.3 |          576.0 |  1,348 KiB |
| Glu     |   3   |     54     |     951.6 |          695.2 |  2,888 KiB |
| Lys     |   4   |     73     |   1,592.1 |        1,163.0 |  5,075 KiB |
| Arg     |   4   |     75     |   1,596.5 |        1,166.2 |  5,214 KiB |
| Gln     |   3   |    108     |   1,834.8 |        1,340.4 |  5,776 KiB |

## Test Environment

- **CPU**: Intel® Core™ i7-13620H (Raptor Lake)
  - 10 Cores (6P + 4E), 16 Threads
  - Max Turbo Frequency: 4.90 GHz
  - Instruction Set: AVX2, FMA3
- **OS**: Linux (Arch Linux, Kernel 6.12.63-1)
- **Profile**: `opt-level=3, lto=true, codegen-units=1`
- **Date**: February 2026

## Performance Analysis

### 1. Linear Scaling with Rotamer Count

Query time scales linearly with `N_ROTAMERS`, as expected from the O(R × N) complexity:

```text
Val  (R=3):     31 ns  →  10.3 ns/rotamer
Leu  (R=9):    143 ns  →  15.9 ns/rotamer
Asn  (R=36):   552 ns  →  15.3 ns/rotamer
Arg  (R=75):  1143 ns  →  15.2 ns/rotamer
Gln  (R=108): 1317 ns  →  12.2 ns/rotamer
```

The per-rotamer cost is remarkably consistent at ~12–16 ns, with smaller residues showing slightly higher overhead due to fixed setup costs.

### 2. Precomputed Sin/Cos Eliminates 8N Trig Calls

By storing `(sin χᵢ, cos χᵢ)` pairs in `GridEntry` instead of raw angles, we reduce per-query trigonometric calls from 9N (4 sin + 4 cos + 1 atan2 per χ angle) to just N (1 atan2 per χ angle). For Arg (N=4):

- **Before optimization**: 36 trig calls/query
- **After optimization**: 4 atan2 calls/query
- **Speedup factor**: ~3× on trig-heavy residues

### 3. Custom `atan2f` Delivers Sub-Microsecond Latency

The custom Taylor-series `atan2f` implementation:

- Uses 7th-degree polynomial with exact rational coefficients
- Achieves ±0.002° maximum error (25× better than the 0.05° precision requirement)
- Contains **zero branches** (fully branchless via `copysign`, `min`, `max`)
- Eliminates libm dependency, enabling `#![no_std]`

### 4. Compact Memory Layout

`GridEntry<N>` stores `prob + chi_sin + chi_cos + chi_sigma` — all `f32`, size `4 + 12N`, **zero padding**. The static table is a contiguous `f32` array with no slack. The const-generic `N` monomorphizes each residue to its exact field count; no worst-case allocation.

The static table holds raw `GridEntry` data only — no pre-built rotamers or iterators. `rotamers(phi, psi)` eagerly interpolates all `R` rotamers into a stack array on every call (one-time O(R×N) cost); `next()` is then a branch-free index increment with zero FP work.

`Rotamer<N>` carries `r: [u8; N]` before `prob: f32`, so `#[repr(C)]` inserts up to 3 bytes of padding — but this only affects the stack output, not the static table.

```text
Query working set (4 corners × R entries): 192 B (Val) … 17,280 B (Gln)
Sweep table size (37 × 37 × R entries): 64 KiB (Val) … 5,776 KiB (Gln)
```

### 5. Zero Heap Allocation

All rotamer results are returned via `RotamerIter<N, R>`, a stack-allocated array of `Rotamer<N>` structs. The largest iterator (Gln) occupies:

```text
[Rotamer<3>; 108] = 108 × 32 bytes = 3,456 bytes
```

This fits comfortably on the stack, avoiding heap allocation overhead entirely.

### 6. Compile-Time Data Embedding

The entire 28 MB rotamer database is embedded in `.rodata` at compile time:

- **Startup latency**: 0 ms (no runtime parsing)
- **Memory mapping**: OS can share read-only pages across processes
- **Cache warmth**: First query may incur page faults; subsequent queries are hot

## Reproducibility

All benchmarks can be reproduced by running:

```bash
cargo bench --bench rotamers
```

Full HTML reports are generated in `target/criterion/`.
