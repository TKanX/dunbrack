# BBDEP2010 Data License & Attribution

The data file (`dunbrack-2010.lib.csv`) and its compiled binary derivatives in this project are sourced from the **2010 Backbone-Dependent Rotamer Library (BBDEP2010)** developed by the Dunbrack Lab at Fox Chase Cancer Center.

## License Declaration

As of July 25, 2019, the BBDEP2010 database is distributed under the **Open Data Commons Attribution License (ODC-By)**.

This dataset is free for both **academic and commercial** use. You are legally free to copy, distribute, and produce derivative works from this database, provided you comply with the attribution requirement below.

## User Obligations

If you integrate this Rust crate into your software, or publish research results based on calculations using this library, you **MUST** explicitly acknowledge the source and cite the following paper:

> Shapovalov, M.V., and Dunbrack, R.L., Jr. (2011). A smoothed backbone-dependent rotamer library for proteins derived from adaptive kernel density estimates and regressions. _Structure_, 19, 844-858.

**Implementation note for downstream developers:**
If your software exposes calculation outputs derived from this crate to end-users, please include the above citation in your software's documentation, `README`, or runtime initialization logs.

## Code vs. Data

- **Data:** The statistical rotamer data in this directory (and its binary equivalents) are subject to the **ODC-By** license.
- **Code:** The Rust source code, algorithms, and binary parsers in this repository are licensed under the **MIT License** (see the root `LICENSE` file).
