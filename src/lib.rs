#![no_std]

mod grid;
mod interp;
mod residue;
mod rotamer;
mod sealed;

// Generated static tables and trait implementations.
include!(concat!(env!("OUT_DIR"), "/tables.rs"));

pub use interp::RotamerIter;
pub use residue::Residue;
pub use residue::{
    Arg, Asn, Asp, Cpr, Cyd, Cyh, Cys, Gln, Glu, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Tpr,
    Trp, Tyr, Val,
};
pub use rotamer::Rotamer;
