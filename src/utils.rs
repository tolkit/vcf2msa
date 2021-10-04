// take the genotype
use rust_htslib::bcf::record::{Genotype, GenotypeAllele};
// the ref and alternate allele,
// and make a call. 0 is ref; 1 is alt.
// these should be the only two calls we care about.
// if 0/0 -> return REF
// if 1/1 -> return ALT
// if 0/1 or 1/0 -> return the appropriate IUPAC code

pub fn to_iupac(reference: &str, alternate: &str) -> String {
    match (reference, alternate) {
        ("A", "G") => "R".to_string(),
        ("G", "A") => "R".to_string(),
        ("C", "T") => "Y".to_string(),
        ("T", "C") => "Y".to_string(),
        ("G", "C") => "S".to_string(),
        ("C", "G") => "S".to_string(),
        ("A", "T") => "W".to_string(),
        ("T", "A") => "W".to_string(),
        ("G", "T") => "K".to_string(),
        ("T", "G") => "K".to_string(),
        ("A", "C") => "M".to_string(),
        ("C", "A") => "M".to_string(),
        _ => "N".to_string(),
    }
}

pub fn return_base(gt: Genotype, reference: &str, alternate: &str) -> String {
    let first = match gt.iter().next() {
        Some(allele) => match allele {
            GenotypeAllele::Unphased(x) => Some(*x),
            GenotypeAllele::Phased(x) => Some(*x),
            GenotypeAllele::UnphasedMissing => None,
            GenotypeAllele::PhasedMissing => None,
        },
        None => panic!("Should be an allele here."),
    };
    let second = match gt.iter().next() {
        Some(allele) => match allele {
            GenotypeAllele::Unphased(x) => Some(*x),
            GenotypeAllele::Phased(x) => Some(*x),
            GenotypeAllele::UnphasedMissing => None,
            GenotypeAllele::PhasedMissing => None,
        },
        None => panic!("Should be an allele here for diploid calls."),
    };

    if first.unwrap() == 0 && second.unwrap() == 0 {
        reference.to_string()
    } else if first.unwrap() == 1 && second.unwrap() == 1 {
        alternate.to_string()
    } else {
        to_iupac(reference, alternate)
    }
}
