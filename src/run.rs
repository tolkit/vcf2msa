use bio::io::fasta;
use rust_htslib::bcf::{Read, Reader};
use std::fs::OpenOptions;
use std::fs::{create_dir_all, File};
use std::io::Read as IoRead;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::utils;

// in desperate need of some optimisations

pub fn run(matches: &clap::ArgMatches) {
    println!("[+]\tBegin vcf2msa run.");
    // input paths
    let path = matches.value_of("vcf").unwrap();
    let fasta_path = matches.value_of("fasta").unwrap();
    // outdir
    let outdir = matches.value_of("outdir").unwrap();

    // create directory for output
    if let Err(e) = create_dir_all(outdir) {
        eprintln!("[-]\tCreate directory error: {}", e.to_string());
    }

    let fasta_path = Path::new(fasta_path).canonicalize().unwrap();
    let fasta_reader = fasta::Reader::from_file(&fasta_path).unwrap().records();

    let bcf = Reader::from_path(path).expect("Error opening file.");
    let header = bcf.header();
    // get sample names
    let mut s_nms = String::new();
    let mut s_mns_vec = Vec::new();
    for mut x in header.samples().into_iter() {
        x.read_to_string(&mut s_nms).expect("?");
        s_mns_vec.push(s_nms.clone());
        s_nms.clear();
    }

    // open N files
    // save paths to vec
    let mut path_vec = Vec::new();
    for file in &s_mns_vec {
        // remove all dots and slashes...
        let s = file.replace(&['.', '/'][..], "");
        let path = format!("{}/{}.fasta", outdir, s);
        File::create(&path).expect("Unable to create file");
        path_vec.push(path);
    }

    // init position here?
    let mut pos = 0usize;
    let mut bcf_record_postion = 0usize;

    // iterate over fasta file
    for record in fasta_reader {
        let fasta_record = record.expect("[-]\tError during fasta record parsing.");
        let cur_header = fasta_record.id();

        // add headers for the next record
        for file in &path_vec {
            let f = OpenOptions::new()
                .append(true)
                .open(file)
                .expect("Could not open file.");
            let mut f = BufWriter::new(f);
            let header = format!(">{}\n", cur_header);
            f.write(header.as_bytes()).expect("Unable to write data");

            f.flush().expect("Could not flush.");
        }

        let mut bcf = Reader::from_path(path).expect("Error opening file.");
        // iterate through each row of the vcf body.
        for record_result in bcf.records() {
            // as we iterate over the bcf records,
            // we also want to populate the fasta files dynamically.
            // write first chromosome to each of N files

            // skip indels
            // skip records with more than one alternative allele
            // keep only SNP's.

            let bcf_record = record_result.expect("Fail to read record");
            let contig = match bcf_record.header().rid2name(match bcf_record.rid() {
                Some(rid) => rid,
                None => panic!("[-]\tRecord ID not found."),
            }) {
                Ok(v) => match std::str::from_utf8(v) {
                    Ok(v) => v,
                    Err(e) => panic!("[-]\tInvalid UTF-8 sequence: {}", e),
                },
                Err(e) => panic!("[-]\tInvalid name: {}", e),
            };

            // header and contig name must match.
            if cur_header != contig {
                continue;
            }

            // if ref/alt > 1, skip.
            // these are the indels/polyploid... etc...
            let alleles = bcf_record.alleles();
            if alleles.len() > 2 {
                continue;
            }
            // get alleles
            let ref_allele = std::str::from_utf8(&alleles[0]).unwrap();
            let alt_allele = std::str::from_utf8(&alleles[1]).unwrap();
            // get genotypes
            let genotypes = bcf_record.genotypes().expect("Error reading genotypes");

            for (index, path) in (0..s_mns_vec.len()).zip(path_vec.iter()) {
                let gt = genotypes.get(index);
                // write up to the variant position and the variant itself
                // bcf record, first base has position 1 (so minus one?)
                bcf_record_postion = (bcf_record.pos() as usize) - 1;
                let sequence = match fasta_record.seq().get(pos..bcf_record_postion - 1) {
                    Some(seq) => seq,
                    None => &[],
                };

                // get the file
                let f = OpenOptions::new()
                    .append(true)
                    .open(path)
                    .expect("Could not open file.");
                let mut f = BufWriter::new(f);

                // write the sequence
                f.write(sequence).expect("Unable to write data");
                // write the variant
                let base = utils::return_base(gt, ref_allele, alt_allele);
                f.write(base.as_bytes()).expect("Unable to write data");
                f.flush().expect("Could not flush.");
                // TODO: currently not writing the sequence from the last variant to the chromosome end...
                // annoying but not end-of-the-worldy
            }
            // is this right?
            pos = bcf_record_postion + 1;
        }
        // add newlines to each of the fasta ends.
        for path in path_vec.iter() {
            let f = OpenOptions::new()
                .append(true)
                .open(path)
                .expect("Could not open file.");
            let mut f = BufWriter::new(f);
            // write the newlines
            f.write("\n".as_bytes()).expect("Unable to write data");
        }

        println!("[+]\tFasta ID: {} processed.", cur_header);
        // set position back to zero
        pos = 0;
    }
}
