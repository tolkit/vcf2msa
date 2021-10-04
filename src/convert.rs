// convert the `run` output to MSA
// take matching chromosome names from different files
// and place into one file with the name of the file == chromosome name
// clap can take multiple arguments (e.g. *.fasta)

// take the first fasta, get the header names
// make the files named with the header names.

// iterate over each fasta
// if entry header matches contig name
// write to file

use bio::io::fasta;
use std::io::Write;
use std::path::Path;
use std::{
    fs::{create_dir_all, File},
    io::BufWriter,
};

pub fn convert(matches: &clap::ArgMatches) {
    println!("[+]\tBegin vcf2msa convert.");
    // parse command line args
    let outdir = matches.value_of("outdir").unwrap();
    // create directory for output
    match outdir {
        "." => {}
        other => {
            if let Err(e) = create_dir_all(other) {
                eprintln!("[-]\tCreate directory error: {}", e.to_string());
            }
        }
    }

    let mut paths = Vec::new();

    if let Some(i) = matches.values_of("fastas") {
        for el in i {
            // save paths
            paths.push(el);
        }
    }

    // now let's do something with the paths
    // iterate over the first fasta
    let paths_clone = paths.clone();
    let fasta_reader = fasta::Reader::from_file(paths_clone[0]).unwrap().records();

    // assuming the ID's of the first fasta are in all the others
    let mut id_vec = Vec::new();

    for record in fasta_reader {
        let record = record.expect("[-]\tError during fasta record parsing.");
        // copy to String
        id_vec.push(record.id().to_string());
    }

    // make a file for each output contig
    for id in &id_vec {
        let _path = format!("{}/{}.fasta", outdir, id);
        let current_file = File::create(&_path).expect("[-]\tUnable to create file");
        let mut current_file = BufWriter::new(current_file);

        // now iterate over each fasta in turn
        for path in &paths {
            let fasta_reader = fasta::Reader::from_file(path).unwrap().records();
            for record in fasta_reader {
                let record = record.expect("[-]\tError during fasta record parsing.");
                if !record.seq().len() > 0 {
                    continue;
                }
                // if the iteration of the fasta == header name
                if record.id() == id {
                    writeln!(
                        current_file,
                        ">{}\n{}",
                        // janky as hell.
                        Path::new(path).file_stem().unwrap().to_str().unwrap(),
                        std::str::from_utf8(record.seq()).unwrap(),
                    )
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
                }
            }
            println!(
                "[+]\t\tEntry: {} finished.",
                Path::new(path).file_stem().unwrap().to_str().unwrap()
            );
        }
        println!("[+]\tID: {} finished.", id);
    }
}
