// Max Brown, Wellcome Sanger Institute 2021
// Convert a VCF to multiple sequence alignment.

use clap::{App, Arg};
use std::process;
use vcf2msa::run;

fn main() {
    let matches = App::new("vcf2msa")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Convert a VCF to multiple sequence alignment.")
        .subcommand(
            clap::SubCommand::with_name("run")
                .about("Main program; convert VCF to multiple sequence alignment.")
                .arg(
                    Arg::with_name("vcf")
                        .short("v")
                        .long("vcf")
                        .takes_value(true)
                        .help("The input VCF file."),
                )
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required_unless("print")
                        .help("The input fasta file."),
                ),
        )
        .get_matches();

    let subcommand = matches.subcommand();
    match subcommand.0 {
        "run" => {
            let matches = subcommand.1.unwrap();
            run::run(matches);
        }
        _ => {
            println!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            process::exit(1);
        }
    }
}
