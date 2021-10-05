# vcf2msa

Given an input reference fasta file and a VCF, make a multiple sequence alignment.

No testing, no guarantees. Usual rust installation and build.

Currently only concerns diploid SNPs. Indels/polyploid calls are ignored.

## Installation

```bash
git clone https://github.com/tolkit/vcf2msa && cd vcf2msa && cargo build --release
```

It's a quick rust solution, though the code itself is not optimised, and pretty slow.

## Usage

`./vcf2msa run`

```bash
vcf2msa-run 
Main program; convert VCF to multiple sequence alignment.

USAGE:
    vcf2msa run [OPTIONS] --fasta <fasta>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>      The input fasta file.
    -o, --outdir <outdir>    The name of the output directory. [default: vcf2msa]
    -v, --vcf <vcf>          The input VCF file.
```

Currently outputs a bunch of fastas in the executed dir, one for each sample in the VCF. Working on combining these into actual MSA's.

The below should do this. Again yet to test properly. E.g:

`vcf2msa convert -f *.fasta -o outdir`

```bash
vcf2msa-convert 
Convert output fastas to one fasta per chromosome.

USAGE:
    vcf2msa convert [OPTIONS] --fastas <fastas>...

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fastas <fastas>...    The input fasta files.
    -o, --outdir <outdir>       The name of the output directory. [default: .]
```
