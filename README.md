# vcf2tsvpy: genomic VCF to tab-separated values (TSV)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/vcf2tsvpy/badges/version.svg)](https://conda.anaconda.org/bioconda) &nbsp;[![Anaconda-Server Badge](https://anaconda.org/bioconda/vcf2tsvpy/badges/latest_release_date.svg)](https://anaconda.org/bioconda/vcf2tsvpy)

A small Python program that converts genomic variant data encoded in [VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) into a tab-separated values (TSV) file.

The program utilizes the [cyvcf2](https://github.com/brentp/cyvcf2) library to parse the VCF file. By default, the program prints the fixed VCF columns, all `INFO` tag values (as defined in the VCF header, `INFO` tags not present in a given record is appended with a `'.'`), and all genotype data (FORMAT columns) for heterozygotes and homozygotes. If genotype data is present, it prints one line per sample, and a column denoted `VCF_SAMPLE_ID` indicates data for a given sample. Importantly, the program has optional arguments to

-   skip sample genotype data (i.e. `FORMAT` colums) - print only variant information
-   keep rejected genotypes (i.e. `FILTER` != `'PASS'` / `GT` == `'./.'`)
-   skip `INFO` data.
-   compress output TSV
-   print data types of VCF columns as a header line

<br>

**IMPORTANT**: If you run *vcf2tsvpy* with a large multi-sample VCF file, the file size of the output TSV will quickly grow fairly large, since there is, by default, one line per sample genotype in the output. Turn on `--skip_genotype_data` if you are primarily interested in the variant INFO elements, file size of output TSV will also be considerably smaller.

## News

* March 9th 2023: **0.6.1 release**
  - Handling of cases where a tag is found __both__ in `INFO` and `FORMAT` columns of VCF (e.g. `DP`). For such cases, the `INFO` tag name is now prepended with a *INFO_* string (e.g. `INFO_DP`), ensuring non-duplicate columns in the final output TSV file.

## Installation

The software can be installed with the [Conda](https://docs.conda.io/en/latest/) package manager, using the following command:

`conda install -c bioconda vcf2tsvpy`

## Usage

    vcf2tsvpy --input_vcf <INPUT_VCF> --out_tsv <OUTPUT_TSV>
           -h [options]

    vcf2tsvpy:  Convert a VCF (Variant Call Format) file with genomic variants to a file with
            tab-separated values (TSV). One entry (TSV line) per sample genotype.

    Required arguments:
        --input_vcf INPUT_VCF   Bgzipped input VCF file with input variants (SNVs/InDels)
        --out_tsv OUT_TSV       Output TSV file with one line per non-rejected sample genotype
                            (variant, genotype and annotation data as tab-separated values)

    Optional arguments:
        --skip_info_data        Skip output of data in INFO column
        --skip_genotype_data    Skip output of genotype_data (FORMAT columns)
        --keep_rejected_calls   Output data also for rejected (non-PASS) calls
        --print_data_type_header    Output a header line with data types of VCF annotations
        --compress              Compress output TSV file with gzip
        --version               Show program's version number and exit
