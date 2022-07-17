# Genomic VCF to tab-separated values

A small script that converts genomic variant data encoded in VCF format into a tab-separated values file. The script utilizes the [brentp/cyvcf2](https://github.com/brentp/cyvcf2) library to parse the VCF file. By default, the program prints the fixed VCF columns, all INFO tag values (as defined in the VCF header, INFO tags not present in a given record is appended with a '.'), and all genotype data (FORMAT columns) for heterozygotes and homozygotes. If genotype data is present, it prints one line per sample, and a column denoted VCF\_SAMPLE_ID indicates data for a given sample. The script has optional arguments to

* skip sample genotype data (i.e. FORMAT colums)
* keep rejected genotypes (i.e. FILTER != 'PASS' / GT == './.')
* skip INFO data.
* compress output TSV
* print data types of VCF columns as a header line

__IMPORTANT__: If you run vcf2tsv with a large multi-sample VCF file, the size of the output TSV will quickly grow large, since there is one line per sample genotype in the output by default. Turn on `--skip_genotype_data` if you are primarily interested in the variant INFO elements, filesize of output will also be considerably smaller.

## Installation

**IN PROGRESS**:
`conda install vcf2tsv`

## Usage:

	usage:
		vcf2tsv
		--input_vcf <INPUT_VCF>
		--out_tsv <OUTPUT_TSV>
		-h [options]

	vcf2tsv: Convert a VCF (Variant Call Format) file with genomic variants to a file with tab-separated values (TSV). One entry (TSV line) per sample genotype.

	Required arguments:
		--input_vcf INPUT_VCF	Bgzipped input VCF file with input variants (SNVs/InDels)
		--out_tsv OUT_TSV     	Output TSV file with one line per non-rejected sample genotype (variant, genotype and annotation data as tab-separated values)

	Optional arguments:
		--skip_info_data      	Skip output of data in INFO column
		--skip_genotype_data  	Skip output of genotype_data (FORMAT columns)
		--keep_rejected_calls	Output data also for rejected (non-PASS) calls
		--print_data_type_header	Output a header line with data types of VCF annotations
		--compress            	Compress output TSV file with gzip
		--version             	Show program's version number and exit
