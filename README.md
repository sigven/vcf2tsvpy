# Genomic VCF to tab-separated values
__Python script for conversion of VCF data to tab-separated values (TSV)__

A small script that converts genomic variant data encoded in VCF format into a tab-separated values file. The script utilizes [brentp/cyvcf2](https://github.com/brentp/cyvcf2) to parse the VCF file. By default, the program prints the fixed VCF columns, all INFO tag values (as defined in the VCF header, INFO tags not present in a given record is appended with a '.'), and all genotype data (FORMAT columns) for heterozygotes and homozygotes. If genotype data is present, it prints one line per sample, and a column denoted VCF\_SAMPLE_ID indicates data for a given sample. The script has optional arguments to

* skip sample genotype data (i.e. FORMAT colums)
* keep rejected genotypes (i.e. FILTER != 'PASS' / GT == './.')
* skip INFO data.
* compress output TSV
* print data types of VCF columns as a header line

## Installation

Running vcf2tsv requires Python 3. It also requires that you have [cyvcf2](https://github.com/brentp/cyvcf2) and [numpy](https://scipy.org/install.html) installed. This can be achieved through the use of [pip](https://pip.pypa.io/en/stable).

## Usage:

	usage: vcf2tsv.py [-h] [--skip_info_data] [--skip_genotype_data]
		   [--keep_rejected_calls] [--print_data_type_header]
		   [--compress]
		   query_vcf out_tsv

	Convert a VCF file with genomic variants to a file with tab-separated values
	(TSV). One entry (TSV line) per sample genotype

	positional arguments:
	query_vcf             Bgzipped input VCF file with query variants
			    (SNVs/InDels)
	out_tsv               Output TSV file with one line pr non-rejected sample
			    genotype (Variant, genotype and annotation data as
			    tab-separated values)

	optional arguments:
	-h, --help            show this help message and exit
	--skip_info_data      Skip printing of data in INFO column (default: False)
	--skip_genotype_data  Skip printing of genotype_data (FORMAT columns)
			    (default: False)
	--keep_rejected_calls
			    Print data for rejected calls (default: False)
	--print_data_type_header
			    Print a header line with data types of VCF annotations
			    (default: False)
	--compress            Compress TSV file with gzip (default: False)
	--version             show program's version number and exit
