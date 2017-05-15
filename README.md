# vcf2tsv
Python script for conversion of VCF data to tab-separated values (TSV)

A small script that converts genomic variant data encoded in VCF format into a tab-separated values. The script utilizes [brentp/cyvcf2](https://github.com/brentp/cyvcf2) to parse the VCF file. By default, the program prints the fixed VCF columns, all INFO tag values (as defined in the VCF header), and all genotype data (FORMAT columns) for heterozygotes/homozygotes. If genotype data is present, it prints one line per sample. It has optional arguments to 

* skip sample genotype data (i.e. FORMAT colums)
* keep rejected genotypes (i.e. genotypes with GT = './.') 
* skip INFO data.

## Usage:

		usage: vcf2tsv.py [-h] [--skip_info_data] [--skip_genotype_data]
                  [--keep_rejected_calls]
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
 
	

