#!/usr/bin/env python

from cyvcf2 import VCF
import numpy as np
import re
import logging
import math
import subprocess
import os
import sys
import errno
from vcf2tsvpy._version import __version__
import argparse
from argparse import RawTextHelpFormatter


def cli():

    program_description = (
        f'vcf2tsvpy: Convert a VCF (Variant Call Format) file with genomic variants to a file with tab-separated values (TSV). One entry (TSV line) per sample genotype.')
    program_options = "\n\t--input_vcf <INPUT_VCF>\n\t--out_tsv <OUTPUT_TSV>\n\t"

    parser = argparse.ArgumentParser(description=program_description,
                                     formatter_class=RawTextHelpFormatter,
                                     usage=f'\n\t%(prog)s {program_options}-h [options] \n\n')

    parser._action_groups.pop()
    required_args = parser.add_argument_group("Required arguments")
    optional_args = parser.add_argument_group("Optional arguments")
    required_args.add_argument(
        '--input_vcf', help='Bgzipped input VCF file with input variants (SNVs/InDels)')
    required_args.add_argument(
        '--out_tsv', help='Output TSV file with one line per non-rejected sample genotype (variant, genotype and annotation data as tab-separated values)')
    optional_args.add_argument(
        "--skip_info_data", action="store_true", help="Skip output of data in INFO column")
    optional_args.add_argument("--skip_genotype_data", action="store_true",
                               help="Skip output of genotype_data (FORMAT columns)")
    optional_args.add_argument("--keep_rejected_calls", action="store_true",
                               help="Output data also for rejected (non-PASS) calls")
    optional_args.add_argument("--print_data_type_header", action="store_true",
                               help="Output a header line with data types of VCF annotations")
    optional_args.add_argument(
        "--compress", action="store_true", help="Compress output TSV file with gzip")
    optional_args.add_argument(
        '--version', action='version', version='%(prog)s ' + str(__version__))
    args = parser.parse_args()
    arg_dict = vars(args)

    logger = getlogger("vcf2tsvpy")

    check_args(arg_dict, logger)
    run_vcf2tsv(arg_dict, logger)


def getlogger(logger_name):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    # add ch to logger
    logger.addHandler(ch)
    # create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")
    # add formatter to ch
    ch.setFormatter(formatter)
    return logger


def error_message(message, logger):
    logger.error("")
    logger.error(message)
    logger.error("")
    sys.exit(1)


def warn_message(message, logger):
    logger.warning(message)

def info_message(message, logger):
    logger.info(message)

def check_args(arg_dict, logger):

    # Check the existence of required arguments
    if arg_dict['input_vcf'] is None or not os.path.exists(arg_dict['input_vcf']):
        err_msg = f"Required argument '--input_vcf' does not exist ({arg_dict['input_vcf']}). Type 'vcf2tsvpy -h' to see required arguments and optional ones."
        error_message(err_msg, logger)

    if arg_dict['out_tsv'] is None:
        err_msg = f"Required argument '--out_tsv' has no/undefined value ({arg_dict['out_tsv']}). Type 'vcf2tsvpy -h' to see required arguments and optional ones."
        error_message(err_msg, logger)


def check_subprocess(command):
    try:
        output = subprocess.check_output(
            str(command), stderr=subprocess.STDOUT, shell=True)
        if len(output) > 0:
            print(str(output.decode()).rstrip())
    except subprocess.CalledProcessError as e:
        print(e.output)
        exit(0)


def run_vcf2tsv(arg_dict, logger):

    vcf = VCF(arg_dict['input_vcf'], gts012=True)

    ofile = None
    try:
        with open(arg_dict['out_tsv'], 'w') as ofile:
            ofile.write(
                '#https://github.com/sigven/vcf2tsvpy version=' + str(__version__) + '\n')

            fixed_columns_header = ['CHROM', 'POS',
                                    'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
            fixed_columns_header_type = ['String', 'Integer',
                                         'String', 'String', 'String', 'Float', 'String']
            samples = vcf.samples           
            vcf_header_tagtypes = {}
            vcf_header_tagtypes['INFO'] = {}
            vcf_header_tagtypes['FORMAT'] = {}
            vcf_header_tags = {}
            vcf_header_tags['INFO'] = []
            vcf_header_tags['FORMAT'] = []
            vcf_header_tags['SAMPLE'] = []
            gt_present_header = 0
            duplicate_tag = {}

            dup_tag_prefix = "INFO"

            if len(samples) > 0:
                vcf_header_tags['SAMPLE'].append('VCF_SAMPLE_ID')

            for e in vcf.header_iter():
                header_element = e.info()

                ## record INFO/FORMAT types (String, Float etc), and append all INFO/FORMAT elements to 
                if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
                    if header_element['HeaderType'] == 'INFO' or header_element['HeaderType'] == 'FORMAT':
                        vcf_header_tagtypes[header_element['HeaderType']][header_element['ID']] = header_element['Type']
                    if header_element['HeaderType'] == 'INFO':
                        if not arg_dict['skip_info_data']:
                            vcf_header_tags['INFO'].append(header_element['ID'])
                    if header_element['HeaderType'] == 'FORMAT':
                        if len(vcf_header_tags['SAMPLE']) > 0 and not arg_dict['skip_genotype_data']:
                            if header_element['ID'] != 'GT':
                                vcf_header_tags['FORMAT'].append(header_element['ID'])
                            else:
                                gt_present_header = 1

            ## check for identical INFO and FORMAT tags, 
            ## ensure there are no duplicate column names in output TSV by prepending any such INFO tags with 'INFO_'
            for info_tag in vcf_header_tags['INFO']:
                for format_tag in vcf_header_tags['FORMAT']:
                    if info_tag == format_tag and info_tag in vcf_header_tagtypes['INFO']: 

                        #duplicate_tags[info_tag] = 1
                        revised_info_tag = f'{dup_tag_prefix}_{info_tag}'

                        ## if the 'INFO_DP' already exists as a tag in the input VCF, change to INFO2_DP
                        if revised_info_tag in vcf_header_tagtypes['INFO']:
                            revised_info_tag = f'{dup_tag_prefix}2_{info_tag}'
                            dup_tag_prefix = "INFO2"

                        duplicate_tag[revised_info_tag] = 1

                        warn_message(message = f'Found "{info_tag}" both as INFO and FORMAT tag - renaming INFO tag to "{revised_info_tag}"', logger = logger)                      
                        ## add new header element and record type
                        vcf_header_tags['INFO'].append(f'{revised_info_tag}')
                        vcf_header_tagtypes['INFO'][f'{revised_info_tag}'] = vcf_header_tagtypes['INFO'][info_tag]
                            
                        ##remove old tag from info column headers and column type dictionary
                        vcf_header_tags['INFO'].remove(info_tag)
                        vcf_header_tagtypes['INFO'].pop(info_tag, None)

            
            ## make String with header tags (FIXED + INFO + FORMAT)
            header_tags = fixed_columns_header
            if not arg_dict['skip_info_data']:
                header_tags = fixed_columns_header + \
                    sorted(vcf_header_tags['INFO'])
                if len(vcf_header_tags['SAMPLE']) > 0:
                    if not arg_dict['skip_genotype_data']:
                        header_tags = fixed_columns_header + \
                            sorted(vcf_header_tags['INFO']) + vcf_header_tags['SAMPLE'] + \
                            sorted(vcf_header_tags['FORMAT']) + ['GT']
                    else:
                        header_tags = fixed_columns_header + \
                            sorted(vcf_header_tags['INFO'])
            else:
                if len(vcf_header_tags['SAMPLE']) > 0:
                    if not arg_dict['skip_genotype_data']:
                        header_tags = fixed_columns_header + vcf_header_tags['SAMPLE'] + \
                            sorted(vcf_header_tags['FORMAT']) + ['GT']
                    else:
                        header_tags = fixed_columns_header
            header_line = '\t'.join(header_tags)

            if arg_dict['print_data_type_header']:
                header_types = []
                for h in header_tags:
                    if h in vcf_header_tagtypes['INFO']:
                        header_types.append(str(vcf_header_tagtypes['INFO'][h]))
                    if h in vcf_header_tagtypes['FORMAT']:
                        header_types.append(str(vcf_header_tagtypes['FORMAT'][h]))
                header_line_type = '\t'.join(
                    fixed_columns_header_type + header_types)
                ofile.write('#' + str(header_line_type) + '\n')
                ofile.write(str(header_line) + '\n')
            else:
                ofile.write(str(header_line) + '\n')

            for rec in vcf:
                rec_id = '.'
                rec_qual = '.'
                rec_filter = '.'
                alt = ",".join(str(n) for n in rec.ALT)
                if not rec.ID is None:
                    rec_id = str(rec.ID)
                if not rec.QUAL is None:
                    rec_qual = str("{0:.2f}".format(rec.QUAL))
                rec_filter = str(rec.FILTER)
                if rec.FILTER is None:
                    rec_filter = 'PASS'

                pos = int(rec.start) + 1
                fixed_fields_string = f'{rec.CHROM}\t{pos}\t{rec_id}\t{rec.REF}\t{alt}\t{rec_qual}\t{rec_filter}'

                if not 'PASS' in rec_filter and not arg_dict['keep_rejected_calls']:
                    continue

                variant_info = rec.INFO
                vcf_info_data = []
                if not arg_dict['skip_info_data']:
                    for info_field in sorted(vcf_header_tags['INFO']):

                        info_field_get = info_field
                        ## change revised duplicate INFO tag back to its origin in order to get variant data stored in rec.INFO
                        if info_field in duplicate_tag:
                            info_field_get = re.sub(r'^' + dup_tag_prefix + '_', '', info_field)

                        if vcf_header_tagtypes['INFO'][info_field] == 'Flag':
                            if variant_info.get(info_field_get) is None:
                                vcf_info_data.append('False')
                            else:
                                vcf_info_data.append('True')
                        elif vcf_header_tagtypes['INFO'][info_field] == 'Float' or vcf_header_tagtypes['INFO'][info_field] == 'Integer' or  \
                            vcf_header_tagtypes['INFO'][info_field] == 'String' or vcf_header_tagtypes['INFO'][info_field] == 'Character':
                            if type(variant_info.get(info_field_get)) is list or type(variant_info.get(info_field_get)) is tuple:
                                vcf_info_data.append(",".join(str(n) for n in variant_info.get(info_field)))
                            else:
                                if variant_info.get(info_field_get) is None:
                                    vcf_info_data.append('.')
                                else:
                                    if vcf_header_tagtypes['INFO'][info_field] == 'Float':
                                        if not isinstance(variant_info.get(info_field_get), float):
                                            warn_msg = f'INFO tag {info_field} is defined in the VCF header as type \'Float\', yet parsed as other type: {type(variant_info.get(info_field_get))}'
                                            warn_message(warn_msg, logger)
                                            if not ',' in str(alt):
                                                warn_msg = f'Multiple values in INFO tag for single ALT allele (VCF multiallelic sites not decomposed properly?): {fixed_fields_string}\t{info_field}\t{variant_info.get(info_field)}'
                                                warn_message(warn_msg, logger)
                                            vcf_info_data.append('.')
                                        else:
                                            val = str("{0:.7f}".format(
                                                variant_info.get(info_field_get)))
                                            vcf_info_data.append(val)
                                    else:
                                        if vcf_header_tagtypes['INFO'][info_field] == 'String' or vcf_header_tagtypes['INFO'][info_field] == 'Character':
                                            if isinstance(variant_info.get(info_field_get), str):
                                                vcf_info_data.append(variant_info.get(info_field_get).encode(
                                                    'ascii', 'ignore').decode('ascii'))
                                            else:
                                                vcf_info_data.append('.')
                                                if vcf_header_tagtypes['INFO'][info_field] == 'String':
                                                    warn_msg = f'INFO tag {info_field} is defined in the VCF header as type \'String\', yet parsed as other type: {type(variant_info.get(info_field_get))}'
                                                    warn_message(
                                                        warn_msg, logger)
                                                if vcf_header_tagtypes['INFO'][info_field] == 'Character':
                                                    warn_msg = f'INFO tag {info_field} is defined in the VCF header as type \'Character\', yet parsed as other type: {type(variant_info.get(info_field_get))}'
                                                    warn_message(
                                                        warn_msg, logger)
                                        else:
                                            if isinstance(variant_info.get(info_field_get), int):
                                                vcf_info_data.append(
                                                    str(variant_info.get(info_field_get)))
                                            else:
                                                warn_msg = f'INFO tag {info_field} is defined in the VCF header as type \'Integer\', yet parsed as other type: {type(variant_info.get(info_field_get))}'
                                                warn_message(warn_msg, logger)
                                                vcf_info_data.append(re.sub(r'\(|\)', '', variant_info.get(info_field_get).encode('ascii', 'ignore').decode('ascii')))

                # dictionary, with sample names as keys, values being genotype data (dictionary with format tags as keys)
                vcf_sample_genotype_data = {}

                if len(samples) > 0 and not arg_dict['skip_genotype_data']:
                    gt_cyvcf = rec.gt_types
                    i = 0
                    while i < len(samples):
                        vcf_sample_genotype_data[samples[i]] = {}
                        gt = './.'
                        if gt_present_header == 1:
                            if gt_cyvcf[i] == 0:
                                gt = '0/0'
                            if gt_cyvcf[i] == 1:
                                gt = '0/1'
                            if gt_cyvcf[i] == 2:
                                gt = '1/1'
                        vcf_sample_genotype_data[samples[i]]['GT'] = gt
                        i = i + 1

                for format_tag in sorted(vcf_header_tags['FORMAT']):
                    if len(samples) > 0 and not arg_dict['skip_genotype_data']:
                        sample_dat = rec.format(format_tag)
                        if sample_dat is None:
                            k = 0
                            while k < len(samples):
                                if samples[k] in vcf_sample_genotype_data:
                                    vcf_sample_genotype_data[samples[k]][format_tag] = '.'
                                k = k + 1
                            continue
                        dim = sample_dat.shape
                        j = 0
                        # sample-wise
                        while j < dim[0]:
                            if sample_dat[j].size > 1:

                                d = ','.join(str(e) for e in np.ndarray.tolist(sample_dat[j]))
                                if vcf_header_tagtypes['FORMAT'][format_tag] == 'Float':
                                    d = ','.join(str(round(e, 4)) for e in np.ndarray.tolist(sample_dat[j]))
                                # undefined/missing value
                                if '-2147483648' in d:
                                    d = d.replace('-2147483648', '.')
                                if 'nan' in d.casefold():
                                    d = d.casefold().replace('nan', '.')
                                if samples[j] in vcf_sample_genotype_data:
                                    vcf_sample_genotype_data[samples[j]
                                                             ][format_tag] = d
                            else:
                                d = '.'
                                if vcf_header_tagtypes['FORMAT'][format_tag] == 'Float':
                                    if not math.isnan(sample_dat[j]):
                                        d = str(sample_dat[j][0])
                                if vcf_header_tagtypes['FORMAT'][format_tag] == 'String':
                                    d = str(sample_dat[j])
                                if vcf_header_tagtypes['FORMAT'][format_tag] == 'Integer':
                                    d = str(sample_dat[j][0])
                                # undefined/missing value
                                if d == '-2147483648' or d.casefold() == 'nan':
                                    d = '.'
                                if samples[j] in vcf_sample_genotype_data:
                                    vcf_sample_genotype_data[samples[j]
                                                             ][format_tag] = d
                            j = j + 1

                tsv_elements = []
                tsv_elements.append(fixed_fields_string)
                if not arg_dict['skip_info_data']:
                    if not arg_dict['skip_genotype_data']:
                        if len(vcf_header_tags['SAMPLE']) > 0:
                            tsv_elements.append("\t".join(str(n)
                                                for n in vcf_info_data))
                            # one line per sample variant
                            for s in sorted(vcf_sample_genotype_data.keys()):
                                sample = s
                                line_elements = []
                                line_elements.extend(tsv_elements)
                                line_elements.append(sample)
                                gt_tag = '.'
                                for tag in sorted(vcf_sample_genotype_data[sample].keys()):
                                    if tag != 'GT':
                                        line_elements.append(vcf_sample_genotype_data[sample][tag].encode(
                                            'ascii', 'ignore').decode('ascii'))
                                    else:
                                        gt_tag = vcf_sample_genotype_data[sample][tag].encode(
                                            'ascii', 'ignore').decode('ascii')
                                line_elements.append(gt_tag)
                                if gt_tag == './.' or gt_tag == '.':
                                    if arg_dict['keep_rejected_calls']:
                                        ofile.write(
                                            '\t'.join(line_elements) + '\n')
                                else:
                                    ofile.write("\t".join(str(n)
                                                for n in line_elements) + '\n')

                        else:
                            tsv_elements.append("\t".join(str(n)
                                                for n in vcf_info_data))
                            line_elements = []
                            line_elements.extend(tsv_elements)
                            ofile.write('\t'.join(line_elements) + '\n')
                    else:
                        tsv_elements.append("\t".join(str(n)
                                            for n in vcf_info_data))
                        line_elements = []
                        line_elements.extend(tsv_elements)
                        ofile.write('\t'.join(line_elements) + '\n')
                else:
                    if not arg_dict['skip_genotype_data']:
                        if len(vcf_header_tags['SAMPLE']) > 0:
                            # one line per sample variant
                            for s in sorted(vcf_sample_genotype_data.keys()):
                                sample = s
                                line_elements = []
                                line_elements.extend(tsv_elements)
                                line_elements.append(sample)
                                gt_tag = '.'
                                for tag in sorted(vcf_sample_genotype_data[sample].keys()):
                                    if tag != 'GT':
                                        line_elements.append(
                                            vcf_sample_genotype_data[sample][tag])
                                    else:
                                        gt_tag = vcf_sample_genotype_data[sample][tag]
                                line_elements.append(gt_tag)
                                if gt_tag == './.' or gt_tag == '.':
                                    if arg_dict['keep_rejected_calls']:
                                        ofile.write(
                                            '\t'.join(line_elements) + '\n')
                                else:
                                    ofile.write(
                                        '\t'.join(line_elements) + '\n')
                    else:
                        line_elements = []
                        line_elements.extend(tsv_elements)
                        line_elements = tsv_elements
                        ofile.write('\t'.join(line_elements) + '\n')

            if not ofile is None:
                ofile.close()

                if arg_dict['compress']:
                    gzip_command = f'gzip -f {arg_dict["out_tsv"]}'
                    check_subprocess(gzip_command)

    except IOError as error:
        print(error)
        # The error itself may include information about why the file write failed
        if error.errno == errno.EACCES:
            print('No permission to write to file ' + str(arg_dict['out_tsv']))
        elif error.errno == errno.EISDIR:
            print('Cannot write to file - a directory with that name exists')

if __name__ == "__main__":
    cli()
