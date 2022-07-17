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
from vcf2tsv._version import __version__
import argparse
from argparse import RawTextHelpFormatter

def cli():

   program_description = (f'vcf2tsv: Convert a VCF (Variant Call Format) file with genomic variants to a file with tab-separated values (TSV). One entry (TSV line) per sample genotype.')
   program_options = "\n\t--input_vcf <INPUT_VCF>\n\t--out_tsv <OUTPUT_TSV>\n\t"

   parser = argparse.ArgumentParser(description=program_description,
                                     formatter_class=RawTextHelpFormatter,
                                     usage=f'\n\t%(prog)s {program_options}-h [options] \n\n')
    
   parser._action_groups.pop()
   required_args = parser.add_argument_group("Required arguments")
   optional_args = parser.add_argument_group("Optional arguments")
   parser.add_argument_group(required_args)
   required_args.add_argument('--input_vcf', help='Bgzipped input VCF file with input variants (SNVs/InDels)')
   required_args.add_argument('--out_tsv', help='Output TSV file with one line per non-rejected sample genotype (variant, genotype and annotation data as tab-separated values)')
   optional_args.add_argument("--skip_info_data",action = "store_true", help="Skip output of data in INFO column")
   optional_args.add_argument("--skip_genotype_data", action="store_true", help="Skip output of genotype_data (FORMAT columns)")
   optional_args.add_argument("--keep_rejected_calls", action="store_true", help="Output data also for rejected (non-PASS) calls")
   optional_args.add_argument("--print_data_type_header", action="store_true", help="Output a header line with data types of VCF annotations")
   optional_args.add_argument("--compress", action="store_true", help="Compress output TSV file with gzip")
   optional_args.add_argument('--version', action='version', version='%(prog)s ' + str(__version__))
   args = parser.parse_args()
   arg_dict = vars(args)
   #print(arg_dict)
   
   logger = getlogger("vcf2tsv")

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

def check_args(arg_dict, logger):

    # Check the existence of required arguments
    if arg_dict['input_vcf'] is None or not os.path.exists(arg_dict['input_vcf']):
        err_msg = f"Required argument '--input_vcf' does not exist ({arg_dict['input_vcf']}). Type 'vcf2tsv -h' to see required arguments and optional ones."
        error_message(err_msg, logger)

    if arg_dict['out_tsv'] is None:
        err_msg = f"Required argument '--out_tsv' has no/undefined value ({arg_dict['out_tsv']}). Type 'vcf2tsv -h' to see required arguments and optional ones."
        error_message(err_msg, logger)


def check_subprocess(command):
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print (str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print (e.output)
      exit(0)


def run_vcf2tsv(arg_dict, logger):
   
   vcf = VCF(arg_dict['query_vcf'], gts012 = True)
   
   try:
      with open(arg_dict['out_tsv'], 'w') as ofile:
        ofile.write('#https://github.com/sigven/vcf2tsv version=' + str(__version__) + '\n')
        #file.close()
   except IOError as error: 
      # You could also catch Exception instead of IOError to check for problems but this may be casting too-wide a net
      # The file was not accessible for some reason
      # Do something here to alert the user or perform another action instead
      print(error)
     # The error itself may include information about why the file write failed
      if error.errno == errno.EACCES:
         print('No permission to write to file')
      elif error.errno == errno.EISDIR:
         print('Cannot write to file - a directory with that name exists')

   fixed_columns_header = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER']
   fixed_columns_header_type = ['String','Integer','String','String','String','Float','String']
   samples = vcf.samples
   info_columns_header = []
   format_columns_header = []
   sample_columns_header = []
   column_types = {}
   gt_present_header = 0
   
   if len(samples) > 0:
      sample_columns_header.append('VCF_SAMPLE_ID')
   
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO' or header_element['HeaderType'] == 'FORMAT':
            column_types[header_element['ID']] = header_element['Type']
         if header_element['HeaderType'] == 'INFO':
            if arg_dict['skip_info_data'] is False:
               info_columns_header.append(header_element['ID'])
         if header_element['HeaderType'] == 'FORMAT':
            if len(sample_columns_header) > 0 and arg_dict['skip_genotype_data'] is False:
               if header_element['ID'] != 'GT':
                  format_columns_header.append(header_element['ID'])
               else:
                  gt_present_header = 1

   header_tags = fixed_columns_header
   if arg_dict['skip_info_data'] is False:
      header_tags = fixed_columns_header + sorted(info_columns_header)
      if len(sample_columns_header) > 0:
         if arg_dict['skip_genotype_data'] is False:
            header_tags = fixed_columns_header + sorted(info_columns_header) + sample_columns_header + sorted(format_columns_header) + ['GT']
         else:
            header_tags = fixed_columns_header + sorted(info_columns_header)
   else:
      if len(sample_columns_header) > 0:
         if arg_dict['skip_genotype_data'] is False:
            header_tags = fixed_columns_header + sample_columns_header + sorted(format_columns_header) + ['GT']
         else:
            header_tags = fixed_columns_header
   header_line = '\t'.join(header_tags)
   
   if arg_dict['print_data_type_header'] is True:
      #header_tags = header_line.rstrip().split('\t')
      header_types = []
      for h in header_tags:
         if h in column_types:
            header_types.append(str(column_types[h]))
      header_line_type = '\t'.join(fixed_columns_header_type + header_types)
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
      #fixed_fields_string = str(rec.CHROM) + '\t' + str(pos) + '\t' + str(rec_id) + '\t' + str(rec.REF) + '\t' + str(alt) + '\t' + str(rec_qual) + '\t' + str(rec_filter)
            
      if not 'PASS' in rec_filter and not arg_dict['keep_rejected_calls']:
         continue
      
      variant_info = rec.INFO
      vcf_info_data = []
      if arg_dict['skip_info_data'] is False:
         for info_field in sorted(info_columns_header):
            if column_types[info_field] == 'Flag':
               if variant_info.get(info_field) is None:
                  vcf_info_data.append('False')
               else:
                  vcf_info_data.append('True')
            elif column_types[info_field] == 'Float' or column_types[info_field] == 'Integer' or column_types[info_field] == 'String' or column_types[info_field] == 'Character':
               if type(variant_info.get(info_field)) is list or type(variant_info.get(info_field)) is tuple:
                  vcf_info_data.append(",".join(str(n) for n in variant_info.get(info_field)))
               else:
                  if variant_info.get(info_field) is None:
                     vcf_info_data.append('.')
                  else:
                     if column_types[info_field] == 'Float':
                        if not isinstance(variant_info.get(info_field), float):
                           warn_msg = f'INFO tag {info_field} is defined in the VCF header as type \'Float\', yet parsed as other type: {type(variant_info.get(info_field))}'
                           warn_message(warn_msg, logger)
                           if not ',' in str(alt):
                              warn_msg = f'Multiple values in INFO tag for single ALT allele (VCF multiallelic sites not decomposed properly?): {fixed_fields_string}\t{info_field}\t{variant_info.get(info_field)}'
                              warn_message(warn_msg, logger)
                           vcf_info_data.append('.')
                        else:
                           val = str("{0:.7f}".format(variant_info.get(info_field)))
                           vcf_info_data.append(val)
                     else:
                        if column_types[info_field] == 'String' or column_types[info_field] == 'Character':
                           if isinstance(variant_info.get(info_field), str):
                              vcf_info_data.append(variant_info.get(info_field).encode('ascii','ignore').decode('ascii'))
                           else:
                              vcf_info_data.append('.')
                              if column_types[info_field] == 'String':
                                    warn_msg = f'INFO tag {info_field} is defined in the VCF header as type \'String\', yet parsed as other type: {type(variant_info.get(info_field))}'
                                    warn_message(warn_msg, logger)
                              if column_types[info_field] == 'Character':
                                    warn_msg = f'INFO tag {info_field} is defined in the VCF header as type \'Character\', yet parsed as other type: {type(variant_info.get(info_field))}'
                                    warn_message(warn_msg, logger)
                        else:
                           if isinstance(variant_info.get(info_field), int):
                              vcf_info_data.append(str(variant_info.get(info_field)))
                           else:
                              warn_msg = f'INFO tag {info_field} is defined in the VCF header as type \'Integer\', yet parsed as other type: {type(variant_info.get(info_field))}'
                              warn_message(warn_msg, logger)
                              vcf_info_data.append(re.sub(r'\(|\)', '', variant_info.get(info_field).encode('ascii','ignore').decode('ascii')))

      
      #dictionary, with sample names as keys, values being genotype data (dictionary with format tags as keys)
      vcf_sample_genotype_data = {}

      if len(samples) > 0 and arg_dict['skip_genotype_data'] is False:
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
               
      for format_tag in sorted(format_columns_header):
         if len(samples) > 0 and arg_dict['skip_genotype_data'] is False:
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
            ## sample-wise
            while j < dim[0]:
               if sample_dat[j].size > 1:

                  d = ','.join(str(e) for e in np.ndarray.tolist(sample_dat[j]))
                  if column_types[format_tag] == 'Float':
                     d = ','.join(str(round(e, 4)) for e in np.ndarray.tolist(sample_dat[j]))
                  ## undefined/missing value
                  if '-2147483648' in d:
                     d = d.replace('-2147483648', '.')
                  if 'nan' in d.casefold():
                     d = d.casefold().replace('nan', '.')
                  if samples[j] in vcf_sample_genotype_data:
                     vcf_sample_genotype_data[samples[j]][format_tag] = d
               else:
                  d = '.'
                  if column_types[format_tag] == 'Float':
                     if not math.isnan(sample_dat[j]):
                        d = str(sample_dat[j][0])
                  if column_types[format_tag] == 'String':
                     d = str(sample_dat[j])
                  if column_types[format_tag] == 'Integer':
                     d = str(sample_dat[j][0])
                  ## undefined/missing value
                  if d == '-2147483648' or d.casefold() == 'nan':
                     d = '.'
                  if samples[j] in vcf_sample_genotype_data:
                     vcf_sample_genotype_data[samples[j]][format_tag] = d
               j = j + 1
      
      tsv_elements = []
      tsv_elements.append(fixed_fields_string)
      if arg_dict['skip_info_data'] is False:
         if arg_dict['skip_genotype_data'] is False:
            if len(sample_columns_header) > 0:
               tsv_elements.append("\t".join(str(n) for n in vcf_info_data))
               ## one line per sample variant
               for s in sorted(vcf_sample_genotype_data.keys()):
                  sample = s
                  line_elements = []
                  line_elements.extend(tsv_elements)
                  line_elements.append(sample)
                  gt_tag = '.'
                  for tag in sorted(vcf_sample_genotype_data[sample].keys()):
                     if tag != 'GT':
                        line_elements.append(vcf_sample_genotype_data[sample][tag].encode('ascii','ignore').decode('ascii'))
                     else:
                        gt_tag = vcf_sample_genotype_data[sample][tag].encode('ascii','ignore').decode('ascii')
                  line_elements.append(gt_tag)
                  if gt_tag == './.' or gt_tag == '.':
                     if arg_dict['keep_rejected_calls']:
                        ofile.write('\t'.join(line_elements) + '\n')
                  else:
                     ofile.write("\t".join(str(n) for n in line_elements) + '\n')
                    
            else:
               tsv_elements.append("\t".join(str(n) for n in vcf_info_data))
               line_elements = []
               line_elements.extend(tsv_elements)
               ofile.write('\t'.join(line_elements) + '\n')
         else:
            tsv_elements.append("\t".join(str(n) for n in vcf_info_data))
            line_elements = []
            line_elements.extend(tsv_elements)
            ofile.write('\t'.join(line_elements) + '\n')
      else:
         if arg_dict['skip_genotype_data'] is False:
            if len(sample_columns_header) > 0:
               ## one line per sample variant
               for s in sorted(vcf_sample_genotype_data.keys()):
                  sample = s
                  line_elements = []
                  line_elements.extend(tsv_elements)
                  line_elements.append(sample)
                  gt_tag = '.'
                  for tag in sorted(vcf_sample_genotype_data[sample].keys()):
                     if tag != 'GT':
                        line_elements.append(vcf_sample_genotype_data[sample][tag])
                     else:
                        gt_tag = vcf_sample_genotype_data[sample][tag]
                  line_elements.append(gt_tag)
                  if gt_tag == './.' or gt_tag == '.':
                     if arg_dict['keep_rejected_calls']:
                        ofile.write('\t'.join(line_elements) + '\n')
                  else:
                     ofile.write('\t'.join(line_elements) + '\n')
         else:
            line_elements = []
            line_elements.extend(tsv_elements)
            line_elements = tsv_elements
            ofile.write('\t'.join(line_elements) + '\n')
       
   ofile.close()
   
   if arg_dict['compress'] is True:
      command = 'gzip -f ' + str(arg_dict['out_tsv'])
      check_subprocess(command)

if __name__ == "__main__":
    cli()