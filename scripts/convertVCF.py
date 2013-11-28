#!/usr/bin/env python

"""Convert VCF files to count format.

This script converts a VCF file containing SNP data to a file in count
format. A reference genome needs to be provided in fasta format as
input. The SNP data will be merged with the reference genome to create
a complete sequence that contains SNPs as well as unchanged bases
(which are also informative for phylogenetic analysis).

"""

import argparse
import import_libPoMo  # noqa
import libPoMo.fasta as fa  # noqa
import libPoMo.vcf as vcf  # noqa

ver = "1.0"

parser = argparse.ArgumentParser(prog='convertVCF',
                                 description="convertVCF script version" + ver)
parser.add_argument('reference',
                    help="path to reference genome in fasta format")
parser.add_argument('vcfFile',
                    help="path to vcf file with SNP information")
args = parser.parse_args()

seq = vcf.open_vcf(args.vcfFile)
searched_base = seq.get_nuc_base('chr1', 824126)
searched_base.print_info()
