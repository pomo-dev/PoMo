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
import libPoMo.countsformat as cf  # noqa

descr = """convertVCF script version 1.0

This script converts a VCF file containing SNP data to a file in count
format. A reference genome needs to be provided in fasta format as
input. The SNP data will be merged with the reference genome to create
a complete sequence that contains SNPs as well as unchanged bases
(which are also informative for phylogenetic analysis).

The sequence names of the reference fasta file have to match the
chromosome names in the VCF file.

The filenames without extensions have to match if option --force is
not given.

"""

parser = argparse.ArgumentParser(prog="convertVCF",
                                 description=descr)
parser.add_argument("reference",
                    help="path to reference genome in fasta format")
parser.add_argument("VCFFile",
                    help="path to vcf file with SNP information")
parser.add_argument("output",
                    help="name of outputfile in count format")
parser.add_argument("-f", "--force", action="count",
                    help="do not check if filenames match")
args = parser.parse_args()

if args.force is None:
    vcfS = vcf.open_seq(args.VCFFile)
    faR = fa.open_seq(args.reference)
else:
    vcfS = vcf.open_seq(args.VCFFile, name='f')
    faR = fa.open_seq(args.reference, name='f')
cf.save_as_countsformat(vcfS, faR, args.output)
