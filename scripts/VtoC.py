#!/usr/bin/env python

"""Convert VCF files to count format.

This script converts VCF files containing SNP data to a file in count
format. A reference genome needs to be provided in fasta format as
input. The SNP data will be merged with the reference genome to create
a complete sequence that contains SNPs as well as unchanged bases
(which are also informative for phylogenetic analysis).

"""

import argparse
import import_libPoMo  # noqa
import libPoMo.fasta as fa  # noqa
import libPoMo.vcf as vcf  # noqa
import libPoMo.cf as cf  # noqa

descr = """convertVCF script version 1.0

This script converts VCF files containing SNP data to a file in count
format. A reference genome needs to be provided in fasta format as
input. The SNP data will be merged with the reference genome to create
a complete sequence that contains SNPs as well as unchanged bases
(which are also informative for phylogenetic analysis).

The filenames of the reference and the VCF files have to match. This
is a check, if the files are really ment to be compared (the referene
of the VCF file has to match the given reference). If you do not want
to rename your files, run this script with the `-f` option.

The SNPs in the VCF files have to be sorted according to their
chromosomes AND according to their positions on the
chromosomes. Furhtermore chromosome names in the VCF file have to
match the sequence names given in the reference.

Please check these important requirements before you use any data
created by this script.

The input as well as the output files can additionally be gzipped
(indicated by a .gz file ending).

"""

parser = argparse.ArgumentParser(
    prog="convertVCF",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=descr)

parser.add_argument("reference",
                    help="path to (gzipped) reference genome in fasta format")
parser.add_argument("VCFFiles", metavar="VCFFileN", nargs='+',
                    help="path to (gzipped) vcf files with SNP information")
parser.add_argument("output",
                    help="name of (gzipped) outputfile in counts format")
parser.add_argument("-f", "--force", action="count",
                    help="do not check if filenames match")
parser.add_argument("-m", "--merge", action="count",
                    help="merge individuals within all given VCFFiles")
args = parser.parse_args()

vcfStrL = []
sName = None
if args.force is not None:
    sName = 'f'
mergeL = None
nameL = None
if args.merge is not None:
    mergeL = []
    nameL = []
    for fn in args.VCFFiles:
        mergeL.append(True)
        nameL.append(fn.split('.', maxsplit=1)[0])
for fn in args.VCFFiles:
    vcfStrL.append(vcf.init_seq(fn, name=sName))
faR = fa.init_seq(args.reference, name=sName)

cf.save_as_cf(vcfStrL, faR, args.output, mergeL=mergeL, nameL=nameL)

for stream in vcfStrL:
    stream.close_fo()
