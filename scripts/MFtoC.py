#!/usr/bin/env python

"""Convert multiple sequence alignments to counts format.

This script converts the alignments given in a multiple sequence
alignment file to counts format.  It uses the SNP data from variant
call format files.  The SNP data will be merged with the reference
genome to create a complete sequence that contains SNPs as well as
unchanged bases (which are also informative for phylogenetic
analysis).

"""

import argparse
import import_libPoMo  # noqa
import libPoMo.fasta as fa
import libPoMo.cf as cf  # noqa

descr = """convertVCF script version 1.0

This script converts the alignments given in a multiple sequence
alignment file to counts format.  It uses the SNP data from variant
call format files.  The SNP data will be merged with the reference
genome to create a complete sequence that contains SNPs as well as
unchanged bases (which are also informative for phylogenetic
analysis).

Tabix index files need to be provided for all VCF files. They can be
created from the terminal with $(tabix -p vcf "vcf-file.vcf.gz") if
tabix is installed.

The chromosome names given in the multiple alignment reference
(cf. CCDS alignments from UCSC) have to match a chromosome name of the
VCF file.

The command line argument -m collapses all individuals of all VCF
files so that they are represented by a single population for each VCF
file.

The input as well as the output files can additionally be gzipped
(indicated by a .gz file ending).

"""

parser = argparse.ArgumentParser(
    prog="convertVCF",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=descr)

parser.add_argument("reference",
                    help="path to (gzipped) multiple sequence alignment file")
parser.add_argument("VCFFiles", metavar="VCFFileN", nargs='+',
                    help="path to (gzipped) vcf files with SNP information")
parser.add_argument("output",
                    help="name of (gzipped) outputfile in counts format")
parser.add_argument("-m", "--merge", action="count",
                    help="merge individuals within all given VCFFiles")
args = parser.parse_args()

MFaRefFN = args.reference
vcfFL = args.VCFFiles
output = args.output

if args.merge is None:
    cfw = cf.CFWriter(vcfFL, output)
else:
    l = len(vcfFL)
    mergeList = [True for v in vcfFL]
    cfw = cf.CFWriter(vcfFL, output, mergeL=mergeList)

mFaStr = fa.MFaStream(MFaRefFN)

cfw.write_HLn()
cf.write_cf_from_MFaStream(mFaStr, cfw)

cfw.close()
