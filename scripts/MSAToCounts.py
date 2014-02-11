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
import os
import libPoMo.fasta as fa
import libPoMo.cf as cf  # noqa

descr = """Convert multiple sequence alignments to counts format.

This script converts the alignments given in a multiple sequence
alignment file to counts format.  It uses the SNP data from bgzipped
variant call format files that need to be aligned to the reference
used in the multiple sequence alignment file.  The SNP data will be
merged with the reference genome to create a complete sequence that
contains SNPs as well as unchanged bases (which are also informative
for phylogenetic analysis).

To speed up data conversion, tabix index files need to be provided for
all VCF files. They can be created from the terminal with $(tabix -p
vcf "vcf-file.vcf.gz") if tabix is installed.  The files to be indexed
have to be zipped with bgzip.

In order to correctly match the SNPs to the reference, the chromosome
names given in the multiple alignment reference (cf. CCDS alignments
from UCSC) have to match the chromosome names of the VCF file.

The command line argument -m collapses all individuals of all VCF
files so that they are represented by a single population for each VCF
file.

The input as well as the output files can additionally be gzipped
(indicated by a .gz file ending).

"""

parser = argparse.ArgumentParser(
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
parser.add_argument("-s", "--synonymous", action="count",
                    help="only print position with a synonymous base")
parser.add_argument("-v", "--verbosity", action="count",
                    help="turn on verbosity")
args = parser.parse_args()

MFaRefFN = args.reference
vcfFnL = args.VCFFiles
output = args.output
vb = args.verbosity

if args.merge is None:
    cfw = cf.CFWriter(vcfFnL, output, verb=vb)
else:
    mergeList = []
    nameList = []
    for fn in vcfFnL:
        mergeList.append(True)
        strippedFn = os.path.basename(fn)
        nameList.append(strippedFn.split('.', maxsplit=1)[0])
    cfw = cf.CFWriter(vcfFnL, output, mergeL=mergeList,
                      nameL=nameList, verb=vb)

if args.synonymous is not None:
    cfw.onlySynonymous = True

mFaStr = fa.MFaStream(MFaRefFN)

cfw.write_HLn()
cf.write_cf_from_MFaStream(mFaStr, cfw)

cfw.close()
