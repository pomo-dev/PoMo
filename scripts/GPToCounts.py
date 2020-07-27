#!/usr/bin/env python

"""Convert gene prediction file with reference to counts format.

EXPERIMENTAL.  This script converts the gene predictions given in a
gene prediction file together with a reference into to counts format.
It uses SNP data from variant call format files.  The SNP data will be
merged with the reference genome to create a complete sequence that
contains SNPs as well as unchanged bases (which are also informative
for phylogenetic analysis).

"""

import argparse
import logging
import libPoMo.gp as gp
import libPoMo.cf as cf

descr = """Convert gene prediction file with reference to counts format.

EXPERIMENTAL.  This script converts the gene predictions given in a
gene prediction file together with a reference into to counts format.
It uses SNP data from variant call format files.  The SNP data will be
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

The input as well as the output files can additionally be gzipped
(indicated by a .gz file ending).

"""

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=descr)

parser.add_argument("gpfile", help="path to gp file")
parser.add_argument("reference",
                    help="path to (gzipped) multiple sequence alignment file")
parser.add_argument("VCFFiles", metavar="VCFFileN", nargs='+',
                    help="path to (gzipped) vcf files with SNP information")
parser.add_argument("output",
                    help="name of (gzipped) outputfile in counts format")
parser.add_argument("-s", "--synonymous", action="store_true",
                    help="only print position with a synonymous base")
parser.add_argument("-v", "--verbose", action="count",
                    help="turn on verbosity (-v or -vv)")
parser.add_argument("-i", "--one-indiv", action="store_true",
                    help="randomly choose one indivual per population")
args = parser.parse_args()

gp_fn = args.gpfile
rf_fn = args.reference
vcfFnL = args.VCFFiles
output = args.output
vb = args.verbose
oneI = args.one_indiv

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()
if args.verbose == 0:
    logger.setLevel(logging.WARN)
elif args.verbose == 1:
    logger.setLevel(logging.INFO)
elif args.verbose == 2:
    logger.setLevel(logging.DEBUG)

cfw = cf.CFWriter(vcfFnL, output, oneIndividual=oneI)

if args.synonymous is True:
    cfw.onlySynonymous = True

gp_stream = gp.GPStream(gp_fn, rf_fn)

cfw.write_HLn()
cf.write_cf_from_gp_stream(gp_stream, cfw)

cfw.close()
