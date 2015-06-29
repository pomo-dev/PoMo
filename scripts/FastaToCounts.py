#!/usr/bin/env python

"""Convert fasta to counts format.

The (aligned) sequences in the fasta file are read in and the data is
written to a counts format file.

"""

import argparse
import os
import logging
import libPoMo.fasta as fa
import libPoMo.cf as cf  # noqa

descr = """Convert fasta to counts format.

The (aligned) sequences in the fasta file are read in and the data is
written to a counts format file.

Sequence names are stripped at the first dash.  If the stripped
sequence name coincide, individuals are put into the same population.

E.g., homo_sapiens-XXX and homo_sapiens-YYY will be in the same
population homo_sapiens.

Take care with large files, this uses a lot of memory.

The input as well as the output files can additionally be gzipped
(indicated by a .gz file ending).

"""

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=descr)

parser.add_argument("fastaFile",
                    help="path to (gzipped) fasta file")
parser.add_argument("output",
                    help="name of (gzipped) outputfile in counts format")
parser.add_argument("-v", "--verbose", action="count",
                    help="turn on verbosity (-v or -vv)")
# TODO
# parser.add_argument("-i", "--one-indiv", action="store_true",
#                     help="randomly choose one indivual per population")
args = parser.parse_args()

FaRefFN = args.fastaFile
output = args.output
vb = args.verbose
# oneI = args.one_indiv

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()
if args.verbose == 0:
    logger.setLevel(logging.WARN)
elif args.verbose == 1:
    logger.setLevel(logging.INFO)
elif args.verbose == 2:
    logger.setLevel(logging.DEBUG)

cf.fasta_to_cf(FaRefFN, output)
