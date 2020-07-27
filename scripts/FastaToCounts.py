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

If heterozygotes are encoded with IUPAC codes (e.g., 'r' for A or G),
homozygotes need to be counted twice so that the level of polymorphism
stays correct.  This can be done with the `--iupac` flag.

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
parser.add_argument("--iupac", action="store_true",
                    help="heteorzygotes are encoded with IUPAC codes")
# TODO
# parser.add_argument("-i", "--one-indiv", action="store_true",
#                     help="randomly choose one indivual per population")
args = parser.parse_args()

FaRefFN = args.fastaFile
output = args.output
vb = args.verbose
iupac_flag = args.iupac
# oneI = args.one_indiv

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()
if args.verbose == 0:
    logger.setLevel(logging.WARN)
elif args.verbose == 1:
    logger.setLevel(logging.INFO)
elif args.verbose == 2:
    logger.setLevel(logging.DEBUG)

cf.fasta_to_cf(FaRefFN, output, double_fixed_sites=iupac_flag)
