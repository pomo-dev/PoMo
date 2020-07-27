#!/usr/bin/env python

"""Convert counts to fasta format.

Extracts the sequences of a counts file.  If more than one base is
present at a single site, one base is sampled out of all present ones
according to its abundance.

The consensus sequence can be extracted (e.g., no sampling but the
bases with highest counts for each individual or population are
chosen) as well (see command line options).

"""

import argparse
import logging
import libPoMo

descr = """Convert counts to fasta format.

Extracts the sequences of a counts file.  If more than one base is
present at a single site, one base is sampled out of all present ones
according to its abundance.

The consensus sequence can be extracted (e.g., no sampling but the
bases with highest counts for each individual or population are
chosen) as well (see command line options).

"""

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=descr)

parser.add_argument("countsInFile",
                    help="path to (gzipped) counts file")
parser.add_argument("fastaOutFile",
                    help="name of (gzipped) fasta outputfile")
parser.add_argument("-v", "--verbose", action="count",
                    help="turn on verbosity (-v or -vv)")
parser.add_argument("-c", "--consensus", action="store_true",
                    help="extract consensus sequence")
args = parser.parse_args()

cfFN = args.countsInFile
output = args.fastaOutFile
vb = args.verbose
consFl = args.consensus

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()
if args.verbose == 0:
    logger.setLevel(logging.WARN)
elif args.verbose == 1:
    logger.setLevel(logging.INFO)
elif args.verbose == 2:
    logger.setLevel(logging.DEBUG)

print("Initialize Counts File Stream.")
cfStream = libPoMo.cf.CFStream(cfFN)

print("Convert to fasta.")
libPoMo.cf.cf_to_fasta(cfStream, output, consensus=consFl)

print("Done!")
