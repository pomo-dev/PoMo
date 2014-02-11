#!/usr/bin/env python

"""Filter a multiple sequence alignment file.

This script filters a multiple fasta sequence alignment file.
Cf. UCSC table browser:
http://genome.ucsc.edu/goldenPath/help/hgTablesHelp.html

The filter options cannot be set with this script, but a file can be
easily filtered using personal filter options with
libPoMo.fasta.filter_mfa_str().  See this function for further
details.

"""

import argparse
import libPoMo.seqbase as sb
import libPoMo.fasta as fa
import libPoMo.cf as cf  # noqa

descr = """Filter a multiple sequence alignment file.

This script filters a multiple fasta sequence alignment file.
Cf. UCSC table browser:
http://genome.ucsc.edu/goldenPath/help/hgTablesHelp.html

The filter options cannot be set with this script, but a file can be
easily filtered using personal filter options with
libPoMo.fasta.filter_mfa_str().  See this function for further
details.

"""

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=descr)

parser.add_argument("msaFile",
                    help="path to (gzipped) multiple sequence alignment file")
parser.add_argument("nSpecies",
                    help="number of aligned species in the given alignment")
parser.add_argument("output",
                    help="name of (gzipped) msa output file")
parser.add_argument('-v', "--verbosity", action="count",
                    help="turn on verbosity")
args = parser.parse_args()

mfaFN = args.msaFile
nSpecies = int(args.nSpecies)
output = args.output
vb = args.verbosity

mfa = fa.MFaStream(mfaFN)
fp = fa.MFaStrFilterProps(nSpecies)

oF = sb.gz_open(output, mode='w')

while True:
    if fa.filter_mfa_str(mfa, fp, vb) is True:
        mfa.print_msa(fo=oF)
    if mfa.read_next_align() is None:
        break

oF.close()
mfa.close()
