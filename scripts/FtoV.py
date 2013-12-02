#!/usr/bin/env python

"""Convert fasta files to VCF files.

This script converts a fasta file and a given fasta reference file to
a file in VCF. A reference genome can be provided in fasta format as
input. If no reference is given, the first sequence in the fasta file
will be used as reference.

"""

import argparse
import import_libPoMo  # noqa
import libPoMo.fasta as fa  # noqa
import libPoMo.vcf as vcf  # noqa

ver = "1.0"

parser = argparse.ArgumentParser(prog="FtoV",
                                 description="convert fast to VCF;" +
                                 " script version" + ver)
parser.add_argument("fastafile",
                    help="path to fasta file")
parser.add_argument("output",
                    help="name of VCF output file")
parser.add_argument("-r", "--reference",
                    help="path to reference genome in fasta format")
args = parser.parse_args()

faSeq = fa.open_seq(args.fastafile)
if args.reference is not None:
    faRef = fa.open_seq(args.reference)
    refSeq = faRef.get_seq_by_id(0)
else:
    refSeq = faSeq.get_seq_by_id(0)
fa.save_as_vcf(faSeq, refSeq, args.output)
