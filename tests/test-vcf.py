#!/usr/bin/env python

"""Test libPoMo.fasta module."""

import import_libPoMo  # noqa
import libPoMo.vcf as vcf  # noqa

test_sequence = 'data/vcf-homo.dat'
print("Testing libPoMo/vcf module.")
######################################################################
print("Read in test sequence ", test_sequence, '.', sep='')
seq = vcf.open_seq(test_sequence)

print("Print info of first base.")
base = seq.get_nuc_base(seq.bases[0].chrom, seq.bases[0].pos)
base.print_info()

print("\n")
print("Print sequence info.")
seq.print_info()

print("Print header line.")
indiv = ['testA', 'testB', 'testC']
vcf.print_header_line(indiv)
