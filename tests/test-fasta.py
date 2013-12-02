#!/usr/bin/env python

"""Test libPoMo.fasta module."""

import os
import import_libPoMo  # noqa
import libPoMo.fasta as fa  # noqa


test_sequence = 'data/fasta-sample1.dat'
print("Testing libPoMo/fasta module.")
######################################################################
# print("Read in test sequence ", test_sequence, '.', sep='')
# seq = fa.open_seq(test_sequence)

# print("Print sequence information.")
# seq.print_info()

test_sequence = "data/fasta-sample-wolfs.dat"
ref_sequence = "data/fasta-reference-wolf.dat"
fn = "vcf-test-tmp.dat"
print("Compare ", test_sequence, " to ", ref_sequence, '.', sep='')
faSeq = fa.open_seq(test_sequence)
faRef = fa.open_seq(ref_sequence)
refSeq = faRef.get_seq_by_id(0)
fa.save_as_vcf(faSeq, refSeq, fn)
with open(fn) as file:
    for line in file:
        print(line, end='')
os.remove(fn)
