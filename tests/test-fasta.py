#!/usr/bin/env python

"""Test libPoMo.fasta module."""

import os
import import_libPoMo  # noqa
import libPoMo.fasta as fa  # noqa


test_sequence = 'data/fasta-sample2.dat'
print("Testing libPoMo/fasta module.")
######################################################################
# print("\n##################################################")
# print("Read in test sequence ", test_sequence, '.', sep='')
# seq = fa.open_seq(test_sequence)

# print("Print sequence information.")
# seq.print_info()

######################################################################
print("\n##################################################")
print("Test FaStream object.")
faStr = fa.init_seq(test_sequence)
faStr.print_info(maxB=None)
while faStr.read_next_seq() is not None:
    faStr.print_info()
faStr.close_fo()


######################################################################
# print("\n##################################################")
# test_sequence = "data/fasta-sample-wolfs.dat"
# ref_sequence = "data/fasta-reference-wolf.dat"
# fn = "vcf-test-tmp.dat"
# print("Compare ", test_sequence, " to ", ref_sequence, '.', sep='')
# faSeq = fa.open_seq(test_sequence)
# faRef = fa.open_seq(ref_sequence)
# refSeq = faRef.get_seq_by_id(0)
# fa.save_as_vcf(faSeq, refSeq, fn)
# with open(fn) as file:
#     for line in file:
#         print(line, end='')
# os.remove(fn)
