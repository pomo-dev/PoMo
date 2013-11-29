#!/usr/bin/env python

"""Test libPoMo.fasta module."""

import import_libPoMo  # noqa
import libPoMo.fasta as fa  # noqa


test_sequence = 'data/fasta-sample1.dat'
print("Testing libPoMo/fasta module.")
######################################################################
print("Read in test sequence ", test_sequence, '.', sep='')
seq = fa.open_seq(test_sequence)

print("Print sequence information.")
seq.print_info()
