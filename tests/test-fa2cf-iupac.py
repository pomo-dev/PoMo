#!/usr/bin/env python

"""Test libPoMo.cf.fasta_to_cf with IUPAC codes."""

import os
import libPoMo.cf as cf

######################################################################
print("##################################################")
test_sequence = "data/fasta-iupac.dat"
fn = "vcf-test-tmp.dat"
print("Convert", test_sequence, "to Counts File.")
cf.fasta_to_cf(test_sequence, fn)
with open(fn) as file:
    for line in file:
        print(line, end='')
os.remove(fn)
