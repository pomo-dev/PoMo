#!/usr/bin/env python

"""Test libPoMo.countsformat module."""

import os
import import_libPoMo  # noqa
import libPoMo.countsformat as cf  # noqa
import libPoMo.vcf as vcf
import libPoMo.fasta as fa


print("Testing libPoMo/countsformat module.")
######################################################################
vcf_sequence = "data/vcf-wolfs.dat"
ref_sequence = "data/fasta-reference-wolf.dat"
fn = "cf-test-tmp.dat"
print("Try to save VCF file as counts format.")
vcfSeq = vcf.open_seq(vcf_sequence, name="wolfs")
faRef = fa.open_seq(ref_sequence, name="wolfs")
cf.save_as_countsformat(vcfSeq, faRef, fn)
print("Output:\n")
with open(fn) as file:
    for line in file:
        print(line, end='')
os.remove(fn)
