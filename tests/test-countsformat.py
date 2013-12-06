#!/usr/bin/env python

"""Test libPoMo.countsformat module."""

import os
import import_libPoMo  # noqa
import libPoMo.countsformat as cf  # noqa
import libPoMo.vcf as vcf
import libPoMo.fasta as fa


saveAsCfSingleChrom = False
saveAsCfMultChrom = False
vcfWithError = True
gzipped = False
largeFile = False

print("Testing libPoMo/countsformat module.")
######################################################################
if saveAsCfSingleChrom is True:
    print("\n##################################################")
    vcf_sequence = "data/vcf-wolfs.dat"
    ref_sequence = "data/fasta-reference-wolf.dat"
    fn = "cf-test-tmp.dat"
    print("Try to save VCF file as counts format with a single chromosome.")
    vcfStr = vcf.init_seq(vcf_sequence, name="wolfs")
    vcfStr.print_info()
    refFaStr = fa.init_seq(ref_sequence, name="wolfs")
    cf.save_as_cf(vcfStr, refFaStr, fn)
    print("\nOutput:")
    with open(fn) as file:
        for line in file:
            print(line, end='')
    vcfStr.close_fo()
    os.remove(fn)
######################################################################
if saveAsCfMultChrom is True:
    print("\n##################################################")
    vcf_sequence = "data/vcf-chroms.dat"
    ref_sequence = "data/fasta-chroms-ref.dat"
    fnStr = "cf-test-tmp.dat"
    print("Try to save VCF file as counts format with a single chromosome.")
    vcfStr = vcf.init_seq(vcf_sequence, name="homo")
    refFaStr = fa.init_seq(ref_sequence, name="homo")
    cf.save_as_cf(vcfStr, refFaStr, fnStr, add=True, name='homo')
    print("\nOutput:")
    with open(fnStr) as file:
        for line in file:
            print(line, end='')
    vcfStr.close_fo()
    refFaStr.close_fo()
    os.remove(fnStr)
######################################################################
if vcfWithError is True:
    print("\n##################################################")
    vcf_sequence = "data/vcf-chroms-with-errors.dat"
    ref_sequence = "data/fasta-chroms-ref.dat"
    fnStr = "cf-test-tmp.dat"
    print("Try to save VCF file as counts format with a single chromosome.")
    vcfStr = vcf.init_seq(vcf_sequence, name="homo")
    refFaStr = fa.init_seq(ref_sequence, name="homo")
    cf.save_as_cf(vcfStr, refFaStr, fnStr, add=True, name='homo')
    print("\nOutput:")
    with open(fnStr) as file:
        for line in file:
            print(line, end='')
    vcfStr.close_fo()
    refFaStr.close_fo()
    os.remove(fnStr)
######################################################################
if gzipped is True:
    print("\n##################################################")
    vcf_sequence = "data/vcf-chroms.dat.gz"
    ref_sequence = "data/fasta-chroms-ref.dat"
    fnStr = "cf-test-tmp.dat"
    print("Save gzipped VCF file as counts format with a single chromosome.")
    vcfStr = vcf.init_seq(vcf_sequence, name="homo")
    refFaStr = fa.init_seq(ref_sequence, name="homo")
    cf.save_as_cf(vcfStr, refFaStr, fnStr, add=True, name='homo')
    print("\nOutput:")
    with open(fnStr) as file:
        for line in file:
            print(line, end='')
    vcfStr.close_fo()
    refFaStr.close_fo()
    os.remove(fnStr)
######################################################################
if largeFile is True:
    print("\n##################################################")
    vcf_sequence = "data/vcf-homo.dat"
    ref_sequence = "data/hg18-chr1.fa.gz"
    fnStr = "cf-test-tmp.dat"
    print("Test save_as_cf on real data.")
    vcfStr = vcf.init_seq(vcf_sequence, name="homoChr1")
    refFaStr = fa.init_seq(ref_sequence, name="homoChr1")
    cf.save_as_cf(vcfStr, refFaStr, fnStr, add=True, name='homoChr1', diploid=True)
    vcfStr.close_fo()
    refFaStr.close_fo()
