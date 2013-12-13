#!/usr/bin/env python

"""Test libPoMo.cf.CFWriter object."""

import argparse
import import_libPoMo  # noqa
import libPoMo.seqbase as sb
import libPoMo.cf as cf  # noqa

parser = argparse.ArgumentParser(
    prog="test-CFWriter.py",
    description="Test `libPoMo.cf.CFWriter` object")
parser.add_argument("start", type=int,
                    help="start position on chr1")
parser.add_argument("end", type=int,
                    help="end position on chr1")
args = parser.parse_args()

start = args.start
end = args.end

vcfFL = ["/home/dominik//PopGen/Data/Great-Ape-Genome-Project/SNPs/Gorilla.vcf.gz",          # noqa
         "/home/dominik//PopGen/Data/Great-Ape-Genome-Project/SNPs/Homo.vcf.gz",             # noqa
         "/home/dominik//PopGen/Data/Great-Ape-Genome-Project/SNPs/Pan_paniscus.vcf.gz",     # noqa
         "/home/dominik//PopGen/Data/Great-Ape-Genome-Project/SNPs/Pan_troglodytes.vcf.gz",  # noqa
         "/home/dominik//PopGen/Data/Great-Ape-Genome-Project/SNPs/Pongo_abelii.vcf.gz",     # noqa
         "/home/dominik//PopGen/Data/Great-Ape-Genome-Project/SNPs/Pongo_pygmaeus.vcf.gz"]   # noqa

rg = sb.Region("chr1", start, end)

cfw = cf.CFWriter("data/hg18-chr1.fa.gz", vcfFL, "test-out.gz")

cfw.write_HLn()
cfw.write_Rn(rg)

cfw.close()
