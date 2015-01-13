#!/usr/bin/env python

"""Test libPoMo.cf.CFWriter object."""

import os
import libPoMo.fasta as fa
import libPoMo.cf as cf  # noqa
import libPoMo.seqbase as sb

vcfFL = ["/home/dominik/Population-Genetics/Data/Great-Ape-Genome-Project/SNPs/Gorilla.vcf.gz",          # noqa
         "/home/dominik/Population-Genetics/Data/Great-Ape-Genome-Project/SNPs/Homo.vcf.gz",             # noqa
         "/home/dominik/Population-Genetics/Data/Great-Ape-Genome-Project/SNPs/Pan_paniscus.vcf.gz",     # noqa
         "/home/dominik/Population-Genetics/Data/Great-Ape-Genome-Project/SNPs/Pan_troglodytes.vcf.gz",  # noqa
         "/home/dominik/Population-Genetics/Data/Great-Ape-Genome-Project/SNPs/Pongo_abelii.vcf.gz",     # noqa
         "/home/dominik/Population-Genetics/Data/Great-Ape-Genome-Project/SNPs/Pongo_pygmaeus.vcf.gz"]   # noqa

outfile = ".test-out.gz"

cfw = cf.CFWriter(vcfFL, outfile)
mFaStr = fa.MFaStream("/home/dominik/Population-Genetics/Data/hg18-CCDS-Alignments-Exons-Human-Chimp-Gorilla-Orang/hg18-chrY-exons.fa.gz")  # noqa

cfw.write_HLn()
cf.write_cf_from_MFaStream(mFaStr, cfw)

cfw.close()

fo = sb.gz_open(outfile)
i = 0
for line in fo:
    print(line, end='')
    i += 1
    if i > 100:
        break
fo.close()
os.remove(outfile)
