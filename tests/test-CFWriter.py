#!/usr/bin/env python

"""Test libPoMo.cf.CFWriter object."""

import libPoMo.fasta as fa
import libPoMo.cf as cf  # noqa

vcfFL = ["/home/dominik/PopGen/Data/Great-Ape-Genome-Project/SNPs/Gorilla.vcf.gz",          # noqa
         "/home/dominik/PopGen/Data/Great-Ape-Genome-Project/SNPs/Homo.vcf.gz",             # noqa
         "/home/dominik/PopGen/Data/Great-Ape-Genome-Project/SNPs/Pan_paniscus.vcf.gz",     # noqa
         "/home/dominik/PopGen/Data/Great-Ape-Genome-Project/SNPs/Pan_troglodytes.vcf.gz",  # noqa
         "/home/dominik/PopGen/Data/Great-Ape-Genome-Project/SNPs/Pongo_abelii.vcf.gz",     # noqa
         "/home/dominik/PopGen/Data/Great-Ape-Genome-Project/SNPs/Pongo_pygmaeus.vcf.gz"]   # noqa

cfw = cf.CFWriter(vcfFL, "test-out.gz")
mFaStr = fa.MFaStream("/home/dominik/PopGen/Data/CCDC-Alignments-Exons-Human-Chimp-Gorilla-Orang/chrY-exons.fa.gz")  # noqa

cfw.write_HLn()
cf.write_cf_from_MFaStream(mFaStr, cfw)

cfw.close()
