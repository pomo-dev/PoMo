#!/usr/bin/env python

"""Convert a fasta file with VCF files to counts format.

This script converts bgzipped VCF files containing SNP data to a file
in count format.  The VCF files need to be aligned to a reference
genome that needs to be provided in fasta format as input.  The SNP
data will be merged with the reference genome to create a complete
sequence that contains SNPs as well as unchanged bases (which are also
informative for phylogenetic analysis).

"""

import argparse
import libPoMo.fasta as fa  # noqa
import libPoMo.vcf as vcf  # noqa
import libPoMo.cf as cf  # noqa

descr = """Convert a fasta file with VCF files to counts format.

This script converts bgzipped VCF files containing SNP data to a file
in count format.  The VCF files need to be aligned to a reference
genome that needs to be provided in fasta format as input.  The SNP
data will be merged with the reference genome to create a complete
sequence that contains SNPs as well as unchanged bases (which are also
informative for phylogenetic analysis).

One of the chromosome names in the VCF file has to match the sequence
name given in the reference.

To speed up data conversion, tabix index files need to be provided for
all VCF files. They can be created from the terminal with $(tabix -p
vcf "vcf-file.vcf.gz") if tabix is installed.  The files to be indexed
have to be zipped with bgzip.

If the fasta reference sequence does not start at position 1 but at
position n+1, use the optional command line argument `--offset n`.

In VCF files, usually the bases of all copies of the same chromosomes
are given and separated by '/' or '|'.  If the species is not diploid,
this ploidy has to be set manually with the optional command line
argument `--ploidy n`.

Please check these important requirements before you use any data
created by this script.

The script can read and save standard text files or gzipped files.
This has to be indicated by .gz file endings.

"""

parser = argparse.ArgumentParser(
    prog="FastaToCounts.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=descr)

parser.add_argument("reference",
                    help="path to (gzipped) reference genome in fasta format")
parser.add_argument("VCFFiles", metavar="VCFFileN", nargs='+',
                    help="path to bgzipped vcf files with SNP information")
parser.add_argument("output",
                    help="name of (gzipped) outputfile in counts format")
parser.add_argument("-m", "--merge", action="count",
                    help="merge individuals within all given VCFFiles")
parser.add_argument("-o", "--offset", nargs=1,
                    help="offset of the sequence")
parser.add_argument("-p", "--ploidy", nargs=1,
                    help="ploidy of the sample")
parser.add_argument("-v", "--verbosity", action="count",
                    help="turn on verbosity")
args = parser.parse_args()

fastaRef = args.reference
vcfFnL = args.VCFFiles
output = args.output
offset = args.offset
vb = args.verbosity
if args.ploidy is not None:
    ploidy = args.ploidy[0]
else:
    ploidy = None

if args.merge is None:
    cfw = cf.CFWriter(vcfFnL, output, verb=vb)
else:
    mergeList = []
    nameList = []
    for fn in vcfFnL:
        mergeList.append(True)
        nameList.append(fn.split('.', maxsplit=1)[0])
    cfw = cf.CFWriter(vcfFnL, output, mergeL=mergeList,
                      nameL=nameList, verb=vb)

if ploidy is not None:
    cfw.set_ploidy(int(ploidy))

faR = fa.init_seq(fastaRef)
if offset is None:
    rg = faR.seq.get_region_no_description()
else:
    rg = faR.seq.get_region_no_description(offset)

cfw.set_seq(faR.seq)
cfw.write_HLn()
cfw.write_Rn(rg)

cfw.close()
faR.close()
