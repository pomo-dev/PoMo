PoMo 1.1.0
==========

Implementation of a polymorphism aware phylogenetic model using HYPHY.

Created by:
- Nicola De Maio
Contributors:
- Dominik Schrempf
- Carolin Kosiol

For a reference, please see and cite: De Maio, Schlotterer, Kosiol
(MBE, 2013), and/or: De Maio, Schrempf and Kosiol (ms. in
preparation).  You can use this software for non-commercial purposes
but please always acknowledge the authors.

Feel free to post any suggestions, doubts and bugs.

libPoMo
=======

PoMo comes with a small python package `libPoMo`. It contains several
modules that ease the handling of data files in fasta, variant call
(vcf), or counts format (the file format PoMo works with).

For a detailed and comprehensive documentation visit
[libPoMo documentation](http://pomo.readthedocs.org).

However, it is highly probable, that you only need to consider the
[file conversion scripts](#file-conversion) to prepare counts files to
be analyzed with PoMo.

Installation Requirements
=========================

Before installing, check the following requirements:
- `git` (https://github.com/)
- `python3` (https://www.python.org/)
PoMo also uses different python libraries that need to be installed
separately:
- [scipy](http://www.scipy.org/),
- [numpy](http://www.numpy.org/) and
- [pysam](http://code.google.com/p/pysam/)

Installation
============

Download PoMo with:
```sh
git clone git://github.com/pomo-dev/PoMo
```

This will create a folder `PoMo`.

PoMo heavily relies on HyPhy.  It comes with its own, patched version
that can be found in the folder created above.  To build this version
of HyPhy, extract the `HYPHY.tar.gz` (this archive is part of the git
repo; PoMo will not work if you use another version of HyPhy) and do
```sh
cd path_to_HYPHY/ cmake ./
make MP2
```

Running PoMo
============
To run PoMo, execute with `python3`
```sh
python path_to_PoMo/PoMo.py  path_to_HYPHY/HYPHYMP path_to_data/data.txt
```

The data file must either be in fasta format, or in counts format (see
example files in `./sample-data/`).

If your dataset is called "data.txt" the final output of PoMo will be
called "data_PoMo_output.txt". It will contain the log-likelihood, the
estimated parameters, and the estimated species tree. The output will
be placed in the folder that you run PoMo from.

Command Line Arguments
======================
PoMo supports various command line arguments. Run `python PoMo.py
--help` or visit [PoMo-help](PoMo-help.txt "PoMo-help") for a detailed
help message about the different flags that are available.

File Conversion
===============
Various scripts for file conversion are provided in the folder
`./scripts/`. Run `python scriptXY.py --help` for detailed help
messages. These scripts use `libPoMo` that has to be in the python
path.

* [CountsToFasta.py](./scripts/CountsToFasta.py): Convert a counts
  file to a fasta file.
* [FastaToCounts.py](./scripts/FastaToCounts.py): Convert a fasta file
  to counts format.
* [FastaToVCF.py](./scripts/FastaToVCF.py): Convert a fasta file to
  variant call format.
* [FastaVCFToCounts.py](./scripts/FastaVCFToCounts.py): Convert a
  fasta reference with VCF files to counts format.
* [FilterMSA.py](./scripts/FilterMSA.py): Filter a multiple sequence
  alignment file (apply standard filters; cf. `libPoMo`).
* [MSAToCounts.py](./scripts/MSAToCounts.py): Convert multiple
  sequence alignments with VCF files to counts format.

Data
====
Sample data can be found in
[`pomo-dev/data`](https://github.com/pomo-dev/data).  This repository
also includes the data sets of our manuscript that has been accepted
to Systematic Biology.  A proper citation will follow, for now please
refer to [BioRxiv](http://dx.doi.org/10.1101/016360).

Notes
=====
PoMo works best if data about within-population variation is
provided. So, in case you do only have a single sample per species,
you need to provide PoMo with some meaningful estimate of theta. In
any other case, you do not have this problem.

Virtual population size is 10. This is the only case we have
extensively simulated. I will relax this constraint in the future, but
for now, population samples larger than 10 are down-sampled to 10.

The model assumes that all data from the loci considered is provided:
both variable sites, and sites that are fixed through species and
populations. If you have only a SNP dataset (no fixed sites), there is
probably a way to work around it, let me know.

The publication that is out in MBE about PoMo discusses estimation of
population parameters (mutation rates, fixation biases, etc) and not
of species trees. Soon we are going to submit our manuscript regarding
species tree estimation with PoMo. Here the focus is solely on
estimating species trees. It should be easy for me to extend the
software to perform both analyses simultaneously in the future.

The computational time should be feasible for up to few dozens of
species, no matter of the sample sizes and number of loci.

PoMo is heavily based on HyPhy. You should probably acknowledge
HyPhy's authors (e.g. Sergei L kosakovsky-Pond) if you use PoMo.

In the fasta format all individuals must be named as "species_n" where
"species" is the name of the species and "n" is a number, it does not
matter which number.
