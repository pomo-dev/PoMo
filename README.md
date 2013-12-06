PoMo 1.0.2
====
Implementation of a polymorphism aware phylogenetic model using HYPHY.

Created by: Nicola De Maio  
Contributors: Carolin Kosiol and Dominik Schrempf

For a reference, please see and cite: De Maio, Schlotterer, Kosiol
(MBE, 2013), and/or: De Maio, Kosiol (in preparation).  You can use
this software for non-commercial purposes but please always
acknowledge the authors.

Feel free to post any suggestions, doubts and bugs.

Installation
====
PoMo heavily relies on HYPHY. To build HYPHY, extract the
`HYPHY.tar.gz` to a folder and do
```sh
cd path_to_HYPHY/
cmake ./
make MP2
```

Running PoMo
====
To run PoMo, execute
```sh
python path_to_PoMo/PoMo.py  path_to_HYPHY/HYPHYMP path_to_data/data.txt
```

The data file must either be in fasta format, or in allele count format (see
example files in `./sample-data/`).

If your dataset is called "data.txt" the final output of PoMo will be
called "data_PoMo_output.txt", and contains the log-likelihood, the
parameters estimated, and the species tree estimated. The output will
be placed in the folder that you run PoMo from.

Command Line Arguments
====
PoMo supports various command line arguments. Run `python PoMo.py
--help` or visit [PoMo-help](PoMo-help.txt "PoMo-help") for a detailed
help message about the different flags that are available.

File Conversion
====
Various scripts for file conversion are provided in the folder
`./scripts/`. TODO

Notes
====
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
