PoMo 1.0.2
====
Implementation of a polymorphism aware phylogenetic model using HYPHY.

Created by: Nicola De Maio  
Contributor: Dominik Schrempf

For a reference, please see and cite: De Maio, Schlotterer, Kosiol
(MBE, 2013), and/or: De Maio, Kosiol (in preparation).  You can use
this software for non-commercial purposes but please always
acknowledge the authors.

Feel free to post any suggestions, doubts and bugs.

Installation
====
PoMo heavily relies on HYPHY. To build HYPHY, extract the
`HYPHY.tar.gz` to a folder and do
`cd path_to_HYPHY/`  
`cmake ./`  
`make MP2`

Running PoMo
====
To run PoMo, execute

`python path_to_PoMo/PoMo.py  path_to_HYPHY/HYPHYMP path_to_data/data.txt`

The data file must either be in fasta format, or in allele count format (see
example files in `./data/`).

Options
====
There are some further optional parameters that can
be specified from the command line. Type `python path_to_PoMo/Pomo.py
--help` for detailed information.

`--help`
: prints detailed help.

`-m 0` turns off the molecular clock.

`-u GTR` allows to choose a mutation model different from the HKY
(default option).  "GTR" corresponds to the general time reversible,
"F81" to the Felsenstein 1981 (reversible, equal mutation rates), and
"NONREV" to the general nonreversible model (all substitution rates
are independent).

`-d 0.90` keeps at least 90% of the sites when downsampling (or any
other proportion if you specify a different one between 0 and
1). Downsampling is done when sites do not have the same coverage
along the genome. In such a case, all sites with coverage higher than
a new fixed sample size are downsampled, those with lower coverage are
discarded. The new sample size is chosen as the highest possible among
those that leave the number of sites above the specified
threshold. The default threshold is 0.66.

`-v` turns on verbosity.

If your dataset is called "data.txt" the final output of PoMo will be
called "data_PoMo_output.txt", and contains the log-likelihood, the
parameters estimated, and the species tree estimated. The output will
be placed in the folder that you run PoMo from.

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
