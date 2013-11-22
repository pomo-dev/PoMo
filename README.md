PoMo
====

Implementation of a polymorphism aware phylogenetic model using HYPHY.

PoMo version 1.0  
Created by: Nicola De Maio.  
Contributors: Dominik Schrempf

For a reference, please see and cite: De Maio, Schlotterer, Kosiol
(MBE, 2013), and/or: De Maio, Kosiol (in preparation).  You can use
this software for non-commercial purposes but please always
acknowledge the authors.  For suggestions, doubts, bugs, etc., please
contact nicola.de.maio.85@gmail.com.

INSTALLATION
====
Type from terminal:

`cd path_to_PoMo/PoMo/
cmake ./
make MP2`

RUNNING PoMo
====
Type from terminal:

`python path_to_PoMo/PoMo/PoMo.py  path_to_data/Data.txt`

To run PoMo. Data must either be in Fasta format, or in allele count
format (see example files).  There are some further optional
parameters that can be specified from command line:

`-m 0` Turns off the molecular clock.

`-v` Turns on verbosity.

`-d 0.90` At least 90% of the sites are kept when downsampling (or any
other proportion if you specify a different one between 0 and
1). Downsampling is done when sites do not have the same coverage
along the genome. In such a case, all sites with coverage higher than
a new fixed sample size are downsampled, those with lower coverage are
discarded. The new sample size is chosen as the highest possible among
those that leave the number of sites above the specified
threshold. The default threshold is 0.66.

If your dataset is called "Data.txt" the final output of PoMo will be
called "Data_PoMo_output.txt", and contains the log-likelihood, the
parameters estimated, and the species tree estimated. The output will
be placed in the same folder that contains your data.

NOTES
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
