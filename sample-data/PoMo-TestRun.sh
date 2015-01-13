#!/bin/bash

# Do a PoMo test run.

Python=python
PoMo=~/Population-Genetics/PoMo-git-repo/PoMo-Development/PoMo.py
HYPHY=~/Population-Genetics/PoMo-git-repo/HYPHY/HYPHYMP
data=counts_file_sample.cf

$Python $PoMo $HYPHY $data -v
