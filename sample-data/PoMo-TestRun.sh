#!/bin/bash

# Do a PoMo test run.
# Please adjust the file paths.

Python=python
PoMo=../PoMo.py
HYPHY=../../HYPHY/HYPHYMP
data=counts_file_sample.cf

$Python $PoMo $HYPHY $data -v
