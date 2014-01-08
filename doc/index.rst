.. PoMo documentation master file, created by
   sphinx-quickstart on Fri Dec 13 15:06:09 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to PoMo's documentation!
================================
Implementation of a polymorphism aware phylogenetic model using HYPHY.

Created by:
  * Nicola De Maio

Contributors:
  * Carolin Kosiol 
  * Dominik Schrempf

For a reference, please see and cite: De Maio, Schlotterer, Kosiol
(MBE, 2013), and/or: De Maio, Kosiol (in preparation). You can use
this software for non-commercial purposes but please always
acknowledge the authors.

Feel free to post any suggestions, doubts and bugs.


libPoMo
=======

PoMo comes with a small python package :doc:`libPoMo`. It contains
several modules that ease the handling and preparation of data files
in variant call format (vcf), fasta format and counts format (cf).

This documentation mainly focusses on this library. For a general
introduction to PoMo, please visit `PoMo at github.com
<https://github.com/fazky/PoMo>`_.

The *libPoMo* package is split into the following modules:
  * :doc:`main <main>`: Contains functions that are used by PoMo.
  * :doc:`seqbase <seqbase>`: Provides basic functions and classes
    needed to work with sequence data.
  * :doc:`fasta <fasta>`: Provides functions to read, write and access
    fasta files.
  * :doc:`vcf <vcf>`: Provides functions to read, write and access vcf
    files.
  * :doc:`cf <cf>`: Provides functions to read, write and access files
    that are in counts format.

Contents:
=========
.. toctree::
   :maxdepth: 2

   libPoMo


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

