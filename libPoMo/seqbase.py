#!/usr/bin/env python

"""libPomo.seqbase
==================

This module provides basic functions and classes needed to work with
sequence data.

Objects
-------
Classes:
  - :class:`Seq`, stores a single sequence
  - :class:`Region`, region in a genome

Exception Classes:
  - :class:`SequenceDataError`
  - :class:`NotAValidRefBase`

Functions:
  - :func:`stripFName()`, strip filename off its ending

----

"""

__docformat__ = 'restructuredtext'

import os


class SequenceDataError(Exception):
    """General sequence data error exception."""
    pass


class NotAValidRefBase(SequenceDataError):
    """Reference base is not valid."""
    pass


class Region():
    """Region in a genome.

    The start and end points need to be given 1-based and are
    converted to 0-based positions that are used internally to save
    all positional data.
    
    :param str chrom: Chromosome name.
    :param int start: 1-based start position.
    :param int end: 1-based end position.

    :ivar str chrom: Chromosome name.
    :ivar int start: 0-based start position.
    :ivar int end: 0-base end position.

    """
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start - 1
        self.end = end - 1


class Seq:
    """A class that stores sequence data.
    .. _seqbase-seq:

    :ivar str name: Name of the sequence (e.g. species or individual
                    name).
    :ivar str descr: Description of the sequence.
    :ivar str data: String with sequence data.
    :ivar int dataLen: Number of saved bases.

    """
    def __init__(self):
        self.name = ''
        self.descr = ''
        self.data = ''
        self.dataLen = 0

    def print_seq_header(self):
        """Print the sequence header line in fasta format."""
        print('>', self.name, ' ', self.descr, sep='')
        return

    def get_base(self, pos):
        """Returns base at 1-based position `pos`."""
        if pos > self.dataLen:
            raise SequenceDataError("Position out of range.")
        return self.data[pos-1]

    def purge(self):
        """Purge data saved in this sequence."""
        self.name = ''
        self.descr = ''
        self.data = ''
        self.dataLen = 0


def stripFName(fn):
    """Convenience function to strip filename off the ".xyz" ending."""
    filename_without_path = os.path.split(fn)[-1]
    return filename_without_path.rsplit('.', maxsplit=1)[0]
