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
import gzip


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
    :param str name: Optional, region name.

    :ivar str chrom: Chromosome name.
    :ivar int start: 0-based start position.
    :ivar int end: 0-base end position.
    :ivar str name: Region name.
    """
    def __init__(self, chrom, start, end, name=None):
        self.chrom = chrom
        self.start = start - 1
        self.end = end - 1
        self.name = name

    def print_info(self):
        """Print information about the region."""
        if self.name is not None:
            print("Region name:", self.name)
        print("Chromosome name:", self.chrom)
        print("1-based start position:", self.start)
        print("1-based end position:", self.end)


class Seq:
    """A class that stores sequence data.
    .. _seqbase-seq:

    :ivar str name: Name of the sequence (e.g. species or individual
                    name).
    :ivar str descr: Description of the sequence.
    :ivar str data: String with sequence data.
    :ivar int dataLen: Number of saved bases.
    :ivar Boolean rc: True if *self.data* stores the
                      reverse-complement of the real sequence.

    """
    def __init__(self):
        self.name = ''
        self.descr = ''
        self.data = ''
        self.dataLen = 0
        self.rc = False

        self.__lowered = False

    def print_fa_header(self):
        """Print the sequence header line in fasta format."""
        print('>', self.name, ' ', self.descr, sep='')
        return

    def get_base(self, pos):
        """Returns base at 1-based position `pos`."""
        if pos > self.dataLen:
            raise SequenceDataError("Position out of range.")
        return self.data[pos-1]

    def toggle_rc(self):
        """Toggle the state of *self.rc*."""
        self.rc = not self.rc

    def set_rc(self):
        """Set the *self.rc*.

        The instance variable *self.rc* is a Boolean value that is
        true if the saved sequence is reversed and complemented.  This
        function sets this value according to the last character in
        the sequence description.

        :raises: *ValueError()* if state could not be detected.

        """

        if self.descr[-1] == '-':
            self.rc = True
        elif self.descr[-1] == '+':
            self.rc = False
        else:
            raise ValueError("State could not be detected.")

    def get_rc(self):
        """Return True if the sequence is reversed and complemented.

        :rtype: Boolean

        """
        return self.rc

    def rev_comp(self):
        """Reverses and complements the sequence.

        This is rather slow for long sequences.

        """
        compDict = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
        self.data = self.data[::-1]
        if self.__lowered is not True:
            self.data = self.data.lower()
        rcData = []
        for i in range(self.dataLen):
            rcData.append(compDict[self.data[i]])
        self.data = ''.join(rcData)
        self.toggle_rc()

    def get_region(self):
        """Try to find the :class:`Region` that the sequence spans.

        The sequence might not physically start at position 1 but at
        some arbitrary value that is indicated in the sequence
        description.  This function gets this physical
        :class:`Region`, if the description of the sequence is of the
        form (cf. `UCSC Table Browser
        <http://genome.ucsc.edu/goldenPath/help/hgTablesHelp.html#FASTA>`_)::

          918 0 0 chr1:58954-59871+

        :raises: :class:`SequenceDataError`, if format of description
          is invalid.

        """
        if self.descr[-1] not in ['+', '-']:
            raise SequenceDataError("Direction character is missing.")
        descrL = self.descr[:-1].rsplit(maxsplit=1)
        if len(descrL) >= 1:
            rgStr = descrL[1]
            rgStrL = rgStr.split(':', maxsplit=1)
            if len(rgStrL) >= 1:
                chromName = rgStrL[0]
                posStr = rgStrL[1]
                posStrL = posStr.split('-', maxsplit=1)
                if len(posStrL) >= 1:
                    start = int(posStrL[0])
                    end = int(posStrL[1])
            else:
                raise SequenceDataError("Regional information is invalid.")
        else:
            raise SequenceDataError("Description format is invalid.")
        rg = Region(chromName, start, end, name=self.name)
        return rg

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


def gz_open(fn, mode='r'):
    """Open file with io.open() or gzip.open().

    :param str fn: Name of the file to open.
    :param char md: Mode '**r**' | 'w'.

    """
    if fn[-2:] == "gz":
        fo = gzip.open(fn, mode=mode+'t')
    else:
        fo = open(fn, mode=mode)
    return fo
