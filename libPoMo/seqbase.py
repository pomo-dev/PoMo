#!/usr/bin/env python

"""libPoMo.seqbase
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
import sys


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
    def __init__(self, chrom, start, end, name=None, orientation="+"):
        self.chrom = chrom
        self.start = start - 1
        self.end = end - 1
        self.name = name
        self.orientation = orientation

    def print_info(self):
        """Print information about the region."""
        if self.name is not None:
            print("Region name:", self.name)
        print("Chromosome name:", self.chrom)
        print("0-based start position:", self.start)
        print("0-based end position:", self.end)


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
        self.name = None
        self.descr = None
        self.data = None
        self.dataLen = 0
        self.rc = False

        self.__lowered = False

    def print_fa_header(self, fo=sys.stdout):
        """Print the sequence header line in fasta format.

        :ivar fileObject fo: Print to file object fo. Defaults to
          stdout.

        """
        print('>', self.name, ' ', self.descr, sep='', file=fo)
        return

    def print_fa_entry(self, maxB=None, fo=sys.stdout):
        """Print a fasta file entry with header and sequence data.

        :ivar int maxB: Print a maximum of maxB bases. Default: print
          all bases.

        """
        self.print_fa_header(fo)
        if maxB is None:
            print(self.data, file=fo)
        else:
            print("First", maxB, "bases: ", end='', file=fo)
            print(self.data[:maxB], file=fo)
        return

    def print_data(self, fo=sys.stdout):
        """Print the sequence data.

        :ivar fileObject fo: Print to file object fo. Defaults to
          stdout.

        """
        print(self.data, file=fo)
        return

    def get_base(self, pos):
        """Returns base at 1-based position `pos`."""
        if pos > self.dataLen:
            raise SequenceDataError("Position out of range.")
        return self.data[pos-1]

    def print_info(self, maxB=50):
        """Print sequence information.

        Print sequence name, description, the length of the sequence
        and a maximum of `maxB` bases (defaults to 50).

        """
        print("Sequence name:", self.name)
        print("Sequence description:", self.descr)
        print("Sequence length:", self.dataLen)
        print("First", maxB, "bases:", end='')
        print(self.data[0:maxB])
        return

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

    def rev_comp(self, change_sequence_only=False):
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
        if self.descr[-1] == '+':
            tempDescr = self.descr[:-1] + '-'
        elif self.descr[-1] == '-':
            tempDescr = self.descr[:-1] + '+'
        if change_sequence_only is False:
            self.descr = tempDescr
            self.toggle_rc()

    def get_exon_nr(self):
        """Try to find the current and the total exon number of the sequence.

        Extract the exon number and the total number of exons, if the
        name of the sequence is of the form (cf. `UCSC Table Browser
        <http://genome.ucsc.edu/goldenPath/help/hgTablesHelp.html#FASTA>`_)::

          >CCDS3.1_hg18_2_19

        :rtype: (int nEx, int nExTot)

        :raises: :class:`SequenceDataError`, if the format of the
          sequence name is invalid.

        """
        nameL = self.name.rsplit(sep='_', maxsplit=2)
        if len(nameL) >= 3:
            try:
                nEx = int(nameL[-1])
            except ValueError:
                raise SequenceDataError("Exon number not valid.")
            try:
                nExTot = int(nameL[-2])
            except ValueError:
                raise SequenceDataError("Total exon number not valid.")
        else:
            raise SequenceDataError("Exon information not valid.")
        return (nEx, nExTot)

    def get_in_frame(self):
        """Try to find the `inFrame` of the gene.

        `inFrame`: the frame number of the first nucleotide in the
        exon. Frame numbers can be 0, 1, or 2 depending on what
        position that nucleotide takes in the codon which contains it.
        This function gets the `inFrame`, if the description of the
        sequence is of the form (cf. `UCSC Table Browser
        <http://genome.ucsc.edu/goldenPath/help/hgTablesHelp.html#FASTA>`_)::

          918 0 0 chr1:58954-59871+

        :rtype: int

        :raises: :class:`SequenceDataError`, if format of description
          is invalid.

        """
        descrL = self.descr[:-1].split(maxsplit=2)
        if len(descrL) >= 2:
            try:
                inFrame = int(descrL[1])
            except ValueError:
                raise SequenceDataError("Description format is invalid.")
        else:
            raise SequenceDataError("Description format is invalid.")
        return inFrame

    def is_synonymous(self, pos):
        """Return True if the base at `pos` is 4-fold degenerate.

        This function checks if the base at `pos` is a synonymous one.
        The description of the sequence has to be of the form
        (cf. `UCSC Table Browser
        <http://genome.ucsc.edu/goldenPath/help/hgTablesHelp.html#FASTA>`_)::

          918 0 0 chr1:58954-59871+

        :ivar int pos: Position of the base in the sequence (0 to
          self.dataLen).

        :rtype Boolean: True if base is 4-fold degenerate.

        :raises: :class:`SequenceDataError`, if format of description
          is invalid.

        """
        raise ValueError("Reverse complemented genes not handled correctly.")
        if self.rc is True:
            raise ValueError("Examination of reverse complemented sequence.")
        degTriplets = ["tc", "ct", "cc", "cg", "ac", "gt", "gc", "gg"]
        inFr = self.get_in_frame()
        if pos < 2:
            # Degeneracy can not be determined.
            return False
        elif (pos + 1 + inFr) % 3 != 0:
            # Position within a Frame is not the third one.
            return False
        else:
            triplet = self.data[pos-2:pos+1]
            triplet = triplet.lower()
            if triplet[0:2] in degTriplets:
                return True
        return False

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

    def get_region_no_description(self, offset=0):
        """Get the region of the sequence.

        If no regional information is available in the sequence
        description (cf. :func:`get_region`), the position of the
        first base in the reference genome can be given
        manually. E.g., if the first base of the sequence does not
        correspond to the first but to the 11th base of the reference
        sequence, the offset should be 10.

        The name of the chromosome will be set to the name of the
        sequence.

        :param int offset: Optional, offset of the sequence.

        """
        chromName = self.name
        start = offset + 1
        end = offset + self.dataLen
        return Region(chromName, start, end)

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


def gz_open(fn, mode='r', buffering=-1):
    """Open file with io.open() or gzip.open().

    :param str fn: Name of the file to open.
    :param char md: Mode '**r**' | 'w'.

    """
    if fn[-2:] == "gz":
        # Tue Jan 12 12:50:50 CET 2016 Ignore buffering because it
        # leads to an error.
        fo = gzip.open(fn, mode=mode+'t')
        # fo = gzip.open(fn, mode=mode+'t', buffering=buffering)
    else:
        fo = open(fn, mode=mode, buffering=buffering)
        # fo = open(fn, mode=mode, buffering=buffering)
    return fo
