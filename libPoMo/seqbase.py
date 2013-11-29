#!/usr/bin/env python

"""libPomo.seqbase
----------------------------------------------------------------------

This module provides basic functions and classes needed to work with
sequence data.

"""

import os


class SequenceDataError(Exception):
    pass


class Seq:
    """A class that stores sequence data.

    self.name = name of the total sequence (e.g. species name)
    self.names = names of the sequences (e.g. name of the individual
                 or of the chromosome)
    self.data = sequence data
    self.dataLen = length of sequence data
    self.nSpecies = number of species (individuals, chromosomes) saved
                    in the object

    """
    def __init__(self):
        self.name = ""
        self.names = []
        self.data = []
        self.dataLen = []
        self.nSpecies = 0

    def print_seq_header(self, i):
        print('>', self.names[i], sep='')
        return

    def print_info(self, maxB=50):
        """Print sequence information.

        Print species names, the length of the sequence and a maximum
        of `maxB` bases (defaults to 50).

        """
        print("Sequence name:", self.name)
        for i in range(0, self.nSpecies):
            self.print_seq_header(i)
            print("Printing", maxB, "out of a total of",
                  self.dataLen[i], "bases.")
            print(self.data[i][0:maxB])
        return

    def get_base(self, seq, pos):
        """Returns base at position `pos` in sequence with name `seq`."""
        try:
            i = self.names.index(seq)
        except:
            raise SequenceDataError("Sequence name not found.")
        if pos > self.dataLen[i]:
            raise SequenceDataError("Position out of range.")
        return self.data[i][pos-1]


def stripFName(fn):
    """Convenience function to strip filename off the `.xyz` ending."""
    filename_without_path = os.path.split(fn)[-1]
    return filename_without_path.rsplit('.', maxsplit=1)[0]
