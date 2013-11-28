#!/usr/bin/env python

"""libPomo.seqbase
----------------------------------------------------------------------

This module provides basic functions and classes needed to work with
sequence data.

"""


class SequenceDataError(Exception):
    pass


class Seq:
    """A class that stores sequence data."""
    def __init__(self):
        self.names = []
        self.data = []
        self.dataLen = []
        self.nSpecies = 0

    def print_info(self, maxB=50):
        """Print sequence information.

        Print species names, the length of the sequence and a maximum
        of `maxB` bases (defaults to 50).

        """
        for i in range(0, self.nSpecies):
            print('>', self.names[i])
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
