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

    self.name = name of the sequence (e.g. species or individual name)
    self.description = description of the sequence
    self.data = sequence data
    self.dataLen = number of saved bases

    """
    def __init__(self):
        self.name = ''
        self.descr = ''
        self.data = ''
        self.dataLen = 0

    def print_seq_header(self):
        print('>', self.name, ' ', self.descr, sep='')
        return

    def get_base(self, pos):
        """Returns base at position `pos`."""
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
    """Convenience function to strip filename off the `.xyz` ending."""
    filename_without_path = os.path.split(fn)[-1]
    return filename_without_path.rsplit('.', maxsplit=1)[0]
