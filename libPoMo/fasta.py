#!/usr/bin/env python

"""libPomo.fasta
----------------------------------------------------------------------

This module provides functions to read, write and access fasta files.

"""


class NotAFastaFileError(Exception):
    """Exception raised if given fasta file is not valid."""
    pass


class SequenceDataError(Exception):
    # TODO
    pass


class FaSeq:
    """A class that stores data retrieved from a fasta file."""
    def __init__(self):
        self.names = []
        self.descr = []
        self.data = []
        self.dataLen = []
        self.nSpecies = 0

    def print_info(self, maxB=50):
        """Print fasta sequence information.

        Print species names, description, the length of the sequence
        and a maximum of `maxB` bases (defaults to 50).

        """
        for i in range(0, self.nSpecies):
            print('>', self.names[i], self.descr[i])
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


def get_sp_name_and_description(fa_header_line):
    """Extracts species name and description from a fasta file header line."""
    linelist = fa_header_line.split()
    name = linelist[0].replace(">", "")
    description = ""
    if len(linelist) > 1:
        description = linelist[1]
    return (name, description)


def test_sequence(seq):
    """Tests if sequences contain data."""
    l = seq.nSpecies
    if l != len(seq.data):
        raise SequenceDataError("List of sequence names and data "
                                "are not of the same length.")
    for i in range(0, l):
        if seq.names[i] == '' or seq.data[i] == '':
            raise SequenceDataError("Sequence name or data is missing.")
    if l > len(set(seq.names)):
        raise SequenceDataError("Sequence names are not unique.")
    return


def open_fa(faFileName, maxskip=50):
    """Opens a fasta file.

    This function tries to open the given fasta file, checks if it is
    in fasta format and reads the sequence(s). It only looks `maxskip`
    lines for the start of a sequence (defaults to 50). It returns an
    FaSeq class object that contains a list of species names, a list
    of the respective desriptions and a list with the sequences.

    """
    seq = FaSeq()

    flag = False
    with open(faFileName) as faFile:
        # Find the start of the first sequence.
        for i in range(0, maxskip):
            line = faFile.readline()
            if line == '':
                raise NotAFastaFileError("File contains no data.")
            if line[0] == '>':
                # species name found in line
                flag = True
                break
        if flag is False:
            raise NotAFastaFileError("Didn't find a species header within " +
                                     str(i+1) + " lines.")
        (name, desc) = get_sp_name_and_description(line)
        seq.names.append(name)
        seq.descr.append(desc)
        data = ""
        for line in faFile:
            if line[0] == '>':
                # new species found in line
                seq.data.append(data)
                seq.dataLen.append(len(data))
                seq.nSpecies += 1
                (name, desc) = get_sp_name_and_description(line)
                seq.names.append(name)
                seq.descr.append(desc)
                data = ""
                line = faFile.readline()
            data += line.replace("\n", "")
    seq.data.append(data)
    seq.dataLen.append(len(data))
    seq.nSpecies += 1
    test_sequence(seq)
    return seq


# TODO
def get_base(chrom, pos, seqNames, seqData):
    """Gets nucleotide base at position `pos` on chromosome `chrom`."""
    















