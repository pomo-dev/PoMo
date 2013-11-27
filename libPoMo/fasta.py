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


def get_sp_name_and_description(fa_header_line):
    """Extracts species name and description from a fasta file header line."""
    linelist = fa_header_line.split()
    name = linelist[0].replace(">", "")
    description = ""
    if len(linelist) > 1:
        description = linelist[1]
    return (name, description)


def test_sequence(seqNames, seqData):
    """Tests if sequences contain data."""
    l = len(seqNames)
    if l != len(seqData):
        raise 

def open_fa(faFileName, maxskip=50):
    """Opens a fasta file.

    This function tries to open the given fasta file, checks if it is
    in fasta format and reads the sequence(s). It only looks `maxskip`
    lines for the start of a sequence (defaults to 50). It returns a
    tuple that contains a list of species names, a list of the
    respective desriptions and a list with the sequences.

    """
    seqNames = []
    seqDescr = []
    seqData = []

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
            raise NotAFastaFileError()
        (name, desc) = get_sp_name_and_description(line)
        seqNames.append(name)
        seqDescr.append(desc)
        data = ""
        for line in faFile:
            if line[0] == '>':
                # new species found in line
                seqData.append(data)
                (name, desc) = get_sp_name_and_description(line)
                seqNames.append(name)
                seqDescr.append(desc)
                data = ""
                line = faFile.readline()
            data += line.replace("\n", "")
    seqData.append(data)
    test_sequence(seqNames, seqData)
    return (seqNames, seqDescr, seqData)


# TODO
def get_base(chrom, pos):
    """Gets nucleotide base at position `pos` on chromosome `chrom`."""
    pass















