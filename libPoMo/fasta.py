#!/usr/bin/env python

"""libPomo.fasta
----------------------------------------------------------------------

This module provides functions to read, write and access fasta files.

"""

import libPoMo.seqbase as sb
import libPoMo.vcf as vcf


class NotAFastaFileError(sb.SequenceDataError):
    """Exception raised if given fasta file is not valid."""
    pass


class FaSeq(sb.Seq):
    """A class that stores sequence data retrieved from a fasta file."""
    def __init__(self):
        super(FaSeq, self).__init__()
        self.descr = []

    def print_seq_header(self, i):
        print('>', self.names[i], ' ', self.descr[i], sep='')
        return


def get_sp_name_and_description(fa_header_line):
    """Extracts species name and description from a fasta file header line."""
    lineList = fa_header_line.split()
    name = lineList[0].replace(">", "")
    description = ""
    if len(lineList) > 1:
        description = lineList[1]
    return (name, description)


def test_sequence(seq):
    """Tests if sequences contain data."""
    l = seq.nSpecies
    if l != len(seq.data):
        raise sb.SequenceDataError("List of sequence names and data "
                                   "are not of the same length.")
    for i in range(0, l):
        if seq.names[i] == '' or seq.data[i] == '':
            raise sb.SequenceDataError("Sequence name or data is missing.")
    if l > len(set(seq.names)):
        raise sb.SequenceDataError("Sequence names are not unique.")
    return


def open_seq(faFileName, maxskip=50):
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
        seq.name = sb.stripFName(faFileName)
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
                                     maxskip + " lines.")
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


def save_as_vcf(faSeq, faRef, VCFFileName):
    """Saves the given FaSeq in VCF format.

    This function saves the SNPs of `faSeq`, a given FaSeq (fasta
    sequence) object in VCF format to the file `VCFFileName`.  The
    reference genome `faRef` to which `faSeq` is compared to needs
    also be passed as an FaSeq object.

    The function compares all sequences in `faSeq` to the sequence
    given in `faRef`.  The names of the individuals in the saved VCF
    file will be the sequence names of the `faSeq` object.  The name
    of the chromosome (#CHROM) will be the sequence name of the
    reference.

    """
    if faRef.nSpecies != 1:
        raise sb.SequenceDataError('Reference contains more than 1 sequence.')
    for i in range(0, faSeq.nSpecies):
        if faSeq.dataLen[i] != faRef.dataLen[0]:
            raise sb.SequenceDataError(
                'Sequence ' + faSeq.names[i] +
                ' has different length than reference.')
    # initialize VCFFile
    with open(VCFFileName, 'w') as VCFFile:
        vcf.file_init(VCFFile)
        # loop over sequences in faSeq
        for s in range(0, faSeq.nSpecies):
            # loop over bases
            for i in range(0, faRef.dataLen[0]):
                pass




















