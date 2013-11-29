#!/usr/bin/env python

"""libPomo.fasta
----------------------------------------------------------------------

This module provides functions to read, write and access fasta files.

"""

import copy
import libPoMo.seqbase as sb
import libPoMo.vcf as vcf


class NotAFastaFileError(sb.SequenceDataError):
    """Exception raised if given fasta file is not valid."""
    pass


class FaSeq():
    """A class that stores sequence data retrieved from a fasta file.

    self.id = fasta sequence identifier
    self.seq = vector of sb.Seq objects (these store the actual sequence data)
    self.nSpecies = number of species (individuals, chromosomes) saved
                    in the object

    """

    def __init__(self):
        self.id = ""
        self.seq = []
        self.nSpecies = 0

    def print_info(self, maxB=50):
        """Print sequence information.

        Print fasta sequence identifier, species names, the length of
        the sequence and a maximum of `maxB` bases (defaults to 50).

        """
        print("Sequence identifier:", self.id)
        for i in range(0, self.nSpecies):
            self.seq[i].print_seq_header()
            print("Printing", maxB, "out of a total of",
                  self.seq[i].dataLen, "bases.")
            print(self.seq[i].data[0:maxB])
        return

    def get_seq_names(self):
        """Returns a list with sequence names."""
        names = []
        for i in range(0, self.nSpecies):
            names.append(self.seq[i].name)
        return names

    def get_seq_by_id(self, i):
        """Return sequence number `i` as Seq object."""
        seq = sb.Seq()
        seq = self.seq[i]
        return seq

    def get_seq_base(self, seq, pos):
        """Returns base at position `pos` in sequence with name `seq`."""
        try:
            i = self.names.index(seq)
        except:
            raise sb.SequenceDataError("Sequence name not found.")
        if pos > self.dataLen[i]:
            raise sb.SequenceDataError("Position out of range.")
        return self.data[i][pos-1]


def get_sp_name_and_description(fa_header_line):
    """Extracts species name and description from a fasta file header line."""
    lineList = fa_header_line.split()
    name = lineList[0].replace(">", "")
    description = ""
    if len(lineList) > 1:
        description = lineList[1]
    return (name, description)


def test_sequence(faSequence):
    """Tests if sequences contain data."""
    l = faSequence.nSpecies
    names = []
    for i in range(0, l):
        names.append(faSequence.seq[i].name)
        if faSequence.seq[i].name == '' or faSequence.seq[i].data[i] == '':
            raise sb.SequenceDataError("Sequence name or data is missing.")
    if l > len(set(names)):
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
    fastaSeq = FaSeq()
    sequence = sb.Seq()

    flag = False
    with open(faFileName) as faFile:
        fastaSeq.name = sb.stripFName(faFileName)
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

        sequence.name = name
        sequence.descr = desc
        data = ""
        for line in faFile:
            if line[0] == '>':
                # new species found in line
                sequence.data = data
                sequence.dataLen = len(data)
                # cp sequence to fastaSeq and purge sequence
                fastaSeq.seq.append(copy.copy(sequence))
                fastaSeq.nSpecies += 1
                sequence.purge()
                (name, desc) = get_sp_name_and_description(line)
                sequence.name = name
                sequence.descr = desc
                data = ""
                line = faFile.readline()
            data += line.replace("\n", "")
    sequence.data = data
    sequence.dataLen = len(data)
    # cp sequence to fastaSeq and purge sequence
    fastaSeq.seq.append(copy.copy(sequence))
    fastaSeq.nSpecies += 1
    test_sequence(fastaSeq)
    return fastaSeq


def save_as_vcf(faSeq, ref, VCFFileName):
    """Saves the given FaSeq in VCF format.

    This function saves the SNPs of `faSeq`, a given FaSeq (fasta
    sequence) object in VCF format to the file `VCFFileName`.  The
    reference genome `ref` to which `faSeq` is compared to needs
    to be passed as a Seq object.

    The function compares all sequences in `faSeq` to the sequence
    given in `ref`.  The names of the individuals in the saved VCF
    file will be the sequence names of the `faSeq` object.

    #CHROM = sequence name of the reference
    POS    = position relative to reference
    ID     = .
    REF    = base of reference
    ALT    = SNP (e.g. 'C' or 'G,T' if 2 different SNPs are present)
    QUAL   = .
    FILTER = .
    INFO   = .
    FORMAT = GT

    """
    def get_altBases_string(sAltBases):
        """Returns ALT bases string from given `altBases`."""
        l = len(sAltBases)
        if l == 0:
            return ''
        string = str(sAltBases[0])
        if l > 1:
            for i in range(1, l):
                string += ',' + sAltBases[i]
        return string

    def get_indiv_string(indivData, altBases, sAltBases):
        """Returns the string of the individual data.

        Returns the string extracted from the indivudal data
        `indivData` with SNPs `altBases`.

        E.g.:
        REF = A
        ALT = C,G
        individual i1 has A
        individual i2 has C
        individual i3 has G

        Then the string should look like:
        '0\t1\t2'
        -> 0 for REF, 1 for first ALT and 2 for second ALT

        """
        l = len(indivData)
        if not (indivData[0] in altBases):
            string = '0'
        else:
            string = str(sAltBases.index(indivData[0]) + 1)
        if l > 1:
            for i in range(1, len(indivData)):
                if not (indivData[0] in altBases):
                    string += '\t' + '0'
                else:
                    string += '\t' + str(sAltBases.index(indivData[0]) + 1)
        return string

    def print_vcf_line_to_file(chromName, pos,
                               refBase, altBaseSring, indivString):
        """Prints a VCF file line with given data to file `VCFFile`."""
        print(altBaseString, indivString)

    if (not isinstance(faSeq, FaSeq)):
        raise sb.SequenceDataError("`faSeq` is not an FaSeq object.")
    if (not isinstance(ref, sb.Seq)):
        raise sb.SequenceDataError("`ref` is not a Seq object.")
    if faSeq.nSpecies == 0:
        raise sb.SequenceDataError("`faSeq` has no saved sequences.")
    for i in range(0, faSeq.nSpecies):
        if faSeq.seq[i].dataLen != ref.dataLen:
            raise sb.SequenceDataError(
                "Sequence " + faSeq.names[i] +
                " has different length than reference.")
    # initialize VCFFile
    with open(VCFFileName, 'w') as VCFFile:
        print(vcf.get_header_line_string(faSeq.get_seq_names()), file=VCFFile)
        # loop over bases
        refBase = ''
        for i in range(0, ref.dataLen):
            refBase = ref.data[i]
            altBases = set()
            indivData = []
            # loop over sequences in faSeq and check if there is a SNP
            for s in range(0, faSeq.nSpecies):
                indivData.append(faSeq.seq[s].data[i])
                if faSeq.seq[s].data[i] != refBase:
                    altBases.add(faSeq.seq[s].data[i])
            sAltBases = sorted(altBases)
            altBaseString = get_altBases_string(sAltBases)
            indivString = get_indiv_string(indivData, altBases, sAltBases)
            # TODO TODO
            print_vcf_line_to_file(
                ref.name, i, refBase, altBaseString, indivString)
    return
