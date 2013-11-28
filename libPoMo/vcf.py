#!/usr/bin/env python

"""libPomo.vcf
----------------------------------------------------------------------

This module provides functions to read, write and access vcf files.

"""

import libPoMo.seqbase as sb


class NotAVariantCallFormatFileError(sb.SequenceDataError):
    """Exception raised if given VCF file is not valid."""
    pass


class NotANucBaseError(sb.SequenceDataError):
    """Exception raised if given nucleotide base is not valid."""
    pass

hdList = ['#CHROM', 'POS', 'ID', 'REF', 'ALT',
          'QUAL', 'FILTER', 'INFO', 'FORMAT']


class NucBase():
    """Stores a nucleotide base.

    A class that stores a single nucleotide base and related
    information retrieved from a VCF file.

    self.chrom = chromosome name
    self.pos = position on chromosome
    self.id
    self.ref
    self.alt
    self.qual
    self.filter

    """
    def __init__(self):
        self.chrom = ''
        self.pos = 0
        self.id = ''
        self.ref = ''
        self.alt = ''
        self.qual = ''
        self.filter = ''
        # self.info = ''
        # self.format = ''
        # self.indivuals = ''

    def print_header_line(self):
        print(*hdList, sep='\t')
        return

    def print_info(self, printHeader=False):
        """Print nucleotide base information.

        Prints the stored single nucleotide base and related
        information from the VCF file. If `printHeader=True` is
        specified, the header line is printed before the nucleotide
        base.

        """
        if printHeader is True:
            self.print_header_line()
        print(self.chrom, self.pos, self.id, self.ref,
              self.alt, self.qual, self.filter,
              sep='\t')
        return


class VCFSeq():
    """A class that stores data retrieved from a VCF file.

    self.header = header information
    self.bases = list with stored NucBase(s)
    self.nbases = number of bases stored
    """
    def __init__(self):
        self.header = []
        self.bases = []
        self.nbases = 0

    def print_info(self, maxB=50, printHeader=False):
        """Print VCF sequence information.

        Print vcf header, the total number of nucleotides and a
        maximum of `maxB` bases (defaults to 50). Only prints header
        if `printHeader=True` is given.

        """
        if printHeader is True:
            print(self.header)
        if self.nbases > 0:
            self.bases[0].print_header_line()
        if self.nbases < maxB:
            maxB = self.nbases
        for i in range(0, maxB):
            self.bases[i].print_info()
        return

    def append_nuc_base(self, base):
        """Appends a given NucBase to the VCFSeq object."""
        self.bases.append(base)
        self.nbases += 1
        return

    def get_nuc_base(self, chrom, pos):
        """Returns base at position `pos` of chromosome `chrom`."""
        startIndex = 0
        while True:
            try:
                i = self.pos[startIndex:].index(pos)
            except:
                raise sb.SequenceDataError("Base position not found.")
            if self.chrom[i+startIndex] == chrom:
                # Base has been found.
                break
            else:
                # Start next search right after this pos
                startIndex = startIndex + i + 1
        return self.ref[i+startIndex]


def check_fixed_field_header(ln):
    """Checks if the given line is the header of the fixed fields.

    Sample header line:
    #CHROM\t POS\t ID\t REF\t ALT\t QUAL\t FILTER\t INFO\t FORMAT\t Individuals

    """
    lnList = ln.split('\t', maxsplit=9)
    if lnList[0:9] != hdList:
        raise NotAVariantCallFormatFileError('Header line is invalid.')
    return


def get_nuc_base_from_line(ln):
    """Retrieves base data from a VCF file line.

    Splits a given VCF file line and returns a NucBase object.

    """
    base = NucBase()
    lnList = ln.split('\t', maxsplit=9)
    if len(lnList) >= 7:
        base.chrom = lnList[0]
        base.pos = lnList[1]
        base.id = lnList[2]
        base.ref = lnList[3]
        base.alt = lnList[4]
        base.qual = lnList[5]
        base.filter = lnList[6]
    else:
        raise NotANucBaseError('Line ' + ln + ' is not a NucBase.')
    return base


def test_sequence(seq):
    """Tests a given VCF sequence."""
    pass                        # TODO


def open_vcf(VCFFileName, maxskip=100):
    """Opens a VCF4.2 file.

    This function tries to open the given VCF file, checks if it is in
    VCF format and reads the bases(s). It only looks `maxskip` lines
    for the start of the bases (defaults to 80). It returns an VCFSeq
    class object that contains all the information. For help, refer to
    `VCFSeq.__doc__`.

    """
    seq = VCFSeq()
    seq.header = ""

    flag = False
    with open(VCFFileName) as VCFFile:
        # Find the start of the first base
        for i in range(0, maxskip):
            line = VCFFile.readline()
            if line[0:2] == '##':
                seq.header += line
            if line[0:6] == '#CHROM':
                # Here starts the data.
                check_fixed_field_header(line)
                flag = True
                break
            if line == '':
                raise NotAVariantCallFormatFileError("File contains no data.")
        if flag is False:
            raise NotAVariantCallFormatFileError(
                "Didn't find any data within " + str(maxskip) + " lines.")
        for line in VCFFile:
            base = get_nuc_base_from_line(line)
            seq.append_nuc_base(base)

        test_sequence(seq)
    return seq


















