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
    information retrieved from a VCF file.  Please see
    http://www.1000genomes.org/ for a detailed description of the vcf
    format.

    self.chrom = chromosome name
    self.pos = position on chromosome
    self.id
    self.ref
    self.alt
    self.qual
    self.filter
    self.info
    self.format
    self.dataIndividuals

    """
    def __init__(self):
        self.chrom = ''
        self.pos = 0
        self.id = ''
        self.ref = ''
        self.alt = ''
        self.qual = ''
        self.filter = ''
        self.info = ''
        self.format = ''
        self.dataIndividuals = []

    def print_info(self):
        """Print nucleotide base information.

        Prints the stored single nucleotide base and related
        information from the VCF file.

        """
        print(self.chrom, self.pos, self.id, self.ref,
              self.alt, self.qual, self.filter,
              self.info, self.format,
              sep='\t', end='\t')
        print(self.dataIndividuals, sep='\t', end='')
        return


def get_header_line_string(indiv):
    """Returns a standard VCF File header string with individuals `indiv`."""
    string = ''
    for s in hdList:
        string += s + '\t'
    for i in indiv:
        string += i + '\t'
    # we added one tab at the end that we do not need
    return string[:-1]


def print_header_line(indiv):
    """Prints a standard VCF File header with individuals `indiv`."""
    print(get_header_line_string(indiv), end='')


class VCFSeq():
    """A class that stores data retrieved from a VCF file.

    self.header = header information
    self.bases = list with stored NucBase(s)
    self.nbases = number of bases stored
    """
    def __init__(self):
        self.header = []
        self.individuals = []
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
        print_header_line(self.individuals)
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
        for i in range(0, self.nbases):
            if pos == self.bases[i].pos \
               and chrom == self.bases[i].chrom:
                return self.bases[i]
        raise sb.SequenceDataError('Base at position ' + str(pos) +
                                   ' on chromosome ' + str(chrom) +
                                   ' not found.')


def check_fixed_field_header(ln):
    """Checks if the given line is the header of the fixed fields.

    Sample header line:
    #CHROM\t POS\t ID\t REF\t ALT\t QUAL\t FILTER\t INFO\t FORMAT\t Individuals

    """
    lnList = ln.split('\t', maxsplit=9)
    if lnList[0:9] != hdList:
        raise NotAVariantCallFormatFileError('Header line is invalid.')
    return


def get_indiv_from_field_header(ln):
    """Returns individuals from a fixed field header.

    Sample header line:
    #CHROM\t POS\t ID\t REF\t ALT\t QUAL\t FILTER\t INFO\t FORMAT\t Individuals

    """
    individuals = []
    lnList = ln.split('\t', maxsplit=9)
    if len(lnList) == 10:
        individuals = lnList[9].split('\t')
    else:
        raise NotAVariantCallFormatFileError('No individuals in header line.')
    return individuals


def get_nuc_base_from_line(ln):
    """Retrieves base data from a VCF file line.

    Splits a given VCF file line and returns a NucBase object.

    """
    base = NucBase()
    lnList = ln.split('\t', maxsplit=9)
    if len(lnList) >= 10:
        base.chrom = lnList[0]
        base.pos = int(lnList[1])
        base.id = lnList[2]
        base.ref = lnList[3]
        base.alt = lnList[4]
        base.qual = lnList[5]
        base.filter = lnList[6]
        base.info = lnList[7]
        base.format = lnList[8]
        base.dataIndividuals = lnList[9]
    else:
        raise NotANucBaseError('Line ' + ln + ' is not a NucBase.')
    return base


def test_sequence(seq):
    """Tests a given VCF sequence."""
    pass                        # TODO


def open_seq(VCFFileName, maxskip=100):
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
                seq.individuals = get_indiv_from_field_header(line)
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
