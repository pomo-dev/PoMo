#!/usr/bin/env python

"""libPomo.vcf
----------------------------------------------------------------------

This module provides functions to read, write and access vcf files.

Classes:
- `NucBase`, store a nucleotide base
- `VCFStream`, a variant call format (VCF) stream object
- `VCFSeq`, a VCF file sequence object

Exception Classes:
- `NotAVariantCallFormatFileError`
- `NotANucBaseError`

Functions:
- `update_base()`: Read a line into a base.
- `get_nuc_base_from_line()`: Create a new `NucBase` from a line.
- `check_fixed_field_header()`: Check a VCF fixed field header string.
- `get_indiv_from_field_header()`: Extract list of individuals from header.

"""

__docformat__ = 'restructuredtext'

import libPoMo.seqbase as sb
import gzip


class NotAVariantCallFormatFileError(sb.SequenceDataError):
    """Exception raised if given VCF file is not valid."""
    pass


class NotANucBaseError(sb.SequenceDataError):
    """Exception raised if given nucleotide base is not valid."""
    pass

hdList = ['#CHROM', 'POS', 'ID', 'REF', 'ALT',
          'QUAL', 'FILTER', 'INFO', 'FORMAT']


def update_base(ln, base, info=True):
    """Read line `ln` into base `base`.

    Split a given VCF file line and returns a NucBase object. If
    `info` is set to False, only #CHROM, REF, ALT and speciesData will
    be read.

    """
    lnList = ln.split('\t', maxsplit=9)
    if len(lnList) >= 10:
        base.chrom = lnList[0]
        base.pos = int(lnList[1])
        base.ref = lnList[3]
        base.alt = lnList[4]
        base.speciesData = lnList[9].rstrip().split('\t')
        # base.speciesData = [s.split(':', maxsplit=1)[0]
        #                    for s in lnList[9].rstrip().split('\t')]
        if info is True:
            base.id = lnList[2]
            base.qual = lnList[5]
            base.filter = lnList[6]
            base.info = lnList[7]
            base.format = lnList[8]
    else:
        raise NotANucBaseError('Line ' + ln + ' is not a NucBase.')
    return base


def get_nuc_base_from_line(ln, info=False, ploidy=None):
    """Retrieve base data from a VCF file line `ln`.

    Split a given VCF file line and returns a NucBase object. If
    `info` is set to False, only #CHROM, POS, REF, ALT and speciesData will
    be read.
    
    - `info`: Boolean; determines if info is retrieved from `ln`.
    - `ploidy`: If ploidy is known and given, it is set.

    """
    base = NucBase()
    update_base(ln, base, info)
    if ploidy is not None:
        base.ploidy = ploidy
    return base


class NucBase():
    # TODO SPECIES DATA WITH / and |
    """Stores a nucleotide base.

    A class that stores a single nucleotide base and related
    information retrieved from a VCF file.  Please see
    http://www.1000genomes.org/ for a detailed description of the vcf
    format.

    self.haploid: True if self.speciesData contains haploid data. Can
    be set with self.set_haploid()

    """
    def __init__(self):
        self.chrom = ''
        """Chromosome name."""
        self.pos = 0
        """Position on the chromosome."""
        self.id = ''
        self.ref = ''
        self.alt = ''
        self.qual = ''
        self.filter = ''
        self.info = ''
        self.format = ''
        self.speciesData = []
        self.ploidy = None
        """Number of sets of chromosomes. Set by `self.set_ploidy()`."""

    def print_info(self):
        """Print nucleotide base information.

        Print the stored single nucleotide base and related
        information from the VCF file.

        """
        print(self.chrom, self.pos, self.id, self.ref,
              self.alt, self.qual, self.filter,
              self.info, self.format,
              sep='\t', end='\t')
        print('\t'.join(self.speciesData))
        return

    def get_alt_base_list(self):
        """Return alternative bases as a list."""
        altBases = [b.lower() for b in self.alt.split(',')]
        return altBases

    def set_ploidy(self):
        """Set self.ploidy."""
        baseInfo = self.speciesData[0].split(':')[0]
        self.ploidy = len(baseInfo.split('/'))
        return self.ploidy

    def get_speciesData(self):
        """Return species data as a list.
        
        data[0][0] = data of first species/individual on chromatide A
        data[0][1] = only set for non-haploids; data of first
        species/individual on chromatide B

        Return None if bases saved are not valid ("./."). This means,
        that the information from all other individuals or species is
        lost too. Maybe this should be changed in the future.

        """
        data = []
        for i in range(0, len(self.speciesData)):
            if self.ploidy == 1:
                # Haploid.
                baseInfo = self.speciesData[i].split(':')[0]
                try:
                    baseInfo = int(baseInfo)
                except ValueError:
                    # Invalid Base.
                    baseInfo = None
                data.append([baseInfo])
            else:
                # Diploid or even more
                baseInfo = self.speciesData[i].split(':')[0]
                baseInfoL = baseInfo.split('/')
                for j in range(len(baseInfoL)):
                    try:
                        baseInfoL[j] = int(baseInfoL[j])
                    except ValueError:
                        # Invalid Base.
                        baseInfoL[j] = None
                data.append(baseInfoL)
        return data

    def purge(self):
        """Purge the data associated with this `NucBase`."""
        self.__init__()


class VCFStream():
    """Store base data from a VCF file line per line.

    This class stores a single base retrieved from a VCF file and the
    file itself.  It is used to parse through a VCF file line by line
    processing the bases without having to read the whole file at one.

    It can be initialized with init_seq().

    """

    def __init__(self, seqName, vcfFileObject, speciesList, firstBase):
        """Initialize a `VCFStream` object; add state objects.

        Parameters:

        - `seqName`: Name of the stream.
        - `vcfFileObject`: File object associated with the stream.
        - `speciesList`: List with species/individuals.
        - `firstBase`: First base to be saved.

        """
        self.name = seqName
        """Sequence name."""
        self.fo = vcfFileObject
        """Stored VCF file object."""
        self.speciesL = speciesList
        """List with species/individuals."""
        self.nSpecies = len(speciesList)
        """Number of species/individuals."""
        self.base = firstBase
        """Stored `NucBase`."""

    def print_info(self):
        """Prints VCFStream information."""
        print("Name:", self.name)
        print("File object:", self.fo)
        print("List of species/individuals:", self.speciesL)
        print("Number of species/individuals:", self.nSpecies)
        print("Saved base:")
        self.base.print_info()

    def read_next_base(self):
        """Read the next base.

        Return position of next base.
        Raise a `ValueError` if no next base is found.

        """
        line = self.fo.readline()
        if line != '':
            update_base(line, self.base)
            return self.base.pos
        else:
            raise ValueError("End of VCFStream.")
            return None

    def close_fo(self):
        """Closes the linked file."""
        self.fo.close()


class VCFSeq():
    """Store data retrieved from a VCF file."""
    def __init__(self):
        """Initialize a `VCFSeq` object; add state objects."""
        self.name = ''
        """Sequence name."""
        self.header = []
        """Header information."""
        self.speciesL = []
        """List with species/individuals."""
        self.nSpecies = 0
        """Number of species/individuals."""
        self.baseL = []
        """List with stored `NucBase` objects."""
        self.nBases = 0
        """Number of `NucBase` objects stored."""

    def get_header_line_string(self, indiv):
        """Return a standard VCF File header string with individuals `indiv`.

        """
        string = ''
        for s in hdList:
            string += s + '\t'
        for i in indiv:
            string += i + '\t'
        # we added one tab at the end that we do not need
        return string[:-1]

    def print_header_line(self, indiv):
        """Print a standard VCF File header with individuals `indiv`."""
        print(self.get_header_line_string(indiv))

    def print_info(self, maxB=50, printHeader=False):
        """Print VCF sequence information.

        Print vcf header, the total number of nucleotides and a
        maximum of `maxB` bases (defaults to 50). Only prints header
        if `printHeader=True` is given.

        """
        if printHeader is True:
            print(self.header)
        self.print_header_line(self.speciesL)
        if self.nBases < maxB:
            maxB = self.nBases
        for i in range(0, maxB):
            self.baseL[i].print_info()
        return

    def append_nuc_base(self, base):
        """Append `base`, a given `NucBase`, to the VCFSeq object."""
        self.baseL.append(base)
        self.nBases += 1
        return

    def has_base(self, chrom, pos):
        """Return True (False) if base is (not) found.

        - `chrom`: Chromosome name.
        - `pos`: Position on `chrom`.

        """
        for i in range(0, self.nBases):
            if pos == self.baseL[i].pos \
               and chrom == self.baseL[i].chrom:
                return True
        return False

    def get_nuc_base(self, chrom, pos):
        """Return base at position `pos` of chromosome `chrom`."""
        for i in range(0, self.nBases):
            if pos == self.baseL[i].pos \
               and chrom == self.baseL[i].chrom:
                return self.baseL[i]
        raise sb.SequenceDataError('Base at position ' + str(pos) +
                                   ' on chromosome ' + str(chrom) +
                                   ' not found.')


def check_fixed_field_header(ln):
    """Check if the given line `ln` is the header of the fixed fields.

    Sample header line:
    #CHROM\t POS\t ID\t REF\t ALT\t QUAL\t FILTER\t INFO\t FORMAT\t SpeciesL

    """
    lnList = ln.split('\t', maxsplit=9)
    if lnList[0:9] != hdList:
        raise NotAVariantCallFormatFileError('Header line is invalid.')
    return


def get_indiv_from_field_header(ln):
    """Return species from a fixed field header line `ln`.

    Sample header line:
    #CHROM\t POS\t ID\t REF\t ALT\t QUAL\t FILTER\t INFO\t FORMAT\t SpeciesL

    """
    speciesL = []
    lnList = ln.split('\t', maxsplit=9)
    if len(lnList) == 10:
        speciesL = lnList[9].rstrip().split('\t')
    else:
        raise NotAVariantCallFormatFileError('No species in header line.')
    return speciesL


def init_seq(VCFFileName, maxskip=100, name=None):
    """Open a (gzipped) VCF4.2 file.

    Try to open the given VCF file, checks if it is in VCF format.
    Initialize a `VCFStream` object that contains the first base.  For
    help, refer to `VCFSeq.__doc__`.

    Please close the associated file object with
    yourVCFStream.close_fo() when you don't need it anymore.

    - `maxskip`: Only look `maxskip` lines for the start of the bases
    (defaults to 80).

    - `name`: Set the name of the sequence to `name`, otherwise set it
    to the filename.

    """
    flag = False
    if VCFFileName[-2:] == "gz":
        VCFFile = gzip.open(VCFFileName, mode='rt')
    else:
        VCFFile = open(VCFFileName)
    # set the vcf sequence name
    if name is None:
        name = sb.stripFName(VCFFileName)
    # Find the start of the first base
    for i in range(0, maxskip):
        line = VCFFile.readline()
        if line == '':
            raise NotAVariantCallFormatFileError("File contains no data.")
        if line[0:6] == '#CHROM':
            # Here starts the data.
            check_fixed_field_header(line)
            speciesL = get_indiv_from_field_header(line)
            flag = True
            break
    if flag is False:
        raise NotAVariantCallFormatFileError(
            "Didn't find any data within " + str(maxskip) + " lines.")
    line = VCFFile.readline()
    base = get_nuc_base_from_line(line, info=False)
    base.set_ploidy()
    return VCFStream(name, VCFFile, speciesL, base)


def open_seq(VCFFileName, maxskip=100, name=None):
    """Open a VCF4.2 file.

    Try to open the given VCF file, checks if it is in VCF format and
    reads the bases(s).  It returns an `VCFSeq` object that contains
    all the information. For help, refer to `VCFSeq.__doc__`.

    - `maxskip`: Only look `maxskip` lines for the start of the bases
    (defaults to 80).

    - `name`: Set the name of the sequence to `name`, otherwise set it
    to the filename.

    """
    def test_sequence(seq):
        """Test a given VCF sequence. TODO."""
        pass

    seq = VCFSeq()
    seq.header = ""

    flag = False
    if VCFFileName[-2:] == "gz":
        VCFFile = gzip.open(VCFFileName, mode='rt')
    else:
        VCFFile = open(VCFFileName)
    # set the vcf sequence name
    if name is not None:
        seq.name = name
    else:
        seq.name = sb.stripFName(VCFFileName)
    # Find the start of the first base
    for i in range(0, maxskip):
        line = VCFFile.readline()
        if line == '':
            raise NotAVariantCallFormatFileError("File contains no data.")
        if line[0:2] == '##':
            seq.header += line
        if line[0:6] == '#CHROM':
            # Here starts the data.
            check_fixed_field_header(line)
            seq.speciesL = get_indiv_from_field_header(line)
            seq.nSpecies = len(seq.speciesL)
            flag = True
            break
    if flag is False:
        raise NotAVariantCallFormatFileError(
            "Didn't find any data within " + str(maxskip) + " lines.")
    for line in VCFFile:
        base = get_nuc_base_from_line(line)
        seq.append_nuc_base(base)

    VCFFile.close()
    test_sequence(seq)
    return seq
