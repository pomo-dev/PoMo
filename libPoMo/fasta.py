#!/usr/bin/env python

"""libPoMo.fasta
================

This module provides functions to read, write and access fasta files.

Objects
-------
Classes:
  - :class:`FaStream`, fasta file sequence stream object
  - :class:`MFaStream`, multiple alignment fasta file sequence stream object
  - :class:`FaSeq`, fasta file sequence object
  - :class:`MFaStrFilterProps`, define multiple fasta file filter preferences

Exception Classes:
  - :class:`NotAFastaFileError`

Functions:
  - :func:`filter_mfa_str()`, filter a given :class:`MFaStream`
    according to the filters defined in :class:`MFaStrFilterProps`
  - :func:`init_seq()`, initialize fasta sequence stream from file
  - :func:`open_seq()`, open fasta file
  - :func:`save_as_vcf()`, save a given :class:`FaSeq` in variant call
    format (VCF)
  - :func:`read_seq_from_fo()`, read a single sequence from file object
  - :func:`read_align_from_fo()`, read an alignment from file object

----

"""
__docformat__ = 'restructuredtext'

import libPoMo.seqbase as sb
import libPoMo.vcf as vcf
import sys
import re


class NotAFastaFileError(sb.SequenceDataError):
    """Exception raised if given fasta file is not valid."""
    pass


def read_seq_from_fo(line, fo, getAlignEndFlag=False):
    """Read a single fasta sequence.

    Read a single fasta sequence from file object *fo* and save it to
    a new :class:`Seq <libPoMo.seqbase.Seq>` sequence object. Return
    the header line of the next fasta sequence and the newly created
    sequence. If no new sequence is found, the next header line will
    be set to None.

    :param str line: Header line of the sequence.
    :param fo fo: File object of the fasta file.

    :param Boolean getAlignFlag: If set to true, an additional Boolean
      value that specifies if a multiple sequence alignment ends, is
      returned.

    :rtype: (str, Seq) | (str, Seq, Boolean)

    """
    def get_sp_name_and_description(fa_header_line):
        """Extract species name and description.

        Extract species name and description from a fasta file header
        line `fa_header_line`.

        """
        lineList = fa_header_line.rstrip().split(maxsplit=1)
        name = lineList[0][1:]
        description = ""
        if len(lineList) > 1:
            description = lineList[1]
        return (name, description)

    def fill_seq_from_fo(line, fo, seq):
        """Read a single fasta sequence.

        Read a single fasta sequence from file object `fo` and save it
        to `seq`. Returns the next header line and a flag that is set
        to true if the end of an alignment is reached (a line only
        contains a newline character).  If no new sequence is found,
        the next header line will be set to None.

        :param str line: Header line of the sequence.
        :param fo for: File object of the fasta file.
        :param Seq seq: The sequence that will be filled.

        """
        (name, descr) = get_sp_name_and_description(line)
        seq.name = name
        seq.descr = descr
        data = ""
        alignEndFl = False
        for line in fo:
            if line == '\n':
                # Newline found, end of alignment.
                alignEndFl = True
            elif line[0] == '>':
                # New species found in line.
                break
            else:
                data += line.rstrip()
        seq.data = data
        seq.dataLen = len(data)
        if line[0] != '>':
            # We reached the end of file.
            line = None
        return (line, alignEndFl)

    seq = sb.Seq()
    (newHeaderLine, alignEndFl) = fill_seq_from_fo(line, fo, seq)
    if getAlignEndFlag is False:
        return (newHeaderLine, seq)
    else:
        return (newHeaderLine, seq, alignEndFl)


class FaStream():
    """A class that stores a fasta file sequence stream.

    The sequence of one species / individual / chromosome is saved and
    functions are provided to read in the next sequence in the file,
    if there is any. This saves memory if files are huge and doesn't
    increase runtime.

    This object is usually initialized with :func:`init_seq`.

    :param str name: Name of the stream.
    :param Seq firstSeq: First sequence (:class:`Seq
                         <libPoMo.seqbase.Seq>` object) to be saved.
    :param str nextHL: Next header line.
    :param fo faFileObject: File object associated with the stream.

    :ivar str name: Stream name.
    :ivar Seq seq: Saved sequence (:class:`Seq
                   <libPoMo.seqbase.Seq>` object)
    :ivar str nextHeaderLine: Next header line.
    :ivar fo fo: File object that points to the start of the data of
                 the next sequence.

    """

    def __init__(self, name, firstSeq, nextHL, faFileObject):
        """Initialize an `FaStream` object; add state objects."""
        self.name = name
        self.seq = firstSeq
        self.nextHeaderLine = nextHL
        self.fo = faFileObject

    def print_info(self, maxB=50):
        """Print sequence information.

        Print information about this FaStream object, the fasta
        sequence stored at the moment the length of the sequence and a
        maximum of `maxB` bases (defaults to 50).

        """
        print("Associated file object:", self.fo)
        print("Next header line:", self.nextHeaderLine)
        print("Saved Sequence:")
        self.seq.print_fa_header()
        print("Printing", maxB, "out of a total of",
              self.seq.dataLen, "bases.")
        print(self.seq.data[0:maxB])
        return

    def read_next_seq(self):
        """Read next fasta sequence in file.

        The return value is the name of the next sequence or None if
        no next sequence is found.

        """
        if self.nextHeaderLine is None:
            return None
        else:
            self.seq.purge()
            (nextHL, self.seq) = read_seq_from_fo(self.nextHeaderLine, self.fo)
            self.nextHeaderLine = nextHL
            return self.seq.name

    def close(self):
        """Close the linked file."""
        self.fo.close()


def read_align_from_fo(line, fo):
    """Read a single fasta alignment.

    Read a single fasta alignment from file object *fo* and save it to
    new :class:`Seq <libPoMo.seqbase.Seq>` sequence objects.  Return
    the header line of the next fasta alignment and the newly created
    sequences in a list.  If no new alignment is found, the next header
    line will be set to None.

    :param str line: Header line of the sequence.
    :param fo fo: File object of the fasta file.
    :rtype: (str, [Seq])

    """
    alignEndFl = False
    seqL = []
    headerLn = line
    while (alignEndFl is not True) and (headerLn is not None):
        (newHeaderLn, seq, alignEndFl) = read_seq_from_fo(headerLn, fo,
                                                          getAlignEndFlag=True)
        headerLn = newHeaderLn
        seq.set_rc()
        seqL.append(seq)
    return (newHeaderLn, seqL)


class MFaStream():
    """Store a multiple alignment fasta file sequence stream.

    The sequences of one gene / alignment are saved for all species /
    individuals / chromosomes.  Functions are provided to read in the
    next gene / alignment in the file that fulfills the given
    criteria, if there is any.  This saves memory if files are huge and
    doesn't increase runtime.

    Initialization of an :class:`MFaStream` opens the given fasta
    file, checks if it is in fasta format and reads the first
    alignment.  The end of an alignment is reached when a line only
    contains the newline character.  This object can later be used to
    parse the whole multiple alignment fasta file.

    Alignments can be filtered with :func:`filter_mfa_str()`.

    :param str faFileName: File name of the multiple alignment fasta file.
    :param int maxskip: Only look *maxskip* lines for the start of a
      sequence (defaults to 50).
    :param str name: Set the name of the stream to *name*, otherwise
      set it to the stripped filename.

    :ivar str name: Stream name.
    :ivar [Seq] seqL: Saved sequences (:class:`Seq
                      <libPoMo.seqbase.Seq>` objects) in a list.
    :ivar int nSpecies: Number of saved sequences / species in the alignment.
    :ivar str nextHeaderLine: Next header line.
    :ivar fo fo: File object that points to the start of the data of
                 the next sequence.

    Please close the associated file object with
    :func:`FaStream.close` when you don't need it anymore.

    """

    def __init__(self, faFileName, maxskip=50, name=None):
        """Open a fasta file and initialize :class:`MFaStream`."""
        def add_instance_variables(name, firstSeqL, nextHL, faFileObject):
            """Add state objects."""
            self.name = name
            self.seqL = firstSeqL
            self.nSpecies = len(self.seqL)
            self.nextHeaderLine = nextHL
            self.fo = faFileObject

        flag = False
        faFile = sb.gz_open(faFileName)
        if name is None:
            name = sb.stripFName(faFileName)
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
        (nextHL, seqL) = read_align_from_fo(line, faFile)
        try:
            nextHL = nextHL.rstrip()
        except:
            pass
        add_instance_variables(name, seqL, nextHL, faFile)

    def print_info(self, maxB=50):
        """Print sequence information.

        Print information about this MFaStream object, the fasta
        sequence stored at the moment the length of the sequence and a
        maximum of `maxB` bases (defaults to 50).

        """
        print("Associated file object:", self.fo)
        print("Next header line:", self.nextHeaderLine)
        print("Saved Sequences:")
        for i in range(self.nSpecies):
            self.seqL[i].print_fa_header()
            if self.seqL[i].get_rc() is True:
                print("Sequence is reversed and complemented.")
            print("Printing", maxB, "out of a total of",
                  self.seqL[i].dataLen, "bases.")
            print(self.seqL[i].data[0:maxB])
        return

    def read_next_align(self):
        """Read next alignment in fasta file.

        The return value is the name of the newly saved alignment or
        None if no next alignment is found.

        """
        if self.nextHeaderLine is None:
            return None
        else:
            (nextHL, self.seqL) = read_align_from_fo(self.nextHeaderLine,
                                                     self.fo)
            self.nextHeaderLine = nextHL
            self.nSpecies = len(self.seqL)
            return self.seqL[0].name

    def orient(self, firstOnly=False):
        """Orient all sequences of the alignment to be in forward direction.

        This is rather slow for long sequences.

        :param Boolean firstOnly: If true, orient the first sequence only.

        """
        if firstOnly is False:
            l = self.nSpecies
        elif firstOnly is True:
            l = 1
        else:
            raise ValueError()

        for i in range(l):
            if self.seqL[i].get_rc() is True:
                self.seqL[i].rev_comp()

    def print_msa(self, fo=sys.stdout):
        """Print multiple sequence alignment at point.

        :ivar fileObject fo: Print to file object fo. Defaults to
          stdout.

        """
        for s in self.seqL:
            pass
            s.print_fa_header(fo=fo)
            s.print_data(fo=fo)
        print('\n', file=fo)
        return

    def close(self):
        """Close the linked file object."""
        self.fo.close()


class MFaStrFilterProps():
    """Define filter preferences for multiple fasta alignments.


    Define the properties of the filter to be applied to an
    :class:`MFaStream`.

    By default, all filters are applied (all variables are set to
    True).

    :param int nSpecies: Number of species that are aligned.

    :ivar Boolean check_all_aligned: Check if all treated species
      are available in the alignment (`nSpecies` gives the number
      of species, given to the object upon initialization).
    :ivar Boolean check_divergence: Check if the divergence of the
      reference genome (the first sequence in the alignment) is lower
      than `maxDiv` (defaults to 10 percent).
    :ivar Boolean check_start_codons: Check if all start codons
      are conserved.
    :ivar Boolean check_stop_codons: Check if all stop codons are
      conserved.
    :ivar Boolean check_frame_shifting_gaps: Check, that there
      are no frame-shifting gaps.
    :ivar Boolean check_for_long_gaps: Check if no gap is longer
      than `maxGapLength` (defaults to 30) bases.
    :ivar Boolean check_nonsense_codon: Check if there is no
      premature stop codon).
    :ivar Boolean check_exon_length: Check that the exon is
      longer than `minExonLen` (defaults to 21).
    :ivar Boolean check_exon_numbers: Check if exon number match
      for all sequences in the alignment.

    """

    def __init__(self, nSpecies):
        self.check_all_aligned = True
        self.nSpecies = nSpecies
        self.check_divergence = True
        self.maxDiv = 0.1
        self.check_start_codons = True
        self.check_stop_codons = True
        self.check_frame_shifting_gaps = True
        self.check_for_long_gaps = True
        self.maxGapLength = 30
        self.check_nonsense_codon = True
        self.check_exon_length = True
        self.minExonLen = 21
        self.check_exon_numbers = True


def filter_mfa_str(mfaStr, fp, verb=None):
    """Check multiple sequence alignment of an MFaStream.

    Multiple sequence alignments usually include alignments that are
    not apt for analysis.  These low quality alignments need to be
    filtered out of the original multiple sequence alignment fasta
    file.  If `verb` is unset from None, information about any
    possible rejection is printed to the standard output.

    :ivar MFaStream mfaStr: :class:`MFaStream` object to check.
    :ivar MFaStrFilterProps fp: :class:`MFaStrFilterProps`; Properties
      of the filter to be applied.
    :ivar Boolean verb: Verbosity.

    :rtype: Boolean, True if all filters have been passed.

    """
    # Define start and stop codon regex strings
    startCodon = r"(atg)"
    stopCodons = r"(tag|taa|tga)"
    # Define regex pattern for indel
    indel = r'-'

    def check_all_aligned():
        if len(mfaStr.seqL) == fp.nSpecies:
            return True
        else:
            if verb is not None:
                print(mfaStr.seqL[0].name, "rejection;",
                      "Not all species are aligned.")
            return False

    def check_divergence():
        s0Data = mfaStr.seqL[0].data
        for s in mfaStr.seqL[1:]:
            sData = s.data
            counts = 0
            for i in range(len(s0Data)):
                if s0Data[i].lower() != sData[i].lower():
                    counts += 1
            if (counts / len(s0Data)) > fp.maxDiv:
                if verb is not None:
                    print(mfaStr.seqL[0].name, "rejection;",
                          "Sequences are too diverged.")
                    return False
        return True

    def check_start_codons():
        pattern = r'^' + startCodon
        for s in mfaStr.seqL:
            (nEx, nExTot) = s.get_exon_nr()
            if nEx == 1:
                dataString = s.data
                m = re.search(pattern, dataString, re.I)
                if m is None:
                    if verb is not None:
                        print(s.name, "rejection;",
                              "Start codons are not preserved.")
                    return False
        return True

    def check_stop_codons():
        pattern = stopCodons + r'$'
        for s in mfaStr.seqL:
            (nEx, nExTot) = s.get_exon_nr()
            if nEx == nExTot:
                dataString = s.data
                m = re.search(pattern, dataString, re.I)
                if m is None:
                    if verb is not None:
                        print(s.name, "rejection;",
                              "Stop codons are not preserved.")
                    return False
        return True

    def check_frame_shifting_gaps():
        pattern = indel + r'+'
        for s in mfaStr.seqL:
            dataString = s.data
            i = re.finditer(pattern, dataString, re.I)
            for m in i:
                # A gap has been found.  Check for frame shift.
                if ((m.end() - m.start()) % 3) != 0:
                    if verb is not None:
                        print(s.name, "rejection;",
                              "Frame-shifting gap.")
                    return False
        return True

    def check_for_long_gaps():
        pattern = indel + r'{' + repr(fp.maxGapLength + 1) + r',}'
        for s in mfaStr.seqL:
            dataString = s.data
            m = re.search(pattern, dataString, re.I)
            if m is not None:
                if verb is not None:
                    print(s.name, "rejection;",
                          "A gap is too long.")
                return False
        return True

    def check_nonsense_codon():
        pattern = r'(' + stopCodons + r')' + r'(?!$)'
        for s in mfaStr.seqL:
            dataString = s.data
            m = re.search(pattern, dataString, re.I)
            if m is not None:
                # Stop codon pattern has been found.  Check if frame
                # is not shifted.
                inFr = s.get_in_frame()
                if (m.start() + inFr) % 3 is 0:
                    if verb is not None:
                        print(s.name, "rejection;",
                              "A nonsense codon has been found.")
                    return False
        return True

    def check_exon_length():
        dataStr = mfaStr.seqL[0].data
        if len(dataStr) < fp.minExonLen:
            if verb is not None:
                print(mfaStr.seqL[0].name, "rejection;",
                      "Exon is too short.")
            return False
        return True

    def check_exon_numbers():
        nExTotL = []
        for s in mfaStr.seqL:
            (nEx, nExTot) = s.get_exon_nr()
            nExTotL.append(nExTot)
        for i in range(len(nExTotL)):
            if nExTotL[0] != nExTotL[i]:
                if verb is not None:
                    print(mfaStr.seqL[0].name, "rejection;",
                          "Exon numbers do not match.")
                return False
        return True

    if fp.check_all_aligned:
        if not check_all_aligned():
            return False
    if fp.check_for_long_gaps:
        if not check_for_long_gaps():
            return False
    if (fp.nSpecies > 1) and fp.check_divergence:
        if not check_divergence():
            return False
    if fp.check_start_codons:
        if not check_start_codons():
            return False
    if fp.check_stop_codons:
        if not check_stop_codons():
            return False
    if fp.check_frame_shifting_gaps:
        if not check_frame_shifting_gaps():
            return False
    if fp.check_nonsense_codon:
        if not check_nonsense_codon():
            return False
    if fp.check_exon_length:
        if not check_exon_length():
            return False
    if fp.check_exon_numbers:
        if not check_exon_numbers():
            return False
    return True


class FaSeq():
    """Store sequence data retrieved from a fasta file.

    :ivar str name: Name of the `FaSeq` object.
    :ivar [Seq] seqL: List of :class:`Seq <libPoMo.seqbase.Seq>`
                      objects that store the actual sequence data.
    :ivar int nSepcies: Number of saved species / individuals /
                        chromosomes.

    """
    def __init__(self):
        self.name = ""
        self.seqL = []
        self.nSpecies = 0

    def print_info(self, maxB=50):
        """Print fasta sequence information.

        Print fasta sequence identifier, species names, the length of
        the sequence and a maximum of `maxB` bases (defaults to 50).

        """
        print("Sequence identifier:", self.name)
        for i in range(0, self.nSpecies):
            self.seqL[i].print_fa_header()
            print("Printing", maxB, "out of a total of",
                  self.seqL[i].dataLen, "bases.")
            print(self.seqL[i].data[0:maxB])
        return

    def get_seq_names(self):
        """Return a list with sequence names."""
        names = []
        for i in range(0, self.nSpecies):
            names.append(self.seqL[i].name)
        return names

    def get_seq_by_id(self, i):
        """Return sequence number `i` as `Seq` object."""
        seq = sb.Seq()
        seq = self.seqL[i]
        return seq

    def get_seq_base(self, seq, pos):
        """Return base at 1-based position `pos` in sequence with name
        `seq`."""
        names = self.get_seq_names()
        try:
            i = names.index(seq)
        except:
            raise sb.SequenceDataError("Sequence name not found.")
        if pos > self.seqL[i].dataLen:
            raise sb.SequenceDataError("Position out of range.")
        return self.seqL[i].get_base(pos)

    def get_distance(self):
        """Number of segregating bases.
        """
        count = 0
        for i in range(self.seqL[0].dataLen):
            base = self.seqL[0].get_base(i)
            for s in range(self.nSpecies):
                if base != self.seqL[s].get_base(i):
                    count += 1
                    break
        return count

def init_seq(faFileName, maxskip=50, name=None):
    """Open a fasta file and initialize an :class:`FaStream`.

    This function tries to open the given fasta file, checks if it is
    in fasta format and reads the first sequence.  It returns an
    :class:`FaStream` object. This object can later be used to parse
    the whole fasta file.

    Please close the associated file object with
    :func:`FaStream.close` when you don't need it anymore.

    :param str faFileName: File name of the fasta file.
    :param int maxskip: Only look *maxskip* lines for the start of a
                        sequence (defaults to 50).
    :param str name: Set the name of the sequence to *name*, otherwise
                     set it to the stripped filename.

    """
    flag = False
    faFile = sb.gz_open(faFileName)
    if name is None:
        name = sb.stripFName(faFileName)
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
    (nextHL, seq) = read_seq_from_fo(line, faFile)
    try:
        nextHL = nextHL.rstrip()
    except:
        pass
    faStr = FaStream(name, seq, nextHL, faFile)
    return faStr


def open_seq(faFileName, maxskip=50, name=None):
    """Open and read a fasta file.

    This function tries to open the given fasta file, checks if it is
    in fasta format and reads the sequence(s).  It returns an
    :class:`FaSeq` object that contains a list of species names, a
    list of the respective desriptions and a list with the sequences.

    :param str faFileName: Name of the fasta file.
    :param int maxskip: Only look *maxskip* lines for the start of a sequence
                        (defaults to 50).
    :param str name: Set the name of the sequence to *name* otherwise
                     set it to the stripped filename.

    """
    def test_sequence(faSequence):
        """Tests if sequences contain data."""
        l = faSequence.nSpecies
        names = []
        for i in range(l):
            names.append(faSequence.seqL[i].name)
            if faSequence.seqL[i].name == '' or faSequence.seqL[i].data == '':
                raise sb.SequenceDataError("Sequence name or data is missing.")
        if l > len(set(names)):
            raise sb.SequenceDataError("Sequence names are not unique.")
        return

    fastaSeq = FaSeq()

    flag = False
    faFile = sb.gz_open(faFileName)
    if name is not None:
        fastaSeq.name = name
    else:
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
    while line is not None:
        (nextLine, seq) = read_seq_from_fo(line, faFile)
        line = nextLine
        fastaSeq.seqL.append(seq)
        fastaSeq.nSpecies += 1
    faFile.close()
    test_sequence(fastaSeq)
    return fastaSeq


def save_as_vcf(faSeq, ref, VCFFileName):
    """Save the given :classL`FaSeq` in VCF format.

    In general, we want to convert a fasta file with various
    individuals with the help of a reference that contains one
    sequence to a VCF file that contains all the SNPs.  This can be
    done with this function.  Until now it is not possible to do this
    conversion for several chromosomes for each individual in one run.
    Still, the conversion can be done chromosome by chromosome.

    This function saves the SNPs of *faSeq*, a given :class:`FaSeq`
    (fasta sequence) object in VCF format to the file *VCFFileName*.
    The reference genome *ref*, to which *faSeq* is compared to, needs
    to be passed as a :class:`Seq <libPoMo.seqbase.Seq>` object.

    The function compares all sequences in *faSeq* to the sequence
    given in *ref*.  The names of the individuals in the saved VCF
    file will be the sequence names of the *faSeq* object.

    ::

      #CHROM = sequence name of the reference
      POS    = position relative to reference
      ID     = .
      REF    = base of reference
      ALT    = SNP (e.g. 'C' or 'G,T' if 2 different SNPs are present)
      QUAL   = .
      FILTER = .
      INFO   = .
      FORMAT = GT

    :param FaSeq faSeq: :class:`FaSeq` object to be converted.
    :param Seq ref: :class:`Seq <libPoMo.seqbase.Seq>` object of the
                    reference sequence.
    :param str VCFFileName: Name of the VCF output file.

    """
    def get_altBases_string(sAltBases):
        """Return ALT bases string from given `sAltBases`."""
        l = len(sAltBases)
        if l == 0:
            return ''
        string = str(sAltBases[0])
        if l > 1:
            for i in range(1, l):
                string += ',' + sAltBases[i]
        return string

    def get_indiv_string(indivData, altBases, sAltBases):
        """Return the string of the individual data.

        Return the string extracted from the indivudal data
        `indivData` with SNPs `altBases`. `sAltBases` is the string
        with the alternative bases.

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
                if not (indivData[i] in altBases):
                    string += '\t' + '0'
                else:
                    string += '\t' + str(sAltBases.index(indivData[i]) + 1)
        return string

    def get_vcf_line(chromName, pos,
                     refBase, altBaseString, indivString):
        """Print a VCF file line with given data to file `VCFFile`."""
        string = chromName + '\t'
        string += str(pos) + '\t'
        string += '.' + '\t'    # id
        string += refBase + '\t'
        string += altBaseString + '\t'
        string += '.' + '\t'    # qual
        string += '.' + '\t'    # filter
        string += '.' + '\t'    # info
        string += "GT" + '\t'   # format
        string += indivString
        return string

    if (not isinstance(faSeq, FaSeq)):
        raise sb.SequenceDataError("`faSeq` is not an FaSeq object.")
    if (not isinstance(ref, sb.Seq)):
        raise sb.SequenceDataError("`ref` is not a Seq object.")
    if faSeq.nSpecies == 0:
        raise sb.SequenceDataError("`faSeq` has no saved sequences.")
    for i in range(0, faSeq.nSpecies):
        if faSeq.seqL[i].dataLen != ref.dataLen:
            raise sb.SequenceDataError(
                "Sequence " + faSeq.seqL[i].name +
                " has different length than reference.")
    VCFFile = sb.gz_open(VCFFileName, mode='w')
    print(vcf.get_header_line_string(faSeq.get_seq_names()), file=VCFFile)
    # loop over bases
    refBase = ''
    for i in range(0, ref.dataLen):
        refBase = ref.data[i]
        altBases = set()
        indivData = []
        # loop over sequences in faSeq and check if there is a SNP
        for s in range(0, faSeq.nSpecies):
            indivData.append(faSeq.seqL[s].data[i])
            if faSeq.seqL[s].data[i] != refBase:
                altBases.add(faSeq.seqL[s].data[i])
        sAltBases = sorted(altBases)
        altBaseString = get_altBases_string(sAltBases)
        indivString = get_indiv_string(indivData, altBases, sAltBases)
        if altBases != set():
            print(
                get_vcf_line(ref.name, i+1, refBase,
                             altBaseString, indivString),
                file=VCFFile)
    VCFFile.close()
    return
