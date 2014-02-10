#!/usr/bin/env python

"""libPomo.cf
=============

This model provides functions to read, write and access files that are
in counts format.

The Counts Format
-----------------
This file format is used by PoMo and lists the base
counts for every position.

It contains:
  - 1 headerline with tab separated sequence names
  - N lines with counts of A, C, G and T bases at position n

::

  CHROM \t Pos   \t Sheep   \t BlackSheep \t RedSheep \t Wolf    \t RedWolf
  chrA  \t s     \t 0,0,1,0 \t 0,0,1,0    \t 0,0,1,0  \t 0,0,5,0 \t 0,0,0,1
  chrA  \t s + 1 \t 0,0,0,1 \t 0,0,0,1    \t 0,0,0,1  \t 0,0,0,5 \t 0,0,0,1
  .
  .
  .
  chrF  \t 8373  \t 0,0,0,1 \t 1,0,0,0    \t 0,1,0,0  \t 0,1,4,0 \t 0,0,1,0
  .
  .
  .
  chrE  \t end   \t 0,0,0,1 \t 0,1,0,0    \t 0,1,0,0  \t 0,5,0,0 \t 0,0,1,0

Convert to Counts Format
------------------------

To convert a fasta reference file with SNP information from a variant
call format (VCF) to counts format use the :class:`CFWriter`. If you
want to convert a multiple alignment fasta file, use the
:class:`CFWriter` together with the convenience function
:func:`write_cf_from_MFaStream()`.

Tabix index files need to be provided for all VCF files. They can be
created from the terminal with $(tabix -p vcf "vcf-file.vcf.gz") if
tabix is installed.

A code example is::

  import import_libPoMo
  import libPoMo.fasta as fa
  import libPoMo.cf as cf

  vcfFL = ["/path/to/vcf/file1", "/path/to/vcf/file2", "..."]

  cfw = cf.CFWriter(vcfFL, "name-of-outfile")
  mFaStr = fa.MFaStream("/path/to/fasta/reference")

  cfw.write_HLn()
  cf.write_cf_from_MFaStream(mFaStr, cfw)

  cfw.close()

Objects
-------
Classes:
  - :class:`CFWriter`, write a counts format file

Exception Classes:
  - :class:`CountsFormatWriterError`

Functions:
  - :func:`write_cf_from_MFaStream()`, write counts file using the
    given MFaStream and CFWriter

----

"""

__docformat__ = 'restructuredtext'

import pysam as ps

import libPoMo.seqbase as sb
import libPoMo.vcf as vcf

dna = {'a': 0, 'c': 1, 'g': 2, 't': 3}


class CountsFormatWriterError(sb.SequenceDataError):
    """General `CFWriter` object error."""
    pass


class CFWriter():
    """Write a counts format file.

    Save information that is needed to write a CF file and use this
    information to write a CF file.  Initialize with a list of vcf
    file names and an output file name::

      CFWriter([vcfFileNames], "output")

    Tabix index files need to be provided for all VCF files. They can
    be created from the terminal with $(tabix -p vcf
    "vcf-file.vcf.gz") if tabix is installed.

    Before the count file can be written, a reference sequence has to
    be specified.  A single reference sequence can be set with
    :func:`set_seq`.

    Write a header line to output::

       self.write_HLn()

    Write lines in counts format from 1-based positions *start* to
    *end* on chromosome *chrom* to output::

       rg = sb.Region("chrom", start, end)
       self.write_Rn(rg)

    If you want to compare the SNPs of the VCF files to a multiple
    alingment fasta stream (:class:`MFaStream
    <libPoMo.fasta.MFaStream>`) consider the very convenient function
    :func:`write_cf_from_MFaStream`.

    To determine the different populations present in the VCF files,
    the names of the individuals will be cropped at a specific char
    that can be set at initialization (standard value = '-'). It is
    also possible to collapse all individuals of determined VCF files
    to a single population (cf. mergeL and nameL).

    Remember to close the attached file objectsL with :func:`close`.

    :param [str] vcfFileNameL: List with names of vcf files.
    :param str outFileName: Output file name.
    :param int verb: Optional; verbosity level.
    :param char splitChar: Optional; set the split character so that
      the individuals get sorted into the correct populations.
    :param [Boolean] mergeL: Optional; a list of truth values.  If
      *mL[i]* is True, all individuals of *self.vcfL[i]* are treated as
      one population orspecies independent of their name.  The
      respective counts are summed up.  If *self.nL[i]* is given, the
      name of the summed sequence will be *self.nL[i]*.  If not, the
      name of the first individual in *vcfL[i]* will be used.
    :param [str] nameL: Optional; a list of names. Cf. *self.mL*.

    :ivar str refFN: Name of reference fasta file.
    :ivar [str] vcfL: List with names of vcf files.
    :ivar str outFN: Output file name.
    :ivar int v: Verbosity.
    :ivar [Boolean] mL: A list of truth values.  If *mL[i]* is True,
        all individuals of *self.vcfL[i]* are treated as one
        population orspecies independent of their name.  The
        respective counts are summed up.  If *self.nL[i]* is given,
        the name of the summed sequence will be *self.nL[i]*.  If not,
        the name of the first individual in *vcfL[i]* will be used.
    :ivar [str] nL: A list of names. Cf. *self.mL*.
    :ivar int nV: Number of vcf files.
    :ivar [fo] vcfTfL: List with *pysam.Tabixfile* objects. Filled by
        *self.__init_vcfTfL()* during initialization.
    :ivar fo outFO: File object of the outfile. Filled by
        *self.__init_outFO()* during initialization.
    :ivar cD: List with allele or base counts. The alleles of
        individuals from the same population are summed up.  Hence,
        *self.cD[p]* gives the base counts of population *p* in the
        form: [0, 0, 0, 0].  Population *p`*does not need to be the
        one from *self.vcfL[p]* because several populations might be
        present in one vcf file.  *self.assM* connects the individual
        j from *self.vcfL[i]* such that *self.assM[i][j]* is *p*.
    :ivar str chrom: Name of the current chromosome. Set and updated
        by :func:`write_Rn`.
    :ivar int pos: Current position on chromosome. Set and updated by
        :func:`write_Rn`.
    :ivar int offset: Value that can be set with :func:`set_offset`,
                      if the reference sequence does not start at the
                      1-based position 1 but at the 1-based position
                      *offset*.
    :ivar indM: Matrix with individuals from vcf files. *self.indM[i]*
        is the list of individuals found in *self.vcfL[i]*.
    :ivar [int] nIndL: List with number of individuals in
        *self.vcfL[i]*.
    :ivar assM: Assignment matrix that connects the individuals from
        the vcf files to the correct *self.cd* index.  Cf. *self.cD*
    :ivar int nPop: Number of different populations in count format
        output file (e.g. number of populations).  Filled by
        *self.__init_assM()* during initialization.
    :ivar Seq refSeq: :class:`Seq <libPoMo.seqbase.Seq>` object of the
        reference Sequence. This has to be set with :class:`set_seq`.
    :ivar int ploidy: Ploidy of individuals in vcf files.  This has to
        be set manually to the correct value for non-diploids!
    :ivar char splitCh: Character that is used to split the
        individual names.
    :ivar Boolean __force: If set to true, skip name checks.

    """
    def __init__(self, vcfFileNameL, outFileName, verb=None,
                 splitChar='-', mergeL=None, nameL=None):
        # Passed variables.
        self.vcfL = vcfFileNameL
        self.outFN = outFileName
        self.v = verb
        self.mL = mergeL
        self.nL = nameL
        # Variables that are filled during initialization.
        self.nV = len(self.vcfL)
        self.vcfTfL = []
        self.outFO = None
        self.cD = []
        self.chrom = None
        self.pos = None
        self.offset = 0
        self.indM = []
        self.nIndL = []
        self.assM = []
        self.nPop = 0
        # Variables that have to be set manually.
        self.refSeq = None
        self.ploidy = 2
        self.splitCh = splitChar

        self.__force = False

        self.__init_vcfTfL()
        self.__init_outFO()
        self.__init_indM()
        self.__init_nIndL()
        self.__init_assM()
        self.__init_nL()
        self.__init_cD()

    def __init_vcfTfL(self):
        """Open vcf files given in *self.vcfL*.

        Tabix index files need to be provided. They can be created
        from the terminal with $(tabix -p vcf "vcf-file.vcf.gz"). The
        tabix file objects are stored in *self.vcfTfL*. They need to
        be closed with :func:`close`.

        """
        for fn in self.vcfL:
            self.vcfTfL.append(ps.Tabixfile(fn))

    def __init_outFO(self):
        """Open *self.outFN*.

        If the file name ends with "gz", the outfile will be
        compressed and is opened with gzip.open().

        """
        self.outFO = sb.gz_open(self.outFN, mode='w')

    def __init_indM(self):

        """Extract individuals from the vcf files."""
        # Get individuals from the vcf files.
        for tf in self.vcfTfL:
            for ln in tf.header:
                hLn = ln
            self.indM.append(
                vcf.get_indiv_from_field_header(hLn.decode("utf-8")))

    def __init_nIndL(self):
        """Count individuals in each vcf file."""
        for indL in self.indM:
            self.nIndL.append(len(indL))

    def __init_assM(self):
        """Fill assignment matrix *self.assM*."""
        def collapse_and_append(n, dN):
            """Collapse individual names of *self.vcfL[n]*.

            Appends the collapsed individual names of *self.vcfL[n]*
            to *self.assM*.

            :param int n: Index.
            :param int dN: Offset in assL.
            :rtype: int

            Returns new offset.

            """
            l = [e.rsplit(self.splitCh, maxsplit=1)[0]
                 for e in self.indM[n]]
            aL = []
            cL = [l[0]]
            ddN = 0
            for s in l:
                try:
                    index = cL.index(s)
                    aL.append(n+dN+index)
                except ValueError:
                    ddN += 1
                    cL.append(s)
                    aL.append(n+dN+ddN)
            self.assM.append(aL)
            return dN + ddN

        i = 0
        dI = 0
        if self.mL is None:
            for i in range(self.nV):
                dI = collapse_and_append(i, dI)
        elif len(self.mL) == self.nV:
            for i in range(self.nV):
                if self.mL[i] is True:
                    self.assM.append([i+dI]*self.nIndL[i+dI])
                elif self.mL[i] is False:
                    dI = collapse_and_append(i, dI)
                else:
                    raise CountsFormatWriterError("Merge list is not " +
                                                  "a list of boolean values.")
        else:
            raise CountsFormatWriterError("`mergeL` is not valid.")
        self.nPop = i + dI + 1

    def __init_nL(self):
        """Fill *self.nL*."""
        def append_to_nL(i):
            for j in range(len(self.assM[i])):
                try:
                    self.nL[self.assM[i][j]]
                except IndexError:
                    self.nL.append(
                        self.indM[i][j].rsplit(
                            self.splitCh, maxsplit=1)[0])

        if self.nL is None:
            self.nL = []
            for i in range(self.nV):
                append_to_nL(i)
        elif len(self.nL) != self.nPop:
            raise CountsFormatWriterError("`nameL` is not valid.")

    def __init_cD(self):
        """Initialize the list with counts data."""
        self.cD = [[0, 0, 0, 0] for i in range(self.nPop)]

    def __snp(self, rg):
        """Generate SNPs in region *rg* out of *self.vcfL*.

        Generator that returns the next SNP in region *rg*
        (cf. :class:`Region <libPoMo.seqbase.Region>`) as a :class:`NucBase`
        object.  To loop over all SNPs in region *rg*:

        >>> rg = sb.Region("chr1", 500000, 1000000)
        >>> for s in self.snp(rg):
        ....:   s.print_info()

        """
        snpL = []
        snpIterL = []
        for i in range(self.nV):
            snpIterL.append(self.vcfTfL[i].fetch(reference=rg.chrom,
                                                 start=rg.start, end=rg.end))
        for i in range(self.nV):
            try:
                snpL.append(vcf.get_nuc_base_from_line(next(snpIterL[i]),
                                                       ploidy=self.ploidy))
            except StopIteration:
                snpL.append(None)
        while True:
            if snpL == [None] * self.nV:
                raise StopIteration()
            for i in range(self.nV):
                if snpL[i] is not None:
                    minPos = snpL[i].pos
                    minI = i
                    break
            for j in range(i+1, self.nV):
                if snpL[j] is not None and \
                   snpL[j].pos < minPos:
                    minPos = snpL[j].pos
                    minI = j
            yield (minI, snpL[minI])
            try:
                snpL[minI] = vcf.get_nuc_base_from_line(next(snpIterL[minI]),
                                                        ploidy=self.ploidy)
            except StopIteration:
                snpL[minI] = None

    def __purge_cD(self):
        self.__init_cD()

    def __fill_cD(self, iL=None, snpL=None):
        """Fill *self.cF*.

        Fill *self.cF* with data from reference at chromosome
        *self.chrom* and position *self.pos*. Possible SNPs in
        *slef.vcfL* at this position are considered.

        :param [int] iL: List with vcf indices of the SNPs in *snpL*,
            must be sorted.
        :param [NucBase] snpL: List with :class:`NucBase
            <libPoMo.vcf.NucBase>` SNPs at this position. None, if
            there is no SNP.
        :raises: :class:`NotAValidRefBase
            <libPoMo.seqbase.NotAValidRefBase>`,
            :class:`SequenceDataError
            <libPoMo.seqbase.SequenceDataError>`

        :class:`NotAValidRefBae <libPoMo.seqbase.NotAValidRefBase>` is
        raised if the reference base is not valid (e.g. N).

        :class:`SequenceDataError <libPoMo.seqbase.SequenceDataError>`
        is raised if the chromosome names do not match.

        """
        if self.v is not None:
            if snpL is not None:
                print("Next SNP(s):")
                for s in snpL:
                    s.print_info()

        def get_refBase():
            """Get reference base on *chrom* at *pos*."""
            return self.refSeq.data[self.pos].lower()

        self.__purge_cD()
        refBase = get_refBase()
        try:
            r = dna[refBase]
        except KeyError:
            raise sb.NotAValidRefBase()
        # If there are no SNPS, fill *self.cD* with data from reference.
        if iL is None:
            for i in range(self.nV):
                for sp in self.assM[i]:
                    self.cD[sp][r] += self.ploidy
        elif (snpL is not None) and (len(iL) == len(snpL)):
            # Else, only fill *self.cD* where the individual has no SNP.
            for i in range(self.nV):
                if i not in iL:
                    for sp in self.assM[i]:
                        self.cD[sp][r] += self.ploidy
            # Now traverse the SNPs.
            for sI in range(len(iL)):
                # Check if the reference bases match.
                vcfRefBase = snpL[sI].get_ref_base().lower()
                if dna[vcfRefBase] != r:
                    print("Error at NucBase:")
                    snpL[sI].print_info()
                    print("The reference base at position", self.pos,
                          "on chromosome", self.chrom,
                          "is", refBase, end=".\n")
                    print("The reference base of the VCF file is",
                          vcfRefBase, end=".\n")
                    raise sb.SequenceDataError("Reference bases do not match.")
                altBases = snpL[sI].get_alt_base_list()
                spData = snpL[sI].get_speciesData()
                vI = iL[sI]
                # Loop over individuals.
                for i in range(0, len(spData)):
                # Loop over chromatides (e.g. diploid).
                    for d in range(0, self.ploidy):
                        if spData[i][d] is None:
                            pass
                        elif spData[i][d] == 0:
                            bI = r
                            self.cD[self.assM[vI][i]][bI] += 1
                        else:
                            bI = dna[altBases[spData[i][d]-1]]
                            self.cD[self.assM[vI][i]][bI] += 1
        else:
            raise sb.SequenceDataError("SNP information is not correct.")

    def __get_Ln(self):
        """Return string with a line in counts format. Positional information
        is written 1-based.

        """
        stringL = [self.chrom, str(self.pos + 1 + self.offset)]
        for data in self.cD:
            stringL.append(','.join(map(str, data)))
        return '\t'.join(stringL)

    def __get_HLn(self):
        """Return a string containing the headerline in counts format."""
        strL = ["CHROM", "POS"]
        strL.extend(self.nL)
        return '\t'.join(strL)

    def set_force(self, val):
        """Sets *self.__force* to *val*.

        :param Boolean val:

        """
        self.__force = val

    def set_seq(self, seq):
        "Set the reference sequence."""
        if (not isinstance(seq, sb.Seq)):
            raise sb.SequenceDataError("`seq` is not a Seq object.")
        self.refSeq = seq

    def set_ploidy(self, ploidy):
        """Set the ploidy.

        In VCF files, usually the bases of all copies of the same
        chromosomes are given and separated by '/' or '|'.  If the
        species is not diploid, this ploidy has to be set manually
        with this function.

        """
        self.ploidy = ploidy

    def set_offset(self, offset):
        """Set the offset of the sequence.

        :param int offset: Value that can be set, if the reference
                           sequence does not start at the 1-based
                           position 1 but at the 1-based position
                           *offset*.

        """
        self.offset = offset
        if self.v is not None:
            print('Offset in CFWriter:', self.offset)

    def __write_Ln(self):
        """Write a line in counts format to *self.outFN*."""
        print(self.__get_Ln(), file=self.outFO)

    def write_HLn(self):
        """Write the counts format header line to *self.outFN*."""
        print(self.__get_HLn(), file=self.outFO)

    def write_Rn(self, rg):
        """Write lines in counts format to *self.outFN*.

        :param Region rg: :class:`Region <libPoMo.seqbase.Region>`
                          object that determines the region that is
                          covered.

        """
        self.set_offset(rg.start)
        snpsG = self.__snp(rg)
        try:
            (nI, nSNP) = next(snpsG)
        except StopIteration:
            nI = None
            nSNP = None

        for rPos in range(rg.start, rg.end + 1):
            snpL = None
            iL = None
            while (nI is not None) or (nSNP is not None):
                if nSNP.pos - 1 == rPos:
                    if (snpL is None) and (iL is None):
                        snpL = []
                        iL = []
                    snpL.append(nSNP)
                    iL.append(nI)
                    try:
                        (nI, nSNP) = next(snpsG)
                    except StopIteration:
                        nI = None
                        nSNP = None
                else:
                    break
            self.chrom = rg.chrom
            self.pos = rPos - self.offset
            try:
                self.__fill_cD(iL, snpL)
            except sb.NotAValidRefBase:
                pass
            else:
                self.__write_Ln()

    def close(self):
        """Close fileobjects."""
        for tf in self.vcfTfL:
            tf.close()
        self.outFO.close()


def write_cf_from_MFaStream(refMFaStr, cfWr):
    """Write counts file using the given MFaStream and CFWriter.

    Write the counts format file using the first sequences of all
    alignments in the MFaStream.  The sequences are automatically
    reversed and complemented if this is needed (indicated in the
    header line).  This is very useful if you e.g. want to compare the
    VCF files to a CCDC alignment.

    :param FMaStream refMFaStr: The reference :class:`MFaStream
      <libPoMo.fasta.MFaStream>`.
    :param CFWriter cfWf: The :class:`CFWriter` object that contains
      the VCF files.

    """
    while True:
        refMFaStr.orient(firstOnly=True)
        rg = refMFaStr.seqL[0].get_region()
        cfWr.set_seq(refMFaStr.seqL[0])
        cfWr.write_Rn(rg)
        if refMFaStr.read_next_align() is None:
            break
