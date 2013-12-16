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

Objects
-------
Classes:
  - :class:`CFWriter`, write a counts format file

Exception Classes:
  - :class:`CountsFormatWriterError`

Functions:
  - :func:`save_as_cf()`, deprecated; save given sequences to a counts
    format file

----

"""

__docformat__ = 'restructuredtext'

import numpy as np
import gzip
import pysam as ps

import libPoMo.seqbase as sb
import libPoMo.vcf as vcf
import libPoMo.fasta as fa

dna = {'a': 0, 'c': 1, 'g': 2, 't': 3}


class CountsFormatWriterError(sb.SequenceDataError):
    """General `CFWriter` object error."""
    pass


class CFWriter():
    """Write a counts format file.

    Save information that is needed to write a CF file and use this
    information to write a CF file.  Initialize with a reference fasta
    file name, a list of vcf file names and an output file name::

      CFWriter("refFasta", [vcfFileNames], "output")

    Write a header line to output::

       self.write_HLn()

    Write lines in counts format from 1-based positions *start* to
    *end* on chromosome *chrom* to output::

       rg = sb.Region("chrom", start, end)
       self.write_Rn(rg)

    Remember to close the attached file objectsL with :func:`close`.

    :param str refFileName: Name of reference fasta file.
    :param [str] vcfFileNameL: List with names of vcf files.
    :param str outFileName: Output file name.
    :param int verb: Optional; verbosity.
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
    :ivar FaStr refFaStr: :class:`FaStream <libPoMo.fasta.FaStream>`
        object of the reference genome. This might be changed to an
        *MFaStream* object in the future.
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
    :ivar indM: Matrix with individuals from vcf files. *self.indM[i]*
        is the list of individuals found in *self.vcfL[i]*.
    :ivar [int] nIndL: List with number of individuals in
        *self.vcfL[i]*.
    :ivar assM: Assignment matrix that connects the individuals from
        the vcf files to the correct *self.cd* index.  Cf. *self.cD*
    :ivar int nPop: Number of different populations in count format
        output file (e.g. number of populations).  Filled by
        *self.__init_assM()* during initialization.
    :ivar int ploidy: Ploidy of individuals in vcf files.  This has to
        be set manually to the correct value for non-diploids!
    :ivar char __splitCh: Character that is used to split the
        individual names.

    """
    def __init__(self, refFileName, vcfFileNameL, outFileName,
                 verb=None, mergeL=None, nameL=None):
        # Passed variables.
        self.refFN = refFileName
        self.vcfL = vcfFileNameL
        self.outFN = outFileName
        self.v = verb
        self.mL = mergeL
        self.nL = nameL
        # Variables that are filled during initialization.
        self.nV = len(self.vcfL)
        self.vcfTfL = []
        self.refFaStr = fa.init_seq(self.refFN)
        self.outFO = None
        self.cD = []
        self.chrom = None
        self.pos = None
        self.indM = []
        self.nIndL = []
        self.assM = []
        self.nPop = 0
        self.ploidy = 2

        self.__splitCh = '-'

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
        if self.outFN[-2:] == "gz":
            self.outFO = gzip.open(self.outFN, mode="wt")
        else:
            self.outFO = open(self.outFN, mode='w')

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
            l = [e.rsplit(self.__splitCh, maxsplit=1)[0]
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
                            self.__splitCh, maxsplit=1)[0])

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

        Fill *self.cF* with data from reference at chromosome *chrom*
        and position *pos*. Possible SNPs in *slef.vcfL* at this
        position are considered.

        :param [int] iL: List with vcf indices of the SNPs in *snpL*,
            must be sorted.

        :param [NucBase] snpL: List with :class:`NucBase
            <libPoMo.vcf.NucBase>` SNPs at this position. None, if
            there is no SNP.
        :raises: :class:`NotAValidRefBase
            <libPoMo.seqbase.NotAValidRefBase>`.

        Raises the excpetion :class:`NotAValidRefBae
        <libPoMo.seqbase.NotAValidRefBase>` if the reference base is
        not valid (e.g. N).

        """
        def get_refBase():
            """Get reference base on *chrom* at *pos*."""
            if self.chrom == self.refFaStr.seq.name:
                return self.refFaStr.seq.data[self.pos]
            else:
                raise sb.SequenceDataError("Chromosome name invalid.")
        self.__purge_cD()
        refBase = get_refBase()
        try:
            r = dna[refBase.lower()]
        except KeyError:
            raise sb.NotAValidRefBase()
        # If there are no SNPS, fill *self.cD* with data from reference.
        if iL is None:
            for i in range(self.nV):
                for sp in self.assM[i]:
                    self.cD[sp][r] += self.ploidy
        elif (snpL is not None) and (len(iL) == len(snpL)):
            # Else, only fill *self.cD* where the individual has no SNP
            for i in range(self.nV):
                if i not in iL:
                    for sp in self.assM[i]:
                        self.cD[sp][r] += self.ploidy
            # Now traverse the SNPs
            for sI in range(len(iL)):
                altBases = snpL[sI].get_alt_base_list()
                spData = snpL[sI].get_speciesData()
                vI = iL[sI]
                # Loop over individuals.
                for i in range(0, len(spData)):
                    # ipdb.set_trace()
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
        stringL = [self.chrom, str(self.pos + 1)]
        for data in self.cD:
            stringL.append(','.join(map(str, data)))
        return '\t'.join(stringL)

    def __get_HLn(self):
        """Return a string containing the headerline in counts format."""
        strL = ["CHROM", "POS"]
        strL.extend(self.nL)
        return '\t'.join(strL)

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
        snpsG = self.__snp(rg)
        try:
            (nI, nSNP) = next(snpsG)
        except StopIteration:
            nI = None
            nSNP = None

        for pos in range(rg.start, rg.end + 1):
            snpL = None
            iL = None
            while (nI is not None) or (nSNP is not None):
                if nSNP.pos - 1 == pos:
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
            self.pos = pos
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


###############################################################################
## Deprecated stuff follows.


def get_cf_headerline(species):
    """Deprecated. Return a string containing the headerline in counts
    format."""
    stringL = [["CHROM", "POS"]]
    stringL.extend(species)
    return '\t'.join([e for lst in stringL for e in lst])


def collapse(speciesL, merge=False, name=None):
    """Deprecated. Collapse the species names.

    Collapse the species names using the naming rules given in
    save_as_countsformat. Returns a dictionary with collapsed
    species names as keys and an assignment list.

    :param Boolean merge: if set to True, all species/individuals will be
      collapsed to a single one with name *name* (if *name* is given).

    """
    l = len(speciesL)
    if merge is True:
        if name is None:
            name = speciesL[0].rsplit('_', maxsplit=1)[0]
        assList = [name for i in range(l)]
        collapsedSp = [name]
    else:
        assList = [speciesL[i].rsplit('_', maxsplit=1)[0]
                   for i in range(l)]
        collapsedSp = sorted(set(assList))
    return (collapsedSp, assList)


def find_next_SNP_pos(nSNPChromL, nSNPPosL, ref):
    """Deprecated. Find the position of next SNP.

    Return the position of the next SNP in *ref* and a list of
    indices of these SNPs in *nSNPChromL* and *nSNPPosL*. If no
    next SNP is found on this sequence (this might happen if all
    next SNPs are on the next chromosome) a ValueError is raised.

    :param nSNPChromL: List with chromosome names of the next SNPs of
                       the *vcfStrL*.
    :param nSNPPosL: List with positions of the next SNPs of the
                     *vcfStrL*.
    :param ref: Seq object of the sequence of the reference.

    """
    ind = -1
    indL = []
    posL = []
    for i in range(len(nSNPChromL)):
        ind += 1
        if ref.name == nSNPChromL[i]:
            posL.append(nSNPPosL[i])
            indL.append(ind)
    if len(posL) > 0:
        npPosL = np.array(posL)
        minIndL = np.where(npPosL == npPosL.min())[0].tolist()
        fIndL = [indL[i] for i in minIndL]
        return (posL[minIndL[0]], fIndL)
    else:
        raise ValueError("No next SNP found in `ref`.")


def fill_species_dict_ref(spDi, assL, refBase, ploidy):
    """Deprecated. Fills the species dictionary if no SNP is present."""
    # reset species dictionary to 0 counts per base
    for key in spDi.keys():
        spDi[key] = [0, 0, 0, 0]
    # if no altBase is given, count species
    try:
        r = dna[refBase.lower()]
    except KeyError:
        # Base is not valid (reference is masked on this position).
        return None
    for sp in assL:
        spDi[sp][r] += ploidy
    return


def fill_species_dict_alt(spDi, assL, refBase, altBases, spData):
    """Deprecated. Fill the species dictionary when a SNP is present.

    Raise *SequenceDataError*, if *spDi* could not be filled.
    """
    # reset species dictionary to 0 counts per base
    for key in spDi.keys():
        spDi[key] = [0, 0, 0, 0]
    if (altBases is not None) \
       and (spData is not None):
        bases = [refBase.lower()]
        for b in altBases:
            bases.append(b.lower())
        # loop over species
        l = len(spData[0])
        for i in range(0, len(spData)):
            # loop over chromatides (e.g. diploid)
            for d in range(0, l):
                if spData[i][d] is not None:
                    bI = dna[bases[spData[i][d]]]
                    spDi[assL[i]][bI] += 1
        return
    raise sb.SequenceDataError("Could not fill species dictionary.")


def get_counts_line(chrom, pos, spDi, speciesL):
    """Deprecated. Return line string in counts format."""
    stringL = [chrom, str(pos+1)]
    for i in range(len(spDi)):
        for j in range(len(speciesL[i])):
            stringL.append(','.join(map(str, spDi[i][speciesL[i][j]])))
    return '\t'.join(stringL)


def save_as_cf(vcfStrL, refFaStr, CFFileName, verb=False,
               mergeL=None, nameL=None):
    """Deprecated. Save the given sequence in counts format.

    This function saves the SNPs from *vcfStrL*, a given list of
    *VCFStream* (variant call format sequence stream) objects in
    counts format to the file *CFFileName*.  The reference genome
    *refFaStr*, to which *VCFSeqStr* is compared to, needs to be
    passed as an *FaStream* object.

    The name of all streams in *vcfStrL* should be the same as the
    name of *faRef*.  The names of the sequences in the reference
    should be the names of the chromosomes found in the *vcfStr*
    object, otherwise we do not know where to compare the sequences
    to.  They must also be in the same order!

    Individuals with the same name and suffix "_n", where n is a
    number, will be saved in one column without the suffix.

    :param vcfStrL: List with VCF Streams containing SNPs.
    :param refFaStr: Reference fasta sequence stream object.
    :param CFFileName: Name of output file (counts format).
    :param verb: If *verb* is set to True, additional information is
                 printed to the output file.
    :param mergeL: A list of truth values. If *mergeL[i]* is True, all
                  individuals of *vcfStrL[i]* are treated as one species
                  independent of their name. The respective counts are
                  summed up.  If *nameL[i]* is given, the name of the
                  summed sequence will be *nameL[i]*. If not, the name of
                  the first individual in *vcfStrL[i]* will be used.
    :param nameL: A list of names. Cf. *mergeL*.

    """
    lVcfStrL = len(vcfStrL)
    if (not isinstance(refFaStr, fa.FaStream)):
        raise sb.SequenceDataError("`faRef` is not an FaStream object.")
    for i in range(lVcfStrL):
        if (not isinstance(vcfStrL[i], vcf.VCFStream)):
            raise sb.SequenceDataError("`vcfStr` " + vcfStrL[i] +
                                       " is not a VCFStream object.")
        if vcfStrL[i].name != refFaStr.name:
            raise sb.SequenceDataError("VCF sequence name " + vcfStrL[i].name +
                                       " and reference name " + refFaStr.name +
                                       " do not match.")
        if vcfStrL[i].nSpecies == 0:
            raise sb.SequenceDataError("`VCFSeq` has no saved data.")

    # Set mergeL and nameL if they are not given.
    if mergeL is None:
        mergeL = []
        for i in range(lVcfStrL):
            mergeL.append(False)
    if nameL is None:
        nameL = []
        for i in range(lVcfStrL):
            nameL.append(None)

    # ploidy = number of haploid sets given in the VCF file
    # allSpeciesL = list of vcfStr species
    # assL = assignment list; allSpeciesL[i]=Wolf_n => assL[i]=Wolf
    # collSpeciesL = collapsed species list (_n removed)
    # spDiL = list of dictionaries with speciesL as keys and list of counts
    # hence, spDi[assL[i]] is the list of counts for Wolf
    ploidy = vcfStrL[0].base.set_ploidy()
    allSpeciesL = []
    collSpeciesL = []
    assL = []
    spDiL = []
    for i in range(lVcfStrL):
        allSpeciesL.extend(vcfStrL[i].speciesL)
        (collSpecies, ass) = collapse(vcfStrL[i].speciesL, mergeL[i], nameL[i])
        collSpeciesL.append(collSpecies)
        assL.append(ass)
        spDiL.append(dict.fromkeys(collSpecies, None))

    if CFFileName[-2:] == "gz":
        fo = gzip.open(CFFileName, mode='wt')
    else:
        fo = open(CFFileName, mode='w')

    if verb is True:
        print("#Sequence name =", refFaStr.name, file=fo)
    print(get_cf_headerline(collSpeciesL), file=fo)
    # get chromosomes and positions of next SNPs from the vcfStrL
    nSNPChromL = []
    nSNPPosL = []
    for i in range(lVcfStrL):
        nSNPChromL.append(vcfStrL[i].base.chrom)
        nSNPPosL.append(vcfStrL[i].base.pos - 1)
    # Loop over sequences (chromosomes) in refFaStr.
    while True:
        # Initialize the saved old position of the previous SNP
        # and the sequence of the reference genome.
        oldPos = 0
        ref = refFaStr.seq
        if verb is True:
            print("#Chromosome name =", ref.name, file=fo)
        # Find next SNP position and the corresponding indices
        try:
            (nSNPPos, nSNPIndL) = find_next_SNP_pos(nSNPChromL,
                                                    nSNPPosL, ref)
        except ValueError:
            nSNPPos = ref.dataLen
            nSNPIndL = None
        # Loop over positions in sequence *ref* of the reference genome.
        while True:
            # Loop from previous SNP to one position before this one.
            for pos in range(oldPos, nSNPPos):
                try:
                    dna[ref.data[pos].lower()]
                except KeyError:
                    # Base is not valid (reference is masked on this position).
                    continue
                for i in range(lVcfStrL):
                    fill_species_dict_ref(spDiL[i], assL[i],
                                          ref.data[pos], ploidy)
                print(get_counts_line(
                    ref.name, pos, spDiL, collSpeciesL), file=fo)
            # Process the SNP at position pos on chromosome ref.
            if nSNPIndL is None:
                # There is no more SNP on this chromosome.
                break
            # Set spDiL to reference
            for i in range(lVcfStrL):
                fill_species_dict_ref(spDiL[i], assL[i],
                                      ref.data[nSNPPos], ploidy)
            # Loop over SNPs at position nSNPPos.
            for s in range(len(nSNPIndL)):
                altBases = vcfStrL[nSNPIndL[s]].base.get_alt_base_list()
                spData = vcfStrL[nSNPIndL[s]].base.get_speciesData()
                fill_species_dict_alt(spDiL[nSNPIndL[s]], assL[nSNPIndL[s]],
                                      ref.data[nSNPPos], altBases,
                                      spData)
                # Read next SNP.
                try:
                    vcfStrL[nSNPIndL[s]].read_next_base()
                    nSNPChromL[nSNPIndL[s]]\
                        = vcfStrL[nSNPIndL[s]].base.chrom
                    nSNPPosL[nSNPIndL[s]] = vcfStrL[nSNPIndL[s]].base.pos-1
                except ValueError:
                    # VCFStream ends here.
                    nSNPChromL[nSNPIndL[s]] = None
                    nSNPPosL[nSNPIndL[s]] = None
            # Save old SNP position.
            oldPos = nSNPPos + 1
            print(get_counts_line(
                ref.name, nSNPPos, spDiL, collSpeciesL), file=fo)
            # Find next SNP position and the corresponding indices
            try:
                (nSNPPos, nSNPIndL) = find_next_SNP_pos(nSNPChromL,
                                                        nSNPPosL, ref)
            except ValueError:
                nSNPPos = ref.dataLen
                nSNPIndL = None
        # Read next sequence in reference and break if none is found.
        if refFaStr.read_next_seq() is None:
            break
    fo.close()
    return
