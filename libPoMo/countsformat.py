#!/usr/bin/env python

"""libPomo.countsformat
----------------------------------------------------------------------

This model provides function to read, write and access files that are
in counts format. This file format is used by PoMo and lists the base
counts for every position. It contains:

1 headerline with tab separated sequence names
N lines with counts of A, C, G and T bases at position n
#TODO positional information

Sheep   \t BlackSheep \t RedSheep \t Wolf    \t RedWolf
0,0,1,0 \t 0,0,1,0    \t 0,0,1,0  \t 0,0,5,0 \t 0,0,0,1
0,0,0,1 \t 0,0,0,1    \t 0,0,0,1  \t 0,0,0,5 \t 0,0,0,1
.
.
.

"""

import libPoMo.seqbase as sb
import libPoMo.vcf as vcf
import libPoMo.fasta as fa


def get_counts_format_headerline(species):
    """Returns a string containing the headerline in counts format."""
    string = species[0]
    l = len(species)
    if l > 1:
        for i in range(1, l):
            string += '\t' + species[i]
    return string


def save_as_countsformat(VCFSeq, faRef, CFFileName):
    """Saves the given sequence in counts format.

    This function saves the SNPs of `VCFSeq`, a given VCFSeq (variant
    call format sequence) object in counts format to the file
    `CFFileName`.  The reference genome `ref`, to which `faSeq` is
    compared to, needs to be passed as an FaSeq object.

    The name of VCFSeq should be the same as the name of faRef.  The
    names of the sequences should be the names of the chromosomes
    found in the VCFSeq object, otherwise we do not know where to
    compare the sequences to.

    Individuals with the same name and suffix "_n", where n is a
    number, will be saved in one column without the suffix.

    """
    
    dna = {'a': 0, 'c': 1, 'g': 2, 't': 3}

    def collapse(speciesL):
        """Collapse the species names.
        
        Collapse the species names using the naming rules given in
        save_as_countsformat. Returns a dictionary with collapsed
        species names as keys and an assignment list.

        """
        l = len(speciesL)
        assList = []
        for i in range(0, l):
            assList.append(speciesL[i].rsplit('_', maxsplit=1)[0])
        collapsedSp = sorted(set(assList))
        return (collapsedSp, assList)

    def fill_species_dict(spDi, assList, refBase, altBases=None, spData=None):
        """Fills the species dictionary."""
        # reset species dictionary to 0 counts per base
        for key in spDi.keys():
            spDi[key] = [0, 0, 0, 0]
        r = dna[refBase.lower()]
        # if no altBase is given, count species
        if altBases is None:
            for sp in assList:
                spDi[sp][r] += 1
            return
        if (altBases is not None) \
           and (spData is not None):
            bases = [refBase.lower()]
            for b in altBases:
                bases.append(b.lower())
            # loop over species
            for i in range(0, len(spData)):
                bI = dna[bases[spData[i]]]
                spDi[assList[i]][bI] += 1
            return
        raise sb.SequenceDataError("Could not fill species dictionary.")
        return

    def get_counts_line(spDi, speciesL):
        """Returns line in counts format."""
        string = ','.join(map(str, spDi[speciesL[0]]))
        l = len(speciesL)
        if l > 1:
            for i in range(1, l):
                string += '\t' + ','.join(map(str, spDi[speciesL[i]]))
        return string

    if (not isinstance(VCFSeq, vcf.VCFSeq)):
        raise sb.SequenceDataError("`VCFSeq` is not a VCFSeq object.")
    if (not isinstance(faRef, fa.FaSeq)):
        raise sb.SequenceDataError("`faRef` is not an FaSeq object.")
    if VCFSeq.name != faRef.name:
        print(VCFSeq.name)
        print(faRef.name)
        raise sb.SequenceDataError("Sequence names do not match")
    if VCFSeq.nBases == 0:
        raise sb.SequenceDataError("`VCFSeq` has no saved bases.")

    allSpeciesL = VCFSeq.speciesL
    # speciesL = sorted list of unique species names
    # assList = assignment list; allSpeciesL[i]=Wolf_n => assList[i]=Wolf
    # spDi = dictionary with speciesL as keys and list of counts
    # hence, spDi[assList[i]] is the list of counts for Wolf
    (speciesL, assList) = collapse(allSpeciesL)
    spDi = dict.fromkeys(speciesL, None)
    i = 0                       # count SNPs in VCFSeq
    j = 0                       # count position in ref
    pos = -1                    # save SNP position
    oldPos = 0                  # save previous SNP position
    # TODO FIX
    ref = faRef.get_seq_by_id(0)

    # TODO TODO implement faRef so that I can merge with more than one
    # chromosome

    with open(CFFileName, 'w') as fo:
        print("#Sequence name =", faRef.name, file=fo)
        print(get_counts_format_headerline(speciesL), file=fo)
        # loop over SNPs
        for i in range(0, VCFSeq.nBases):
            oldPos = pos + 1
            pos = VCFSeq.baseL[i].pos
            # loop from previous SNP to one position before this one
            for j in range(oldPos, pos):
                fill_species_dict(spDi, assList, ref.data[j])
                print(get_counts_line(spDi, speciesL), file=fo)
            # process the SNP at pos
            altBases = VCFSeq.baseL[i].get_alt_base_list()
            spData = VCFSeq.baseL[i].get_speciesData()
            fill_species_dict(spDi, assList, ref.data[pos], altBases, spData)
            print(get_counts_line(spDi, speciesL), file=fo)
        # finish until the end of file
        if ref.dataLen > pos+1:
            for i in range(pos+1, ref.dataLen):
                fill_species_dict(spDi, assList, ref.data[i])
                print(get_counts_line(spDi, speciesL), file=fo)









