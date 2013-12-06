#!/usr/bin/env python

"""libPomo.countsformat
----------------------------------------------------------------------

This model provides functions to read, write and access files that are
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


def save_as_cf(vcfStr, refFaStr, CFFileName, verb=False,
               add=False, name=None, diploid=False):
    """Save the given sequence in counts format.

    This function saves the SNPs of `vcfStr`, a given VCFStream
    (variant call format sequence stream) object in counts format to
    the file `CFFileName`.  The reference genome `refFaStr`, to which
    `VCFSeqStr` is compared to, needs to be passed as an FaStream
    object.

    The name of `vcfStr` should be the same as the name of `faRef`.
    The names of the sequences in the reference should be the names of
    the chromosomes found in the `vcfStr` object, otherwise we do not
    know where to compare the sequences to.  They must also be in the
    same order!

    Individuals with the same name and suffix "_n", where n is a
    number, will be saved in one column without the suffix.

    `verb` - If `verb` is True, additional information is printed to
    the output file.

    `add` - If `add` is True, all individuals are treated as one
    species independent of their name and counts are summed up.  If
    `name` is given, the name of the summed sequence will be
    `name`. If not, the name of the first individual will be used.

    `diploid` - set to true if vcfStr contains diploid data ("1/2")

    """

    dna = {'a': 0, 'c': 1, 'g': 2, 't': 3}

    def get_cf_headerline(species):
        """Returns a string containing the headerline in counts format."""
        return '\t'.join(species)

    def collapse(speciesL, add, name):
        """Collapse the species names.

        Collapse the species names using the naming rules given in
        save_as_countsformat. Returns a dictionary with collapsed
        species names as keys and an assignment list.

        `add`: if set to true, all species/individuals will be
        collapsed to a single one with name `name` (if `name` is
        given).

        """
        l = len(speciesL)
        if add is True:
            if name is None:
                name = speciesL[0].rsplit('_', maxsplit=1)[0]
            assList = [name for i in range(l)]
            collapsedSp = [name]
        else:
            assList = [speciesL[i].rsplit('_', maxsplit=1)[0]
                       for i in range(l)]
            collapsedSp = sorted(set(assList))
        return (collapsedSp, assList)

    def fill_species_dict(spDi, assList, refBase, altBases=None, spData=None):
        """Fills the species dictionary.

        Return True if all went well.
        Return None if a base is not valid.
        """
        # reset species dictionary to 0 counts per base
        for key in spDi.keys():
            spDi[key] = [0, 0, 0, 0]
        # if no altBase is given, count species
        if altBases is None:
            try:
                r = dna[refBase.lower()]
            except KeyError:
                # Base is not valid (reference is masked on this position).
                return None
            for sp in assList:
                spDi[sp][r] += 1
            return True
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
                        spDi[assList[i]][bI] += 1
            return True
        # raise sb.SequenceDataError("Could not fill species dictionary.")
        return None

    def get_counts_line(spDi, speciesL):
        """Returns line string in counts format."""
        string = ','.join(map(str, spDi[speciesL[0]]))
        l = len(speciesL)
        if l > 1:
            for i in range(1, l):
                string += '\t' + ','.join(map(str, spDi[speciesL[i]]))
        return string

    if (not isinstance(vcfStr, vcf.VCFStream)):
        raise sb.SequenceDataError("`vcfStr` is not a VCFStream object.")
    if (not isinstance(refFaStr, fa.FaStream)):
        raise sb.SequenceDataError("`faRef` is not an FaStream object.")
    if vcfStr.name != refFaStr.name:
        raise sb.SequenceDataError("VCF sequence name " + vcfStr.name +
                                   " and reference name " + refFaStr.name +
                                   " do not match.")
    if vcfStr.nSpecies == 0:
        raise sb.SequenceDataError("`VCFSeq` has no saved data.")

    allSpeciesL = vcfStr.speciesL
    # speciesL = sorted list of unique species names
    # assList = assignment list; allSpeciesL[i]=Wolf_n => assList[i]=Wolf
    # spDi = dictionary with speciesL as keys and list of counts
    # hence, spDi[assList[i]] is the list of counts for Wolf
    (collSpeciesL, assList) = collapse(allSpeciesL, add, name)
    spDi = dict.fromkeys(collSpeciesL, None)
    # The following
    j = 0                       # count position in ref

    with open(CFFileName, 'w') as fo:
        if verb is True:
            print("#Sequence name =", refFaStr.name, file=fo)
        print(get_cf_headerline(collSpeciesL), file=fo)
        # Loop over chromosomes in refFaStr.
        while True:
            pos = -1                    # initialize current SNP position
            oldPos = 0                  # initialize previous SNP position
            ref = refFaStr.seq
            if verb is True:
                print("#Chromosome name =", ref.name, file=fo)
            # Loop over SNPs.
            while ref.name == vcfStr.base.chrom:
                oldPos = pos + 1
                pos = vcfStr.base.pos - 1
                # Loop from previous SNP to one position before this one.
                for j in range(oldPos, pos):
                    if fill_species_dict(spDi, assList,
                                         ref.data[j]) is True:
                        print(get_counts_line(spDi, collSpeciesL), file=fo)
                # Process the SNP at position pos on chromosome ref.
                altBases = vcfStr.base.get_alt_base_list()
                spData = vcfStr.base.get_speciesData(diploid)
                if fill_species_dict(spDi, assList,
                                     ref.data[pos], altBases,
                                     spData) is True:
                    print(get_counts_line(spDi, collSpeciesL), file=fo)
                if vcfStr.read_next_base() is None:
                    break
            # Finish until the end of the chromosome.
            if ref.dataLen > pos+1:
                for j in range(pos+1, ref.dataLen):
                    if fill_species_dict(spDi, assList,
                                         ref.data[j]) is True:
                        print(get_counts_line(spDi, collSpeciesL), file=fo)
            # Read next sequence in reference and break if none is found.
            if refFaStr.read_next_seq() is None:
                break
