#!/usr/bin/env python

"""libPoMo.gp
=============

A short library to read gene prediction files (.gp).

Two sample lines:

  ensGene.ENST00000391546.1     AHZZ01158546   -              9167           9392           9167           9392           1              9167           9392
  knownGene.uc010gkw.1.1.inc    chr3           +              32554111       32570175       32554111       32570175       4              32554111,32557251,32560104,32570077          32554339,32557360,32560206,325701

At the moment, I am not quite sure why the start and end positions are
given two times.

"""

import libPoMo.fasta as fasta
import libPoMo.seqbase as sb

class Exon():
    """An exon with start and end position."""
    def __init__(self, start, end):
        self.start = start
        self.end = end


class Gene():
    """A gene stored in a GP file line."""
    def __init__(self, ln):
        lnl = ln.split()
        self.name = lnl[0]
        self.chrom = lnl[1]
        self.orientation = lnl[2]
        self.start = int(lnl[3])
        self.end = int(lnl[4])
        # lnl[5]?
        # lnl[6]?
        self.nr_exons = int(lnl[7])
        starts_str = lnl[8]
        ends_str = lnl[9]
        exons_starts = starts_str.split(",")
        exons_ends = ends_str.split(",")
        self.exons = [Exon(int(exons_starts[i]), int(exons_ends[i]))
                      for i in range(self.nr_exons)]

        if self.orientation == "-":
            self.is_rc = True
        elif self.orientation == "+":
            self.is_rc = False
        else:
            raise ValueError("Invalid orientation.")

        # TODO: Correctly handle in and out frame for reverse
        # complemented genes.
        if self.is_rc:
            self.exons.reverse()

        # self.regions = self.__get_regions__()

    # def __get_regions__(self):
    #     if self.orientation not in ["+", "-"]:
    #         raise sb.SequenceDataError("Direction character is missing.")
    #     rgs = []
    #     for i in range(self.nr_exons):
    #         rgs.append(sb.Region(self.chrom, self.exons[i].start,
    #                              self.exons[i].end,
    #                              orientation=self.orientation))
    #     return rgs


class GPStream():
    """Read a GP file line per line.

    In order to interpret the data, a reference fasta file name is
    needed.

    """

    def __init__(self, fn, rf_fn):
        self.fn = fn
        self.fo = open(fn, mode="r")
        self.rf = fasta.open_seq(rf_fn)
        self.read_next_gene()

    def read_next_gene(self):
        """Get next gene."""
        ln = self.fo.readline()
        if ln != "":
            self.gene = Gene(ln)
            inframe = 0
            self.seqs = []
            for e in self.gene.exons:
                self.seqs.append(convert_exon_to_seq(self.gene, e,
                                                     inframe, self.rf))
                inframe = self.seqs[-1].get_out_frame()
        else:
            raise ValueError("End of CFStream.")

    def close(self):
        self.fo.close()


def convert_exon_to_seq(gene, exon, inframe, rf):
    """Extract the sequence information of gene from a reference.

    The `Gene()` only contains the positional information.  To get a
    valid sequence, information from a reference genome `rf`
    `fasta.FaSeq()` is needed.  The `inframe` should give the position
    in the triplet of the first base of the exon (0, 1 or 2).

    For the baboon GP file, the `inframe` is 0 for the first exon.
    The next exon always continues with the outframe of the previous
    one (i.e., the exons are directly connected).

    A `seqbase.Seq()` object is returned.

    """
    seq = sb.Seq()
    seq.name = gene.name
    # TODO: How are exon start and end defined?  Here: Indexing starts
    # with 0 but the exon end is not included in the sequence, this
    # seems to be most coherent with start codons in the data.
    seq.data = rf.seqD[gene.chrom].data[exon.start:exon.end]
    seq.dataLen = exon.end - exon.start
    # import pdb; pdb.set_trace()
    # Get orientation.
    # if gene.orientation == "+":
    #     rc = False
    # elif gene.orientation == "-":
    #     rc = True
    # else:
    #     raise ValueError("Invalid orientation.")
    seq.rc = False

    # in_frame = (exon.start - gene.start) % 3
    in_frame = inframe
    out_frame = (seq.dataLen + in_frame) % 3
    # FIXME: For the baboon data the chromosome names in the VCF file
    # do not include 'chr'.  That's why I remove this here, but this
    # should be made more general in the future.
    if gene.chrom[:3] == "chr":
        chrom = gene.chrom[3:]
        # print("Fixme; chromosome name changed for baboon data.")
    else:
        chrom = gene.chrom
    # Also set appropriate description so that the helper functions of
    # `sb.Seq()` can be used.
    seq.descr = str(seq.dataLen) + " " + str(in_frame) + " "
    seq.descr += str(out_frame) + " " + chrom + ":"
    seq.descr += str(exon.start+1) + "-" + str(exon.end) + gene.orientation
    # if gene.orientation == "-":
    #     seq.rev_comp(change_sequence_only=True)
    seq.set_gene_is_rc_from_descr()
    return seq
