#!/usr/bin/env python

"""libPoMo.main
===============

This library contains functions that are used by PoMo.

"""

import argparse
import random
from scipy.misc import comb as choose
import libPoMo as lp
import pdb


# define PoMo10 states
codons = ["aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act",
          "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att",
          "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct",
          "cga", "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt",
          "gaa", "gac", "gag", "gat", "gca", "gcc", "gcg", "gct",
          "gga", "ggc", "ggg", "ggt", "gta", "gtc", "gtg", "gtt",
          "taa", "tac", "tag", "tat", "tca", "tcc", "tcg", "tct",
          "tga", "tgc"]
nucs = ["A", "C", "G", "T"]


# Define mutation models.
mutmod = {}
mutmod["F81"] = ["global mu=0.01;\n", "mac:=mu;\n", "mag:=mu;\n",
                 "mat:=mu;\n", "mca:=mu;\n", "mct:=mu;\n",
                 "mcg:=mu;\n", "mgc:=mu;\n", "mga:=mu;\n",
                 "mgt:=mu;\n", "mta:=mu;\n", "mtc:=mu;\n",
                 "mtg:=mu;\n"]
mutmod["HKY"] = ["global kappa=0.01;\n", "global mu=0.01;\n",
                 "mac:=mu;\n", "mag:=kappa;\n", "mat:=mu;\n",
                 "mca:=mu;\n", "mct:=kappa;\n", "mcg:=mu;\n",
                 "mgc:=mu;\n", "mga:=kappa;\n", "mgt:=mu;\n",
                 "mta:=mu;\n", "mtc:=kappa;\n", "mtg:=mu;\n"]
mutmod["GTR"] = ["global muac=0.01;\n", "global muag=0.01;\n",
                 "global muat=0.01;\n", "global mucg=0.01;\n",
                 "global muct=0.01;\n", "global mugt=0.01;\n",
                 "mac:=muac;\n", "mag:=muag;\n", "mat:=muat;\n",
                 "mca:=muac;\n", "mct:=muct;\n", "mcg:=mucg;\n",
                 "mgc:=mucg;\n", "mga:=muag;\n", "mgt:=mugt;\n",
                 "mta:=muat;\n", "mtc:=muct;\n", "mtg:=mugt;\n"]
mutmod["NONREV"] = ["global mac=0.01;\n", "global mag=0.01;\n",
                    "global mat=0.01;\n", "global mcg=0.01;\n",
                    "global mct=0.01;\n", "global mgt=0.01;\n",
                    "global mca=0.01;\n", "global mga=0.01;\n",
                    "global mta=0.01;\n", "global mgc=0.01;\n",
                    "global mtc=0.01;\n", "global mtg=0.01;\n"]


# Define selection models.
selmod = {}
selmod["NoSel"] = ["sc := 0.0;\n", "sa := 0.0;\n", "st := 0.0;\n",
                   "sg := 0.0;\n"]
selmod["GCvsAT"] = ["global Sgc=0.0001;\n", "sc := Sgc;\n", "sa := 0.0;\n",
                    "st := 0.0;\n", "sg := Sgc;\n"]
selmod["AllNuc"] = ["global sc=0.0003;\n", "global sg=0.0003;\n",
                    "sa := 0.0;\n", "global st=0.0001;\n"]


def mutModel(mm):
    """Mutation model **type** for argparse."""
    value = str(mm)
    if not (mm in mutmod.keys()):
        msg = "%r is not a valid mutation model" % mm
        raise argparse.ArgumentTypeError(msg)
    return value


def selModel(sm):
    """Selection model **type** for argparse."""
    value = str(sm)
    if not (sm in selmod.keys()):
        msg = "%r is not a valid selection model" % sm
        raise argparse.ArgumentTypeError(msg)
    return value


def dsRatio(dsR):
    """Downsampling ratio **type** for argparse."""
    value = float(dsR)
    if not (0 < value <= 1):
        msg = "%r is not a valid downsampling ratio" % dsR
        raise argparse.ArgumentTypeError(msg)
    return value


def setGM(gm):
    """Set variable mutation rate, if `gm` is given."""
    if gm > 0:
        mutgamma = ["global shape;\n",
                    "category rateCatMut =(" + str(gm) +
                    ", EQUAL, MEAN, GammaDist(_x_,shape,shape), "
                    "CGammaDist(_x_,shape,shape),0,1e25);\n"]
    else:
        mutgamma = ["rateCatMut := 1.0;\n"]
    return mutgamma


def setGS(gs):
    """Set fixation bias, if `gs` is given."""
    if gs > 0:
        selgamma = ["global shape2;\n",
                    "category rateCatSel =(" + str(gs) +
                    ", EQUAL, MEAN, GammaDist(_x_,shape2,shape2), "
                    "CGammaDist(_x_,shape2,shape2),0,1e25);\n"]
    else:
        selgamma = ["rateCatSel := 1.0;\n"]
    return selgamma


def a(n):
    """Calculate the Watterson's Theta coefficient."""
    ret = 0
    for i in range(n-1):
        ret += (float(1.0)/(i+1))
    return ret


def is_number(s):
    """Determine if value is an integer."""
    try:
        int(s)
        return True
    except ValueError:
        return False


def binom(s, p, n):
    """Binomial Distribution

    Calculate the binomial sampling probability (not very efficient,
    but not much effieciency is needed with small samples).

    """
    prob = (choose(n, s) * p**s * (1-p)**(n-s))
    return prob


def probability_matrix(n):
    """Create probability matrices for the HyPhy batch file."""
    o = n-1
    #ignore values below this threshold (keeps the matrix sparse,
    #avoiding increase in computational demands)
    lim = 0.0001
    s = ""
    polys = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]

    #write matrix
    s += "matrixto"+str(o+1)+" ={\n"
    for nucs in range(4):
        s += "{"
        for l in range(58-1):
            if l == nucs:
                s += "1.0,"
            else:
                s += "0.0,"
        s += "0.0}\n"
    for pol in range(6):
        for fre in range(9):
            s += "{"
            for nucs in range(4):
                if nucs == polys[pol][0]:
                    val = binom(o+1, float(9-fre)/10, o+1)
                    if val > lim:
                        s += str(val)+","
                    else:
                        s += "0.0,"
                elif nucs == polys[pol][1]:
                    val = binom(o+1, float(fre+1)/10, o+1)
                    if val > lim:
                        s += str(val)+","
                    else:
                        s += "0.0,"
                else:
                    s += "0.0,"
            for pol2 in range(6):
                for fre2 in range(9):
                    if pol == pol2:
                        if fre2 < o:
                            val = binom(fre2+1, float(fre+1)/10, o+1)
                            if val > lim:
                                s += str(val)
                            else:
                                s += "0.0"
                        else:
                            s += "0.0"
                        if pol2*fre2 != 40:
                            s += ","
                        else:
                            s += "}\n"
                    else:
                        if pol2*fre2 != 40:
                            s += "0.0,"
                        else:
                            s += "0.0}\n"
    s += "};\n\n\n\n"
    s += "Model Mto" + str(o+1) + " = (\"matrixto" + \
         str(o+1) + "\", Freqs, EXPLICIT_FORM_MATRIX_EXPONENTIAL);\n\n"
    return s


def get_species_from_cf_headerline(line):
    """Get the number of species and the names fom a counts format header line.

    :param str line: The header line.

    :rtype: (int n_species, [str] sp_names)

    """
    sp_names = line.split()[2:]
    n_species = len(sp_names)

    if n_species < 2:
        print("Error: Not sufficiently many species (<2).\n")
        raise ValueError()

    return (n_species, sp_names)


def get_data_from_cf_line(l, n_species):
    """Read in the data of a single counts format line.

    The return type is a list with the number of samples and a two
    dimensional array of the form data[species][nucleotide], where
    species is the index of the species and nucleotide is the index of
    the nucleotide (0,1,2 or 3 for a,c,g and t, respectively).

    :param str l: The line that contains the data.
    :param int n_species: Number of species to be expected.

    :rtype: ([int] n_samples, [[int]] data)

    """
    linelist = l.split()[2:]

    if len(linelist) != n_species:
        print("Error: input line \"" + l +
              "\" does not fit number of species.\n")
        raise ValueError()

    n_samples = []
    data = []
    for i in range(n_species):
        p = linelist[i].split(",")
        q = []
        summ = 0
        if len(p) == 1:
            p = p[0].split("/")
        for j in range(4):
            q.append(int(p[j]))
            summ += q[j]

        n_samples.append(summ)
        data.append(q)

    return (n_samples, data)


def read_data_write_HyPhy_input(fn, N, thresh, path_bf,
                                muts, mutgamma,
                                sels, selgamma,
                                PoModatafile, PoModatafile_cons,
                                vb=None):
    """Read the count data and write the HyPhy input file.

    The provided filename has to point to a data file in counts format
    (cf. :doc:`cf <cf>`).  The data will be downsampled if necessary
    and the HyPhy batch and input files will be written.  The number
    of species, the species names, the number of species samples and
    the theta value (usr_def) will be returned in a tuple.

    :param str fn: Counts format file name.
    :param int N: Virtual population size.
    :param float thresh: Trheshold of data discard for downsampling.
    :param str path_bf: Path to the HyPhy batch files
    :param str muts: Mutation model (:func:`mutModel`).
    :param str mutgamma: Gamma of the mutation model (:func:`setGM`).
    :param str sels: Selection model (:func:`selModel`).
    :param str selgamma: Gamma of selection model (:func:`setGS`).
    :param str PoModatafile: Path to HyPhy input file.
    :param str PoModatafile_cons: Path to HyPhy input file.

    :param Boolean vb: Verbosity.

    :rtype: (int n_species, [str] sp_names, [str] sp_samples, Boolean all_one,
             float usr_def)

    """
    # define variables
    # number of species
    n_species = 0
    # species names
    sp_names = []
    # sample size of each species
    sp_samples = []
    # actual data; it is a 3-dimensional array sp_data[species][pos][base]
    sp_data = []

    # Fri Feb 14 13:40:21 CET 2014
    # Depcrecated, see note to fasta file format input below.
    # line = infile.readline()
    # while line[0] != ">":
    #     line = infile.readline()
    #     if line == "":
    #         break

    # if line == "":

    # TODO use np arrays, Problem: fixed sized arrays
    # infile = lp.seqbase.gz_open(in_name)
    # for line in infile:

    if vb is not None:
        print("Starting to read input file.")
    line = ''
    infile = lp.seqbase.gz_open(fn)
    while len(line) > 0 and line[0] == "#":
        line = infile.readline()
    while len(line.split()) == 0:
        line = infile.readline()
        if line == "":
            print("Error: No Data.\n")
            exit()

    # Assign species names (first two columns are Chrom and Pos).
    (n_species, sp_names) = get_species_from_cf_headerline(line)
    # Initialize the number of species samples to 0.
    for i in range(n_species):
        sp_data.append([])
        sp_samples.append(0)

    # Read in the data.
    leng = 0
    for l in infile:
        leng += 1
        (n_samples, data) = get_data_from_cf_line(l, n_species)
        # Update sp_data and the number of samples.
        for i in range(n_species):
            sp_data[i].append(data[i])
            if n_samples[i] > sp_samples[i]:
                sp_samples[i] = n_samples[i]

    if vb is not None:
        print("Count file has been read.")

    # pdb.set_trace()

    # Fri Feb 14 13:38:43 CET 2014
    # Support for fasta file format input has been removed.
    # Reasons: Performance and clarity.
    # Scripts for fasta to counts file format conversion are provided.

    # # In case of Fasta format:
    # while line != "":
    #     if line[0] == ">":
    #         linelist = line.split()
    #         name = linelist[0].split("_")[0].replace(">", "")
    #         if len(linelist) > 1:
    #             init_data = linelist[len(linelist)-1]
    #         else:
    #             init_data = ""
    #         found = 0
    #         for i in range(len(sp_names)):
    #             if name == sp_names[i]:
    #                 found = 1
    #                 sp_samples[i] += 1
    #                 sp_data[i].append("")
    #                 if init_data != "":
    #                     line = init_data
    #                 else:
    #                     line = infile.readline()
    #                 while (len(line) > 0 and line[0] != ">"):
    #                     sp_data[i][len(sp_data[i])-1] += \
    #                         line.replace("\n", "")
    #                     line = infile.readline()
    #                 if sp_data[i][len(sp_data[i])-1] == "":
    #                     print("\n\n\nSpecies " + sp_names[i] + " sample "
    #                           + linelist[0].split("_")[1] +
    #                           " has no data. PoMo is stopping here."
    #                           "Please check your data file.\n\n\n")
    #                     exit()
    #                 break
    #         if found == 0:
    #             sp_names.append(name)
    #             n_species += 1
    #             sp_samples.append(1)
    #             sp_data.append([""])
    #             if init_data != "":
    #                 line = init_data
    #             else:
    #                 line = infile.readline()
    #             while (len(line) > 0 and line[0] != ">"):
    #                 sp_data[len(sp_data)-1][0] += line.replace("\n", "")
    #                 line = infile.readline()
    #             if sp_data[len(sp_data)-1][0] == "":
    #                 print("\n\n\nSpecies " + sp_names[len(sp_data)-1] +
    #                       " sample " + linelist[0].split("_")[1] +
    #                       " has no data. PoMo is stopping here."
    #                       " Please check your data file.\n\n\n")
    #                 exit()
    #     if len(line) > 0 and line[0] != ">":
    #         line = infile.readline()

    # # Put fasta data in counts format
    # DNA = ["A", "C", "G", "T"]
    # DNA2 = ["a", "c", "g", "t"]
    # if VCF == 0:
    #     sp_data2 = sp_data
    #     sp_data = []
    #     for i in range(n_species):
    #         sp_data.append([])
    #     leng = len(sp_data2[0][0])
    #     for i in range(n_species):
    #         for l in range(sp_samples[i]):
    #             if len(sp_data2[i][l]) != leng:
    #                 print("\n\n\nError: individuals have different number "
    #                       "of bases (not a proper alignment).\n\n\n")
    #                 exit()
    #     for l in range(n_species):
    #         for m in range(leng):
    #             count = [0, 0, 0, 0]
    #             for k in range(sp_samples[l]):
    #                 for d in range(4):
    #                     if sp_data2[l][k][m] == DNA[d] or \
    #                        sp_data2[l][k][m] == DNA2[d]:
    #                         count[d] += 1
    #                         break
    #             p = count
    #             sp_data[l].append(p)

    # Sites where some species have coverage 0 are removed
    to_remove = []
    for i in range(leng):
        total = 1
        for s in range(n_species):
            summ = 0
            for d in range(4):
                summ += sp_data[s][i][d]
            if summ == 0:
                total = 0
                break
        if total == 0:
            to_remove.append(i)
    summ = 0
    for i in range(len(to_remove)):
        for s in range(n_species):
            sp_data[s].pop(to_remove[i]-summ)
        summ += 1

    # Debugging point to improve memory.
    # pdb.set_trace()

    # Now, downsample if necessary
    print("Doing downsampling\n")
    sp_samples2 = []
    for i in range(n_species):
        if sp_samples[i] > N:
            sp_samples2.append(N)
        else:
            sp_samples2.append(sp_samples[i])

    advantages = {}
    covered = 0
    for i in range(len(sp_data[0])):
        summs = []
        newlims = []
        cov = 1
        for s in range(n_species):
            summs.append(0)
            newlims.append(sp_samples2[s])
            for d in range(4):
                summs[s] += sp_data[s][i][d]
            if summs[s] < sp_samples2[s]:
                newlims[s] = summs[s]
                cov = 0
        limkey = ""
        for ne in range(len(newlims)):
            limkey += (str(newlims[ne])+":")
        if cov == 1:
            covered += 1
        elif limkey in advantages.keys():
            advantages[limkey] += 1
        else:
            advantages[limkey] = 1
    ke = list(advantages)
    while float(covered)/leng < thresh:
        increments = []
        advs = []
        for s in range(n_species):
            advs.append(0)
            increments.append(1)
            while advs[s] == 0:
                for k in range(len(ke)):
                    kl = ke[k].split(":")
                    valid = 1
                    for s2 in range(n_species):
                        if s2 != s and int(kl[s2]) < sp_samples2[s2]:
                            valid = 0
                    if valid == 1 and int(kl[s]) >= \
                       sp_samples2[s] - increments[s] \
                       and int(kl[s]) < sp_samples2[s]:
                        advs[s] += advantages[ke[k]]
                if advs[s] == 0:
                    if increments[s] < sp_samples2[s] - 1:
                        increments[s] += 1
                    else:
                        break
        max_ad = 0
        max_ind = -1
        for s in range(n_species):
            if advs[s] > max_ad:
                max_ad = advs[s]
                max_ind = s
        if max_ad == 0:
            print("Downsampling with threshold " + str(thresh) +
                  " reached an empasse. "
                  "Please lower the threshold using option "
                  "--DS, change downsampling strategy, "
                  "or ask for assistance!\n")
            exit()
        sp_samples2[max_ind] = sp_samples2[max_ind] - increments[max_ind]
        covered += max_ad
    sp_samples = sp_samples2

    # Sites where some species have not sufficient coverage are removed
    to_remove = []
    for i in range(len(sp_data[0])):
        total = 1
        for s in range(n_species):
            summ = 0
            for d in range(4):
                summ += sp_data[s][i][d]
            if summ < sp_samples[s]:
                total = 0
                break
        if total == 0:
            to_remove.append(i)
    summ = 0
    for i in range(len(to_remove)):
        for s in range(n_species):
            sp_data[s].pop(to_remove[i]-summ)
        summ += 1
    leng = len(sp_data[0])

    print("Number of species: ", str(n_species))
    print("Sample sizes effectively used: ", sp_samples)
    # print(sp_data)
    all_one = True
    for i in range(n_species):
        if sp_samples[i] != 1:
            all_one = False
        if sp_samples[i] > N:
            print("\n\n\nWarning: the number of samples " +
                  str(sp_samples[i]) +
                  " is bigger than the virtual population size " + str(N) +
                  ". The considered species will be downsampled to " + str(N) +
                  ". This is usually not a problem, "
                  "but if you want to avoid this, "
                  "if possible please increase the virtual population size."
                  "\n\n\n")
    if all_one is True:
        usr_def = float(input("""\n\n\nAll species have a sample size of
        1, therefore there is no information at the population level,
        which is required by PoMo. So, please enter a guessed or otherwise
        estimated value for theta (population diversity):\n"""))
    else:
        usr_def = 0.01
    infile.close()

    if n_species < 2:
        print("Error: cannot calculate a tree with fewer than 2 species.")
        exit()

    # default options
    # TODO Why are they not needed
    # sampling = 1  # noqa
    # onlysampling = 1  # noqa
    # mbin = 0  # noqa

    # Writing the HyPhy batch file for PoMo
    newsamfile = open("PoMo10_root_only_sampling_preliminary_used.bf",
                      "w")
    samfile = open(path_bf + "PoMo10_root_only_sampling_preliminary.bf")
    line = "\n"
    while line != "/*Define global parameters*/\n":
        line = samfile.readline()
        linelist = line.split()
        newsamfile.write(line)
    for i in range(23):
        line = samfile.readline()
    for i in range(len(muts)):
        newsamfile.write(muts[i])
    for i in range(len(sels)):
        newsamfile.write(sels[i])
    for i in range(len(mutgamma)):
        newsamfile.write(mutgamma[i])
    for i in range(len(selgamma)):
        newsamfile.write(selgamma[i])
    while line != "/*Find Root*/\n":
        line = samfile.readline()
        linelist = line.split()
        if len(linelist) > 1 and linelist[0] == "fprintf" \
           and linelist[1] == "(stdout," and vb is None:
            newsamfile.write("/*"+line.replace("\n", "")+"*/\n")
        else:
            newsamfile.write(line)
    samples_num = []
    for i in range(n_species):
        if not (sp_samples[i] in samples_num):
            newsamfile.write(lp.main.probability_matrix(sp_samples[i]))
            samples_num.append(sp_samples[i])
            newsamfile.write("\n\n\n")
    line = "\n"
    while line != "":
        line = samfile.readline()
        linelist = line.split()
        if line.split("=")[0] == "\tNsamples":
            newsamfile.write("\tNsamples={{\"")
            for i in range(n_species-1):
                newsamfile.write(str(sp_samples[i])+"\"}{\"")
            newsamfile.write(str(sp_samples[n_species-1])+"\"}};\n")
        elif len(linelist) > 1 and linelist[0] == "fprintf" \
             and linelist[1] == "(stdout," and vb is None:  # noqa
            newsamfile.write("/*"+line.replace("\n", "")+"*/\n")
        else:
            newsamfile.write(line)
    samfile.close()
    newsamfile.close()

    # Writing the HyPhy batch file for PoMo with NNI
    newsamfile = open("PoMo10_NNI_sampling_preliminary_used.bf", "w")
    samfile = open(path_bf + "PoMo10_NNI_sampling.bf")
    line = "\n"
    while line != "/*Define global parameters*/\n":
        line = samfile.readline()
        linelist = line.split()
        newsamfile.write(line)
    for i in range(23):
        line = samfile.readline()
    for i in range(len(muts)):
        newsamfile.write(muts[i])
    for i in range(len(sels)):
        newsamfile.write(sels[i])
    for i in range(len(mutgamma)):
        newsamfile.write(mutgamma[i])
    for i in range(len(selgamma)):
        newsamfile.write(selgamma[i])
    while line != "/*pre-ML*/\n":
        line = samfile.readline()
        linelist = line.split()
        if len(linelist) > 1 and linelist[0] == "fprintf" \
           and linelist[1] == "(stdout," and vb is None:
            newsamfile.write("/*" + line.replace("\n", "") + "*/\n")
        else:
            newsamfile.write(line)
    samples_num = []
    for i in range(n_species):
        if not (sp_samples[i] in samples_num):
            newsamfile.write(lp.main.probability_matrix(sp_samples[i]))
            samples_num.append(sp_samples[i])
            newsamfile.write("\n\n\n")
    line = "\n"
    while line != "":
        line = samfile.readline()
        linelist = line.split()
        if line.split("=")[0] == "\tNsamples":
            newsamfile.write("\tNsamples={{\"")
            for i in range(n_species-1):
                newsamfile.write(str(sp_samples[i])+"\"}{\"")
            newsamfile.write(str(sp_samples[n_species-1])+"\"}};\n")
        elif len(linelist) > 1 and linelist[0] == "fprintf" \
             and linelist[1] == "(stdout," and vb is None:  # noqa
            newsamfile.write("/*"+line.replace("\n", "")+"*/\n")
        else:
            newsamfile.write(line)
    samfile.close()
    newsamfile.close()

    # creating HyPhy input file
    for l in range(n_species):
        PoModatafile.write(">s" + str(l+1) + "\n")
        PoModatafile_cons.write(">s" + str(l+1) + "\n")
        for m in range(leng):
            count = sp_data[l][m]
            p = count
            maxcount = 0
            i2 = -1
            for j2 in range(4):
                if p[j2] > maxcount:
                    i1 = j2
                    maxcount = p[j2]
            refs3 = codons[i1]
            maxcount = 0
            for j2 in range(4):
                if j2 != i1 and p[j2] > maxcount:
                    i2 = j2
                    maxcount = p[j2]
            if i2 == -1:
                refs = codons[i1]
                # refs2 = codons[i1]
            else:
                if p[i1]+p[i2] > sp_samples[l]:
                    count1 = p[i1]
                    count2 = p[i2]
                    newcount1 = 0
                    newcount2 = 0
                    for j5 in range(sp_samples[l]):
                        num = random.random()
                        if num < float(count1)/(count1+count2):
                            newcount1 += 1
                            count1 = count1 - 1
                        else:
                            newcount2 += 1
                            count2 = count2 - 1
                else:
                    newcount1 = p[i1]
                    newcount2 = p[i2]
                if i1 > i2:
                    i3 = i1
                    i1 = i2
                    i2 = i3
                    newcount3 = newcount1
                    newcount1 = newcount2
                    newcount2 = newcount3
                if newcount1 == sp_samples[l]:
                    refs = codons[i1]
                    # refs2 = codons[i1]
                elif newcount2 == sp_samples[l]:
                    refs = codons[i2]
                    # refs2 = codons[i2]
                else:
                    pol = 0
                    if i1 == 1:
                        pol = 3
                    if i1 == 2:
                        pol = 5
                    pol += (i2-(i1+1))
                    p1 = newcount2 - 1
                    pos = 4+pol*(N-1)+p1
                    refs = codons[pos]
            PoModatafile.write(refs)
            PoModatafile_cons.write(refs3)
        PoModatafile.write("\n")
        PoModatafile_cons.write("\n")
    PoModatafile.close()
    PoModatafile_cons.close()

    # Debugging point if necessary.
    # pdb.set_trace()
    return (n_species, sp_names, sp_samples, all_one, usr_def)
