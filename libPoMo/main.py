# This file contains functions that are used by PoMo.py
import argparse
from scipy.misc import comb as choose

# Define mutation models
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


# define mutModel type for argparse
def mutModel(mm):
    value = str(mm)
    if not (mm in mutmod.keys()):
        msg = "%r is not a valid mutation model" % mm
        raise argparse.ArgumentTypeError(msg)
    return value


# define selection models
selmod = {}
selmod["NoSel"] = ["sc := 0.0;\n", "sa := 0.0;\n", "st := 0.0;\n",
                   "sg := 0.0;\n"]
selmod["GCvsAT"] = ["global Sgc=0.0001;\n", "sc := Sgc;\n", "sa := 0.0;\n",
                    "st := 0.0;\n", "sg := Sgc;\n"]
selmod["AllNuc"] = ["global sc=0.0003;\n", "global sg=0.0003;\n",
                    "sa := 0.0;\n", "global st=0.0001;\n"]


# define selModel type for argparse
def selModel(sm):
    value = str(sm)
    if not (sm in selmod.keys()):
        msg = "%r is not a valid selection model" % sm
        raise argparse.ArgumentTypeError(msg)
    return value


# define downsampling ratio type for argparse
def dsRatio(dsR):
    value = float(dsR)
    if not (0 < value <= 1):
        msg = "%r is not a valid downsampling ratio" % dsR
        raise argparse.ArgumentTypeError(msg)
    return value


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
    """Calculate the binomial sampling probability (not very efficient,
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


def process_infile(infile):
    """ Process input file. Checks the input file format.

    Fasta format:
    ----------------------------------------------------------------------
    HEADER
    >species_nameA_N1
    ACGTACGT
    >species_nameA_N2
    ACGTACGA
    >species_nameB_N1
    ACGTTGCA

    Counts format:
    ----------------------------------------------------------------------
    #Comments
    species_nameA_N1    species_nameA_N2    species_nameB_N1
    
    """
    pass
