# This file contains functions that are used by PoMo.py
import argparse
from scipy.misc import comb as choose


#define downsampling ratio type for argparse
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

# def choose(n,m):
#     """Calculate n choose m."""
#     return math.factorial(n) / (math.factorial(m)*math.factorial(n-m))


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
