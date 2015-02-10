"""Execute PoMo10.

This script executes PoMo. Run this script with `--help` to print help
information and exit.

"""
import argparse
import sys
import os
import logging
import re
import libPoMo as lp
import pdb
import time

# PoMo version
ver = '1.1.0'

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='PoMo',
                                 description="PoMo10 script version "+ver)
parser.add_argument('hyphy_bin', help="""Path of the HYPHY binary used
to maximize the likelihood.""")

parser.add_argument('file', help="""
Name of the counts or fasta file containing the alignment.  Each
individual's name must be \"species-n\" where \"species\" is the name
of its species, and \"n\" is a number identifying the individual among
the other samples of the same population/species (it actually does not
matter which number is given to which individual).  No gaps or Ns are
allowed, which means, you have to remove columns where one of the
individuals has missing characters or unknowns. It does not matter
from which individual each base comes, so you can sub-sample and
randomly assign to individuals within a species.""")

# TODO Check with Nicola, if this is still true (no gaps or Ns are
# allowed ...)

parser.add_argument('-m', '--molecular-clock', type=int,
                    choices=[0, 1], default=1, help="""Determines if
                    you want the molecular clock constraint (and
                    therefore also look for a root) or not. Default is
                    yes. Type `-m 0` to specify no molecular clock.""")
parser.add_argument('-u', '--MM', type=lp.main.mutModel, default="HKY",
                    help="""Allows to choose a mutation model
                    different from the HKY (default option). `GTR`
                    corresponds to the general time reversible,
                    `F81` to the Felsenstein 1981 (reversible, equal
                    mutation rates), and `NONREV` to the general
                    nonreversible model (all substitution rates are
                    independent). To change, type for example `--MM
                    GTR`.""")
parser.add_argument('-s', '--SM', type=lp.main.selModel, default="NoSel",
                    help="""Allows to choose fixation rates.  `-s
                    NoSel`: fixation rates are equal for all
                    nucleotides; default.  `-s GCvsAT`: one parameter
                    describes fixation difference of GC versus AT, as
                    is expected from biased gene conversion.  `-s
                    AllNuc`: each nucleotide has a different fitness
                    (3 free parameters since only fitness differences
                    matter).  Warning: estimating fixation biases will
                    be more time consuming.""")
parser.add_argument('-g', '--GM', type=int, default=0, help="""Allows
to set a variable mutation rate over sites, gamma-distributed,
approximated with a number of classes as specified by the
user. Default: uniform mutation rate. For example, `--GM 6` specifies
a gamma distribution of total mutation rate with 6 discrete
cathegories.""")
parser.add_argument('-f', '--GS', type=int, default=0, help="""Allows
to set a variable fixation bias over sites, gamma-distributed,
approximated with a number of classes as specified by the
user. Default: uniform fixation rate (`f 0`).""")
parser.add_argument('-d', '--ds-ratio', type=lp.main.dsRatio, default=0.66,
                    help="""Determines which proportion of the data is
                    kept after downsampling. Downsampling is done when
                    sites do not have the same coverage along the
                    genome. In such a case, all sites with coverage
                    higher than a certain sample size are downsampled,
                    those with lower coverage are discarded. The
                    threshold sample size is chosen as the highest
                    possible that leaves the kept number of sites
                    above the specified threshold.  By default, sample
                    sizes are decreased until at least 2 thirds of the
                    sites are included (`-d 0.66`).""")
parser.add_argument('-t', "--theta", type=float, default=None,
                    help="""Manyally set population diversity theta.
                    This value can only be set, if all species have a
                    sample size of 1.  If no value is specified, PoMo
                    will ask the user for theta on the command line if
                    necessary.""")
parser.add_argument('-v', '--verbose', action='count',
                    help="""turn on verbosity (-v or -vv)""")
parser.add_argument('--version', action='version', version='%(prog)s '+ver)

args = parser.parse_args()

if args.molecular_clock == 1:
    noMC = 0
elif args.molecular_clock == 0:
    noMC = 1

# Mutation model.
muts = lp.main.mutmod[args.MM]

# Variable mutation rate (+Gamma).
mutgamma = lp.main.setGM(args.GM)

# Fixation bias.
selgamma = lp.main.setGS(args.GS)

# Selection model.
sels = lp.main.selmod[args.SM]

# Population diversity theta.  Needs to be multiplied by a(10), so
# that it refers to the expected number of polymorphisms per base.
theta = lp.main.a(10) * args.theta

# Verbosity and logger.
vb = args.verbose
# Verbose HYPHY output only with -vv or more.
if (vb is None) or (vb == 1):
    vbHyphy = None

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()
if args.verbose == 0:
    logger.setLevel(logging.WARN)
elif args.verbose == 1:
    logger.setLevel(logging.INFO)
elif args.verbose == 2:
    logger.setLevel(logging.DEBUG)

# Threshold of data discard for downsampling
thresh = args.ds_ratio
print("============================================================")
print("PoMo", ver)
print(""" Created by Nicola De Maio; maintained by Dominik Schrempf. For a
reference, please see and cite: De Maio, Schlotterer, Kosiol (MBE,
2013), and/or: De Maio, Schrempf, Kosiol (in preparation). You can use
this software for non-commercial purposes, but please, always
acknowledge the authors. For suggestions, doubts, bugs, etc., please
contact nicola.de.maio.85@gmail.com""")

print("===========================================================")
print("Start Time:", lp.main.timeStr())
start_time = time.time()
if (vb is not None):
    print("Verbose mode.")
    print("===========================================================")
# Define paths to files.
in_name = str(args.file)
in_name_no_extension = in_name.rsplit(".", maxsplit=1)[0]
in_basename_no_extension = os.path.basename(in_name)
in_basename_no_extension = in_basename_no_extension.rsplit(".", maxsplit=1)[0]
out_name = in_basename_no_extension + "_PoMo_output.txt"

# Define the names of the PoMo data files; they are created in the
# current working directory.
PoModata_name = in_basename_no_extension + "_PoMo_HyPhy.txt"
PoModata_name_cons = in_basename_no_extension + "_consensus_HyPhy.txt"

# Create file descriptors.
PoModatafile = open(PoModata_name, "w")
PoModatafile_cons = open(PoModata_name_cons, "w")

# Get path of data file.
path_data = os.path.abspath(os.path.dirname(in_name))
path_data = path_data + "/"

# Gget currnt working directory.
path_cwd = os.getcwd()
path_cwd = path_cwd + "/"

# Get path of PoMo.py.
try:
    path_PoMo = os.path.abspath(os.path.dirname(__file__))
except:
    path_PoMo = os.path.abspath(os.path.dirname(sys.argv[0]))
path_PoMo = path_PoMo + "/"

# Define path of batchfiles.
path_bf = path_PoMo + "batchfiles/"

# Get path of HyPhy.
HyPhy_bin = str(args.hyphy_bin)
path_HyPhy = os.path.abspath(os.path.dirname(HyPhy_bin))
path_HyPhy = path_HyPhy + "/"

# Define virtual population size.
N = 10

# Read in the data and write the HyPhy batch and input files.
(n_species, sp_names, sp_samples, all_one, usr_def) \
    = lp.main.read_data_write_HyPhy_input(in_name, N, thresh, path_bf,
                                          muts, mutgamma,
                                          sels, selgamma,
                                          PoModatafile, PoModatafile_cons,
                                          theta, vb)

# Debugging point if necessary.
# print(n_species, sp_names, sp_samples)
# pdb.set_trace()

print("============================================================")
print("\nRunning 1: NJ consensus\n")
# Run HyPhy concatenation, NJ and root positioning, on consensus data
HPfile = open(path_bf + "Nuc_NJandRoot.bf")
HPfile2 = open("Nuc_NJandRoot_used.bf", "w")
line = HPfile.readline()
line = HPfile.readline()
HPfile2.write("inp=\"" + PoModata_name_cons + "\";\n")
HPfile2.write("out2=\"" + in_basename_no_extension +
              "_consensus_NJandRoot_out.txt\";\n")
while line != "":
    line = HPfile.readline()
    linelist = line.split()
    if len(linelist) > 0 and linelist[0] == "ExecuteAFile":
        HPfile2.write(line.replace("pairwise", path_bf + "pairwise"))
    elif len(linelist) > 1 and linelist[0] == "fprintf" \
         and linelist[1] == "(stdout," and vbHyphy is None:  # noqa
        HPfile2.write("/*" + line.replace("\n", "") + "*/\n")
    else:
        HPfile2.write(line)
HPfile2.close()
HPfile.close()
os.system(HyPhy_bin + " Nuc_NJandRoot_used.bf")

file = open(in_basename_no_extension + "_consensus_NJandRoot_out.txt")
line = file.readline()
linel = line.split()
while line != "":
    line = file.readline()
    linel = line.split()
    if (len(linel) > 1 and linel[0] == "Tree"):
        lasttree = linel[1].replace("givenTree=", "").replace("testTree=", "")
        for i in range((2 * n_species) - 3):
            lasttree = lasttree.replace("Node" + str(n_species-i), "")
NucNJtree_cons = lasttree
file.close()

if n_species > 3:
    print("============================================================")
    print("\nRunning 2: NNI consensus\n")
    # Running HyPhy concatenation, finding topology and root altogether
    # with NNI and rooting, on consensus data
    HPfile = open(path_bf + "Nuc_NNIwithRoot.bf")
    HPfile3 = open("Nuc_NNIwithRoot_used.bf", "w")
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    HPfile3.write("inp = \"" + PoModata_name_cons + "\";\n")
    HPfile3.write("out2=\"" + in_basename_no_extension +
                  "_consensus_NNIwithRoot_out.txt\";\n")
    HPfile3.write("treeString=\"" + NucNJtree_cons + "\";\n")
    while line != "":
        line = HPfile.readline()
        linelist = line.split()
        if len(linelist) > 0 and linelist[0] == "ExecuteAFile":
            HPfile3.write(line.replace("pairwise", path_bf + "pairwise"))
        elif len(linelist) > 0 and linelist[0] == "#include":
            HPfile3.write(line.replace("heuristic", path_bf + "heuristic"))
        elif len(linelist) > 1 and linelist[0] == "fprintf" \
             and linelist[1] == "(stdout," and vbHyphy is None:  # noqa
            HPfile3.write("/*" + line.replace("\n", "") + "*/\n")
        else:
            HPfile3.write(line)
    HPfile3.close()
    HPfile.close()
    os.system(HyPhy_bin + " Nuc_NNIwithRoot_used.bf")
    HPofile = open(in_basename_no_extension + "_consensus_NNIwithRoot_out.txt")
    line = HPofile.readline()
    while line != "":
        line = HPofile.readline()
        linel = line.split()
        if (len(linel) > 1 and linel[0] == "Tree"):
            lasttree = linel[1].replace("givenTree=", "")
            lasttree = lasttree.replace("testTree=", "")
            for i in range((2*n_species) - 3):
                lasttree = lasttree.replace("Node" + str(n_species-i), "")
    consetree = lasttree
    HPofile.close()

    print("============================================================")
    print("\nRunning 3: NNI PoMo\n")
    # Running PoMo10, finding the topology with NNI
    HPfile = open("PoMo10_NNI_sampling_preliminary_used.bf")
    HPfile2 = open("PoMo10_NNI_sampling_used.bf", "w")
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    HPfile2.write("inp = \"" + PoModata_name + "\";\n")
    HPfile2.write("out2=\"PoMo10_NNI_sampling_out.txt\";\n")
    if all_one is True:
        HPfile2.write("user_defining=1;\n")
    else:
        HPfile2.write("user_defining=0;\n")
    HPfile2.write("user_defined_Ppol="+str(usr_def)+";\n")
    if all_one is True:
        HPfile2.write("scale_Ppol:=1.0;\n")
    else:
        a_total = 0.0
        for i in range(n_species):
            a_total += lp.main.a(sp_samples[i])
        a_total = a_total/n_species
        HPfile2.write("scale_Ppol:=" + str(lp.main.a(N)/a_total)+";\n")
    HPfile2.write("sample:=0;\n")
    NJtree2 = consetree
    NJtree2Samp = NJtree2
    for ss in range(n_species):
        NJtree2Samp = NJtree2Samp.replace("s"+str(n_species-ss),
                                          "(z"+str(n_species-ss) +
                                          "{Mto" +
                                          str(sp_samples[(n_species-1)-ss]) +
                                          "}:1.0)u" +
                                          str(n_species-ss)+"{M1}")
    for ss in range(n_species):
        NJtree2Samp = NJtree2Samp.replace("z", "s")
    HPfile2.write("treeString=\""+NJtree2Samp+"\";\n")
    HPfile2.write("NoSampTree=\""+consetree+"\";\n")
    while line != "":
        line = HPfile.readline()
        linelist = line.split()
        if len(linelist) > 0 and linelist[0] == "ExecuteAFile":
            HPfile2.write(line.replace("pairwise", path_bf + "pairwise"))
        elif len(linelist) > 0 and linelist[0] == "#include":
            HPfile2.write(line.replace("heuristic", path_bf + "heuristic"))
        else:
            HPfile2.write(line)
    HPfile2.close()
    HPfile.close()
    os.system(HyPhy_bin + " PoMo10_NNI_sampling_used.bf")
    HPofile = open("PoMo10_NNI_sampling_out.txt")
    line = HPofile.readline()
    while line != "":
        line = HPofile.readline()
        linel = line.split()
        if (len(linel) > 1 and linel[0] == "Tree"):
            lasttree = linel[1].replace("givenTree=", "")
            lasttree = lasttree.replace("testTree=", "")
            for i in range((2*n_species) - 3):
                lasttree = lasttree.replace("Node" + str(n_species-i), "")
    NNItreesamp = lasttree
    HPofile.close()
else:
    NNItreesamp = NucNJtree_cons

# What happens when there is no molecular clock?
if n_species > 3 and noMC == 1:
    # If PoMo NNI has been done, output outcome
    HPofile = open("PoMo10_NNI_sampling_out.txt")
    line = HPofile.readline()
    out_wr = line
    out_wr2 = ""
    while line != "":
        line = HPofile.readline()
        linel = line.split()
        if (len(linel) > 1 and linel[0] == "Tree"):
            lasttree = linel[1].replace("givenTree=", "")
            lasttree = lasttree.replace("testTree=", "")
            for i in range((2*n_species)):
                lasttree = lasttree.replace("Node" + str(n_species-i), "")
            out_wr2 = out_wr
            out_wr = ""
        else:
            out_wr += line
    swap_fast_samp_tree = lasttree
    HPofile.close()


elif n_species <= 3 and noMC == 1:
    # If no PoMo NNI has been done, do a single ML run without looking
    # for best topology.

    # Write the HyPhy batch file for PoMo without NNI
    newsamfile = open("PoMo10_NoMolClock_preliminary.bf", "w")
    samfile = open(path_bf + "PoMo10_NoMolClock.bf")
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
           and linelist[1] == "(stdout," and vbHyphy is None:
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
        if len(linelist) > 1 and linelist[0] == "fprintf" \
           and linelist[1] == "(stdout," and vbHyphy is None:
            newsamfile.write("/*" + line.replace("\n", "") + "*/\n")
        else:
            newsamfile.write(line)
    samfile.close()
    newsamfile.close()

    print("============================================================")
    print("\nRunning PoMo without Molecular clock\n")
    # Running PoMo10, finding the topology with NNI
    HPfile = open("PoMo10_NoMolClock_preliminary.bf")
    HPfile2 = open("PoMo10_NoMolClock_used.bf", "w")
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    HPfile2.write("inp = \"" + PoModata_name + "\";\n")
    HPfile2.write("out2=\"PoMo10_NoMolClock_out.txt\";\n")
    if all_one == 1:
        HPfile2.write("user_defining=1;\n")
    else:
        HPfile2.write("user_defining=0;\n")
    HPfile2.write("user_defined_Ppol="+str(usr_def)+";\n")
    if all_one == 1:
        HPfile2.write("scale_Ppol:=1.0;\n")
    else:
        a_total = 0.0
        for i in range(n_species):
            a_total += lp.main.a(sp_samples[i])
        a_total = a_total/n_species
        HPfile2.write("scale_Ppol:="+str(lp.main.a(N)/a_total)+";\n")
    HPfile2.write("sample:=0;\n")
    NJtree2 = NucNJtree_cons
    NJtree2Samp = NJtree2
    for ss in range(n_species):
        NJtree2Samp = NJtree2Samp.replace("s" + str(n_species-ss), "(z" +
                                          str(n_species-ss) + "{Mto" +
                                          str(sp_samples[(n_species-1)-ss]) +
                                          "}:1.0)u" + str(n_species-ss) +
                                          "{M1}")
    for ss in range(n_species):
        NJtree2Samp = NJtree2Samp.replace("z", "s")
    HPfile2.write("treeString=\""+NJtree2Samp+"\";\n")
    HPfile2.write("NoSampTree=\""+NucNJtree_cons+"\";\n")
    while line != "":
        line = HPfile.readline()
        linelist = line.split()
        if len(linelist) > 0 and linelist[0] == "ExecuteAFile":
            HPfile2.write(line.replace("pairwise", path_bf + "pairwise"))
        elif len(linelist) > 0 and linelist[0] == "#include":
            HPfile2.write(line.replace("heuristic", path_bf + "heuristic"))
        else:
            HPfile2.write(line)
    HPfile2.close()
    HPfile.close()
    os.system(HyPhy_bin + " PoMo10_NoMolClock_used.bf \n")
    HPofile = open("PoMo10_NoMolClock_out.txt")
    line = HPofile.readline()
    out_wr = line
    out_wr2 = ""
    while line != "":
        line = HPofile.readline()
        linel = line.split()
        if (len(linel) > 1 and linel[0] == "Tree"):
            lasttree = linel[1].replace("givenTree=", "")
            lasttree = lasttree.replace("testTree=", "")
            for i in range((2*n_species)):
                lasttree = lasttree.replace("Node" + str(n_species-i), "")
            out_wr2 = out_wr
            out_wr = ""
        else:
            out_wr += line
    swap_fast_samp_tree = lasttree
    HPofile.close()

else:
    print("============================================================")
    print("\nRunning 4: Rooting PoMo\n")
    # Running PoMo10, finding root from the topology estimated with NNI
    # and PoMo10
    HPfile = open("PoMo10_root_only_sampling_preliminary_used.bf")
    HPfile2 = open("PoMo10_root_only_sampling_used.bf", "w")
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    line = HPfile.readline()
    HPfile2.write("inp = \"" + PoModata_name + "\";\n")
    HPfile2.write("out2=\"PoMo10_NNI_sampling_rooted_out.txt\";\n")
    if all_one == 1:
        HPfile2.write("user_defining=1;\n")
    else:
        HPfile2.write("user_defining=0;\n")
    HPfile2.write("user_defined_Ppol=" + str(usr_def) + ";\n")
    if all_one == 1:
        HPfile2.write("scale_Ppol:=1.0;\n")
    else:
        a_total = 0.0
        for i in range(n_species):
            a_total += lp.main.a(sp_samples[i])
        a_total = a_total/n_species
        HPfile2.write("scale_Ppol:="+str(lp.main.a(N)/a_total)+";\n")

    NJtree2 = NNItreesamp
    while line != "":
        line = HPfile.readline()
        if len(line.split("=")) > 1 and line.split("=")[0] == "treeString" \
           and line.split("=")[1].replace("\n", "") == "\"\";":
            lasttree = NJtree2
            if n_species > 3:
                for spes in range(n_species):
                    pap = re.compile("\\(s" + str(spes+1) + ":" +
                                     "\\d+(\\.\\d+)?" + "(e-\\d\\d)?\\)u" +
                                     str(spes+1))
                    mam = pap.search(lasttree)
                    lasttree = lasttree.replace(mam.group(), "s" + str(spes+1))
            HPfile2.write("treeString=\""+lasttree+"\";\n")
        else:
            HPfile2.write(line)
    HPfile2.close()
    HPfile.close()
    os.system(HyPhy_bin + " PoMo10_root_only_sampling_used.bf")

    HPofile = open("PoMo10_NNI_sampling_rooted_out.txt")
    line = HPofile.readline()
    out_wr = line
    out_wr2 = ""
    while line != "":
        line = HPofile.readline()
        linel = line.split()
        if (len(linel) > 1 and linel[0] == "Tree"):
            lasttree = linel[1].replace("givenTree=", "")
            lasttree = lasttree.replace("testTree=", "")
            for i in range((2*n_species)):
                lasttree = lasttree.replace("Node" + str(n_species-i), "")
            out_wr2 = out_wr
            out_wr = ""
        else:
            out_wr += line
    swap_fast_samp_tree = lasttree
    HPofile.close()

# Write final output to file
for spes in range(n_species):
    pap = re.compile("\\(s" + str(n_species-spes) + ":" + "\\d+(\\.\\d+)?" +
                     "(e-\\d\\d)?\\)u" + str(n_species-spes))
    mam = pap.search(swap_fast_samp_tree)
    swap_fast_samp_tree = swap_fast_samp_tree.replace(mam.group(), sp_names[(n_species - 1) - spes])  # noqa
out_wr2 += swap_fast_samp_tree
outfile = open(out_name, "w")
outfile.write(out_wr2)
outfile.write("\n")
outfile.close

os.system("rm -f Nuc_NJandRoot_used.bf")
os.system("rm -f " + in_basename_no_extension + "_consensus_NJandRoot_out.txt")
os.system("rm -f Nuc_NNIwithRoot_used.bf")
os.system("rm -f " + in_basename_no_extension +
          "_consensus_NNIwithRoot_out.txt")
os.system("rm -f PoMo10_NNI_sampling_used.bf")
os.system("rm -f PoMo10_NNI_sampling_out.txt")
os.system("rm -f PoMo10_root_only_sampling_used.bf")
os.system("rm -f PoMo10_NNI_sampling_rooted_out.txt")
os.system("rm -f PoMo10_root_only_sampling_preliminary_used.bf")
os.system("rm -f PoMo10_NNI_sampling_preliminary_used.bf")
os.system("rm -f " + PoModata_name_cons)
os.system("rm -f " + PoModata_name)
os.system("rm -f PoMo10_NoMolClock_out.txt")
os.system("rm -f PoMo10_NoMolClock_preliminary.bf")
os.system("rm -f PoMo10_NoMolClock_used.bf")
print("Done!")
print("===========================================================")
print("End Time:", lp.main.timeStr())
end_time = time.time()
print("Runtime in seconds:", end_time - start_time)
print("============================================================")
exit()
