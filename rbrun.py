
import os
from riverbuilder.core import river
import matplotlib
import matplotlib.pyplot as plt

fname = ".\samples\S1\S1.txt"
outfolder = ".\samples\S1\S1_output"
log = ".\samples\S1\S1_output\S1_log.txt"
#############################################################################################
# Run RiverBuilder
river.buildRiver(fname, outfolder, log)
