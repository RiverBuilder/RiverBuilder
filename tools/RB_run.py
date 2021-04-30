#############################################################################################
# RB_run.py
#   Written by Anzy Lee, Postdoctoral Scholar, Utah State University
#   Date: 10/21/2020
#############################################################################################
import os
from riverbuilder.core import river
import matplotlib
import matplotlib.pyplot as plt

#############################################################################################
# Input: dir, case_name
case_name = "S1_noSlope"   # Case Name
dir = "../samples/S1"      # Folder Name

#############################################################################################
# Optional inputs: fname, outfolder, log
os.chdir(dir)
fname = case_name + ".txt"
outfolder = case_name
log = case_name + "_log.txt"
#############################################################################################
# Run RiverBuilder
river.buildRiver(fname, outfolder, log)

#plt.close('all')
plt.close(2)
plt.close(3)
plt.close(4)
plt.close(5)
plt.axis('scaled')