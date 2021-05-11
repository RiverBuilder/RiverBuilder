from riverbuilder.core import river

fname = ".\samples\S1\S1.txt"
outfolder = ".\samples\S1\S1"
log = outfolder + "\S1_log.txt"
#############################################################################################
# Run RiverBuilder
river.buildRiver(fname, outfolder, log)

