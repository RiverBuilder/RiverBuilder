from riverbuilder.core import river

fname = ".\\samples\\rosgenC4\\vanilla.txt"
outfolder = ".\\samples\\rosgenC4\\vanilla"
log = outfolder + "\\log.txt"
#############################################################################################
# Run RiverBuilder
river.buildRiver(fname, outfolder, log)

