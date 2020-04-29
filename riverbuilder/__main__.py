from riverbuilder.core import river
import sys

if __name__ == "__main__":
    fname = ""
    log = ''
    if len(sys.argv) == 1:
        log += "No input file. An input file should be provided.\n"
        print(log)
        sys.exit()
    elif len(sys.argv) == 2:
        fname = sys.argv[1]
        outfolder = "riverbuilder_output"
    elif len(sys.argv) == 3:
        fname = sys.argv[1]
        outfolder = sys.argv[2]
    else:
        fname = sys.argv[1]
        outfolder = sys.argv[2]
        log = sys.argv[3]+'\n'

    log += 'Given input file is: '+fname+'\n'
    log += 'Given output folder is: '+fname+'\n'
    river.buildRiver(fname, outfolder, log)

