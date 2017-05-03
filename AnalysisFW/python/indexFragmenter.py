import argparse
parser = argparse.ArgumentParser()
parser.add_argument('filename', metavar="FILENAME", type=str, nargs=1)

indexFileName = parser.parse_args().filename[0]
isData = True

indexFile = open(indexFileName)
indexLines = indexFile.readlines()

if not isData:
	outFileName = indexFileName.replace("CMS_MonteCarlo2011_Summer11LegDR_","")
	outFileName = outFileName.replace("_AODSIM_PU_S13_START53_LV6-v1","")
	outFileName = outFileName.replace("file_index.txt","")
else:
	outFileName = indexFileName.replace("CMS_Run2011A_","")
	outFileName = outFileName.replace("file_index.txt","")

fragmentInd = 1
for line in indexLines:
	outFile = open(outFileName+"{0:d}".format(fragmentInd)+"_fragment.txt","w")
	outFile.write(line)
	fragmentInd += 1
