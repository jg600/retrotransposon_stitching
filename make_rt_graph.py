import sys, re, os, subprocess, json, linecache, argparse
from math import *
		
def getNames(bedFileName):
	outFileName = re.sub(r'\.bed$', r'.names', bedFileName)
	cmdString = "awk '{print $4}' %s | sort | uniq > %s" % (bedFileName, outFileName)
	subprocess.Popen(cmdString, shell = True)
	return outFileName
	
	
parser = argparse.ArgumentParser(description = "Create a graph representing overlaps of genomic regions from a BED file.")
parser.add_argument("--regionBed", help = "BED file of genomic regions to be analysed", dest = "regionBedFileName", type = str, required = True)
parser.add_argument("--namesFile", help = "File containing the names of the different kinds of region. If not given, a file will be created using the --regionBED file.", type = str, dest = "uniqNamesFileName", default = "NONE")
parser.add_argument("--maxGap", help = "Largest gap between regions that is still counted as an overlap.", default = 100, type = int, dest = "maxGap")
parser.add_argument("--outJSON", help = "File to dump JSON string of results matrix to.", default = "json_dump", type = str, dest = "outJson")

argDict = vars(parser.parse_args())
regionBedFileName = argDict["regionBedFileName"]
uniqNamesFileName = argDict["uniqNamesFileName"]
maxGap = argDict["maxGap"]
outJson = argDict["outJson"]

if uniqNamesFileName != "NONE":
	uniqNamesFile = open(uniqNamesFileName, 'r')
else:
	print("No names file provided. Creating file ...")
	uniqNamesFileName = getNames(regionBedFileName)
	print("Opening file ...")
	uniqNamesFile = open(uniqNamesFileName, 'r')
	print("Done")


#Initialise the matrix using the unique names
uniqNames = []

for line in uniqNamesFile:
	uniqNames.append(line.strip())
	
countMatrix = {row:{col:0 for col in uniqNames} for row in uniqNames}

#Go through the BED file of regions and add to the count matrix whenever an overlap is found

print("Counting lines in BED file ..."),
numLines = int(subprocess.Popen("wc -l %s | awk '{print $1}'" % regionBedFileName, shell = True, stdout = subprocess.PIPE).stdout.read())
print("Done")

n = 1
print("Starting to get overlaps ...")
sys.stdout.flush()

print('['+' '*50+']'),
print('\r'),
percent = 0.0
i = 0.0

for n in range(1, numLines+1):

	percent = 100.0*(n/numLines)
	numBars = int(floor(percent/2))
	print('\r'),
	print('['+'|'*numBars+' '*(50-numBars)+']'),
	sys.stdout.flush()

	baseLineList = re.split(re.compile('\s+'), linecache.getline(regionBedFileName, n).strip())
	baseDict = {'chrom':baseLineList[0], 'start':int(baseLineList[1]), 'end':int(baseLineList[2]), 'name':baseLineList[3]}
	
	ahead = n
	keepGoing = True
	
	while keepGoing:
		ahead += 1
		if ahead <= numLines:
			newLineList = re.split(re.compile('\s+'), linecache.getline(regionBedFileName, ahead).strip())
			newDict = {'chrom':newLineList[0], 'start':int(newLineList[1]), 'end':int(newLineList[2]), 'name':newLineList[3]}
			
			if newDict['chrom'] != baseDict['chrom']:
				keepGoing = False
			else:
				if newDict['start'] > baseDict['end'] + 100:
					keepGoing = False
				else:
					countMatrix[baseDict['name']][newDict['name']] += 1
		else:
			keepGoing = False

print("\nCount matrix finished.")

print("Normalising counts ..."),
normedCounts = {rowName:{} for rowName in countMatrix.keys()}

for rowName in countMatrix.keys():
	rowSum = float(sum(countMatrix[rowName].values()))
	if rowSum != 0:
		normedCounts[rowName] = {colName:countMatrix[rowName][colName]/rowSum for colName in countMatrix[rowName].keys()}
	else:
		normedCounts[rowName] = {colName:0 for colName in countMatrix[rowName].keys()}
print("Done")

with open(outJson, 'w') as f:
	f.write(json.dumps(normedCounts, sort_keys = True, indent = 4, separators = (',', ':')))














