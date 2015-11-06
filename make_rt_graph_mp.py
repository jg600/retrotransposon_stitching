import sys, re, os, subprocess, json, linecache, argparse
from math import *
import multiprocessing as mp

def getNames(bedFileName):
	outFileName = re.sub(r'\.bed$', r'.names', bedFileName)
	cmdString = "awk '{print $4}' %s | sort | uniq > %s" % (bedFileName, outFileName)
	subprocess.Popen(cmdString, shell = True)
	return outFileName
	
def createMatrix(bedFileName, uniqNames, startLine, endLine):
	countMatrix = {row:{col:0 for col in uniqNames} for row in uniqNames}
	
	for n in range(startLine, endLine+1):
		baseLineList = re.split(re.compile('\s+'), linecache.getline(bedFileName, n).strip())
		baseDict = {'chrom':baseLineList[0], 'start':int(baseLineList[1]), 'end':int(baseLineList[2]), 'name':baseLineList[3]}
		
		ahead = n
		keepGoing = True
		
		while keepGoing:
			ahead += 1
			if ahead <= endLine:
				newLineList = re.split(re.compile('\s+'), linecache.getline(bedFileName, ahead).strip())
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
			
	return countMatrix
	
def combineDicts(dictList):
	allKeys = dictList[0].keys()
	outDict = {k:{k:0 for k in allKeys} for k in allKeys}
	for d in dictList:
		for k in allKeys:
			for j in allKeys:
				outDict[k][j] += d[k][j]
	return outDict
				
def normaliseMatrix(countMatrix):
	normedCounts = {rowName:{} for rowName in countMatrix.keys()}

	for rowName in countMatrix.keys():
		rowSum = float(sum(countMatrix[rowName].values()))
		if rowSum != 0:
			normedCounts[rowName] = {colName:countMatrix[rowName][colName]/rowSum for colName in countMatrix[rowName].keys()}
		else:
			normedCounts[rowName] = {colName:0 for colName in countMatrix[rowName].keys()}
	return normedCounts
	
def getChrLines(bedFileName):
	cmdString = "awk '{print $1}' %s | sort | uniq" % bedFileName
	uniqChrs = (subprocess.Popen(cmdString, shell = True, stdout = subprocess.PIPE).stdout.read()).split()
	chrLinesDict = {chrom:{"first":-1, "last":-1} for chrom in uniqChrs}
	for chrom in uniqChrs:
		print("Processing %s" % chrom)
		#nullFds = [os.open(os.devnull, os.O_RDWR) for x in xrange(2)]
		#save = os.dup(1), os.dup(2)
		getChrFirstLineCmd = "cat -n %s | grep -P \'%s\s+\' | head -n 1 | cut -f 1" % (bedFileName, chrom)
		getChrLastLineCmd = "cat -n %s | grep -P \'%s\s+\' | tail -1 | cut -f 1" % (bedFileName, chrom)
		#os.dup2(nullFds[0],1)
		#os.dup2(nullFds[1],2)
		chrLinesDict[chrom]["first"] = int(subprocess.check_output(getChrFirstLineCmd, shell = True))
		chrLinesDict[chrom]["last"] = int(subprocess.check_output(getChrLastLineCmd, shell = True))
		#print("Done\nFirst:%d\nLast:%d" % (chrLinesDict[chrom]["first"], chrLinesDict[chrom]["last"]))
		#sys.stdout.flush()
	return chrLinesDict

parser = argparse.ArgumentParser(description = "Create a graph representing overlaps of genomic regions from a BED file.")
parser.add_argument("--regionBed", help = "BED file of genomic regions to be analysed", dest = "regionBedFileName", type = str, required = True)
parser.add_argument("--namesFile", help = "File containing the names of the different kinds of region. If not given, a file will be created using the --regionBED file.", type = str, dest = "uniqNamesFileName", default = "NONE")
parser.add_argument("--maxGap", help = "Largest gap between regions that is still counted as an overlap.", default = 100, type = int, dest = "maxGap")
parser.add_argument("--outJSON", help = "File to dump JSON string of results matrix to.", default = "json_dump", type = str, dest = "outJson")
parser.add_argument("--numProcs", help = "Number of parallel processes to run.", default = 2, type = int, dest = "numProcs")

argDict = vars(parser.parse_args())
regionBedFileName = argDict["regionBedFileName"]
uniqNamesFileName = argDict["uniqNamesFileName"]
maxGap = argDict["maxGap"]
outJson = argDict["outJson"]
numProcs = argDict["numProcs"]

if uniqNamesFileName != "NONE":
	uniqNamesFile = open(uniqNamesFileName, 'r')
else:
	print("No names file provided. Creating file ...")
	uniqNamesFileName = getNames(regionBedFileName)
	print("Opening file ...")
	uniqNamesFile = open(uniqNamesFileName, 'r')
	print("Done")

uniqNames = []
for line in uniqNamesFile:
	uniqNames.append(line.strip())

chrLinesDict = getChrLines(regionBedFileName)

pool = mp.Pool(processes = numProcs)
results = [pool.apply_async(createMatrix, args = (regionBedFileName, uniqNames, d["first"], d["last"])) for d in chrLinesDict.values()]
matrices = [p.get() for p in results]
combined = combineDicts(matrices)
normedCounts = normaliseMatrix(combined)

with open(outJson, 'w') as f:
	f.write(json.dumps(normedCounts, sort_keys = True, indent = 4, separators = (',', ':')))






























