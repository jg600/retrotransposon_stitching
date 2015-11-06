import sys, re, subprocess, linecache, json
from math import *

regionBedFileName = sys.argv[1]
uniqNamesFile = open(sys.argv[2], 'r')
maxGap = int(sys.argv[3])

print("Getting unique names ..."),
uniqNames = []

for name in uniqNamesFile:
	uniqNames.append(name.strip())
print("Done")

print("Initialising matrix of counts ..."),
countMatrix = {name:{name:0 for name in (uniqNames+["NONE"])} for name in uniqNames}
print("Done")

print("Counting lines in BED file ..."),
numLines = int(subprocess.Popen("wc -l %s | awk '{print $1}'" % sys.argv[1], shell = True, stdout = subprocess.PIPE).stdout.read())
print("Done")

n = 1
print("Starting to get overlaps ...")
sys.stdout.flush()

print('['+' '*50+']'),
print('\r'),
percent = 0.0
i = 0.0

for n in range(1,numLines+1):

	percent = 100.0*(n/numLines)
	numBars = int(floor(percent/2))
	print('\r'),
	print('['+'|'*numBars+' '*(50-numBars)+']'),
	sys.stdout.flush()

	baseLineList = re.split(re.compile('\s+'), linecache.getline(regionBedFileName, n).strip())
	baseDict = {'chrom':baseLineList[0], 'start':int(baseLineList[1]), 'end':int(baseLineList[2]), 'name':baseLineList[3]}
	
	ahead = n
	total = 0
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
					countMatrix[newDict['name']][baseDict['name']] += 1
					total += 1
		else:
			keepGoing = False

	if total == 0:
		countMatrix[baseDict['name']]['NONE'] += 1

print("\nCount matrix finished")	

#Now we have made the matrix of counts, we want to normalise it by the number of occurrences of
#each classification, so divide each element by the total of its row
print("Normalising counts ..."),
normedCounts = {rowName:{} for rowName in countMatrix.keys()}

for rowName in countMatrix.keys():
	rowSum = float(sum(countMatrix[rowName].values()))
	normedCounts[rowName] = {colName:countMatrix[rowName][colName]/rowSum for colName in countMatrix[rowName].keys()}
print("Done")

with open("json_dump", 'w') as f:
	f.write(json.dumps(normedCounts, sort_keys = True, indent = 4, separators = (',', ':')))

































