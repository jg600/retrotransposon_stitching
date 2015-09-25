import sys, re, subprocess, linecache, json
from math import *
import numpy as np
import matplotlib.pyplot as plt

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
print("Done)
n = 1
print("Starting to get overlaps ..."),
while n <= numLines:
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
					total += 1
		else:
			keepGoing = False
	
	behind = n
	keepGoing = True

	while keepGoing:
		behind -= 1
		if behind >= 1:
			newLineList = re.split(re.compile('\s+'), linecache.getline(regionBedFileName, behind).strip())
			newDict = {'chrom':newLineList[0], 'start':int(newLineList[1]), 'end':int(newLineList[2]), 'name':newLineList[3]}
		
			if newDict['chrom'] != baseDict['chrom']:
				keepGoing = False
			else:
				if newDict['end'] < baseDict['start'] - 100:
					keepGoing = False
				else:
					countMatrix[baseDict['name']][newDict['name']] += 1
					total += 1
		else:
			keepGoing = False

	if total == 0:
		countMatrix[baseDict['name']]['NONE'] += 1
	n += 1#

#print(json.dumps(countMatrix, sort_keys = True, indent = 4, separators = (',', ':')))
print("Count matrix finished")	

'''
rowNames = countMatrix.keys()
rowNames.sort()
colNames = rowNames + ['NONE']

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(len(rowNames)):
	for j in range(len(colNames)):
		count = countMatrix[rowNames[i]][colNames[j]]
		plt.plot(i,j,marker = 's', color = str(1-(count/10.0)), linewidth = 0, markersize = 10)
ax.set_aspect('equal')
plt.show()
'''

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































