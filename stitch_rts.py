import sys, re, json, subprocess, linecache
from math import *
import numpy as np

regionBedFileName = sys.argv[1]
uniqNamesFile = open(sys.argv[2], 'r')
jsonFile = open(sys.argv[3], 'r')
maxGap = int(sys.argv[4])
outputGtf = open(sys.argv[5], 'wa')

print("Loading json file ..."),
probMatrix = json.load(jsonFile)
print("Done")

clusterList = []

print("Counting number of lines in BED file ..."),
numLines = int(subprocess.Popen("wc -l %s | awk '{print $1}'" % sys.argv[1], shell = True, stdout = subprocess.PIPE).stdout.read())
print("Done")

print("Starting to create clusters ...")
n = 1
while n <= numLines:
	print("n=%d" % n)
	startLine = re.split(re.compile('\s+'), linecache.getline(regionBedFileName, n).strip())
	startDict = {'chrom':startLine[0], 'start':int(startLine[1]), 'end':int(startLine[2]), 'name':startLine[3], 'strand':startLine[4]}
	startDict["family"], startDict["type"] = startLine[3].split("/")
	
	cluster = [startDict]
	keepGoing = True
	next = n
	while keepGoing and (next < numLines):
		next += 1
		nextLine = re.split(re.compile('\s+'), linecache.getline(regionBedFileName, next).strip())
		nextDict = {'chrom':nextLine[0], 'start':int(nextLine[1]), 'end':int(nextLine[2]), 'name':nextLine[3], 'strand':startLine[4]}
		nextDict["family"], nextDict["type"] = nextLine[3].split("/")
		
		if nextDict['chrom'] != startDict['chrom']:
			print("Chromosomes don't match")
			keepGoing = False
		elif nextDict['family'] != startDict['family']:
			print("Families don't match")
			keepGoing = False
		elif nextDict['start'] > (startDict['end'] + maxGap):
			print("Neighbouring RTs are too far apart")
			keepGoing = False
		else:
			if probMatrix[startDict['name']][nextDict['name']] > 0:
				print("Types associate with p = %.8f, adding to cluster" % probMatrix[startDict['name']][nextDict['name']])
				cluster.append(nextDict)
			else:
				print("Types do not associate")
				keepGoing = False
	if next == numLines:
		break
	n = next
	clusterList.append(cluster)
print("Clusters finished")

print("Writing results to GTF file ..."),
for cluster, feat_id in zip(clusterList, range(len(clusterList))):
	featureChrom = list(set([rt['chrom'] for rt in cluster]))[0]
	source = "RepeatMasker"
	feature = list(set([rt['family'] for rt in cluster]))[0]
	featureStart = str(min([rt['start'] for rt in cluster]))
	featureEnd = str(min([rt['end'] for rt in cluster]))
	featureAttribute = "stitched_retrotransposon"
	strandList = list(set([rt['strand'] for rt in cluster]))
	
	if len(strandList) == 1:
		featureStrand = strandList[0]
	else:
		featureStrand = '.'
	
	outputGtf.write('\t'.join([featureChrom,source,feature,featureStart,featureEnd,'.',featureStrand,'.','feat_id:"%d"; type:"stitched_retrotransposon"' % feat_id])+"\n")
	for section, sectionNum in zip(cluster, range(len(cluster))):
		outputGtf.write('\t'.join([section['chrom'],source,feature,str(section["start"]),str(section["end"]),'.',section["strand"],'.','feat_id:"%d"; type:"%s"' % (feat_id, section["type"])]))
		if (sectionNum != len(cluster)) and (feat_id != len(clusterList)):
			outputGtf.write("\n")

print("Done")
















