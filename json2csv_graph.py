import json, re, subprocess, sys


jsonFile = open(sys.argv[1], 'r')
outputNodeFile = open(sys.argv[2], 'w')
outputEdgeFile = open(sys.argv[3], 'w')
minWeight = float(sys.argv[4])
try:
	matchNames = bool(sys.argv[5])
except IndexError:
	matchNames = False

print("Loading json file ..."),
probDict = json.load(jsonFile)
print("Done")

nodeNames = probDict.keys()
nodeNames.sort()

outputNodeFile.write("Id\n")
for n in nodeNames:
	outputNodeFile.write(n)
	if n != nodeNames[-1]:
		outputNodeFile.write("\n")
		
outputEdgeFile.write("Source;Target;Type;Weight\n")

for row in nodeNames:
	for col in nodeNames:
		if probDict[row][col] >= minWeight:
			if not matchNames:
				outputEdgeFile.write(row+";"+col+";"+"Directed" + ";" + str(probDict[row][col])+"\n")
			else:
				rowName = re.split(re.compile('/'), row)[0]
				colName = re.split(re.compile('/'), col)[0]
				if rowName != colName:
					outputEdgeFile.write(row+";"+col+";"+"Directed" + ";" + str(probDict[row][col])+"\n")
