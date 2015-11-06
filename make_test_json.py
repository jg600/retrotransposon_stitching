import sys, re, json

namesFile = open(sys.argv[1], 'r')
outputJsonFileName = sys.argv[2]

uniqNames = []
for line in namesFile:
	uniqNames.append(line.strip())
	
namesFile.close()

outputDict = {name:{name:0.0 for name in uniqNames + ["NONE"]} for name in uniqNames}

with open(outputJsonFileName, 'w') as f:
	f.write(json.dumps(outputDict, sort_keys = True, indent = 4, separators = (',', ':')))
