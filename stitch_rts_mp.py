import sys, re, subprocess, linecache, json, argparse
from math import *
import multiprocessing as mp

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
	
class region:
	def __init__(self, chrom, start, end, name):
		self.chrom = chrom
		self.start = int(start)
		self.end = int(end)
		self.name = name
		
class stitchedRegion:
	def __init__(self, name):
		self.regionList = []
		self.chrom = ''
		self.start = -1
		self.end = -1
		self.name = name
		
	def addRegion(self, region):
		self.regionList.append(region)
		
	def calcBounds(self):
		self.start = min(r.start for r in self.regionList)
		self.end= max(r.end for e in self.regionList)
	
def stitchChr(bedFileName, firstLine, lastLine):
	outputList = []
	
	





























