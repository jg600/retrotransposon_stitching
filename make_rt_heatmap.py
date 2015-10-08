import sys, re, subprocess, json
from math import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage, dendrogram

jsonFile = open(sys.argv[1], 'r')

print("Loading json file ..."),
probDict = json.load(jsonFile)
print("Done")

sortedRowNames = sorted(probDict.keys())
sortedColNames = sorted(probDict[sortedRowNames[0]].keys())

probArray = np.array([[probDict[rowName][colName] for colName in sortedColNames] for rowName in sortedRowNames])

probDf = pd.DataFrame(probArray, columns = sortedColNames, index = sortedRowNames)

rowDist = pd.DataFrame(squareform(pdist(probDf, metric='euclidean')), columns=sortedRowNames, index=sortedRowNames)


rowClusters = linkage(pdist(probDf, metric='euclidean'), method='complete')
rowDendr = dendrogram(rowClusters, labels=sortedRowNames)
plt.show()


#http://nbviewer.ipython.org/github/rasbt/pattern_classification/blob/master/clustering/hierarchical/clust_complete_linkage.ipynb

#Heatmaps in R:
#http://sebastianraschka.com/Articles/heatmaps_in_r.html
