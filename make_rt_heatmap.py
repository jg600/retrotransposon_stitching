import sys, re, subprocess, json
from math import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster import hierarchy

jsonFile = open(sys.argv[1], 'r')

print("Loading json file ..."),
probDict = json.load(jsonFile)
print("Done")

sortedRowNames = sorted(probDict.keys())
sortedColNames = sorted(probDict[sortedRowNames[0]].keys())

probArray = np.array([[probDict[rowName][colName] for colName in sortedRowNames] for rowName in sortedRowNames])

probDf = pd.DataFrame(probArray, columns = sortedRowNames, index = sortedRowNames)

'''
fig = plt.figure()

ax = fig.add_subplot(111)

cax = ax.matshow(probDf, interpolation='nearest', cmap='hot_r')
fig.colorbar(cax)

ax.set_xticklabels([''] + list(probDf.columns))
ax.set_yticklabels([''] + list(probDf.index))

plt.show()
'''

rowDist = pd.DataFrame(squareform(pdist(probDf, metric='euclidean')), columns=sortedRowNames, index=sortedRowNames)

rowClusters = linkage(pdist(probDf, metric='euclidean'), method='complete')

hierarchy.set_link_color_palette(['black'])

fig = plt.figure(figsize = (8,8))
axd = fig.add_axes([0.09,0.1,0.2,0.6])

rowDendr = dendrogram(rowClusters, orientation = 'right', color_threshold = np.inf,)
dfRowClust = probDf.ix[rowDendr['leaves'][::-1]]
print(rowDendr['leaves'])

axd.set_xticks([])
axd.set_yticks([])

for i in axd.spines.values():
        i.set_visible(False)


axm = fig.add_axes([0.26,0.1,0.6,0.6]) # x-pos, y-pos, width, height
cax = axm.matshow(probDf, interpolation='nearest')
fig.colorbar(cax)
plt.xticks(range(0, len(list(dfRowClust.columns))))
axm.set_xticklabels(list(dfRowClust.columns), rotation = 90, fontsize = 5, visible = True)
axm.yaxis.tick_right()
axm.set_yticklabels([''] + list(dfRowClust.index), fontsize = 5, visible = True)

plt.show()


#http://nbviewer.ipython.org/github/rasbt/pattern_classification/blob/master/clustering/hierarchical/clust_complete_linkage.ipynb

#Heatmaps in R:
#http://sebastianraschka.com/Articles/heatmaps_in_r.html








