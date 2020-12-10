import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist,squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import logging
#import fpdf
import os
import statistics as stat
import GuiBackground as GB
from tkinter import filedialog



#ask the user for the file input
file = filedialog.askopenfilename()



metab_data = pd.read_excel(file, sheet_name='Medians')

#finding the number of groups in the metabolomics study to allow for proper clustering of the results
num_groups = metab_data.shape[1] - 2

#creating a numpy array that is the size of the data that is being read in.
data = np.zeros((metab_data.shape[0],metab_data.shape[1]-2))

for i in range(num_groups):
    #create the appropriate string to add to the group name
    num = str(i + 1)
    g_name = "M" + num
    #grab the Medians for each group
    medianCur = metab_data[g_name]

    data[:,i] = medianCur

#Standardize the data before clustering the results
for i in range(metab_data.shape[0]):
    data[i,:] = GB.standardize(data[i,:])

#Create the figure that will plot the results
fig = plt.figure(figsize=(8,8))
link ='ward'
dist ='euclidean'
#Create the linkage matrix
linkageMetabOut = linkage(data,link,dist)
print(linkageMetabOut)
#create the dendrogram
metaboliteDend = dendrogram(linkageMetabOut, truncate_mode=None, orientation='top')
metaboliteDendLeaves = metaboliteDend['ivl']
plt.show()

