import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist,squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from tkinter import filedialog
import GuiBackground as GB


file = filedialog.askopenfilename()

dataRaw = pd.read_excel(file,sheet_name='Medians')


#find the number of groups in that data
num_groups = dataRaw.shape[1]-2
#standardize the data before finding the distance matrix
data = np.zeros((dataRaw.shape[0],dataRaw.shape[1]-2))

for i in range(num_groups):
    #create the appropriate string to add to the group name
    num = str(i+1)
    g_name = "M" + num
    #grab the medians  for each group
    
    medianCur = dataRaw[g_name]

    #add the medians data to the array to be clustered
    data[:,i] = medianCur

for i in range(dataRaw.shape[0]):
	data[i,:] = GB.standardize(data[i,:])


#find the distance matrix using the pairwise distance function used in the analysis of Agglomerative clustering. 
pairWise = pdist(data)

pairWise = squareform(pairWise)

#put the distance matrix into the appropriate format for creation of MST using scipy's MST capabilities
mstInput = csr_matrix(pairWise)
mstOutInd = mstInput.nonzero()


# for i in range(mstInput.shape[0]):
	#pull out the values from the third column for the 



# #create the minimum spanning tree
# mstOut = minimum_spanning_tree(mstInput)
# print(mstOut)
# mstOutInd = mstOut.nonzero()

# #create the matrix containing the various connections of the mst
# mstOutMat = mstOut.todense()

# dataMST = np.zeros([data.shape[0]-1,3])

# for i in range(data.shape[0]-1):
#     #input the values of each matrix element to the numpy array for saving to csv file.
#     dataMST[i,0] = mstOutInd[0][i]
#     dataMST[i,1] = mstOutInd[1][i]
#     dataMST[i,2] = mstOut[dataMST[i,0],dataMST[i,1]]






