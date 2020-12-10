import numpy as np
import statistics as stat
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist,squareform
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
import pandas as pd
from tkinter import filedialog
import logging
import time


        
#Standardizing the data that is input to python.
def standardize(data):
    #I would like to add the ability to check for the number of rows or columns that are given in a numpy array.
    #find the mean and standard deviation of the given row(eventually will need to move this to allow for the entire table to input.)
    mean_data = stat.mean(data)
    std_data = stat.stdev(data)
    dataCur = np.zeros(data.shape[0])

    for i in range(data.shape[0]):
        dataCur[i] = (data[i]-mean_data)/std_data

    return dataCur

#standardizing data after it has been normalized to a specific biological profile
def normStandardize(data,leaves):
    ##********** This function has not been tested yet****************
    #The mean for a normalized data set should always be 1 for all of the metabolites.
    mean = data[leaves[0]]

    #initialize the needed arrays
    dataCur = np.zeros(data.shape[0])

    #Calculate the residuals of the data
    residuals = 0
    for j in leaves:
        residuals += (data[j]-mean)**2
    #calculate the standard deviation
    std_data = (residuals/(data.shape[0]-1))**0.5
    #standardize each row of data
    for i in leaves:
        #Input the data into an array of standardized values
        dataCur[i] = (data[i]-mean)/std_data

    return dataCur


#creating function that creates the dendrogram
def create_dendrogram(data, norm=1,link='ward',dist='euclidean'):
    #Create the figure that will plot the results
    fig = plt.figure(figsize=(8,8))
    metabDendAxes =[0.05,0.05,0.1,0.8]
    metabAxes = fig.add_axes(metabDendAxes)
    #Create the linkage matrix
    linkageMetabOut = linkage(data,link,dist)
    #create the dendrogram
    metaboliteDend = dendrogram(linkageMetabOut, orientation='left',  no_labels=True)
    metaboliteDendLeaves = metaboliteDend['leaves']

    #tranpose the data and then run an analysis on the groups.
    groupCluster = np.transpose(data)
    #create axes for the group clusters
    groupAxes = [0.15, 0.85, 0.8, 0.1]
    fig.add_axes(groupAxes,xticks=[4,5,6,7,0,1,2,3])

    #Create a linkage matrix for the data
    linkageGroupOut = linkage(groupCluster,link,dist)
    #create the dendrogram of the output
    groupDend = dendrogram(linkageGroupOut,orientation='top')
    #get the leaves of the dendrogram to allow for the appropriate regrouping of the data
    groupDendLeaves = groupDend['leaves']

    #Rework the data to create the clustergram
    dataFinal = np.zeros((data.shape[0],data.shape[1]))
    for i in range(data.shape[1]):
        #rearranging the data for heatmap
        for j in range(data.shape[0]):
            #going through the metabolites
            dataFinal[j,i] = data[metaboliteDendLeaves[j],groupDendLeaves[i]]



    if norm == 0:
        #create the axes in which the heatmap will be mapped upon
        heatmapAxes = [0.15, 0.05, 0.8, 0.8]
        heatmapAxes = fig.add_axes(heatmapAxes)
        heatmapAxes.matshow(dataFinal,aspect ='auto',origin='left')

    elif norm == 1:
        for i in range(dataFinal.shape[0]):
            data[i,:] = normStandardize(data[i,:],groupDendLeaves)
        #create array that will store final version of reorganized data.
        dataFinalNorm = np.zeros((data.shape[0],data.shape[1]))

        for i in range(data.shape[1]):
            #rearranging the data for heatmap
            for j in range(data.shape[0]):
                #going through the metabolites
                dataFinalNorm[j,i] = data[metaboliteDendLeaves[j],groupDendLeaves[i]]
         
        #create the axes in which the heatmap will be mapped upon
        heatmapAxes = [0.15, 0.05, 0.8, 0.8]
        heatmapAxes = fig.add_axes(heatmapAxes)  
        #output the normalized heatmap 
        heatmapAxes.matshow(dataFinalNorm,aspect='auto',origin='left',cmap="hot")
        
def clustConnect(clusters,curConnections,mstCur1Connections,mstCur2Connections):
    print('In progress')
    

#initialize the plot
def plotting():
    plt.xlabel('Clustered Metabolites')
    plt.savefig('Clustergram.png')
    plt.show()


def Validate(data,dists,num_groups):
    #grab the input dictionary size
    clusterings = len(data)
    val_index = np.zeros((2,clusterings))
    initTime = time.time()

    for i in range(clusterings):
        #grab the current set of metabolite clusters
        curClusters = data[i]

        #from the current clusters determine the length in order to determine the next step
        curClustersLength = len(curClusters)

        #sum of intra cluster distances
        sumIntra = 0

        #create a numpy array for that contains the cluster centers for calculation of the inter cluster distance.
        centersNum = np.zeros((curClustersLength,num_groups))
        for j in range(curClustersLength):

            #pull out the current cluster of metabolites
            cluster = curClusters[j]
            
            #check for instances of the clusters imbedded in the dictionary for each clustering outcome
            if isinstance(cluster, list):
                #check the length of the cluster and then pull in the values for each of the clustered metabolites
                lengthList = len(cluster)
                clustCoordinates = np.zeros((lengthList,num_groups))

                for k in range(lengthList):
                    #grab the cluster coordinates for each metabolite clustered with this particular round of clusters
                    clustCoordinates[k,:] = dists[cluster[k]]

                #Create a numpy array for the coordinates of the center of the current cluster
                center = np.zeros((1,num_groups))

                for m in range(num_groups):
                    #find the mean of each group in the cluster
                    center[0,m] = stat.mean(clustCoordinates[:,m])
                #update dictionary of the cluster centers for later calculation of inter-cluster distances
                centersNum[j,:] = center

                #initialize an array that stores the values that will be sent to the pdist function
                curDistIntra = np.zeros((2,num_groups))

                #calculate the intra_cluster distance for each list of metabolites in the 
                for k in range(lengthList):
                    #grab the first value from the list find its distance and add it to the sum of the Intra cluster distances
                    curMetab = clustCoordinates[k,:]
                    curDistIntra[0,:] = curMetab
                    curDistIntra[1,:] = center
                    sumIntra += pdist(curDistIntra)

            elif isinstance(cluster, np.integer) or isinstance(cluster,int):
                #find the center and put it into the dictionary
                center = np.zeros((1,num_groups))
                center[0,:] = dists[cluster]
                centersNum[j,:] = center

        #calculate the average compactness of the clusters
        print(sumIntra)
        intraDist = sumIntra/(clusterings+1)
        #print(intraDist)
        
        # find the distance between the centers
        centerDists = pdist(centersNum)
        lenCenters = len(centerDists)
        centerDists = squareform(centerDists)
        centerDists = csr_matrix(centerDists)
        centerDistsInd = centerDists.nonzero()

        dataMST = np.zeros([lenCenters,1])

        for k in range(lenCenters):
            #input the values of each matrix element to the numpy array for saving to csv file.
            curVal0 = centerDistsInd[0][k]
            curVal1 = centerDistsInd[1][k]
            dataMST[k,0] = centerDists[curVal0,curVal1]
        

        #calculate the inter-cluster distance for the current cluster set-up
        if len(dataMST[:,0]) > 0:
            interDist = np.min(dataMST[:,0])
            #calculate the validation index
            val_index[0,i] = intraDist/interDist
            val_index[1,i] = clusterings - (i)

        elif len(dataMST[:,0]) == 0:
            interDist = 0
            #set the validation index to a large number since denominator would be zero in current config.
            val_index[0,i] = 1000
            val_index[1,i] = clusterings - (i)

        if i == 0:
            firstTime = time.time() -initTime
            # print(firstTime)
        if i%10 == 0 or i ==clusterings-1:
            curTime = time.time()
            runTime = (curTime - initTime)/3600
            # print(i)
            # print(runTime)

            
    return val_index













