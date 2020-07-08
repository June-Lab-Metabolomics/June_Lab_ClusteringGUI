import numpy as np
import statistics as stat
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from matplotlib import pyplot as plt
import pandas as pd


        
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
def normStandardize(data):
    ##********** This function has not been tested yet****************
    #The mean for a normalized data set should always be 1 for all of the metabolites.
    mean = 1

    #initialize the needed arrays
    dataCur = np.zeros(data.shape[0])
    #stdCur = np.zeros(data.shape[0])

    #standardize each row of data
    for i in range(data.shape[0]):
        residuals = 0
        #calculate the residuals
        for j in range(data.shape[0]):
            #sum up the residuals
            residuals += (data[i]-mean)**2
        #calculate the standard deviation
        std_data = (residuals/(data.shape[0]-1))**0.5

        #Input the data into an array of standardized values
        dataCur[i] = (data[i]-mean)/std_data

    return dataCur


#creating function that creates the dendrogram
def create_dendrogram(data, **kwargs):
    #Create the linkage matrix
    linkageOut = linkage(data,'ward')
    #create the dendrogram
    dendrogram(linkageOut)

#initialize the plot
def plotting():
    #create a figure window
    plt.title('Medians Clustering')
    plt.xlabel('Clustered Metabolites')
    plt.show()

    