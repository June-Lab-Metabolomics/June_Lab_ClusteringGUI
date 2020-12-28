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

def fileCheck(file=''):
    '''
    Check that the selected file is of the appropriate file extension and is able to be read in. 
    '''
    #log that the user called the Create Clustergram function
    logging.warning(': User checking file.')
    #ask the user to select the file that they would like to create a clustergram for.
    if file =='':
        file = filedialog.askopenfilename()
    
    if file == '': 
        logging.error(': Failed to select a file')
        return

    #Open the excel file that the user to looking to use for their clustering analysis
    try:
        data = pd.read_excel(file, sheet_name='Medians')
        del(file)
    except:
        logging.error(': Failed to read in the excel sheet check the sheet name is Medians')
        return
    logging.info(': User file opened and submitted to the function.')
    return data
        
#Standardizing the data that is input to python.
def standardize(data):
    '''
    Standardize the input data for best clustering results.
    '''
    logging.info(': Standardizing data.')
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
    '''
    Normalize input data that has been standardized for normalized clustergrams.
    '''
    logging.info(': Normalizing strandardized data.')
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

#dendrogram function
def create_dendrogram(data, norm=1,link='ward',dist='euclidean',func='clustergram'):
    '''
    Create dendrogram for either the ensemble or clustergram functions
    '''
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
    if func == 'clustergram':
        groupCluster = np.transpose(data)
        print('clustergram work')
    elif func == 'ensemble':
        groupCluster = data
        print('ensemble work')
    else:
        logging.error(': Inappropriate entry for the create_dendrogram function. Check docs.')

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
        heatmapAxes.matshow(dataFinal,aspect ='auto',origin='upper')

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
        heatmapAxes.matshow(dataFinalNorm,aspect='auto',origin='upper',cmap="hot")

def cooccurrence(data):
    '''
    Creation of the cooccurrence matrix for the determination of the number times each set of metabolites is clustered together
    in a set of N clusterings. 
    '''
    logging.info(': Creating the co-occurrence matrix.')

    #find the number of metabolites
    numMetabolites = data.shape[0]

    #create co-occurrence matrix of zeros
    coOcc = np.zeros((numMetabolites,numMetabolites))

    for i in range(numMetabolites):
        #populate the diagonal matrix with ones
        coOcc[i,i] = 1

    return coOcc

def popCooccurrence(clusters,coOcc,numClusterings):
    '''
    Populate the cooccurrence matrix with the connections for each matrix based upon the number of clusterings that occur. 
    '''
    logging.info(': Populating the co-occurrence matrix.')
    dictKeys = list(clusters.keys())
    dictKeys = len(dictKeys)

    #weight to be added to the matrix
    weight = 1/numClusterings

    #populate the coOccurence matrix with the appropriate values based upon the keys using the binary operation +=
    for i in range(dictKeys):
        #get the first entry from the dictionary 
        curCluster = clusters[i]

        #check for a list of metabolites clustered together.
        if isinstance(curCluster, list):
            #determine the length of the list
            lenList = len(curCluster)
            
            for j in range(lenList):
                #populate the Cooccurrence matrix
                remList = lenList-j
                if remList > 1:
                    for k in range(j+1,lenList):
                        #add the weight to the coOcc matrix. 
                        coOcc[curCluster[j],curCluster[k]] += weight
                        coOcc[curCluster[k],curCluster[j]] += weight
    return coOcc

def clustConnectLink(linkageCheck):
    '''
    Determine the connections from the scipy linkage function output.
    '''
    logging.info(': Starting to determine metabolite clusters.')
    #determine the clusters using a new alogorithm based on linkage method.
    #create a dictonary to store the linkage functions. 
    metabolites = np.zeros((linkageCheck.shape[0]+1,2))
    #fill the array with the metabolite identifiers
    for i in range(linkageCheck.shape[0]+1):
        metabolites[i,0] = i
        metabolites[i,1] = i

    metabolites = metabolites.astype(int)
    linkageCheck = linkageCheck.astype(int)
    #set the limit for the metabolites such that if the metabolites match up with the 
    metabLimit = linkageCheck.shape[0]
    limit = linkageCheck.shape[0]

    #create dictionary to store the clusters as they are created.
    clusters = {}

    #create an empty dictionary to store all of the cluster configurations.
    validationClusters = {}
    for i in range(linkageCheck.shape[0]):
        #check the linkage connections to determine the appropriate 
        curCon1 = linkageCheck[i,0]
        curCon2 = linkageCheck[i,1]

        curCon1 = int(curCon1)
        curCon2 = int(curCon2)

        if i == 0:
            #place first connection into the list 
            firstConnect = [curCon1, curCon2]
            metabolites[curCon1,1] = limit + 1
            metabolites[curCon2,1] = limit + 1
            limit += 2

            for j in range(linkageCheck.shape[0]-(i+1)):
                clusters.update({j:j})

            #populate the first part of the dictionary with the first list containing the clustered metabolites
            clusters[0] = firstConnect 

            #search the array for values of the array which are equal to the one plus the number of metabolites studied
            unClustered = np.where(metabolites[:,1] != (metabLimit+1))
            for j in range(1,linkageCheck.shape[0]):
                clusters[j] = unClustered[0][j-1]
  
        else:
            #save the previous clusters dictionary prior to deleting it. 
            clusterPrevious = clusters
            del(clusters)
            clusters = {}

            for j in range(linkageCheck.shape[0]-(i+1)):
                clusters.update({j:j})

            curCon1 = linkageCheck[i,0]
            curCon2 = linkageCheck[i,1]

            curCon1 = int(curCon1)
            curCon2 = int(curCon2)

            if curCon1 > metabLimit and curCon2 <= metabLimit:
                #search the second column of the list for the appropriate value.
                oneCon = np.where(metabolites[:,1] == curCon1)
                connect1 = oneCon[0][0]
                curCon1 = int(connect1)
                curCon2 = int(curCon2)
            elif curCon1 > metabLimit and curCon2 > metabLimit:
                #search the second column of the reference list for the appropriate value.
                oneCon = np.where(metabolites[:,1] == curCon1)
                connect1 = oneCon[0][0]
                curCon1 = int(connect1)
                twoCon = np.where(metabolites[:,1] == curCon2)
                connect2 = twoCon[0][0]
                curCon2 = int(connect2)
            elif curCon1 <= metabLimit and curCon2 > metabLimit:
                #search the second column of the list for the appropriate value. 
                twoCon = np.where(metabolites[:,1] == curCon2)
                connect2 = twoCon[0][0]
                curCon2 = int(connect2)
                curCon1 = int(curCon1)

            curCon1Connect = 0;
            curCon2Connect = 0;
            unchanged = []

            previousKey = list(clusterPrevious.keys())
            previousKey = len(previousKey)

            for k in range(previousKey):
                #Determine matching metabolites in the dictionary.
                curCheck = clusterPrevious[k]

                if isinstance(curCheck, list):
                    #check list for the first connection
                    curCon1Check = curCon1 in curCheck

                    #check list for the second connection
                    curCon2Check = curCon2 in curCheck

                    if curCon1Check == True and curCon2Check == False:
                        curCon1Connect = 1
                        curCon1Location = k
                    elif curCon1Check == False and curCon2Check == True:
                        curCon2Connect = 1
                        curCon2Location = k
                    elif curCon1Check == False and curCon2Check == False:
                        unchanged.append(k)
                    elif curCon1Check == True and curCon2Check == True:
                        logging.warning(': Issue clustering the data a duplication has been discovered.')


                elif isinstance(curCheck, np.integer) or isinstance(curCheck, int):
                    #check list for the first connection
                    curCon1Check = curCon1 == curCheck

                    #check list for the second connection
                    curCon2Check = curCon2 == curCheck

                    if curCon1Check == True and curCon2Check == False:
                        curCon1Connect = 1
                        curCon1Location = k
                    elif curCon1Check == False and curCon2Check == True:
                        curCon2Connect = 1
                        curCon2Location = k
                    elif curCon1Check == False and curCon2Check == False:
                        unchanged.append(k)
                    elif curCon1Check == True and curCon2Check == True:
                        loggging.warning(': Issue clustering the data a duplication has been discovered.')


                else:
                    logCheck = type(curCheck)
                    logging.error(logCheck)
                    logging.error(': Inappropriate data type for the for clustering ID algorithm.')
                    return

                if curCon1Connect == 1 and curCon2Connect == 1 and k == previousKey-1:
                    for m in range(1,len(unchanged)+1):
                        #cluster the appropriate remaining clusters together
                        clusters[m] = clusterPrevious[unchanged[m-1]]

                    newConnect1 = clusterPrevious[curCon1Location]
                    newConnect2 = clusterPrevious[curCon2Location]
                    
                    if isinstance(newConnect1,list) and isinstance(newConnect2, list):
                        # print('passed')
                        newConnect = newConnect1 + newConnect2
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,list) and isinstance(newConnect2, np.integer):
                        # print('passed')
                        newConnect = newConnect1
                        intList = newConnect[:] + [newConnect2]
                        clusters[0] = intList
                        for m in intList:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,list) and isinstance(newConnect2, int):
                        # print('passed')
                        newConnect = newConnect1
                        intList = newConnect[:] + [newConnect2]
                        clusters[0] = intList
                        for m in intList:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2, list):
                        # print('passed')
                        newConnect = newConnect2
                        intList = newConnect[:] + [newConnect1]
                        clusters[0] = intList
                        for m in intList:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,int) and isinstance(newConnect2, list):
                        # print('passed')
                        newConnect = newConnect2
                        intList = newConnect[:] + [newConnect1]
                        clusters[0] = intList
                        for m in intList:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2,np.integer):
                        # print('passed')
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,int) and isinstance(newConnect2,np.integer):
                        # print('passed')
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2,np.integer):
                        # print('passed')
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,int) and isinstance(newConnect2,int):
                        # print('passed')
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1

        validationClusters.update({i:clusters})
    logging.info(': Success! Metabolite clusters determined.')
    return validationClusters

def clustConnect(dataMST,mstOutNp):
    '''
    Clustering connections determination from the minimum spanning tree output. 
    '''
    logging.info(': Determining the metabolite cluster for MST.')
    #Create an empty dictionary that will contain the clusters
    clusters = {}

    #create an empty dictionary that allows us to save the clusters for validation.
    validationClusters = {}

    # create initial list of metabolites that serves as the initial list of metabolites that will be clustered.
    metabolites = np.ones((dataMST.shape[0]+1,1))

    for i in range(dataMST.shape[0]):
        #pull out the connections for the current iteration
        curCon1 = mstOutNp[i,0]
        curCon2 = mstOutNp[i,1]

        curCon1 = int(curCon1)
        curCon2 = int(curCon2)

        if i == 0:
            #convert the metabolites to string for easier comparison
            firstConnect = [curCon1, curCon2]

            #set the metabolites equal to zero in the initial metabolite list
            metabolites[curCon1,0] = 0 
            metabolites[curCon2,0] = 0

            #create dictionary of clusters 
            for j in range(dataMST.shape[0]-(i+1)):
                clusters.update({j:j})
            
            #find the metabolites that are ones and were not clustered
            unClustered = np.where(metabolites == 1)

            #input the connection
            clusters[0] = firstConnect
            

            for j in range(1,dataMST.shape[0]):
                #input the cluster values into the dictionary
                clusters[j] = unClustered[0][j-1]
        else:
            #save the previous dictionary 
            clusterPrevious = clusters
            del(clusters)
            clusters = {}

            #create a new dictionary 
            for j in range(dataMST.shape[0]-(i+1)):
                clusters.update({j:j})

            #grab the latest connections
            curCon1 = mstOutNp[i,0]
            curCon2 = mstOutNp[i,1]

            curCon1 = int(curCon1)
            curCon2 = int(curCon2)

            curCon1Connect = 0;
            curCon2Connect = 0;
            unchanged = []

            previousKey = list(clusterPrevious.keys())
            previousKey = len(previousKey)

            for k in range(previousKey):
                #Determine if any of the new clustered meatbolites 
                curCheck = clusterPrevious[k]

                if isinstance(curCheck, list):
                    #check list for the first connection
                    curCon1Check = curCon1 in curCheck

                    #check list for the second connection
                    curCon2Check = curCon2 in curCheck

                    if curCon1Check == True and curCon2Check == False:
                        curCon1Connect = 1
                        curCon1Location = k
                    elif curCon1Check == False and curCon2Check == True:
                        curCon2Connect = 1
                        curCon2Location = k
                    elif curCon1Check == False and curCon2Check == False:
                        unchanged.append(k)
                    elif curCon1Check == True and curCon2Check == True:
                        logging.warning(': Issue clustering the data a duplication has been discovered.')


                elif isinstance(curCheck, np.integer) or isinstance(curCheck, int):
                    #check list for the first connection
                    curCon1Check = curCon1 == curCheck

                    #check list for the second connection
                    curCon2Check = curCon2 == curCheck

                    if curCon1Check == True and curCon2Check == False:
                        curCon1Connect = 1
                        curCon1Location = k
                    elif curCon1Check == False and curCon2Check == True:
                        curCon2Connect = 1
                        curCon2Location = k
                    elif curCon1Check == False and curCon2Check == False:
                        unchanged.append(k)
                    elif curCon1Check == True and curCon2Check == True:
                        loggging.warning(': Issue clustering the data a duplication has been discovered.')


                else:
                    logCheck = type(curCheck)
                    logging.error(logCheck)
                    logging.error(': Inappropriate data type for the for clustering ID algorithm.')
                    return

                if curCon1Connect == 1 and curCon2Connect == 1 and k == previousKey-1:
                    for m in range(1,len(unchanged)+1):
                        #cluster the appropriate remaining clusters together
                        clusters[m] = clusterPrevious[unchanged[m-1]]

                    newConnect1 = clusterPrevious[curCon1Location]
                    newConnect2 = clusterPrevious[curCon2Location]
                    
                    if isinstance(newConnect1,list) and isinstance(newConnect2, list):
                        newConnect = newConnect1 + newConnect2
                        clusters[0] = newConnect
                    elif isinstance(newConnect1,list) and isinstance(newConnect2, np.integer):
                        newConnect = newConnect1
                        intList = newConnect[:] + [newConnect2]
                        clusters[0] = intList
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2, list):
                        newConnect = newConnect2
                        intList = newConnect[:] + [newConnect1]
                        clusters[0] = intList
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2,np.integer):
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] =newConnect

        validationClusters.update({i:clusters})
    logging.info(': Success! MST clusters determined.')
    return validationClusters
    

#initialize the plot
def plotting():
    '''
    Plotting the clustergram from the above create_dendrogram function. 
    '''
    logging.info(': Plotting Clustergram!')
    plt.xlabel('Clustered Metabolites')
    logging.info(': Saving...')
    plt.savefig('Clustergram.png')
    logging.info(': Success!')
    plt.show()

#convert seconds to HH:MM:SS
def timeConverter(runTime):
    '''
    Convert run time in seconds to a measure of Hours:Minutes:Seconds.
    '''
    #immediately determine the number of full hours that were consumed
    hrs = runTime/3600
    hrs = int(hrs)

    #subtract the number of seconds that hrs took up and then determine the number of full minutes left in the
    #remaining seconds
    runTime -= (hrs*3600)
    mins = runTime/60
    mins = int(mins)
    runTime -= (mins*60)
    secs = int(runTime)
    runTime -= secs

    runTime = str(runTime)
    runTime = runTime.strip("0")
    runTime = runTime[0:3]

    if hrs < 10:
        hrsStr = '0' + str(hrs)
    else:
        hrsStr =str(hrs)

    if mins < 10:
        minsStr = "0" + str(mins)
    else:
        minsStr = str(mins)

    if secs < 10:
        secsStr = "0" + str(secs)
    else:
        secsStr = str(secs)

    runTime = hrsStr +':'+ minsStr +':'+ secsStr

    return runTime

def Validate(data,dists,num_groups):
    '''
    Determine the appropriate number of clusters for the optimal clustering of a input data set. 
    '''
    logging.info(': Starting cluster validation!')
    #grab the input dictionary size
    clusterings = len(data)

    startPoint = 0.6*clusterings 
    startPoint = int(startPoint)
    numIts = clusterings - startPoint
    print(numIts)

    numClusters = clusterings-startPoint
    numClusters = int(numClusters)
    val_index = np.zeros((2,clusterings))
    initTime = time.time()

    for i in range(startPoint,clusterings):
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
        intraDist = sumIntra/(clusterings+1)
        
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
            val_index[0,i] = 1
            val_index[1,i] = clusterings - (i)

        if i == 0:
            logging.info(': 0% completed')

        elif i%100==0:
            #calculate the amount of time the algorithm has been running 
            curTime = time.time()
            runTime = (curTime - initTime)/3600
            percent = 100*(1 - ((clusterings-i)/numIts))
            message1 = round(percent,2) 
            message1 = str(message1) + '% completed'
            logging.info(message1)
            message = str(i)+': ' +str(runTime)
            logging.info(message)   

    runTime = time.time()-initTime
    runTime = timeConverter(runTime)
    logging.info(runTime)
    
    return val_index