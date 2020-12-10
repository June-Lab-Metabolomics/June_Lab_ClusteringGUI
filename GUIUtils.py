#Creating a class containing functions that will be used in GUI
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


class GUIUtils:
    
    def dataIntegrity(file):
        #log that the user called the data integrity function
        logging.warning(': User called the Data Integrity function.')

        #Read in Volcano Plot data
        try:
            volcano = pd.read_excel(file)
        except:
            logging.warning(': Failed to read in the excel file. Please put error in the Github issues tab.')
            return
            
        #grab the first row the volcano data
        check = volcano['Unnamed: 0']

        #create array that can save the fixed data and the data that did not need to be fixed
        correctedArray = np.zeros(check.shape[0])

        #search each of the volcano data rows to determine if they have double decimals.
        for i in range(check.shape[0]):
            #grab the row corresponding to the current run through the loop
            curVal = check[i]

            #reset the number of decimals to 0 before checking the string for decimal points
            decimal = 0

            #creating a string that will contain the corrected string
            corrected = ''

            #Determine if the value is a float to allow for determination of whether or not we need to check it for the appropriate value
            if isinstance(curVal,float) != True:
                #Look through the strings to find data integrity issues. 
                for j in range(len(curVal)):
                    #Check for decimals otherwise add the value to a string
                    value = curVal[j]
                    if value == '.':
                        decimal += 1
                        if decimal == 1:
                            corrected += value
                    else:
                        corrected += value
                        if j == len(curVal)-1:
                            try:
                                correctedArray[i] = float(corrected)
                            except:
                                logging.warning(': Unable to convert values to floats. Make sure all data values are only contain decimals or numberic values')
                                return
            else:
                #save the data that did not need to be corrected to the correctedArray
                correctedArray[i] = curVal

        #Replace the values in the dataframe with the appropriate values
        volcano['Unnamed: 0'] = correctedArray
        del(correctedArray,i,j,curVal,decimal,corrected,value)

        #Replace the file name with the appropriate rename
        file = file[0:len(file)-5]
        file += '_corrected.xlsx'

        #specify the file to write to
        output = pd.ExcelWriter(file)

        #write to excel file
        volcano.to_excel(output,index=False)
        del(file,volcano)
        #save the excel sheet
        output.save()

        #log that the data integrity function has been sucessfully completed. 
        logging.warning(': Data Integrity check sucessfully completed.')

    def createClustergram(norm,linkFunc,distMet):
        #log that the user called the Create Clustergram function
        logging.warning(': User called the Create Clustergram Function.')
        #ask the user to select the file that they would like to create a clustergram for.
        file = filedialog.askopenfilename()
        
        if file == '': 
            logging.warning(': Failed to select a file')
            return

        #Open the excel file that the user to looking to use for their clustering analysis
        try:
            metab_data = pd.read_excel(file, sheet_name='Medians')
            del(file)
        except:
            logging.warning(': Failed to read in the excel sheet check the sheet name is Medians')
            return

        #finding the number of groups in the metabolomics study to allow for proper clustering of the results
        num_groups = metab_data.shape[1] - 2

        #creating a numpy array that is the size of the data that is being read in.
        data = np.zeros((metab_data.shape[0],metab_data.shape[1]-2))

        for i in range(num_groups):
            #create the appropriate string to add to the group name
            num = str(i + 1)
            g_name = "M" + num
            #grab the Medians for each group
            try:
                medianCur = metab_data[g_name]
            except:
                logging.warning(': Improper column header, make sure to follow the convention of M1-MN')
                return

            #add the medians data to the array to be clustered
            data[:,i] = medianCur
        del(num_groups,num,medianCur) 
        #Standardize the data before clustering the results
        for i in range(metab_data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])
        del(metab_data)
        #create dendrogram and plot data
        GB.create_dendrogram(data,norm,link=linkFunc,dist=distMet)
        GB.plotting()
        del(data,norm,linkFunc,distMet)
        #log that you have sucessfully created the clustergram
        logging.warning(': Sucessfully created the wanted Clustergram')

    def groupMedians(file):
        #log that the user called the group medians function
        logging.warning(': User called the Group Medians function.')

        #read in the file containing all of the metabolite values
        try:
            medians = pd.read_excel(file)
        except:
            logging.warning(': Failed to read in the excel file. Please put error in the Github issues tab.')    
            return
        #get the group letters to know which groups to find the medians of
        groups = medians['mz']

        #determine the number of groups and then create a list or array of the appropriate
        #beginning and ending of each group.
        #This assumes that the groups are all of equal size which should be
        #the goal for any and all analysis. Groups with out the same sizes should be considered
        #inappropriate for analysis in this context, additionally it should be noted that statistics
        #with out the same groups sizes can lead to incorrect analysis. 
        num_groups = int(medians['Unnamed: 0'][23][0])
        factor = len(medians['Unnamed: 0'])/num_groups


        #create a numpy array containing 7 columns to allow for input of the m/z values and the groups
        mediansOut = np.zeros((medians.shape[1]-2,num_groups+1))
            

        #populate the first column of the array with the m/z values
        mediansOut[:,0] = medians.columns[2:medians.shape[1]]


        #Get the medians for each group and metabolite
        for i in range(num_groups):
            #calculate the start and end for each set of median calculations
            start = int(factor*i)
            end = int((factor*(i+1)))
            for j in range(mediansOut.shape[0]):
                #find the median for the first groups of values
                curMean = stat.median(medians[medians.columns[j+2]][start:end])
                #put medians into the appropriate table
                mediansOut[j,i+1] = curMean
        del(curMean,start,end,i,j,medians)
        #create list contains the headers for the files
        medianList = ['m/z']
        for i in range(num_groups):
            medianList.append('M' +str(i+1))

        #create dictionary that contains the data with there appropriate headers to input to a dataframe and
        #then be saved to a csv file
        medianDict = {}
        for i in range(num_groups+1):
            #input the appropriate data and key to the dictionary
            medianDict[medianList[i]] = mediansOut[:,i]
        del(mediansOut)
        #create dataframe that prepares the data to be input to a csv file
        mediansCSV = pd.DataFrame(data=medianDict)

        file = file[0:len(file)-5]
        file += '_Medians.csv'

        #specify the file that I want the program to write to.
        mediansCSV.to_csv(file,columns=medianList,index =False)
        del(file,medianList,medianDict,mediansCSV)
        #logging the completion of the group medians function
        logging.warning(': Sucessfully grouped the Medians of each group!')



    def linkageComparison(file,num_comps,linkList):
        #Log that user called linkage comparison function
        logging.warning(': User called the Linkage Comparison function.')

        #Open the excel file that the user to looking to use for their clustering analysis
        try:
            metab_data = pd.read_excel(file, sheet_name='Medians')
        except:
            logging.warning(': Failed to read in the excel file check to make sure that the sheet name is Medians.')
            return

        #finding the number of groups in the metabolomics study to allow for proper clustering of the results
        num_groups = metab_data.shape[1] - 2

        #creating a numpy array that is the size of the data that is being read in.
        data = np.zeros((metab_data.shape[0],metab_data.shape[1]-2))

        for i in range(num_groups):
            #create the appropriate string to add to the group name
            num = str(i + 1)
            g_name = "M" + num
            #grab the Medians for each group
            try:
                medianCur = metab_data[g_name]
            except:
                logging.warning(': Improper column headers, make sure you are following the convention M1-MN')
                return
            #add the medians data to the array to be clustered
            data[:,i] = medianCur
        del(i,num,g_name,medianCur)    
        #Standardize the data before clustering the results
        for i in range(metab_data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])
        del(i)
        if num_comps == 2:
            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0])
            linkageTwo = linkage(data,linkList[1])

            #Create the appropriate plt figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(1,2,figsize=(8,8))

            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0],above_threshold_color='y',orientation='left',no_labels=True)
            dend2 = dendrogram(linkageTwo,ax=axes[1],above_threshold_color='y',orientation='left',no_labels=True)
            del(linkageOne,linkageTwo,num_comps)
            #save the linkage Comparison Figure.
            plt.savefig('Test2.png')
            #show the plot
            plt.show()

            

        elif num_comps == 3:
            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0])
            linkageTwo = linkage(data,linkList[1])
            linkageThree = linkage(data,linkList[2])

            #Create the appropriate plt figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(1,3,figsize=(8,8))

            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0],above_threshold_color='y',orientation='left',no_labels=True)
            dend2 = dendrogram(linkageTwo,ax=axes[1],above_threshold_color='y',orientation='left',no_labels=True)
            dend3 = dendrogram(linkageThree,ax=axes[2],above_threshold_color='y',orientation='left',no_labels=True)
            del(linkageOne,linkageTwo,linkageThree,num_comps)
            #save the linkage Comparison figure.
            plt.savefig('Test3.png')
            #create the plot
            plt.show()
            

        elif num_comps == 4:
            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0])
            linkageTwo = linkage(data,linkList[1])
            linkageThree = linkage(data,linkList[2])
            linkageFour = linkage(data, linkList[3])

            #Create the appropriate figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(2,2,figsize=(8,8))
            plt.title('Linkage Comparison')

            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0,0],above_threshold_color='y',orientation='left',no_labels=True)
            dend2 = dendrogram(linkageTwo,ax=axes[0,1],above_threshold_color='y',orientation='left',no_labels=True)
            dend3 = dendrogram(linkageThree,ax=axes[1,0],above_threshold_color='y',orientation='left',no_labels=True)
            dend4 = dendrogram(linkageFour,ax=axes[1,1],above_threshold_color='y',orientation='left',no_labels=True)
            
            del(linkageOne,linkageTwo,linkageThree,linkageFour,num_comps)
            #save the figure
            plt.savefig('Test4.png')
            plt.show()

        #log the completion of the linkage comparison
        logging.warning(': Sucessfuly completed the comparison of the linkage functions!')
            
    def compoundMatchUp():
        # Reads in our Kegg Compound Dataset (Single Column)
        kegg_data = pd.read_excel("kegg_compound_IDs_3.xlsx")

        # Splits our single column into two more user friendly ones
        kegg_data[["ID","compound"]] = kegg_data["ID"].str.split(" ", 1, expand = True)

        # Pulls in our matched compound data
        my_data = pd.read_csv("mummichog_matched_compound_all.csv")

        # Makes an ID column
        my_data["ID"] = my_data["Matched.Compound"]

        # Deletes the unneeded columns
        gonecolumns = ["Query.Mass", "Matched.Form", "Mass.Diff", "Matched.Compound"]
        my_data = my_data.drop(axis = 1, labels = gonecolumns)

        # Filters the keggs data to only the ids that are
        # in the matched compound dataset
        my_final_data = kegg_data[kegg_data["ID"].isin(my_data["ID"])]

        # Writes this final dataset to a csv
        my_final_data.to_csv(path_or_buf = "filename.csv")

    def ensembleClustering():
        #log that the user called ensemble clustering function
        logging.warning(': User called Ensemble Clustering function.')

        #The distance measures and linkage functions should be consistent but we could also develop
        #a GUI that allows for the users to select various distance measures. The linkage functions 
        #should be consistent for all ensemble clustering techniques
        filename = filedialog.askopenfilename()

        if filename == '':
            logging.warning(': Failed to select a file.')
            return

        #*******LinkageFunctions******************
        #*******Single
        #*******Complete
        #*******Average
        #*******Ward

        #*******Current Distance Measures*****************
        #*******Euclidean
        #*******Squared Euclidean
        #*******Cosine
        #*******Chebyshev
        #*******Other Available Distance Measures*
        #*******Minkowski
        #*******Cityblock
        #*******Standardized Euclidean
        #*******Correlation
        #*******Hamming
        #*******Jaccard
        #*******Canberra
        #*******Braycurtis
        #*******Mahaloanobis
        #*******Yule
        #*******Matching
        #*******Dice
        #*******Kulsinsi
        #*******RogerStanimoto
        #*******RussellRao
        #*******SokalMichener
        #*******SokalSneath
        #*******Wminkowski

        #Read in the data for ensemble clustering
        try:
            metab_data = pd.read_excel(filename, sheet_name='Medians')
        except:
            logging.warning(': Failed to open excel file, make sure that the sheet name is Medians')
            return
        del(filename)
        #List for the use in creating and plotting the clustering results
        linkageList = ['single','complete','average']
        distList = ['euclidean','sqeuclidean','cosine','chebyshev']

        #finding the number of groups in the metabolomics study to allow for proper clustering of the results
        num_groups = metab_data.shape[1] - 2

        #creating a numpy array that is the size of the data that is being read in.
        data = np.zeros((metab_data.shape[0],metab_data.shape[1]-2))

        for i in range(num_groups):
            #create the appropriate string to add to the group name
            num = str(i + 1)
            g_name = "M" + num
            #grab the Medians for each group
            try:
                medianCur = metab_data[g_name]
            except:
                logging.warning(': Improper column headers, make sure to follow the convention of M1-MN')
                return 

            #add the medians data to the array to be clustered
            data[:,i] = medianCur
        del(medianCur,num,g_name,num_groups,metab_data)  
        #Standardize the data before clustering the results
        for i in range(data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])


        #Create a plot with wanted ensemble clusters
        fig, axes = plt.subplots(3,4,figsize=(10,8))
        for i in range(len(linkageList)):
            for j in range(len(distList)):
                linkCur = linkage(data,linkageList[i],distList[j])
                
                dendCur = dendrogram(linkCur,ax=axes[i,j],above_threshold_color='y',orientation='top',no_labels=True)
        del(linkageList,distList)
        plt.savefig('EnsembleTest.png')
        plt.show()

        #Log the completion of the ensemble clustering function
        logging.warning(': Sucessfully completed Ensemble clustering!')

    
    #Function to create a minimum spanning tree
    def MST():
        #log that user called MST
        logging.warning(': User called Minimum Spanning Tree function.')
        #ask user to specify the file they would like to analyze with a minimum spanning tree
        file = filedialog.askopenfilename()

        if file == '':
            logging.warning(': Failed to select a file')
            return

        #read in the file that is needed to calculate the minimum spanning tree
        try:
            dataRaw = pd.read_excel(file,sheet_name='Medians')
        except:
            logging.warning(': Failed to open excel file, make sure sheet name is Medians')
            return
        
        #find the number of groups in that data
        num_groups = dataRaw.shape[1]-2
        #standardize the data before finding the distance matrix
        data = np.zeros((dataRaw.shape[0],dataRaw.shape[1]-2))

        for i in range(num_groups):
            #create the appropriate string to add to the group name
            num = str(i+1)
            g_name = "M" + num
            #grab the medians  for each group
            try:
                medianCur = dataRaw[g_name]
            except:
                logging.warning(': Improper column headers, make sure to follow the convention M1-MN')
                return 

            #add the medians data to the array to be clustered
            data[:,i] = medianCur

        #standardize the data before clustering
        for i in range(dataRaw.shape[0]):
            data[i,:] = GB.standardize(data[i,:])

        #find the distance matrix using the pairwise distance function used in the analysis of Agglomerative clustering. 
        pairWise = pdist(data)

        pairWise = squareform(pairWise)
        
        #put the distance matrix into the appropriate format for creation of MST using scipy's MST capabilities
        mstInput = csr_matrix(pairWise)

        #create the minimum spanning tree
        mstOut = minimum_spanning_tree(mstInput)
        mstOutInd = mstOut.nonzero()

        #create the matrix containing the various connections of the mst
        mstOutMat = mstOut.todense()
        
        dataMST = np.zeros([data.shape[0]-1,3])

        for i in range(data.shape[0]-1):
            #input the values of each matrix element to the numpy array for saving to csv file.
            dataMST[i,0] = mstOutInd[0][i]
            dataMST[i,1] = mstOutInd[1][i]
            dataMST[i,2] = mstOut[dataMST[i,0],dataMST[i,1]]


        #Input the dataMST into the dataframe to save the results of the MST for future use if needed
        mstOut = pd.DataFrame(dataMST, columns=['index1','index2','dist'])
        mstOut = mstOut.sort_values(by='dist')
        mstOutNp = mstOut.to_numpy()

        mstOutMat = pd.DataFrame(mstOutMat)

        #Create an empty dictionary that will contain the clusters
        clusters = {}

        #create an empty dictionary that allows us to save the clusters for validation.
        validationClusters = {}


        # create initial list of metabolites that serves as the initial list of metabolites that will be clustered.
        metabolites = np.ones((data.shape[0],1))

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
                        logging.warning(logCheck)
                        logging.warning(': Inappropriate data type for the for clustering ID algorithm.')
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

        #Validate the number of clusters that should be used in the clustering solutions.
        valIndex = GB.Validate(validationClusters,data,num_groups)
        valIndex = pd.DataFrame(valIndex)

        #save to a csv file
        mstOut.to_csv('MST.csv',index=False)

        #save validation measure to csv file
        valIndex.to_csv('valIndex.csv',index=False)

        #logging the completion of the Minimum spanning tree
        logging.warning(': Sucessfully completed MST!')


    def PDFGenerator():
        #log that the function has been called
        logging.warning(': User called PDF Generator function.')
        print('Commented Out!')

        # #create the pdf and title for each page.
        # pdf = fpdf.FPDF('P','mm','Letter')

        # #Create the title and set the default font
        # os.chdir("C:/Users/bradyhislop/Desktop")


        # #Create the first page
        # pdf.add_page()
        # pdf.set_font('Arial','B',20)
        # pdf.cell(197,10,'Clustering Results',0,0,'C')
        # pdf.ln(10)
        # pdf.cell(197,10,'Clustergram (Single, Euclidean)',0,0,'C')
        # pdf.ln(10)
        # file = filedialog.askopenfilename()
        # pdf.image(file,55,30,100,100)
        # pdf.ln(120)
        # pdf.cell(197,10,'Clustergram (Ward, Euclidean)',0,0,'C')
        # file = filedialog.askopenfilename()
        # pdf.image(file,55,160,100,100)
        

        # #create second page. 
        # pdf.add_page()
        # file = filedialog.askopenfilename()
        # pdf.cell(197,10,'Linkage Comparison',0,0,'C')
        # pdf.ln(10)
        # pdf.cell(197,10,'(Single, Ward, Complete, Average)',0,0,'C')
        # pdf.image(file,55,30,100,100)


        # pdf.output('Example.pdf','F')




        # #log the sucessful creation of the pdf
        # logging.warning(': Sucessfully created a pdf of the results!')


