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
import fpdf
import os
import statistics as stat
import GuiBackground as GB
from tkinter import filedialog


class GUIUtils:
    
    def dataIntegrity(file):
        #log that the user called the data integrity function
        logging.info(': User called the Data Integrity function.')

        #Read in Volcano Plot data
        try:
            volcano = pd.read_excel(file)
        except:
            logging.error(': Failed to read in the excel file. Please put error in the Github issues tab.')
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
                                logging.error(': Unable to convert values to floats. Make sure all data values are only contain decimals or numberic values')
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
        logging.info(': Data Integrity check sucessfully completed.')

    def createClustergram(norm,linkFunc,distMet):
        #log that the user called the Create Clustergram function
        logging.warning(': User called the Create Clustergram Function.')
        #check that the file the user selects is appropriate
        metab_data = GB.fileCheck()
        if metab_data is None:
            #log error message and return for soft exit.
            logging.error(': Error loading in the Excel sheet.')
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
                logging.error(': Improper column header, make sure to follow the convention of M1-MN')
                return

            #add the medians data to the array to be clustered
            data[:,i] = medianCur
        del(num_groups,num,medianCur) 
        #Standardize the data before clustering the results
        for i in range(metab_data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])
        del(metab_data)
        #create dendrogram and plot data
        GB.create_dendrogram(data,norm,link=linkFunc,dist=distMet,func='clustergram')
        GB.plotting()
        del(data,norm,linkFunc,distMet)
        #log that you have sucessfully created the clustergram
        logging.info(': Sucessfully created the wanted Clustergram')

    def groupMedians():
        '''
        Determine the number of groups and then create a list or array of the appropriate
        beginning and ending of each group. This assumes that the groups are all of equal size which should be
        the goal for any and all analysis. Groups with out the same sizes should be considered
        inappropriate for analysis in this context, additionally it should be noted that statistics
        with out the same groups sizes can lead to incorrect analysis.
        '''
        #log that the user called the group medians function
        logging.info(': User called the Group Medians function.')
        file = filedialog.askopenfilename()
        medians = GB.fileCheck(file =file)
        if medians is None:
            #log error message and return the function for soft exit.
            logging.error(': Error reading in the excel sheet')
            return

        #get the group letters to know which groups to find the medians of
        groups = medians['mz']
        numObs = int(medians.shape[0]-1)

        num_groups = int(medians['Unnamed: 0'][numObs][0])
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
        logging.info(': Sucessfully grouped the Medians of each group!')

    def linkageComparison(file,num_comps,linkList):
        #Log that user called linkage comparison function
        logging.info(': User called the Linkage Comparison function.')
        #check that the file is appropriate for our data set
        metab_data = GB.fileCheck(file)    
        if metab_data is None:
            #Logs error and returns function to ensure soft exit.
            logging.error(': Error loading in excel file check log file!')
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
                logging.error(': Improper column headers, make sure you are following the convention M1-MN')
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
            distMeasure = pdist(data)
            distMeasure = squareform(distMeasure)
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
        logging.info(': Sucessfuly completed the comparison of the linkage functions!')
            
    def compoundMatchUp():
        logging.info(': Compound Match-Up function called!')
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
        logging.info(': User called Ensemble Clustering function.')

        #optimum number of clusters from validation index.
        optNum = 3

        #The distance measures and linkage functions should be consistent but we could also develop
        #a GUI that allows for the users to select various distance measures. The linkage functions 
        #should be consistent for all ensemble clustering techniques
        metab_data = GB.fileCheck()

        #List for the use in creating and plotting the clustering results
        linkageList = ['single','complete','average']
        distList = ['euclidean','sqeuclidean','cosine','chebyshev']

        #calculate the number of clusterings based upon the size of the lists
        numClusterings = len(linkageList)*len(distList)

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
                logging.error(': Improper column headers, make sure to follow the convention of M1-MN')
                return 

            #add the medians data to the array to be clustered
            data[:,i] = medianCur
        del(medianCur,num,g_name,num_groups,metab_data)  
        #Standardize the data before clustering the results
        for i in range(data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])

        #determine the the number of clusters and the dictionary location that needs to be called. 
        numMetabs = data.shape[0]
        dictLoc = numMetabs-optNum-1

        #create co-occurrence matrix.
        coOcc = GB.cooccurrence(data)


        for i in range(len(linkageList)):
            for j in range(len(distList)):
                linkCur = linkage(data,linkageList[i],distList[j])
                valid = GB.clustConnectLink(linkCur)
                coOcc = GB.popCooccurrence(valid[dictLoc],coOcc,numClusterings)
                print('passed')
        del(linkageList,distList)
        link = 'ward'
        dist = 'euclidean'
        GB.create_dendrogram(coOcc,norm=0,link=link,dist=dist,func='ensemble')
        GB.plotting()
        #Log the completion of the ensemble clustering function
        logging.info(': Sucessfully completed Ensemble clustering!')

    #Function to create a minimum spanning tree
    def MST():
        #log that user called MST
        logging.info(': User called Minimum Spanning Tree function.')
        #check the stability of the file
        dataRaw = GB.fileCheck()
        if dataRaw is None:
            #log error and return function to ensure a soft closing of the class
            logging.error(': Error loading the Excel sheet.')
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
                logging.error(': Improper column headers, make sure to follow the convention M1-MN')
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
 
        #determine how the minimum spanning tree was created for validation of clusters
        validationClusters = GB.clustConnect(dataMST,mstOutNp)

        #Validate the number of clusters that should be used in the clustering solutions.
        valIndex = GB.Validate(validationClusters,data,num_groups)
        x = valIndex[1,:]
        y = valIndex[0,:]
        #outputs the validity index
        nonZeros = x.nonzero()
        nonZerosLen = len(nonZeros[0])

        xUpdate = np.zeros(nonZerosLen)
        yUpdate = np.zeros(nonZerosLen)
        for i in range(nonZerosLen):
            xUpdate[i] = x[nonZeros[0][i]]
            yUpdate[i] = y[nonZeros[0][i]]

        valIndex = pd.DataFrame(valIndex)

        valMin = min(yUpdate)
        valMin = np.where(yUpdate==valMin)
        print(valMin)

        #save to a csv file
        mstOut.to_csv('MST.csv',index=False)

        #save validation measure to csv file
        valIndex.to_csv('valIndex.csv',index=False)
        #logging the completion of the Minimum spanning tree
        logging.info(': Sucessfully completed MST and clustering validation!')
        plt.plot(xUpdate,yUpdate)
        plt.xlabel('Clusters')
        plt.ylabel('Validation Index')
        plt.title('Cluster Validation')
        plt.show()

    def PDFGenerator():
        #log that the function has been called
        logging.info(': User called PDF Generator function.')

        #create the pdf and title for each page.
        pdf = fpdf.FPDF('P','mm','Letter')

        #Create the title and set the default font
        directory = filedialog.askdirectory()
        os.chdir(directory)

        #Create the first page
        pdf.add_page()
        pdf.set_font('Arial','B',24)
        pdf.cell(197,10,'Metabolanalyst Results',0,0,'C')
        pdf.set_font('Arial','B',14)
        pdf.set_font('')
        pdf.ln(10)
        pdf.cell(197,10,'Normalization',0,0,'L')
        pdf.ln(10)
        file = 'norm_dpi300.png'
        pdf.image(file,55,30,100,100)
        pdf.ln(120)
        pdf.cell(197,10,'Sample Normalization',0,0,'L')
        file = 'snorm_dpi300.png'
        pdf.image(file,55,160,100,100)
        
        #create second page. 
        pdf.add_page()
        file = 'pca_pair_dpi300.png'
        pdf.cell(197,10,'PCA pairs',0,0,'L')
        pdf.ln(10)
        pdf.image(file,55,30,100,100)
        pdf.ln(10)
        pdf.ln(120)
        pdf.cell(197,10,'PCA Scree plots',0,0,'L')
        file = 'pca_scree_dpi300.png'
        pdf.image(file,55,160,100,100)

        #create third page.
        pdf.add_page()
        file = 'pca_loading_dpi300.png'
        pdf.cell(197,10,'PCA Loading',0,0,'L')
        pdf.ln(10)
        pdf.image(file,55,30,100,100)
        pdf.ln(10)
        pdf.ln(120)
        pdf.cell(197,10,'PCA biplot',0,0,'L')
        file = 'pca_biplot_dpi300.png'
        pdf.image(file,55,160,100,100)

        #create fourth page
        pdf.add_page()
        pdf.cell(197,10,'PLS-DA Pairs',0,0,'L')
        pdf.ln(10)
        file = 'pls_pair_dpi300.png'
        pdf.image(file,55,30,100,100)
        pdf.ln(10)
        pdf.ln(120)
        pdf.cell(197,10,'PLS-DA 1 v.2',0,0,'L')
        file ='pls_score2d_dpi300.png'
        pdf.image(file,55,160,100,100)

        #create fifth page
        pdf.add_page()
        pdf.cell(197,10,'PLS-DA 3D plot',0,0,'L')
        pdf.ln(10)
        file = 'pls_score3d_dpi300.png'
        pdf.image(file,55,30,100,100)
        pdf.ln(10)
        pdf.ln(120)
        pdf.cell(197,10,'PLS-DA Loading',0,0,'L')
        file ='pls_loading_dpi300.png'
        pdf.image(file,55,160,100,100)

        #create sixth page
        pdf.add_page()
        pdf.cell(197,10,'PLS-DA CV',0,0,'L')
        pdf.ln(10)
        file = 'pls_cv_dpi300.png'
        pdf.image(file,55,30,100,100)
        pdf.ln(10)
        pdf.ln(120)
        pdf.cell(197,10,'PLS-DA Imp',0,0,'L')
        file ='pls_imp_dpi300.png'
        pdf.image(file,55,160,100,100)

        #create seventh page
        pdf.add_page()
        pdf.cell(197,10,'Dendrogram',0,0,'L')
        pdf.ln(10)
        file = 'tree_dpi300.png'
        pdf.image(file,55,30,100,100)
        pdf.ln(10)
        pdf.ln(120)
        pdf.cell(197,10,'Fold Change',0,0,'L')
        file = 'fc_dpi300.png'
        pdf.image(file,55,160,100,100)

        #create eigth page
        pdf.add_page()
        pdf.cell(197,10,'T-test',0,0,'L')
        file = 'tt_dpi300.png'
        pdf.image(file,55,30,100,100)
        pdf.ln(10)
        pdf.ln(120)
        pdf.cell(197,10,'Volcano Plot',0,0,'L')
        file = 'volcano_dpi300.png'
        pdf.image(file,55,160,100,100)

        #create the pdf of the results
        pdf.output('Example.pdf','F')
        #log the sucessful creation of the pdf
        logging.info(': Sucessfully created a pdf of the results!')