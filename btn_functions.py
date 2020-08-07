#Creating a module that contains the the functions that call out class
# which contains the functions that create all of the needed GUI information
from tkinter import filedialog
from GUIUtils import GUIUtils as GU 
import clustergramGUI as CG
import linkageGUI as LG
import logging
import time

def cluster(*args):
    #ask the user to select a file before the program sends the file to another function to create the clustergram
    #call function as needed it should be a simple callout of the
    try:
        CG.clustergramGUI()
    except:
        logging.info('Failed to Create a Clustergram')
        return

def Medians(*args):
    #ask the user to select a file that will be used to create a medians file.
    filename = filedialog.askopenfilename()
    try:
        #Try to run the group medians function
        GU.groupMedians(filename)
    except:
        error_time = time.asctime()
        error_time += ':  Failed to run the group medians function'
        logging.info(error_time)
        return

def Linkages():
    #ask the user to select a clustergram file to compare linkage functions.
    try:
        LG.linkageGUI()
    except:
        logging.info('Failed to run the Linkage Comparison function')
        return

def Valid():
    #ask the user to select a clustergram file to run through a validition study.
    filename = filedialog.askopenfilename()

def P2P():
    #ask the user to select a file containing the selected clusters to send to a function that will create the peaks to pathways file.
    filename = filedialog.askopenfilename()

def integrity():
    #ask the user to select a volcano plot file to check the integrity of the data against. 
    filename = filedialog.askopenfilename()
    #VC.dataIntegrity(filename)
    try:
        GU.dataIntegrity(filename)
    except:
        error_time = time.asctime()
        error_time += ':   Failed to run the Data Integrity function'
        logging.info(error_time)
        return

def ensemble():
    #ask the user to select the file in which they would like to perfrom ensemble clustering on.
    try:
        GU.ensembleClustering()
    except:
        logging.info('Failed to run the Ensemble Clustering Function')
        return

def mst():
    #ask the user to select the file in which will be used to create a minimum spanning tree.
    try:
        GU.MST()
    except:
        logging.info('Failed to run the Minimum Spanning Tree Function')
        return

def generate():
    #ask the user to select the file in which will be used to create a minimum spanning tree.
    try:
        GU.PDFGenerator()
    except:
        logging.info('Failed to run the PDF Generator Function')
        return
