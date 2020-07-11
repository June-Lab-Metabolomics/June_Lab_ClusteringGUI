#Creating a module that contains the the functions that call out class
# which contains the functions that create all of the needed GUI information
from tkinter import filedialog
from GUIUtils import GUIUtils as GU 
import clustergramGUI as CG
import linkageGUI as LG

def cluster(*args):
    #ask the user to select a file before the program sends the file to another function to create the clustergram
    #call function as needed it should be a simple callout of the
    CG.clustergramGUI()

def Medians(*args):
    #ask the user to select a file that will be used to create a medians file.
    filename = filedialog.askopenfilename()
    GU.groupMedians(filename)
    
def Linkages():
    #ask the user to select a clustergram file to compare linkage functions.
    #filename = filedialog.askopenfilename()
    LG.linkageGUI()
    #GU.linkageComparison(filename,4)

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
    GU.dataIntegrity(filename)

def ensemble():
    #ask the user to select the file in which they would like to perfrom ensemble clustering on.
    GU.ensembleClustering()

def mst():
    #ask the user to select the file in which will be used to create a minimum spanning tree.
    filename = filedialog.askopenfilename()

def generate():
    #ask the user to select the file in which will be used to create a minimum spanning tree.
    filename = filedialog.askopenfilename()
