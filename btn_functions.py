#Creating a module that contains the the functions that call out class
# which contains the functions that create all of the needed GUI information
from tkinter import filedialog
from GUIUtils import GUIUtils as GU 
import clustergramGUI as CG
import linkageGUI as LG

def cluster(*args):
    #call function as needed it should be a simple callout of the
    CG.clustergramGUI()

def Medians(*args):
    #ask the user to select a file that will be used to create a medians file.
    filename = filedialog.askopenfilename()
    GU.groupMedians(filename)

def Linkages():
    #send the user to the linkage GUI function
    LG.linkageGUI()

def Valid():
    #ask the user to select a clustergram file to run through a validition study.
    #Waiting on confirmation...
    filename = filedialog.askopenfilename()

def P2P():
    #Waiting until the MST functionality is complete. 
    filename = filedialog.askopenfilename()

def integrity():
    #ask the user to select a volcano plot file to check the integrity of the data against. 
    filename = filedialog.askopenfilename()
    GU.dataIntegrity(filename)


def ensemble():
    #send the user to the ensemble clustering function
    GU.ensembleClustering()

def mst():
    #send the user to the minimum spanning tree function.
    GU.MST()

def generate():
    #ask the user to select the file in which will be used to create a minimum spanning tree.
    GU.PDFGenerator()
