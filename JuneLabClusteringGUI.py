from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import tkinter as tk
import GuiBackground as GB
from GUIUtils import GUIUtils as GU 
import logging
import time
import getpass
import os
import reinit as RI


class JuneLabClusteringGUI(ttk.Frame):
	def __init__(self, master=None):
		super().__init__(master)
		self.grid(column=0, row=0, sticky=(N, W, E, S))
		self.rowconfigure(0, weight =2)
		self.columnconfigure(1, weight =2)
		self.columnconfigure(2, weight =2)
		self.rowconfigure(1, weight = 2)
		self.rowconfigure(2, weight = 2)
		self.rowconfigure(3, weight = 2)
		self.columnconfigure(3, weight=2) 
		self.pack()
		self.create_widgets()

	def create_widgets(self):
		self.style = ttk.Style()
		self.style.configure("RW.TLabel", foreground="#e7a53a",font=("TkHeadingFont",36))
		self.style.configure("RW.TButton", padding=15, borderwidth=15, foreground="black", background="#e7a53a",font=("TkHeadingFont",20))
		self.JuneLab = ttk.Label(self, text="Welcome to the June Lab Clustering GUI",style="RW.TLabel").grid(column=0,row=0,columnspan=4)
		self.clust = ttk.Button(self,text="Create Clustergram",style="RW.TButton",command=self.createClustergram).grid(column=1,row=1, sticky=(N,S,E,W))
		#Create a button to allow the user to create a medians file for better clustering results. 
		self.med = ttk.Button(self, text="Group Medians", style="RW.TButton", command=self.medians).grid(column=1, row=3, sticky =(N,S,E,W))
		#Create a button to allow the user to compare the four most common linkage functions.
		self.link = ttk.Button(self, text="Compare Linkage Functions", style="RW.TButton", command=self.linkages).grid(column=2, row=1,sticky =(N,S,E,W))
		#Create a button to allow the user to validate the appropriate number of clusters needed for a given set of metabolites.
		self.val = ttk.Button(self, text="Compound Match-Up", style="RW.TButton", command=self.compound).grid(column=3, row=2, sticky =(N,S,E,W))
		#Create a button to allow the user to create the peaks to pathways files needed to analyze the peaks to pathways in Mummichog
		self.peak = ttk.Button(self, text="Peaks to Pathways", style="RW.TButton", command=self.P2P).grid(column=2, row=2, sticky =(N,S,E,W))
		#Create a button to allow the user to check the integrity of the data downloaded from Metaboanalysts Volcano plot results. 
		self.integrity = ttk.Button(self, text="Data Integrity", style="RW.TButton", command=self.integrity).grid(column=2, row=3, sticky =(N,S,E,W))
		#Create a button to allow the user to do an Ensemble clustering on the data.
		self.ensemble = ttk.Button(self,text="Ensemble Clustering", style="RW.TButton", command=self.ensemble).grid(column=1,row=2,sticky=(N,S,E,W))
		#Create a button to allow the user to create a minimum spanning tree on data
		self.mst = ttk.Button(self,text='Minimum Spanning Tree', style="RW.TButton", command=self.mst).grid(column=3,row=1,sticky=(N,S,E,W))
		#Create a button for the generation of a report
		self.generate = ttk.Button(self,text='Generate PDF Report', style="RW.TButton", command=self.generate).grid(column=3,row=3,sticky=(N,S,E,W))
		# pad each widget with 5 pixels on each side to ensure that the buttons do not stay together. 
		for child in self.winfo_children(): child.grid_configure(padx=5, pady=5)


	def home(self):
		self.master.destroy()
		RI.reinit()
		# objects = self.grid_slaves()
		# for i in objects:
		# 	i.destroy()
		# self.create_widgets()


	def createClustergram(self):
		def linkageOutput(*args):
			#grab the current selection of the list
			selection = distListBox.curselection()
			selectionIndex = binaryList[selection[0]]
			for i in range(int(len(selectionIndex)/2)):
				if i == 0:
					linkLoc = (2*int(selectionIndex[0]) + int(selectionIndex[1]))
					link = linkageList[linkLoc]
				elif i == 1:
					distLoc = (2*int(selectionIndex[2]) + int(selectionIndex[3]))
					dist = distList[distLoc]
					GU.createClustergram(0,link,dist)


		objects = self.grid_slaves()
		for i in objects:
			i.destroy()

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Clustergram Input",font=("TkHeadingFont",36)).grid(column=0,row=0,sticky=(N))
		self.home = ttk.Button(self,text="Return to Home",command=self.home).grid(column=0,row=2, sticky=(N))
		distListBox = Listbox(self,height=5)#.grid(column=1,row=1,rowspan=4)
		
		#Create the lists of available options for selection 
		linkageList = ('single','ward','complete','average')
		distList = ('euclidean','sqeuclidean','cosine','chebyshev')
		linkDistList = ('single====euclidean','single====sqeuclidean','single====cosine','single====chebyshev','ward====euclidean','complete====euclidean',\
					'complete====sqeuclidean','complete====cosine','complete====chebyshev','average====euclidean','average====sqeuclidean','average====cosine','average====chebyshev')
		binaryList = ('0000','0001','0010','0011','0100','1000','1001','1010','1011','1100','1101','1110','1111')
		distNames = StringVar(value=linkDistList)

		for i in range(len(linkDistList)):
			distListBox.insert(i,linkDistList[i])

		distListBox.bind('<Double-1>',linkageOutput)
		self.distListBox = distListBox
		self.distListBox.grid(column=0,row=1,columnspan=4)

	def medians(self):
	    #send the user to the groupMedians function
	    GU.groupMedians()

	def linkages(self):
		def linkageComp(*args):
			#get the file containing the groups
			linkList = []
			file = filedialog.askopenfilename()
			selection = linkListBox.curselection()
			selectionIndex = binaryList[selection[0]]
			for i in range(int(len(selectionIndex))):
				if int(selectionIndex[i]) == 1:
					#put the output into the linkList to be sent to the linkageComparison function
					linkList.append(linkageList[i])

			#get the linkage length to tell the linkageComparison how many comparisons to perform
			num_comps = len(linkList)
			
			#send the parameters for linkage comparison 		
			GU.linkageComparison(file,num_comps,linkList)

		objects = self.grid_slaves()
		for i in objects:
			i.destroy()

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Linkage Comparison Options",font=("TkHeadingFont",36)).grid(column=0,row=0,sticky=(N))
		self.home = ttk.Button(self,text="Return to Home",command=self.home).grid(column=0,row=2, sticky=(N))
		linkListBox = Listbox(self,height=5)

		#Create the list of available options for selection
		linkageList = ('single','ward','complete','average')
		linkageOptions = ('single====ward','single====complete','single====average','ward====complete','ward====average','complete====average',\
					  'single==ward==complete','single==ward==average','single==complete==average','ward==complete==average',\
					  'single=ward=complete=average')
		binaryList = ('1100','1010','1001','0110','0101','0011','1110','1101','1011','0111','1111')
		#Create linkages variables
		linkages = StringVar(value=linkageOptions)

		#
		for i in range(len(linkageOptions)):
			linkListBox.insert(i,linkageOptions[i])

		linkListBox.bind('<Double-1>',linkageComp)
		self.linkListBox = linkListBox
		self.linkListBox.grid(column=0,row=1,columnspan=4)

	def compound(self):
	    #ask the user to select a clustergram file to run through a validition study.
	    #Waiting on confirmation...
	    GU.compoundMatchUp()

	def P2P(self):
	    #Waiting until the MST functionality is complete. 
	    GU.peaksToPathways()
	def integrity(self):
	    #ask the user to select a volcano plot file to check the integrity of the data against. 
	    filename = filedialog.askopenfilename()
	    GU.dataIntegrity(filename)

	def mstF(self):
		objects = self.grid_slaves()
		for i in objects:
			i.destroy()
		self.JuneLab = ttk.Label(self, text="Close MST Validation graph prior to running Ensemble!!",foreground="red",font=("TkHeadingFont",24)).grid(column=0,row=1,sticky=(N))
		self.ok = ttk.Button(self,text="OK",command=self.ensemble).grid(column=0,row=2,sticky=(N))
		GU.MST()


	def ensemble(self):
		def ensembleOutput(*args):
			#grab the current selection of the list
			selection = optClustBox.curselection()
			selection = int(selection[0])+1
			
			GU.ensembleClustering(optNum = selection)

		objects = self.grid_slaves()
		for i in objects:
			i.destroy()
		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Input Optimal Number of Clusters",font=("TkHeadingFont",36)).grid(column=0,row=0,sticky=(N))
		self.mstF = ttk.Button(self,text="Run MST",command=self.mstF).grid(column=0,row=2,sticky=(N))
		self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=0,row=3, sticky=(N))

		#create list box of the ensemble optimal clusters
		optClustBox = Listbox(self,height=5)
		#Create the lists of available options for selection 
		optClusters = ('1','2','3','4','5','6',\
					'7','8','9','10','11','12','13')
		distNames = StringVar(value=optClusters)

		for i in range(len(optClusters)):
			optClustBox.insert(i,optClusters[i])

		#create binding event and 
		optClustBox.bind('<Double-1>',ensembleOutput)
		self.optClustBox = optClustBox
		self.optClustBox.grid(column=0,row=1,columnspan=4)


	def mst(self):
	    #send the user to the minimum spanning tree function.
	    GU.MST()

	def generate(self):
	    #ask the user to select the file in which will be used to create a minimum spanning tree.
	    GU.PDFGenerator()


curUser = getpass.getuser()
curUser = GB.who(curUser)
log_time = time.strftime("%a_%b_%d_%Y_%H_%M_%S")
directory = "C:/Users/" + getpass.getuser() + '/Desktop/LogPOutputFiles'
os.chdir(directory)
log_file = curUser + '_' + str(log_time) + '.log' 
logging.basicConfig(filename=log_file,format='%(asctime)s %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO)
logging.info('Started GUI')

root = tk.Tk()
app = JuneLabClusteringGUI(master=None)
app.master.minsize(200,200)
height = app.master.winfo_screenheight() - 400
width = app.master.winfo_screenwidth()-400
app.master.maxsize(width,height)
width = str(width)
height = str(height)
screenSize = width +'x' + height + '+0+100'
app.master.geometry("")
app.master.resizable(0,0)

app.master.title("Home")
app.mainloop()




