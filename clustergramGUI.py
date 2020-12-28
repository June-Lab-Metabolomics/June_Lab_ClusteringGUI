from tkinter import *
from tkinter import ttk
from GUIUtils import GUIUtils as GU

def clustergramGUI():
	root = Tk()
	root.title("Clustergram Options")

	#Create the lists of available options for selection 
	linkageList = ('single','ward','complete','average')
	distList = ('euclidean','sqeuclidean','cosine','chebyshev')
	linkDistList = ('single====euclidean','single====sqeuclidean','single====cosine','single====chebyshev','ward====euclidean','complete====euclidean',\
					'complete====sqeuclidean','complete====cosine','complete====chebyshev','average====euclidean','average====sqeuclidean','average====cosine','average====chebyshev')
	binaryList = ('0000','0001','0010','0011','0100','1000','1001','1010','1011','1100','1101','1110','1111')
	distNames = StringVar(value=linkDistList)


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

	#Create the frame that will contain the parameters to create a clutergram
	mainframe = ttk.Frame(root, padding="15 15 15 15")
	mainframe.grid(column=0,row=0, sticky=(N,W,E,S))

	#Create resizing parameters for each clustergram dataset.
	root.columnconfigure(0, weight = 1)
	root.rowconfigure(0, weight = 1)
	root.rowconfigure(1,weight =1)



	#Create all of the widgets
	distListBox = Listbox(mainframe,height=5)
	distLabel = ttk.Label(mainframe,text="Clustering Parameters")
	instructLabel = ttk.Label(mainframe,text="Double Click to select the wanted parameters")


	for i in range(len(linkDistList)):
		distListBox.insert(i,linkDistList[i])

	#configure grid to handle the list length
	distListBox.grid(column=1,row=1,rowspan=6,sticky=(N,S,E,W))
	distLabel.grid(column=1,row=0,sticky=(N,S,E,W))
	instructLabel.grid(column=1,row=7,sticky=(N,S,E,W))
	root.grid_columnconfigure(0,weight=1)
	root.grid_rowconfigure(5,weight=1)

	distListBox.bind('<Double-1>',linkageOutput)

	distListBox.selection_set(0)


	root.mainloop()



