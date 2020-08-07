from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from GUIUtils import GUIUtils as GU


def linkageGUI():
	root = Tk()
	root.title("Linkage Comparison")

	#Create the list of available options for selection
	linkageList = ('single','ward','complete','average')
	linkageOptions = ('single====ward','single====complete','single====average','ward====complete','ward====average','complete====average',\
					  'single==ward==complete','single==ward==average','single==complete==average','ward==complete==average',\
					  'single=ward=complete=average')
	binaryList = ('1100','1010','1001','0110','0101','0011','1110','1101','1011','0111','1111')


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

	#Create the frame that will contain the parameters to create a clutergram
	mainframe = ttk.Frame(root, padding="15 15 15 15")
	mainframe.grid(column=0,row=0, sticky=(N,W,E,S))
	#Create resizing parameters for each clustergram dataset.
	root.columnconfigure(0, weight = 1)
	root.rowconfigure(0, weight = 1)
	mainframe.columnconfigure(1, weight=2)
	mainframe.rowconfigure(1, weight = 2)

	#Create linkages variables
	linkages = StringVar(value=linkageOptions)

	#Create the listbox widgets to allows the user to allow for the selection of the linkage functions
	linkListBox = Listbox(mainframe,height=5)
	linkLabel = ttk.Label(mainframe,text="Linkage Comprison Options")
	instruct = ttk.Label(mainframe,text="Double-click to select the link functions you would link to compare!")

	#input the linkage comparison value to the list
	for i in range(len(linkageOptions)):
		linkListBox.insert(i,linkageOptions[i])

	#configure the grid to handle the list
	linkListBox.grid(column=1,row=1,rowspan=6,sticky=(N,S,E,W))
	linkLabel.grid(column=1,row=0,sticky=(N,S,E,W))
	instruct.grid(column=1,row=7,sticky=(N,S,E,W))
	root.grid_columnconfigure(0,weight=1)
	root.grid_rowconfigure(5,weight=1)

	linkListBox.bind('<Double-1>',linkageComp)

	linkListBox.selection_set(0)



	root.mainloop()