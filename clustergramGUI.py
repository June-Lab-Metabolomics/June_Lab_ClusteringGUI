from tkinter import *
from tkinter import ttk
from GUIUtils import GUIUtils as GU

def clustergramGUI():
	root = Tk()
	root.geometry("500x500")
	root.title("Clustergram Options")

	def linkageOutput(link,dist):
		dist = int(dist.get())
		link = int(link.get())
		#determine the linkage function based upon radio button
		if link == 0:
			link = 'single'
		elif link == 1:
			link = 'complete'
		elif link == 2:
			link = 'ward'
		elif link == 3:
			link = 'average'


		#determine the distance measure based upon radio button
		if dist == 0:
			dist = 'euclidean'
			GU.createClustergram(link,dist)
		elif dist == 1:
			dist = 'sqeuclidean'
			GU.createClustergram(link,dist)
		elif dist == 2:
			dist = 'cosine'
			GU.createClustergram(link,dist)
		elif dist == 3:
			dist = 'chegyshev'
			GU.createClustergram(link,dist)


	#Create the frame that will contain the parameters to create a clutergram
	mainframe = ttk.Frame(root, padding="15 15 15 15")
	mainframe.grid(column=0,row=0, sticky=(N,W,E,S))

	#Create resizing parameters for each clustergram dataset.
	root.columnconfigure(0, weight = 1)
	root.rowconfigure(0, weight = 1)
	mainframe.columnconfigure(1, weight=2)
	mainframe.columnconfigure(2, weight =2)
	mainframe.rowconfigure(1, weight = 2)
	mainframe.rowconfigure(2, weight = 2)
	mainframe.rowconfigure(3, weight = 2)
	mainframe.rowconfigure(4, weight= 2)
	mainframe.rowconfigure(5, weight=2)

	#initialize variable for the radiobuttons
	output = IntVar()
	dist = IntVar()
	#Create the label for the linkage function caption
	linkageLabel = ttk.Label(mainframe,text="Linkage Functions",font=("TkHeadingFont",20)).grid(column=1,row=0,sticky=(N,S,E,W))
	#Create radiobuttons to allow for the selection of a linkage function
	linkageSingle =   ttk.Radiobutton(mainframe,text='Single',variable=output,value=0).grid(column=1,row=1,sticky=(N,S,E,W))
	linkageComplete = ttk.Radiobutton(mainframe,text='Complete',variable=output,value=1).grid(column=1,row=2,sticky=(N,S,E,W))
	linkageWard =     ttk.Radiobutton(mainframe,text='Ward',variable=output,value=2).grid(column=1,row=3,sticky=(N,S,E,W))
	linkageAverage =  ttk.Radiobutton(mainframe,text='Average',variable=output,value=3).grid(column=1,row=4,sticky=(N,S,E,W))

	#Create the label for the distance measure
	distanceLabel = ttk.Label(mainframe,text="Distance Measures",font=("TkHeadingFont",20)).grid(column=2,row=0,sticky=(N,S,E,W))
	#Create radiobuttons to allow for the selection of a distance measure
	distanceMeasureEuclidean = ttk.Radiobutton(mainframe,text='Euclidean',variable=dist,value=0).grid(column=2,row=1,sticky=(N,S,E,W))
	distanceMeasureSqEuclidean=ttk.Radiobutton(mainframe,text='Squared Euclidean',variable=dist,value=1).grid(column=2,row=2,sticky=(N,S,E,W))
	distanceMeasureCosine =    ttk.Radiobutton(mainframe,text='Cosine',variable=dist,value=2).grid(column=2,row=3,sticky=(N,S,E,W))
	distanceMeasureChebyshev = ttk.Radiobutton(mainframe,text='Chebyshev',variable=dist,value=3).grid(column=2,row=4,sticky=(N,S,E,W))

	#Create a button in the column that doesnt contain any radiobuttons that allow for the user to properly go to the create clustergram function
	clustergramBtn = ttk.Button(mainframe,text="Create Clustergram",command=lambda: linkageOutput(output,dist)).grid(column=2,row=5,sticky=(N,S,E,W))

	




	root.mainloop()



