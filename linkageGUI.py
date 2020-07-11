from tkinter import *
from tkinter import ttk
from GUIUtils import GUIUtils as GU


def linkageGUI():
	root = Tk()
	root.geometry("500x500")
	root.title("Clustergram Options")

	#Create the frame that will contain the parameters to create a clutergram
	mainframe = ttk.Frame(root, padding="15 15 15 15")
	mainframe.grid(column=0,row=0, sticky=(N,W,E,S))
	#Create resizing parameters for each clustergram dataset.
	root.columnconfigure(0, weight = 1)
	root.rowconfigure(0, weight = 1)
	mainframe.columnconfigure(1, weight=2)
	mainframe.rowconfigure(1, weight = 2)
	mainframe.rowconfigure(2, weight = 2)
	mainframe.rowconfigure(3, weight = 2)
	mainframe.rowconfigure(4, weight = 2)

	#Create a check button
	single = StringVar()
	complete = StringVar()
	ward = StringVar()
	average = StringVar()
	singleCheck =   ttk.Checkbutton(mainframe,text='Single',variable=linkages,onvalue='single',offvalue='').grid(column=1,row=1,sticky=(N,S,E,W))
	completeCheck = ttk.Checkbutton(mainframe,text='Complete',variable=linkages,onvalue='complete',offvalue='').grid(column=1,row=2,sticky=(N,S,E,W))
	wardCheck     = ttk.Checkbutton(mainframe,text='Ward',variable=linkages,onvalue='ward',offvalue='').grid(column=1,row=3,sticky=(N,S,E,W))
	averageCheck =  ttk.Checkbutton(mainframe,text='Average',variable=linkages,onvalue='averages',offvalue='').grid(column=1,row=4,sticky=(N,S,E,W))
