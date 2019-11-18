# -*- coding: utf-8 -*-

import math
import sys
import os
import numpy as np
import KMC_test8 as KMC_engine
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QPushButton, QAction, QVBoxLayout, QGraphicsView, QToolBar, QGraphicsScene, QButtonGroup, QHBoxLayout, QGraphicsRectItem, QGraphicsItem, QGraphicsItemGroup, QMenu, QAction, QLabel, QDialog, QLineEdit, QMessageBox, QFileDialog, QListView
from PyQt5.QtGui import QIcon, QPixmap, QPolygon, QColor, QPainter, QPen, QBrush, QTransform, QFont, QFontMetrics, QPolygonF, QPainterPath, QStandardItemModel, QStandardItem
from PyQt5.QtCore import Qt, QPointF, QRectF, QLine, QVariant
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


### Global variables and options
TOOLS= ['species', 'reactAtoB', 'reactABtoC', 'reactAtoBC', 'reactABtoCD', 'reactAtoBCD','reactABtoCDE']

Na=6.02*10**23
nodeFont=QFont("times",16)
speciesFillColor=QColor(79,126,151)
reactionsFillColor=QColor(46,144,66)
plugSideLength=30
plugWidth=plugSideLength*(math.sqrt(3)/2)
lineInteractionRange=20		
speciesCounter=1
speciesList=[]
reactionsCounter=1
reactionsList=[]
connectionsCounter=1
connectionsList=[]
isMoving=False
movingItem=None
isConnecting=False
connectionStart=None
connectionEnd=None
KMCParams=[100000,1,0.01,1,1]#total N of starting particles,repeats,store timestep,max time,volume
fileName=None
lasttVector=[]
lastPVector=[]
lastCVector=[]


######################
### global Methods ###
######################
def calcPopulationVector(totalNParticles):
	global speciesList
	global KMCParams
	Na=6.02*10**23
	concList=[]
	sumConc=0
	if len(speciesList)!=0:
		for species in speciesList:
			concList.append(species.nodeBox.number)
			sumConc+=species.nodeBox.number
	if sumConc !=None:
		V=totalNParticles/(sumConc*Na)
	populationList=np.array([])
	for concentration in concList:
		populationList=np.append(populationList,int(V*concentration*Na))
	if sumConc !=None:
		KMCParams[4]=V
	return populationList
	
#edit KMC parameters
def editKMC():
	global speciesList
	global KMCParams
	getKMCParams=editKMCParams()
#generate output file stream
def generateOutputStream():
	global speciesList
	global reactionsList
	global connectionsList
	global KMCParams
	outputStream=[]
	
	#species list 
	for species in speciesList:
		outputLine='//species '+str(species.pos().x())+" "+str(species.pos().y())+" "+str(species.nodeBox.name)+" "+str(species.nodeBox.number)
		outputStream.append(outputLine)
	
	#reactions list
	for reactions in reactionsList:
		outputLine='//reactions '+str(reactions.pos().x())+" "+str(reactions.pos().y())+" "+str(reactions.nodeBox.name)+" "+str(reactions.nodeBox.number)
		#append reaction type to the end of the line
		if isinstance(reactions,reactionAtoBNode)==True:
			outputLine=outputLine+" AtoB"
		if isinstance(reactions,reactionABtoCNode)==True:
			outputLine=outputLine+" ABtoC"
		if isinstance(reactions,reactionAtoBCNode)==True:
			outputLine=outputLine+" AtoBC"
		outputStream.append(outputLine)
	
	#connections list
	for connection in connectionsList:
		outputLine='//connections '+str(connection.startNode.parentItem().nodeBox.name)+" "+str(connection.startNode.name)+" "+str(connection.endNode.parentItem().nodeBox.name)+" "+str(connection.endNode.name)
		outputStream.append(outputLine)
	#KMC parameters
	outputLine='//KMCparams '+str(KMCParams[0])+" "+str(KMCParams[1])+" "+str(KMCParams[2])+" "+str(KMCParams[3])+" "+str(KMCParams[4])
	outputStream.append(outputLine)
	#population vector
	populationVector=calcPopulationVector(KMCParams[0])
	outputLine='//popVector'
	for item in populationVector:
		outputLine=outputLine+" "+str(item)
	outputStream.append(outputLine)
	#name vector
	outputLine='//nameVector'
	for species in speciesList:
		outputLine=outputLine+" "+str(species.nodeBox.name)
	outputStream.append(outputLine)
	#rate vector
	outputLine='//rateVector'
	for reaction in reactionsList:
		outputLine=outputLine+" "+str(reaction.nodeBox.number)
	outputStream.append(outputLine)
	#connectivity matrix
	connectivityMatrix=np.zeros(shape=(len(reactionsList),len(speciesList)))

	#iterate over all reactions (rows in connectivity matrix)
	i=0
	while i < len(reactionsList):
		#iterate over plugs of reaction
		for reactionChildItem in reactionsList[i].childItems():
			if isinstance(reactionChildItem,plug):
				#iterate over all species (columns in connectivity matrix)
				j=0
				while j< len(speciesList):
					#check all plugs in given species
					for speciesChildItem in speciesList[j].childItems():
						if isinstance(speciesChildItem,plug):
							#iterate over all connections to check if connection exists 
							for connection in connectionsList:
								#check if connection's start and end plugs are identical to current reaction and species plugs
								if connection.startNode==reactionChildItem and connection.endNode==speciesChildItem:
									if reactionChildItem.mode=="in":
										connectivityMatrix[i][j]-=1
									if reactionChildItem.mode=="out":
										connectivityMatrix[i][j]+=1
								if connection.endNode==reactionChildItem and connection.startNode==speciesChildItem:
									if reactionChildItem.mode=="in":
										connectivityMatrix[i][j]-=1
									if reactionChildItem.mode=="out":
										connectivityMatrix[i][j]+=1
					j+=1	

		i+=1
	for line in connectivityMatrix:
		outputLine='//connMatrix'
		for item in line:
			outputLine=outputLine+" "+str(item)
		outputStream.append(outputLine)
	return outputStream
def readInputStream(inputStream):
	global speciesList
	global reactionsList
	global connectionsList
	global speciesCounter
	global reactionsCounter
	global connectionsCounter
	global KMCParams
	speciesCounter=0
	reactionsCounter=0
	connectionsCounter=0

	for line in inputStream:
		lineList=line.split(" ")
		
		if lineList[0] == '//species':
			#generate new species based on this line information
			objectName=lineList[3]
			objectName=speciesNode(QPointF(float(lineList[1]),float(lineList[2])),lineList[3],float(lineList[4]))
			speciesList.append(objectName)
			speciesCounter+=1
			AppWindow.canvas.addItem(objectName)
			AppWindow.canvas.update()
			
		if lineList[0] =='//reactions':
			#generate new reactions based on this line information
			objectName=lineList[3]
			if lineList[5]=='AtoB':
				objectName=reactionAtoBNode(QPointF(float(lineList[1]),float(lineList[2])),lineList[3],float(lineList[4]))
			if lineList[5]=='AtoBC':
				objectName=reactionAtoBCNode(QPointF(float(lineList[1]),float(lineList[2])),lineList[3],float(lineList[4]))
			if lineList[5]=='ABtoC':
				objectName=reactionABtoCNode(QPointF(float(lineList[1]),float(lineList[2])),lineList[3],float(lineList[4]))
			
			reactionsList.append(objectName)
			reactionsCounter+=1
			AppWindow.canvas.addItem(objectName)
			AppWindow.canvas.update()
			
			
		if lineList[0] =='//connections':
			#generate new connections based on this line information
			#get plugs of startNode:
			for species in speciesList+reactionsList:
				for plugItem in species.childItems():
					if isinstance(plugItem,plug):
						if lineList[1] == species.nodeBox.name and lineList[2]==plugItem.name:
							startPlug=plugItem
						if lineList[3] == species.nodeBox.name and lineList[4]==plugItem.name:
							endPlug=plugItem
			objectName='connection'+str(connectionsCounter)
			objectName=connection(startPlug,endPlug)
			connectionsList.append(objectName)
			connectionsCounter+=1
			AppWindow.canvas.addItem(objectName)
		if lineList[0]=='//KMCparams':
			#generate KMC parameters based on this line information
			#KMCParams=[100000,1,0.01,1,1]#total N of starting particles,repeats,store timestep,max time,volume
			KMCParams[0]=int(lineList[1])
			KMCParams[1]=int(lineList[2])
			KMCParams[2]=float(lineList[3])
			KMCParams[3]=float(lineList[4])
			KMCParams[4]=float(lineList[5])
	
		
#calculate the height and width of node box
def getNodeWH(textH,titleTextW,numberTextW):
	h=2.5*textH #2.5 because font "times" has leading -1 (text has one preceeding empty line)
	if numberTextW>titleTextW:
		w=numberTextW+15
	else:
		w=titleTextW+15
	if w<h:
		w=h #if text and number are short, make it a rectangle for aesthetic reasons
	return w,h

# get width of text
def getTextWidth(text):
	global nodeFont
	fontMetrics=QFontMetrics(nodeFont)
	w=fontMetrics.boundingRect(text).width()
	return w

# get height of text
def getTextHeight(text):
	global nodeFont
	fontMetrics=QFontMetrics(nodeFont)
	h=fontMetrics.boundingRect(text).height()
	return h

#create node
def createNode(tool,position):
	global reactionsCounter
	global speciesCounter
	global reactionsList
	global speciesList
	if tool=="unselected":
		AppWindow.statusBar().showMessage("No tool selected",5000)
	if tool=="species":
		objectName='species'+str(speciesCounter)
		objectTitle='S'+str(speciesCounter)
		objectName=speciesNode(position,objectTitle,1.0)#create node; no of molecules=10000
		speciesList.append(objectName)
		speciesCounter+=1
		AppWindow.canvas.addItem(objectName)
		AppWindow.canvas.update()
	if tool=="reactAtoB":
		objectName='reaction'+str(reactionsCounter)
		objectTitle='R'+str(reactionsCounter)
		objectName=reactionAtoBNode(position,objectTitle,10)#create node; no of molecules=10000
		reactionsList.append(objectName)
		reactionsCounter+=1
		AppWindow.canvas.addItem(objectName)
		AppWindow.canvas.update()
	if tool=="reactABtoC":
		objectName='reaction'+str(reactionsCounter)
		objectTitle='R'+str(reactionsCounter)
		objectName=reactionABtoCNode(position,objectTitle,10)#create node; no of molecules=10000
		reactionsList.append(objectName)
		reactionsCounter+=1
		AppWindow.canvas.addItem(objectName)
		AppWindow.canvas.update()
	if tool=="reactAtoBC":
		objectName='reaction'+str(reactionsCounter)
		objectTitle='R'+str(reactionsCounter)
		objectName=reactionAtoBCNode(position,objectTitle,10)#create node; no of molecules=10000
		reactionsList.append(objectName)
		reactionsCounter+=1
		AppWindow.canvas.addItem(objectName)
		AppWindow.canvas.update()
	if tool=="reactABtoCD":
		objectName='reaction'+str(reactionsCounter)
		objectTitle='R'+str(reactionsCounter)
		objectName=reactionABtoCDNode(position,objectTitle,10)#create node; no of molecules=10000
		reactionsList.append(objectName)
		reactionsCounter+=1
		AppWindow.canvas.addItem(objectName)
		AppWindow.canvas.update()
	if tool=="reactAtoBCD":
		objectName='reaction'+str(reactionsCounter)
		objectTitle='R'+str(reactionsCounter)
		objectName=reactionAtoBCDNode(position,objectTitle,10)#create node; no of molecules=10000
		reactionsList.append(objectName)
		reactionsCounter+=1
		AppWindow.canvas.addItem(objectName)
		AppWindow.canvas.update()
	if tool=="reactABtoCDE":
		objectName='reaction'+str(reactionsCounter)
		objectTitle='R'+str(reactionsCounter)
		objectName=reactionABtoCDENode(position,objectTitle,10)#create node; no of molecules=10000
		reactionsList.append(objectName)
		reactionsCounter+=1
		AppWindow.canvas.addItem(objectName)
		AppWindow.canvas.update()


#create connection
def createConnection():
	global isConnecting
	global connectionStart
	global connectionEnd
	global connectionsList
	global connectionsCounter
	if isinstance(connectionStart,plug) and isinstance(connectionEnd,plug):
		legalConnection=True
	#test if valid connection
	if isConnecting==True:
		#put in checks so that only legal connections are allowed
		#1. You cannot connect node to itself
		if connectionStart.parentItem()==connectionEnd.parentItem():
			AppWindow.statusBar().showMessage("You cannot connect node to itself!",5000)
			legalConnection=False
		#2. You cannot create multiple connections between same plug pairs
		for existingConnection in connectionsList:
			if existingConnection.startNode==connectionStart and existingConnection.endNode ==connectionEnd:
				AppWindow.statusBar().showMessage("You cannot create multiple connections between same plug pairs!",5000)
				legalConnection=False
			if existingConnection.endNode==connectionStart and existingConnection.startNode ==connectionEnd:
				AppWindow.statusBar().showMessage("You cannot create multiple connections between same plug pairs!",5000)
				legalConnection=False
		#3. You can only connect different type of plugs
		if connectionStart.mode==connectionEnd.mode:
			AppWindow.statusBar().showMessage("You can only connect different type of plugs!",5000)
			legalConnection=False
		#4. You can only connect reactions to species and vice versa
		if connectionStart.parentItem().nodeType==connectionEnd.parentItem().nodeType:
			AppWindow.statusBar().showMessage("You can only connect reactions to species and vice versa!",5000)
			legalConnection=False
		#5. Only one connection is allowed per reaction plug
		if connectionStart.parentItem().nodeType=='reaction':
			for existingConnection in connectionsList:
				if existingConnection.startNode==connectionStart or existingConnection.endNode ==connectionStart:
					AppWindow.statusBar().showMessage("Only one connection is allowed per reaction plug!",5000)
					legalConnection=False
		if connectionEnd.parentItem().nodeType=='reaction':
			for existingConnection in connectionsList:
				if existingConnection.startNode==connectionEnd or existingConnection.endNode ==connectionEnd:
					AppWindow.statusBar().showMessage("Only one connection is allowed per reaction plug!",5000)
					legalConnection=False

	if legalConnection==True:
		#actually create the connection
		objectName='connection'+str(connectionsCounter)
		objectName=connection(connectionStart,connectionEnd)
		connectionsList.append(objectName)
		connectionsCounter+=1
		AppWindow.canvas.addItem(objectName)
		AppWindow.canvas.update()
		isConnecting=False
		connectionStart=None
		connectionEnd=None

def runKMC():
	global lasttVector
	global lastPVector
	global lastCVector
	global KMCParams
	global speciesList
	global reactionsList
	global connectionsList
	global fileName
	if fileName==None:
		AppWindow.saveFile()
	
	# generate population vector
	populationVector=calcPopulationVector(KMCParams[0])
	# generate rate constants vector
	rateConstantsVector=np.array([])
	for reaction in reactionsList:
		rateConstantsVector=np.append(rateConstantsVector,reaction.nodeBox.number)
	#generate connectivity matrix
	connectivityMatrix=np.zeros(shape=(len(reactionsList),len(speciesList)))

	#iterate over all reactions (rows in connectivity matrix)
	i=0
	while i < len(reactionsList):
		#iterate over plugs of reaction
		for reactionChildItem in reactionsList[i].childItems():
			if isinstance(reactionChildItem,plug):
				#iterate over all species (columns in connectivity matrix)
				j=0
				while j< len(speciesList):
					#check all plugs in given species
					for speciesChildItem in speciesList[j].childItems():
						if isinstance(speciesChildItem,plug):
							#iterate over all connections to check if connection exists 
							for connection in connectionsList:
								#check if connection's start and end plugs are identical to current reaction and species plugs
								if connection.startNode==reactionChildItem and connection.endNode==speciesChildItem:
									if reactionChildItem.mode=="in":
										connectivityMatrix[i][j]-=1
									if reactionChildItem.mode=="out":
										connectivityMatrix[i][j]+=1
								if connection.endNode==reactionChildItem and connection.startNode==speciesChildItem:
									if reactionChildItem.mode=="in":
										connectivityMatrix[i][j]-=1
									if reactionChildItem.mode=="out":
										connectivityMatrix[i][j]+=1
					j+=1	

		i+=1
	
	lasttVector,lastPVector=KMC_engine.runKMC(populationVector,rateConstantsVector,connectivityMatrix,KMCParams[4],KMCParams[1],KMCParams[2],KMCParams[3])
	#calculate concentration vector from population vector and volume
	if len(lasttVector)>0 and len(lastPVector)>0:
		lastCVector=np.empty(shape=lastPVector.shape)
		x=0
		while x< len(lastCVector[:,0]):
			y=0
			while y< len(lastCVector[x,:]):
				lastCVector[x,y]=lastPVector[x,y]/(Na*KMCParams[4])
				y+=1
			x+=1
			
		
	#write population vector output file
	outPopFileName=fileName+'_population.csv'
	
	
	POPOUT=open(outPopFileName, 'w')
	POPOUT.write("t ")
	
	for species in speciesList:
		POPOUT.write(str(species.nodeBox.name)+" ")
	POPOUT.write("\n")
	for i in range(len(lasttVector)):
		POPOUT.write(str(lasttVector[i])+" ")
		for j in range(len(lastPVector[i])):
			POPOUT.write(str(lastPVector[i,j])+" ")
		POPOUT.write("\n")
	POPOUT.close()
	#write concentration vector output file
	outConcFileName=fileName+'_concentration.csv'
	
	CONCOUT=open(outConcFileName, 'w')
	CONCOUT.write("t ")
	
	for species in speciesList:
		CONCOUT.write(str(species.nodeBox.name)+" ")
	CONCOUT.write("\n")
	for i in range(len(lasttVector)):
		CONCOUT.write(str(lasttVector[i])+" ")
		for j in range(len(lastCVector[i])):
			CONCOUT.write(str(lastCVector[i,j])+" ")
		CONCOUT.write("\n")
	CONCOUT.close()
	
	

###############
### Classes ###
###############
class PlotWindow(QMainWindow):
	def __init__(self):
		super(PlotWindow,self).__init__()
		global lasttVector
		global lastPVector
		global speciesList
		global lastCVector
		
		#set central widget
		self.centralWidget=QWidget()
		self.setCentralWidget(self.centralWidget)
		#set layout
		self.HLayout=QHBoxLayout()
		self.centralWidget.setLayout(self.HLayout)
		self.setWindowTitle("Plotting results")
		#generate figure canvas
		self.fig=Figure((10.0,12.0),dpi=100)
		self.canvas=FigureCanvas(self.fig)
		self.canvas.setParent(self)
		self.axes=self.fig.add_subplot(111)
		#add matplotlib standard toolbar
		self.matPlotLibToolbar=NavigationToolbar(self.canvas,self)
		self.plotLayout=QVBoxLayout()
		self.HLayout.addLayout(self.plotLayout)
		self.plotLayout.addWidget(self.canvas)
		self.plotLayout.addWidget(self.matPlotLibToolbar)
		#set layout for buttons and list
		self.VLayout=QVBoxLayout()
		self.HLayout.addLayout(self.VLayout)
		#add listview for selecting data series
		self.listView=QListView()
		self.listModel=QStandardItemModel()
		self.createDataSeries()
		self.listView.setModel(self.listModel)
		self.VLayout.addWidget(self.listView)
		#add button to display graph
		self.testButton=QPushButton("show")
		self.VLayout.addWidget(self.testButton)
		self.testButton.clicked.connect(self.onShow)
	def createDataSeries(self):
		self.listModel.clear()
		for species in speciesList:
			item=QStandardItem(species.nodeBox.name)
			item.setCheckState(Qt.Checked)
			item.setCheckable(True)
			self.listModel.appendRow(item)
	def onShow(self):
		self.axes.clear()
		i=0
		if len(lastCVector)!=0:
			numberOfSpecies=len(lastCVector[0,:])
			for row in range(self.listModel.rowCount()):
				index=self.listModel.index(row,0)
				if self.listModel.data(index,Qt.CheckStateRole)==QVariant(Qt.Checked):
					self.axes.scatter(lasttVector,lastCVector[:,i],label=speciesList[i].nodeBox.name)
				i+=1
		self.axes.legend(loc='best')
		self.canvas.draw()

class confirmWindow(QDialog):
	def __init__(self,text):
		super(confirmWindow,self).__init__()
		self.setGeometry(100,100,50,50)
		#set vertical layout
		self.VLayout=QVBoxLayout()
		self.setLayout(self.VLayout)
		#display text
		self.textDisplay=QLabel()
		self.textDisplay.setText(text)
		self.VLayout.addWidget(self.textDisplay)
		#create horizontal layout
		self.HLayout=QHBoxLayout()
		self.VLayout.addLayout(self.HLayout)
		#create OK button		
		self.OKButton=QPushButton("OK",self)
		self.OKButton.clicked.connect(self.OKPressed)
		self.HLayout.addWidget(self.OKButton)
		#create OK button		
		self.CancelButton=QPushButton("Cancel",self)
		self.CancelButton.clicked.connect(self.CancelPressed)
		self.HLayout.addWidget(self.CancelButton)
		
		#display window
		self.exec()
	def OKPressed(self,pressed):
		self.close()
	def CancelPressed(self,pressed):
		self.close()
# class for editing KMC parameters
class editKMCParams(QDialog):
	def __init__(self):
		super(editKMCParams,self).__init__()
		self.setGeometry(100,100,400,200)
		global KMCParams
		#KMCParams=[1,1,0.01,1]#volume,repeats,timestep,max time
		#set vertical layout
		self.VLayout=QVBoxLayout()
		self.setLayout(self.VLayout)

		#layout for time interval
		self.tIntervalLine=QHBoxLayout()
		self.VLayout.addLayout(self.tIntervalLine)
		self.tIntervalLabel=QLabel()
		self.tIntervalLabel.setText("Timestep for data storage (s):")
		self.tIntervalLine.addWidget(self.tIntervalLabel)
		self.tIntervalEdit=QLineEdit()
		self.tIntervalEdit.setText(str(KMCParams[2]))
		self.tIntervalLine.addWidget(self.tIntervalEdit)
			
		#layout for maximum allowed time
		self.maxTLine=QHBoxLayout()
		self.VLayout.addLayout(self.maxTLine)
		self.maxTLabel=QLabel()
		self.maxTLabel.setText("Max simulation time (s):")
		self.maxTLine.addWidget(self.maxTLabel)
		self.maxTEdit=QLineEdit()
		self.maxTEdit.setText(str(KMCParams[3]))
		self.maxTLine.addWidget(self.maxTEdit)
		
		#layout for total number of starting particles
		self.totalParticlesLine=QHBoxLayout()
		self.VLayout.addLayout(self.totalParticlesLine)
		self.totalParticlesLabel=QLabel()
		self.totalParticlesLabel.setText("Total number of starting molecules:")
		self.totalParticlesLine.addWidget(self.totalParticlesLabel)
		self.totalParticlesEdit=QLineEdit()
		self.totalParticlesEdit.setText(str(KMCParams[0]))
		self.totalParticlesLine.addWidget(self.totalParticlesEdit)

		#layout for repeats
		self.repeatsLine=QHBoxLayout()
		self.VLayout.addLayout(self.repeatsLine)
		self.repeatsLabel=QLabel()
		self.repeatsLabel.setText("Number of repeats:")
		self.repeatsLine.addWidget(self.repeatsLabel)
		self.repeatsEdit=QLineEdit()
		self.repeatsEdit.setText(str(KMCParams[1]))
		self.repeatsLine.addWidget(self.repeatsEdit)
		
		#layout for displaying volume
		self.VolumeLine=QHBoxLayout()
		self.VLayout.addLayout(self.VolumeLine)
		self.VolumeLabel=QLabel()
		self.VolumeLabel.setText("Simulation volume (L):")
		self.VolumeLine.addWidget(self.VolumeLabel)
		self.VolumeValue=QLabel()
		popList=calcPopulationVector(KMCParams[0])
		self.VolumeValue.setText(str(KMCParams[4]))
		self.VolumeLine.addWidget(self.VolumeValue)

		
		#layout for buttons line
		self.ButtonsLine=QHBoxLayout()
		self.VLayout.addLayout(self.ButtonsLine)
		
		#create OK button
		self.OKButton=QPushButton("OK",self)
		self.OKButton.clicked.connect(self.OKPressed)
		self.ButtonsLine.addWidget(self.OKButton)
		#create Cancel button
		self.CancelButton=QPushButton("Cancel",self)
		self.CancelButton.clicked.connect(self.CancelPressed)
		self.ButtonsLine.addWidget(self.CancelButton)
		

		#launch window
		self.exec()
		
	def OKPressed(self,pressed):
		source=self.sender()
		validOutput=True
		global KMCParams
		try:
			float(self.tIntervalEdit.text())
		except:
			invalidWindow=QMessageBox.information(self,"Error","time interval must be a number")
			validOutput=False
		try:
			float(self.maxTEdit.text())
		except:
			invalidWindow=QMessageBox.information(self,"Error","maximum time must be a number")
			validOutput=False
		try:
			int(self.totalParticlesEdit.text())
		except:
			invalidWindow=QMessageBox.information(self,"Error","total number of starting molecules must be an integer")
			validOutput=False
		try:
			int(self.repeatsEdit.text())
		except:
			invalidWindow=QMessageBox.information(self,"Error","number of repats must be an integer")
			validOutput=False
		if validOutput==True:
			#KMCParams=[100000,1,0.01,1,1]#total N of starting particles,repeats,store timestep,max time,volume
			KMCParams[0]=int(self.totalParticlesEdit.text())
			KMCParams[1]=int(self.repeatsEdit.text())
			KMCParams[2]=float(self.tIntervalEdit.text())
			KMCParams[3]=float(self.maxTEdit.text())
			self.close()
			
			
	def CancelPressed(self,pressed):
		#do nothing when cancel is pressed - delete widget and do not save changes
		self.close()

		
# class for editing node objects
class editNodes(QDialog):
	def __init__(self,node,nodeType,name,number):
		super(editNodes,self).__init__()
		self.setGeometry(100,100,200,150)
		
		self.originNode=node
		self.originType=nodeType
		#set vertical layout
		self.VLayout=QVBoxLayout()
		self.setLayout(self.VLayout)

		#layout for node Text and edit
		self.textLine=QHBoxLayout()
		self.VLayout.addLayout(self.textLine)
		self.nameLabel=QLabel()
		self.nameLabel.setText("Name:")
		self.textLine.addWidget(self.nameLabel)
		self.nameEdit=QLineEdit()
		self.nameEdit.setText(name)
		self.textLine.addWidget(self.nameEdit)
			
		#layout for node number and edit
		self.numberLine=QHBoxLayout()
		self.VLayout.addLayout(self.numberLine)
		self.numberLabel=QLabel()
		if nodeType=="species":
			self.numberLabel.setText("Concentration (mol/L):")
		if nodeType=="reaction":
			self.numberLabel.setText("Rate constant:")

		self.numberLine.addWidget(self.numberLabel)
		self.numberEdit=QLineEdit()
		self.numberEdit.setText(str(number))
		self.numberLine.addWidget(self.numberEdit)

		#layout for buttons line
		self.ButtonsLine=QHBoxLayout()
		self.VLayout.addLayout(self.ButtonsLine)
		
		#create OK button
		self.OKButton=QPushButton("OK",self)
		self.OKButton.clicked.connect(self.OKPressed)
		self.ButtonsLine.addWidget(self.OKButton)
		#create Cancel button
		self.CancelButton=QPushButton("Cancel",self)
		self.CancelButton.clicked.connect(self.CancelPressed)
		self.ButtonsLine.addWidget(self.CancelButton)
		

		#launch window
		self.exec()
		
	def OKPressed(self,pressed):
		source=self.sender()
		self.validName=True
		self.validNumber=True
		global reactionsList
		global speciesList
		#check if name field is unique
		for species in speciesList:
			if species.nodeBox==self.originNode:
				pass
			elif species.nodeBox.name==self.nameEdit.text():
				invalidWindow=QMessageBox.information(self,"Error","name in use")
				self.validName=False
		for reactions in reactionsList:
			if reactions.nodeBox==self.originNode:
				pass
			elif reactions.nodeBox.name==self.nameEdit.text():
				invalidWindow=QMessageBox.information(self,"Error","name in use")
				self.validName=False
		#check if number is int or float (for species and reaction, respectively)
		if self.originType=="species":
			try:
				float(self.numberEdit.text())
				#check if concentration is positive or 0
				if float(self.numberEdit.text()) <0:
					self.validNumber=False
					invalidWindow=QMessageBox.information(self,"Error","concentration must be positive or 0")

			except ValueError:
				invalidWindow=QMessageBox.information(self,"Error","concentration must be floating point number")
				self.validNumber=False		
		if self.originType=="reaction":
			try:
				float(self.numberEdit.text())
				#check if rate constant is positive or 0
				if float(self.numberEdit.text()) <0:
					self.validNumber=False
					invalidWindow=QMessageBox.information(self,"Error","rate constant must be positive or 0")

			except ValueError:
				invalidWindow=QMessageBox.information(self,"Error","rate constant must be integer or floating point number")
				self.validNumber=False
				
		#if all values are acceptable, save changes and close widget
		if self.validName==True and self.validNumber==True:
			self.originNode.name=self.nameEdit.text()
			if self.originType=="species":
				self.originNode.number=float(self.numberEdit.text())
			if self.originType=="reaction":
				self.originNode.number=float(self.numberEdit.text())
			self.originNode.updateNode()
			self.close()
			
			
	def CancelPressed(self,pressed):
		#do nothing when cancel is pressed - delete widget and no not save changes
		self.close()

		
# general class for all node objects
class nodeBox(QGraphicsItem):
	def __init__(self,parent,position,objectTitle,number):
		global nodeFont
		global plugSideLength
		global plugWidth
		global reactionsList
		global speciesList
		self.parent=parent
		super(nodeBox,self).__init__()
		self.createNode(self,position,objectTitle,number)
	def createNode(self,parent,position,objectTitle,number):
		#add central box
		self.name=objectTitle
		self.number=number
		self.textH=getTextHeight(self.name)
		self.titleTextW=getTextWidth(self.name)
		self.numberTextW=getTextWidth(str(self.number))
		self.nodeBoxW,self.nodeBoxH=getNodeWH(self.textH,self.titleTextW,self.numberTextW)
		#calculate node width and height
		self.width=self.nodeBoxW+2*plugWidth
		self.height=self.nodeBoxH
		#move node center to cursor location
		self.setPos(0,0) 
	def boundingRect(self):
		return QRectF(0,0,self.width,self.height)
	def paint(self,painter,option,widget):
		global nodeFont
		global plugWidth
		painter.setRenderHint(QPainter.Antialiasing)
		rect=QPainterPath()
		brush=QBrush(self.boxColor)
		rect.addRoundedRect(QRectF(0+plugWidth,0,self.nodeBoxW,self.nodeBoxH),10,10)
		
		painter.setPen(QPen(Qt.SolidLine))
		painter.setFont(nodeFont)
		painter.fillPath(rect,self.boxColor)
		painter.drawPath(rect)
		#painter.fillRect(0+plugWidth,0,self.nodeBoxW,self.nodeBoxH,QBrush(self.boxColor))
		painter.drawText(int(0+plugWidth+self.nodeBoxW*0.5-self.titleTextW*0.5),0+self.textH,self.name)
		painter.drawText(int(0+plugWidth+self.nodeBoxW*0.5-self.numberTextW*0.5),0+2*self.textH,str(self.number))
		
		if self.parent.selected==True:
			painter.setPen(QPen(Qt.DashLine))
			painter.drawRect(self.boundingRect())
		self.update()
	def contextMenuEvent(self,event):
		menu=QMenu()
		editAction=QAction('Edit',None)
		editAction.triggered.connect(self.editNode)
		menu.addAction(editAction)
		deleteAction=QAction('Delete',None)
		deleteAction.triggered.connect(self.deleteNode)
		menu.addAction(deleteAction)
		if self.parentItem().selected==True:
			selectionText='Unselect'
		else:
			selectionText='Select'
		selectAction=QAction(selectionText,None)
		selectAction.triggered.connect(self.selectNode)
		menu.addAction(selectAction)
		menu.exec_(event.screenPos())
	def editNode(self):
		editWidget=editNodes(self,self.parent.nodeType,self.name, self.number)
		
	def selectNode(self):
		if self.parentItem().selected==True:
			self.parentItem().selected=False
		else:
			self.parentItem().selected=True
	def deleteNode(self):
		self.deleteList=[]
		#clean up all connections related to this node
		for connection in connectionsList:
			if connection.startNode.parentItem() == self.parentItem() or connection.endNode.parentItem() == self.parentItem():
				connection.selected=True
				self.deleteList.append(connection)
		for connection in self.deleteList:
			connectionsList.remove(connection)
			AppWindow.canvas.removeItem(connection)
				
				
		#if parent object is species, clean up speciesList
		if isinstance(self.parentItem(),speciesNode):
			for node in speciesList:
				if node==self.parentItem():
					speciesList.remove(node)
		
		#if parent object is reaction, clean up reactionsList		
		if isinstance(self.parentItem(),reactionAtoBNode) or isinstance(self.parentItem(),reactionAtoBCNode) or isinstance(self.parentItem(),reactionABtoCNode):
			for node in reactionsList:
				if node==self.parentItem():
					reactionsList.remove(node)
	
		AppWindow.canvas.removeItem(self.parentItem())
	def updateNode(self):
		self.textH=getTextHeight(self.name)
		self.titleTextW=getTextWidth(self.name)
		self.numberTextW=getTextWidth(str(self.number))
		self.nodeBoxW,self.nodeBoxH=getNodeWH(self.textH,self.titleTextW,self.numberTextW)
		self.width=self.nodeBoxW+2*plugWidth
		self.height=self.nodeBoxH
		#update position of outgoing plugs
		for item in self.parentItem().childItems():
			if isinstance(item,plug) and item.mode=="out":
				item.x=self.width-plugWidth
				item.updateCoords()
			
#species node class
class speciesNode(QGraphicsItem):
	def __init__(self,position,objectTitle,number):
		super(speciesNode,self).__init__()
		global nodeFont
		global plugSideLength
		global plugWidth
		global speciesFillColor
		self.selected=False
		self.nodeType="species"
		self.createNode(position,objectTitle,number)
		self.setZValue(1)
	def createNode(self,position,objectTitle,number):
		self.nodeBox=nodeBox(self,position,objectTitle,number)
		self.nodeBox.setParentItem(self)
		self.nodeBox.boxColor=speciesFillColor
		self.setPos(int(position.x()-self.nodeBox.width/2),int(position.y()-self.nodeBox.height/2))
		self.nodePlugIn=plug(0,self.nodeBox.height/2-plugSideLength/2,"in","in1")
		self.nodePlugIn.setParentItem(self)
		self.nodePlugOut=plug(self.nodeBox.width-plugWidth,self.nodeBox.height/2-plugSideLength/2,"out","out1")
		self.nodePlugOut.setParentItem(self)
	def boundingRect(self):
		return self.nodeBox.boundingRect()
	def updateCoords(self,position):
		self.setPos(position.x()-self.nodeBox.width/2,position.y()-self.nodeBox.height/2)
	def paint(self,painter,option,widget):
		pass
	

#reaction A to B node class
class reactionAtoBNode(QGraphicsItem):
	def __init__(self,position,objectTitle,number):
		super(reactionAtoBNode,self).__init__()
		global nodeFont
		global plugSideLength
		global plugWidth
		global reactionsFillColor
		self.selected=False
		self.nodeType="reaction"
		self.createNode(position,objectTitle,number)
		self.setZValue(1)
	def createNode(self,position,objectTitle,number):
		self.nodeBox=nodeBox(self,position,objectTitle,number)
		self.nodeBox.setParentItem(self)
		self.nodeBox.boxColor=reactionsFillColor
		self.setPos(int(position.x()-self.nodeBox.width/2),int(position.y()-self.nodeBox.height/2))
		self.nodePlugIn=plug(0,self.nodeBox.height/2-plugSideLength/2,"in","in1")
		self.nodePlugIn.setParentItem(self)
		self.nodePlugOut=plug(self.nodeBox.width-plugWidth,self.nodeBox.height/2-plugSideLength/2,"out","out1")
		self.nodePlugOut.setParentItem(self)
	def boundingRect(self):
		return self.nodeBox.boundingRect()
	def updateCoords(self,position):
		self.setPos(position.x()-self.nodeBox.width/2,position.y()-self.nodeBox.height/2)
	def paint(self,painter,option,widget):
		pass

#reaction AB to C node class
class reactionABtoCNode(QGraphicsItem):
	def __init__(self,position,objectTitle,number):
		super(reactionABtoCNode,self).__init__()
		global nodeFont
		global plugSideLength
		global plugWidth
		global reactionsFillColor
		self.selected=False
		self.nodeType="reaction"
		self.createNode(position,objectTitle,number)
		self.setZValue(1)
	def createNode(self,position,objectTitle,number):
		self.nodeBox=nodeBox(self,position,objectTitle,number)
		self.nodeBox.setParentItem(self)
		self.nodeBox.boxColor=reactionsFillColor
		self.setPos(int(position.x()-self.nodeBox.width/2),int(position.y()-self.nodeBox.height/2))
		self.nodePlugIn1=plug(0,(self.nodeBox.height-2*plugSideLength)/3,"in","in1")
		self.nodePlugIn1.setParentItem(self)
		self.nodePlugIn2=plug(0,plugSideLength+2*(self.nodeBox.height-2*plugSideLength)/3,"in","in2")
		self.nodePlugIn2.setParentItem(self)
		self.nodePlugOut=plug(self.nodeBox.width-plugWidth,self.nodeBox.height/2-plugSideLength/2,"out","out1")
		self.nodePlugOut.setParentItem(self)
	def boundingRect(self):
		return self.nodeBox.boundingRect()
	def updateCoords(self,position):
		self.setPos(position.x()-self.nodeBox.width/2,position.y()-self.nodeBox.height/2)
	def paint(self,painter,option,widget):
		pass

#reaction A to BC node class
class reactionAtoBCNode(QGraphicsItem):
	def __init__(self,position,objectTitle,number):
		super(reactionAtoBCNode,self).__init__()
		global nodeFont
		global plugSideLength
		global plugWidth
		global reactionsFillColor
		self.selected=False
		self.nodeType="reaction"
		self.createNode(position,objectTitle,number)
		self.setZValue(1)
	def createNode(self,position,objectTitle,number):
		self.nodeBox=nodeBox(self,position,objectTitle,number)
		self.nodeBox.setParentItem(self)
		self.nodeBox.boxColor=reactionsFillColor
		self.setPos(int(position.x()-self.nodeBox.width/2),int(position.y()-self.nodeBox.height/2))
		self.nodePlugIn=plug(0,self.nodeBox.height/2-plugSideLength/2,"in","in1")
		self.nodePlugIn.setParentItem(self)
		self.nodePlugOut1=plug(self.nodeBox.width-plugWidth,(self.nodeBox.height-2*plugSideLength)/3,"out","out1")
		self.nodePlugOut1.setParentItem(self)
		self.nodePlugOut2=plug(self.nodeBox.width-plugWidth,plugSideLength+2*(self.nodeBox.height-2*plugSideLength)/3,"out","out2")
		self.nodePlugOut2.setParentItem(self)
	def boundingRect(self):
		return self.nodeBox.boundingRect()
	def updateCoords(self,position):
		self.setPos(position.x()-self.nodeBox.width/2,0+position.y()-self.nodeBox.height/2)
	def paint(self,painter,option,widget):
		pass
#reaction AB to CD node class
class reactionABtoCDNode(QGraphicsItem):
	def __init__(self,position,objectTitle,number):
		super(reactionABtoCDNode,self).__init__()
		global nodeFont
		global plugSideLength
		global plugWidth
		global reactionsFillColor
		self.selected=False
		self.nodeType="reaction"
		self.createNode(position,objectTitle,number)
		self.setZValue(1)
	def createNode(self,position,objectTitle,number):
		self.nodeBox=nodeBox(self,position,objectTitle,number)
		self.nodeBox.setParentItem(self)
		self.nodeBox.boxColor=reactionsFillColor
		self.setPos(int(position.x()-self.nodeBox.width/2),int(position.y()-self.nodeBox.height/2))
		self.nodePlugIn1=plug(0,(self.nodeBox.height-2*plugSideLength)/3,"in","in1")
		self.nodePlugIn1.setParentItem(self)
		self.nodePlugIn2=plug(0,plugSideLength+2*(self.nodeBox.height-2*plugSideLength)/3,"in","in2")
		self.nodePlugIn2.setParentItem(self)
		self.nodePlugOut1=plug(self.nodeBox.width-plugWidth,(self.nodeBox.height-2*plugSideLength)/3,"out","out1")
		self.nodePlugOut1.setParentItem(self)
		self.nodePlugOut2=plug(self.nodeBox.width-plugWidth,plugSideLength+2*(self.nodeBox.height-2*plugSideLength)/3,"out","out2")
		self.nodePlugOut2.setParentItem(self)
	def boundingRect(self):
		return self.nodeBox.boundingRect()
	def updateCoords(self,position):
		self.setPos(position.x()-self.nodeBox.width/2,0+position.y()-self.nodeBox.height/2)
	def paint(self,painter,option,widget):
		pass
#reaction A to BCD node class
class reactionAtoBCDNode(QGraphicsItem):
	def __init__(self,position,objectTitle,number):
		super(reactionAtoBCDNode,self).__init__()
		global nodeFont
		global plugSideLength
		global plugWidth
		global reactionsFillColor
		self.selected=False
		self.nodeType="reaction"
		self.createNode(position,objectTitle,number)
		self.setZValue(1)
	def createNode(self,position,objectTitle,number):
		self.nodeBox=nodeBox(self,position,objectTitle,number)
		self.nodeBox.setParentItem(self)
		self.nodeBox.boxColor=reactionsFillColor
		self.setPos(int(position.x()-self.nodeBox.width/2),int(position.y()-self.nodeBox.height/2))
		self.nodePlugIn1=plug(0,(self.nodeBox.height-2*plugSideLength)/3,"in","in1")
		self.nodePlugIn=plug(0,self.nodeBox.height/2-plugSideLength/2,"in","in1")
		self.nodePlugIn.setParentItem(self)
		self.nodePlugOut1=plug(self.nodeBox.width-plugWidth,(self.nodeBox.height-3*plugSideLength)/4,"out","out1")
		self.nodePlugOut1.setParentItem(self)
		self.nodePlugOut2=plug(self.nodeBox.width-plugWidth,plugSideLength+2*(self.nodeBox.height-3*plugSideLength)/4,"out","out2")
		self.nodePlugOut2.setParentItem(self)
		self.nodePlugOut3=plug(self.nodeBox.width-plugWidth,2*plugSideLength+3*(self.nodeBox.height-3*plugSideLength)/4,"out","out3")
		self.nodePlugOut3.setParentItem(self)
		
	def boundingRect(self):
		return self.nodeBox.boundingRect()
	def updateCoords(self,position):
		self.setPos(position.x()-self.nodeBox.width/2,0+position.y()-self.nodeBox.height/2)
	def paint(self,painter,option,widget):
		pass	
#reaction AB to CDE node class
class reactionABtoCDENode(QGraphicsItem):
	def __init__(self,position,objectTitle,number):
		super(reactionABtoCDENode,self).__init__()
		global nodeFont
		global plugSideLength
		global plugWidth
		global reactionsFillColor
		self.selected=False
		self.nodeType="reaction"
		self.createNode(position,objectTitle,number)
		self.setZValue(1)
	def createNode(self,position,objectTitle,number):
		self.nodeBox=nodeBox(self,position,objectTitle,number)
		self.nodeBox.setParentItem(self)
		self.nodeBox.boxColor=reactionsFillColor
		self.setPos(int(position.x()-self.nodeBox.width/2),int(position.y()-self.nodeBox.height/2))
		self.nodePlugIn1=plug(0,(self.nodeBox.height-2*plugSideLength)/3,"in","in1")
		self.nodePlugIn1.setParentItem(self)
		self.nodePlugIn2=plug(0,plugSideLength+2*(self.nodeBox.height-2*plugSideLength)/3,"in","in2")
		self.nodePlugIn2.setParentItem(self)
		self.nodePlugOut1=plug(self.nodeBox.width-plugWidth,(self.nodeBox.height-3*plugSideLength)/4,"out","out1")
		self.nodePlugOut1.setParentItem(self)
		self.nodePlugOut2=plug(self.nodeBox.width-plugWidth,plugSideLength+2*(self.nodeBox.height-3*plugSideLength)/4,"out","out2")
		self.nodePlugOut2.setParentItem(self)
		self.nodePlugOut3=plug(self.nodeBox.width-plugWidth,2*plugSideLength+3*(self.nodeBox.height-3*plugSideLength)/4,"out","out3")
		self.nodePlugOut3.setParentItem(self)
		
	def boundingRect(self):
		return self.nodeBox.boundingRect()
	def updateCoords(self,position):
		self.setPos(position.x()-self.nodeBox.width/2,0+position.y()-self.nodeBox.height/2)
	def paint(self,painter,option,widget):
		pass	
# class for plug items
class plug(QGraphicsItem):
	def __init__(self,x,y,mode,name):
		super(plug,self).__init__()
		self.x=x
		self.y=y
		self.mode=mode
		self.name=name
		global plugSideLength
		self.centre=QPointF(self.x+0.5*plugWidth,self.y+plugSideLength/2)
		self.triangle=QPolygonF()
		self.triangle.append(QPointF(self.x,self.y))
		self.triangle.append(QPointF(self.x,self.y+plugSideLength))
		self.triangle.append(QPointF(self.x+plugSideLength*(math.sqrt(3)/2),self.y+plugSideLength/2))
	
	def boundingRect(self):
		return QRectF(self.x,self.y,plugSideLength*(math.sqrt(3)/2),plugSideLength)
	def paint(self,painter,option,widget):
		painter.setBrush(QBrush(QColor(150,150,150)))
		painter.setPen(QPen(Qt.SolidLine))
		painter.drawPolygon(self.triangle)
		#painter.drawEllipse(self.centre,2,2)
		self.update()
		
	def updateCoords(self):
		self.triangle=QPolygonF()
		self.triangle.append(QPointF(self.x,self.y))
		self.triangle.append(QPointF(self.x,self.y+plugSideLength))
		self.triangle.append(QPointF(self.x+plugSideLength*(math.sqrt(3)/2),self.y+plugSideLength/2))
		self.centre=QPointF(self.x+0.5*plugWidth,self.y+plugSideLength/2)

# class for connection items
class connection(QGraphicsItem):
	def __init__(self,startNode,endNode):
		super(connection,self).__init__()
		self.startNode=startNode
		self.endNode=endNode
		self.selected=False
		self.itemIsMovable=True
		global lineInteractionRange
		global connectionsList
		self.setZValue=0
	def boundingRect(self):
		
		#get the top left corner coordinates and height/width of connection 
		if self.startNode.scenePos().x()+self.startNode.centre.x()<=self.endNode.scenePos().x()+self.endNode.centre.x():
			self.topCornerX=self.startNode.scenePos().x()+self.startNode.centre.x()-lineInteractionRange
			self.bottomCornerX=self.endNode.scenePos().x()+self.endNode.centre.x()+lineInteractionRange
		else:
			self.topCornerX=self.endNode.scenePos().x()+self.endNode.centre.x()-lineInteractionRange
			self.bottomCornerX=self.startNode.scenePos().x()+self.startNode.centre.x()+lineInteractionRange

		self.boundingRectHeight=self.bottomCornerX-self.topCornerX
		if self.startNode.scenePos().y()+self.startNode.centre.y()<=self.endNode.scenePos().y()+self.endNode.centre.y():
			self.topCornerY=self.startNode.scenePos().y()+self.startNode.centre.y()-lineInteractionRange
			self.bottomCornerY=self.endNode.scenePos().y()+self.endNode.centre.y()+lineInteractionRange
		else:
			self.topCornerY=self.endNode.scenePos().y()+self.endNode.centre.y()-lineInteractionRange
			self.bottomCornerY=self.startNode.scenePos().y()+self.startNode.centre.y()+lineInteractionRange
		self.boundingRectWidth=self.bottomCornerY-self.topCornerY
		return QRectF(self.topCornerX,self.topCornerY,self.boundingRectHeight,self.boundingRectWidth)
		
		self.update()
	def paint(self,painter,option,widget):
		#painter.setBrush(QBrush(QColor(200,150,150)))
		painter.setPen(QPen(Qt.SolidLine))
		painter.drawLine(self.startNode.scenePos().x()+self.startNode.centre.x(),self.startNode.scenePos().y()+self.startNode.centre.y(),self.endNode.scenePos().x()+self.endNode.centre.x(),self.endNode.scenePos().y()+self.endNode.centre.y())
		self.selectionArea=self.createSelectionArea()
		if self.selected==True:
			
			painter.setPen(QPen(Qt.DashLine))
			painter.drawPolygon(self.selectionArea)
		self.update()
	def createSelectionArea(self):
		
		if self.endNode.scenePos().y()+self.endNode.centre.y()-self.startNode.scenePos().y()-self.startNode.centre.y() ==0:
			slope =1
		else:
			slope=(self.endNode.scenePos().x()+self.endNode.centre.x()-self.startNode.scenePos().x()-self.startNode.centre.x())/(self.endNode.scenePos().y()+self.endNode.centre.y()-self.startNode.scenePos().y()-self.startNode.centre.y())
		slopeRadians=math.atan(slope)	
		mouseInteractionBox=QPolygonF()
		point1x=self.startNode.scenePos().x()+self.startNode.centre.x()+math.sqrt(lineInteractionRange**2/(1+slope**2))
		point1y=self.startNode.scenePos().y()+self.startNode.centre.y()+(-1)*slope*math.sqrt(lineInteractionRange**2/(1+slope**2))
		mouseInteractionBox.append(QPointF(point1x,point1y))
		point2x=self.startNode.scenePos().x()+self.startNode.centre.x()-math.sqrt(lineInteractionRange**2/(1+slope**2))
		point2y=self.startNode.scenePos().y()+self.startNode.centre.y()-(-1)*slope*math.sqrt(lineInteractionRange**2/(1+slope**2))
		mouseInteractionBox.append(QPointF(point2x,point2y))
		point3x=self.endNode.scenePos().x()+self.endNode.centre.x()-math.sqrt(lineInteractionRange**2/(1+slope**2))
		point3y=self.endNode.scenePos().y()+self.endNode.centre.y()-(-1)*slope*math.sqrt(lineInteractionRange**2/(1+slope**2))	
		mouseInteractionBox.append(QPointF(point3x,point3y))
		point4x=self.endNode.scenePos().x()+self.endNode.centre.x()+math.sqrt(lineInteractionRange**2/(1+slope**2))
		point4y=self.endNode.scenePos().y()+self.endNode.centre.y()+(-1)*slope*math.sqrt(lineInteractionRange**2/(1+slope**2))
		mouseInteractionBox.append(QPointF(point4x,point4y))
		return mouseInteractionBox
	def contextMenuEvent(self,event):
		menu=QMenu()
		deleteAction=QAction('Delete',None)
		deleteAction.triggered.connect(self.deleteConnection)
		menu.addAction(deleteAction)
		if self.selected==True:
			selectionText='Unselect'
		else:
			selectionText='Select'
		selectAction=QAction(selectionText,None)
		selectAction.triggered.connect(self.selectConnection)
		menu.addAction(selectAction)
		if self.selectionArea.containsPoint(event.scenePos(),Qt.OddEvenFill):
			menu.exec_(event.screenPos())
	def selectConnection(self):
		if self.selected==True:
			self.selected=False
		else:
			self.selected=True
	def deleteConnection(self):
		#clean up connectionsList
		for connection in connectionsList:
			if connection==self:
				connectionsList.remove(connection)
		AppWindow.canvas.removeItem(self)
selectableObjects=(speciesNode,reactionAtoBNode,reactionAtoBCNode,reactionABtoCNode,connection)		

# Graphics scene	
class DrawingArea(QGraphicsScene):
	def __init__(self,parent):
		super(DrawingArea,self).__init__(parent)
		self.setSceneRect(0,0,1000,1000)
	def mousePressEvent(self,event):
		global movingItem
		global isMoving
		global isConnecting
		global connectionStart
		global connectionEnd
		self.clickedItem=self.itemAt(event.scenePos(),QTransform())
		#if event.button()==Qt.RightButton:		
		
		if event.button()==Qt.LeftButton:
			self.connectionPresent=False
			self.nodePresent=False
			self.plugPresent=False
			
			#check what items are at mouse press position
			for items in self.items(event.scenePos()):
				if isinstance(items,plug)==True and self.plugPresent==False:
					self.plugPresent=True
				if isinstance(items,connection)==True and items.selectionArea.containsPoint(event.scenePos(),Qt.OddEvenFill)and self.connectionPresent==False:
					self.connectionPresent=True
				if isinstance(items,nodeBox)==True and self.nodePresent==False:
					self.nodePresent=True
					
			
			
			#if clicked on empty space with node creation tool, create that node
			if self.itemAt(event.scenePos(),QTransform()) == None and AppWindow.canvas.currentTool !="unselected":
				createNode(AppWindow.canvas.currentTool,event.scenePos())
				
			#if clicked on plug and not on connection, create connection
			if self.plugPresent == True and self.connectionPresent == False:
				
				if isinstance(self.itemAt(event.scenePos(),QTransform()),plug):
					
					isConnecting=True
					connectionStart=self.itemAt(event.scenePos(),QTransform())
				else:
					print("should create connection but itemAt scenePos is not plug --> zlevel issue")
			#if clicked on node and not on plug or connection, move node
			if self.nodePresent ==True and self.connectionPresent==False and self.plugPresent==False:
				if isinstance(self.itemAt(event.scenePos(),QTransform()),nodeBox)==True:
					isMoving=True
					movingItem=self.itemAt(event.scenePos(),QTransform()).parentItem()
				else:
					print("should be moving item but itemAt scenePos is not nodeBox -->zlevel issue")
			
	def mouseMoveEvent(self,event):
		global movingItem
		global connectionsList
		if movingItem != None:
			movingItem.updateCoords(event.scenePos())
			for connection in connectionsList:
				connection.prepareGeometryChange()
	def mouseReleaseEvent(self,event):
		global movingItem
		global isMoving
		global isConnecting
		global connectionEnd
		if isMoving==True:
			isMoving=False
			movingItem=None
		self.clickedItem=self.itemAt(event.scenePos(),QTransform())
		if isinstance(self.clickedItem,plug) and isConnecting==True:
			connectionEnd=self.clickedItem
			createConnection()
			
		
	def mouseDoubleClickEvent(self,event):
		global selectableObjects
		if event.button()==Qt.LeftButton:
			self.clickedItem=self.itemAt(event.scenePos(),QTransform())
			#if clicked on plug, get parent item
			if isinstance(self.clickedItem,nodeBox):
				self.clickedItem=self.clickedItem.parentItem()
			if isinstance(self.clickedItem,connection) and self.clickedItem.selectionArea.containsPoint(event.scenePos(),Qt.OddEvenFill):
				pass
			#if double click on selectable objects, toggle selection
			if isinstance(self.clickedItem,	selectableObjects):
				if self.clickedItem.selected==False:
					self.clickedItem.selected=True
					
				else:
					self.clickedItem.selected=False
					
				
				

#class for the main application window
class MainWindow(QMainWindow):
	def __init__(self):
		super(MainWindow,self).__init__()
		global TOOLS
		self.initUI()
		
				
				
		#add functionalities not defined in initUI module
		self.canvas.currentTool="unselected"


		#set canvas.currentTool based on the buttons pressed
		for tool in TOOLS:
			btn=getattr(self, '%sButton' % tool)
			btn.pressed.connect(lambda tool=tool: self.setTool(tool))
		
	
		
	def setTool(self,tool): #function for storing current tool
		self.canvas.currentTool=tool

	def clearCanvas(self): #function for clearing all objects - new document
		confirmWindow=QMessageBox.question(self, '', "Clear all objects, are you sure?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
		if confirmWindow == QMessageBox.Yes:
			global speciesCounter
			global speciesList
			global reactionsCounter
			global reactionsList
			global connectionsCounter
			global connectionsList
			global isMoving
			global movingItem
			global isConnecting
			global connectionStart
			global connectionEnd
			speciesCounter=1
			speciesList=[]
			reactionsCounter=1
			reactionsList=[]
			connectionsCounter=1
			connectionsList=[]
			isMoving=False
			movingItem=None
			isConnecting=False
			connectionStart=None
			connectionEnd=None
			self.canvas.clear()
	def saveFile(self): #function for saving all objects to file
		global fileName
		saveDialog=QFileDialog()
		saveDialog.setDefaultSuffix('kmc')
		saveDialog.setAcceptMode(QFileDialog.AcceptSave)
		saveDialog.setNameFilters(['kinetic Monte Carlo (*.kmc)'])
		saveDialog.setOptions(QFileDialog.DontUseNativeDialog)
		if saveDialog.exec_() == QDialog.Accepted:
			filename=saveDialog.selectedFiles()[0].split(".")[0]
			extension=saveDialog.selectedFiles()[0].split(".")[1]
			if extension != saveDialog.defaultSuffix():
				print('wrong extension, "',extension,'", correcting')
				saveFileName=filename+'.'+saveDialog.defaultSuffix()
				
			else:
				saveFileName=saveDialog.selectedFiles()[0]
			fileName=filename	
			#save all reactions
			outputStream=generateOutputStream()
			file=open(saveFileName,"w+")		
			for line in outputStream:
				file.write(str(line)+"\n")
			file.close()
			
			
	def loadFile(self): #function for loading objects from file
		global fileName
		global speciesCounter
		global lastCVector
		global lastPVector
		global lasttVector
		loadDialog=QFileDialog()
		loadDialog.setDefaultSuffix('kmc')
		loadDialog.setAcceptMode(QFileDialog.AcceptOpen)
		loadDialog.setNameFilters(['kinetic Monte Carlo (*.kmc)'])
		loadDialog.setOptions(QFileDialog.DontUseNativeDialog)	
		if loadDialog.exec_() == QDialog.Accepted:
			filename=loadDialog.selectedFiles()[0]
			fileName=loadDialog.selectedFiles()[0].split(".")[0]
			with open(filename,"r") as inputFile:
				inputStream=[]
				line=inputFile.readline()
				while line:
					line=line.rstrip()
					inputStream.append(line)
					line=inputFile.readline()
			self.clearCanvas()
			readInputStream(inputStream)
			#read previous simulation data from file if present
			populationFilename=fileName+"_population.csv"
			
			try:
				with open(populationFilename,"r") as popInput:
					popInStream=[]
					line=popInput.readline()
					while line:
						line=line.rstrip()
						popInStream.append(line)
						line=popInput.readline()
					timeVector=np.array([])
					popVector=np.array([])
					lineCounter=0
					for line in popInStream:
						lineList=line.split(" ")
						if lineCounter>0:
							timeVector=np.append(timeVector,float(lineList[0]))
							lineList.pop(0)
							if lineCounter==1:
								popVector=np.asarray(lineList)
							if lineCounter>1:
								popVector=np.vstack([popVector,np.asarray(lineList)])
							lineCounter+=1
						if lineCounter==0:
							nameVector=lineList
							lineCounter+=1
							if len(nameVector)!=speciesCounter+1:
								print("incompatible population vector file")
					lastPVector=popVector.astype(np.float)
					lasttVector=timeVector
					#calculate concentration vector from population vector and volume
					if len(lasttVector)>0 and len(lastPVector)>0:
						lastCVector=np.empty(shape=lastPVector.shape)
						x=0
						while x< len(lastCVector[:,0]):
							y=0
							while y< len(lastCVector[x,:]):
								lastCVector[x,y]=lastPVector[x,y]/(Na*KMCParams[4])
								y+=1
							x+=1	
			except:
				print("previous simulation data not available")
	def showPlot(self):
		self.plotWindow=PlotWindow()
		self.plotWindow.show()

	def initUI(self):
		
		global speciesList
		global reactionsList
		global connectionsList
		global lasttVector
		global lastPVector

		#super(MainWindow,self).__init__()
		self.resize(800, 800)
		self.centralwidget = QGraphicsView(self)
		self.centralwidget.setObjectName("centralwidget")
		self.setCentralWidget(self.centralwidget)
		#add QGraphicsScene widget to draw on
		self.canvas=DrawingArea(self)
		self.centralwidget.setScene(self.canvas)
		self.setMouseTracking(True)

		# build menubar
		self.mainMenu=QMainWindow.menuBar(self)
		self.mainMenu.setNativeMenuBar(False)
		# build file menu
		self.menuFile = self.mainMenu.addMenu('File')
		self.actionNew = QAction('New',self)
		self.actionSave=QAction('Save',self)
		self.actionOpen = QAction('Open',self)
		self.actionExit = QAction('Exit',self)
		self.menuFile.addAction(self.actionNew)
		self.menuFile.addAction(self.actionSave)
		self.menuFile.addAction(self.actionOpen)
		self.menuFile.addAction(self.actionExit)
		self.actionNew.triggered.connect(self.clearCanvas)
		self.actionSave.triggered.connect(self.saveFile)
		self.actionOpen.triggered.connect(self.loadFile)
		# build edit menu
		self.menuEdit = self.mainMenu.addMenu('Edit')
		self.actionEditReactTable = QAction('Reaction table',self)
		self.actionKMCParams = QAction('KMC parameters',self)
		self.actionRun = QAction('Run',self)
		#self.menuEdit.addAction(self.actionEditReactTable)
		self.menuEdit.addAction(self.actionKMCParams)
		self.menuEdit.addAction(self.actionRun)
		self.actionRun.triggered.connect(runKMC)
		self.actionKMCParams.triggered.connect(editKMC)
		self.plotResults=QAction('Plot results',self)
		self.menuEdit.addAction(self.plotResults)
		
		self.plotResults.triggered.connect(self.showPlot)
		
		#build toolbar
		self.toolBar = QToolBar()
		self.addToolBar(Qt.TopToolBarArea, self.toolBar)
		# add button for A-->B reaction
		self.reactAtoBButton = QPushButton(self)
		self.reactAtoBButton.setObjectName("reactAtoBButton")
		self.reactAtoBButton.setCheckable(True)
		self.AtoBIcon = QIcon()
		self.AtoBIcon.addPixmap(QPixmap("icons/AtoB.png"), QIcon.Normal, QIcon.Off)
		self.reactAtoBButton.setIcon(self.AtoBIcon)
		# add button for A+B-->C reaction
		self.reactABtoCButton = QPushButton(self)
		self.reactABtoCButton.setObjectName("reactABtoCButton")
		self.reactABtoCButton.setCheckable(True)
		self.ABtoCIcon = QIcon()
		self.ABtoCIcon.addPixmap(QPixmap("icons/ABtoC.png"), QIcon.Normal, QIcon.Off)
		self.reactABtoCButton.setIcon(self.ABtoCIcon)
		# add button for A-->B+C reaction
		self.reactAtoBCButton = QPushButton(self)
		self.reactAtoBCButton.setObjectName("reactAtoBCButton")
		self.reactAtoBCButton.setCheckable(True)
		self.AtoBCIcon = QIcon()
		self.AtoBCIcon.addPixmap(QPixmap("icons/AtoBC.png"), QIcon.Normal, QIcon.Off)
		self.reactAtoBCButton.setIcon(self.AtoBCIcon)
		# add button for A+B-->C+D reaction
		self.reactABtoCDButton = QPushButton(self)
		self.reactABtoCDButton.setObjectName("reactABtoCDButton")
		self.reactABtoCDButton.setCheckable(True)
		self.ABtoCDIcon = QIcon()
		self.ABtoCDIcon.addPixmap(QPixmap("icons/ABtoCD.png"), QIcon.Normal, QIcon.Off)
		self.reactABtoCDButton.setIcon(self.ABtoCDIcon)
		# add button for A-->B+C+D reaction
		self.reactAtoBCDButton = QPushButton(self)
		self.reactAtoBCDButton.setObjectName("reactAtoBCDButton")
		self.reactAtoBCDButton.setCheckable(True)
		self.AtoBCDIcon = QIcon()
		self.AtoBCDIcon.addPixmap(QPixmap("icons/AtoBCD.png"), QIcon.Normal, QIcon.Off)
		self.reactAtoBCDButton.setIcon(self.AtoBCDIcon)
		# add button for A+B-->C+D+E reaction
		self.reactABtoCDEButton = QPushButton(self)
		self.reactABtoCDEButton.setObjectName("reactABtoCDEButton")
		self.reactABtoCDEButton.setCheckable(True)
		self.ABtoCDEIcon = QIcon()
		self.ABtoCDEIcon.addPixmap(QPixmap("icons/ABtoCDE.png"), QIcon.Normal, QIcon.Off)
		self.reactABtoCDEButton.setIcon(self.ABtoCDEIcon)
		# add button for species
		self.speciesButton = QPushButton(self)
		self.speciesButton.setObjectName("speciesButton")
		self.speciesButton.setCheckable(True)
		self.speciesIcon = QIcon()
		self.speciesIcon.addPixmap(QPixmap("icons/species.png"), QIcon.Normal, QIcon.Off)
		self.speciesButton.setIcon(self.speciesIcon)
		
		#add buttons to toolbar
		self.toolBar.addWidget(self.reactAtoBButton)
		self.toolBar.addWidget(self.reactABtoCButton)
		self.toolBar.addWidget(self.reactAtoBCButton)
		self.toolBar.addWidget(self.reactABtoCDButton)
		self.toolBar.addWidget(self.reactAtoBCDButton)
		self.toolBar.addWidget(self.reactABtoCDEButton)
		self.toolBar.addWidget(self.speciesButton)
		#set tool buttons as exclusive
		self.toolGroup=QButtonGroup(self)
		self.toolGroup.setExclusive(True)
		self.toolGroup.addButton(self.reactAtoBButton)
		self.toolGroup.addButton(self.reactABtoCButton)
		self.toolGroup.addButton(self.reactAtoBCButton)
		self.toolGroup.addButton(self.reactABtoCDButton)
		self.toolGroup.addButton(self.reactAtoBCDButton)
		self.toolGroup.addButton(self.reactABtoCDEButton)
		self.toolGroup.addButton(self.speciesButton)
		
		
#initialize main app window if program is called
if __name__ == "__main__":
	import sys
	app = QApplication(sys.argv)
	AppWindow=MainWindow()
	AppWindow.show()
	sys.exit(app.exec_())



