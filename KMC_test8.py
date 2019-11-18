#! /usr/bin/python3

import sys
import random
import math
import numpy as np
from decimal import Decimal
import matplotlib.pyplot as plt

Na=6.02*10**23 #Avogadros constant [1/mol]

#function for calculating propensity of reaction A-->B
def propensity1(reaction,no_of_species,k_vec,p_vec,conn_matrix):
	c=0
	propensity=0
	while c<no_of_species:
		if conn_matrix[reaction,c]==-1:
			propensity=k_vec[reaction]*p_vec[c]
		c=c+1
	return propensity
#function for calculating propensity of reaction A+B-->C
def propensity2(reaction,no_of_species,k_vec,p_vec,V,conn_matrix):
	
	c=0
	z=0
	propensity=0
	x=[None,None]
	while c<no_of_species:
		if conn_matrix[reaction,c]==-1:
			x[z]=p_vec[c]
			z=z+1
		c=c+1
	propensity=x[0]*x[1]*k_vec[reaction]/(V*Na)
	return propensity

#function for calculating propensity of reaction 2A-->C
def propensity3(reaction,no_of_species,k_vec,p_vec,V,conn_matrix):
	
	c=0
	propensity=0
	while c<no_of_species:
		if conn_matrix[reaction,c]==-2:
			propensity=(p_vec[c]**2)*k_vec[reaction]/(V*Na)
			
		c=c+1
	return propensity
#function for generating type vector from connectivity matrix
def genTypeVec(conn_matrix):
	no_of_reactions=len(conn_matrix[:,0])
	no_of_species=len(conn_matrix[0,:])
	i=0
	#create empty typeVector
	typeVector=np.array([])
	while i<no_of_reactions:
		sumReactants=0
		minConnectivity=0
		j=0
		while j < no_of_species:
			if conn_matrix[i,j]<0:
				#found reactant
				sumReactants+=conn_matrix[i,j]
				if conn_matrix[i,j]< minConnectivity:
					minConnectivity=conn_matrix[i,j]
			j+=1
		if sumReactants==-1:
			#reaction A-->B, type 1
			typeVector=np.append(typeVector,1)
		if sumReactants==-2 and minConnectivity==-1:
			#reaction A+B-->C, type 2
			typeVector=np.append(typeVector,2)
		if sumReactants==-2 and minConnectivity==-2:
			#reaction A+A-->B, type 3
			typeVector=np.append(typeVector,3)

		i+=1
	return typeVector
def runKMC(p_vec_in,k_vec,conn_matrix,V,repeats,t_interval,tmax):
	no_of_reactions=len(conn_matrix[:,0])
	no_of_species=len(conn_matrix[0,:])
	#create empty vectors for storing sums of all simulation runs
	type_vector=genTypeVec(conn_matrix)
	#create evenly spaced time vector based on max time and time interval
	outputtVector=np.array([])
	z=0
	while z*t_interval<=tmax:
			outputtVector=np.append(outputtVector,z*t_interval)
			z+=1
	del z
	
	#Main program starts here
	#make requested number of repeats
	M=0
	while M < repeats:
		output_counter=0
		t=0.0
		OutputtVector=np.array([]) #create empty vector for storing time
		OutputPVector=np.array([]) #create empty vector for storing population
		#create empty vector for storing propensities
		a_vec=[]
		z=0
		while z< no_of_reactions:
			a_vec.append(0)
			z+=1
		p_vec=p_vec_in
		
		N=0
		#add population vector at t=0 to output time vector and population matrix	
		tVector=np.append(OutputtVector,t) #append t=0 for step 0
		PVector=p_vec #assign vector for storing populations of the run, first element: population at  t=0 for step 0
		output_counter=output_counter+1
		#loop over time
		while True:
			i=0
			k=0	
			#calculate propensity vector
			while i < no_of_reactions: #loop through all reactions
				if type_vector[i]==1:
					
					a_vec[i]=propensity1(i,no_of_species,k_vec,p_vec,conn_matrix)
			
				if type_vector[i]==2:
					
					a_vec[i]=propensity2(i,no_of_species,k_vec,p_vec,V,conn_matrix)
			
				if type_vector[i]==3:
					
					a_vec[i]=propensity3(i,no_of_species,k_vec,p_vec,V,conn_matrix)
			
				i=i+1
			r1=random.random() #generate first random number
			while r1==0: #r1 cannot be 0
				r1=random.random()
	
			r2=random.random() #generate second random number
			while r2==0: #r2 cannot be 0
				r2=random.random()
			q=r2*sum(a_vec)
			if sum(a_vec)==0: #total propensity cannot be 0	
				print("a_vec=0, breaking")
				if t<=tmax:
					print("premature stop of simulation, all reactions finished")
					print("time",t,)
					#continue writing last population vector with interval of t_interval until tmax is reached
					#find next time point corresponding to t_interval
					lastFixedt=t_interval*math.ceil(t/t_interval)
					
					t=lastFixedt+t_interval
					while t<tmax+t_interval:
						tVector=np.append(tVector,t)
						PVector=np.vstack([PVector,p_vec])
						t+=t_interval
						
				break 
			#select which reaction
			sumprop=0
			while k <len(k_vec):
				sumprop=sumprop+a_vec[k]
				k=k+1		
				if sumprop>q:
					reaction=k-1
					break
			#calculate required tau
			tau=(1.0/sum(a_vec))*math.log(1/r1)	
			t=t+tau #increase time
			#change population vector
			p_vec=p_vec+conn_matrix[reaction,:]
			
			N=N+1
			if t>tmax:
				tVector=np.append(tVector,t)
				PVector=np.vstack([PVector,p_vec])
				break
			
			tVector=np.append(tVector,t)
			PVector=np.vstack([PVector,p_vec])
		
		
		#create evenly spaced time series from KMC run data
		#generate time vector for final output vector
		
		#outputLineCounter=1
		outputPVector=np.zeros(shape=(len(outputtVector),no_of_species))
		x=0
		y=0
		while x<len(outputtVector):
			while y<len(tVector):
				if tVector[y]>outputtVector[x]:
					z=0
					while z<no_of_species:
						outputPVector[x][z]=PVector[y-1][z]
						z+=1
					break
				if tVector[y]==outputtVector[x]:
					z=0
					
					while z<no_of_species:
						outputPVector[x][z]=PVector[y][z]
						z+=1
					break
				y+=1
			x+=1
		if M==0:
			FinalPVector=np.zeros(shape=(len(outputPVector),no_of_species))
		#add current simulation run to previous pop vector
		if len(FinalPVector)==len(outputPVector):FinalPVector=FinalPVector+outputPVector #add current run population vector to the vector of all runs
		else:
			print("error, simulation run and output pop vector lengths do not match")
			print("simulation pop vec len",len(outputPVector),"summed pop vec len",len(FinalPVector))		
		M=M+1

	#average output over repeat runs
	FinalPVector=FinalPVector/repeats
	return outputtVector,FinalPVector
	
if __name__=="__main__":
	infile=sys.argv[1]
	with open(infile,"r") as inputFile:
		inputStream=[]
		line=inputFile.readline()
		while line:
			line=line.rstrip()
			inputStream.append(line)
			line=inputFile.readline()
	k=0
	for line in inputStream:
		lineList=line.split(" ")
		if lineList[0] == '//KMCparams':
			repeats=int(lineList[2])
			t_interval=float(lineList[3])
			tmax=float(lineList[4])
			V=float(lineList[5])
		if lineList[0] == '//popVector':
			i=1
			p_vec_in=np.array([])
			while i<len(lineList):
				p_vec_in=np.append(p_vec_in,float(lineList[i]))
				i+=1
		if lineList[0] == '//rateVector':
			i=1
			k_vec=np.array([])
			while i<len(lineList):
				k_vec=np.append(k_vec,float(lineList[i]))
				i+=1
		if lineList[0] == '//nameVector':
			i=1
			name_vec=np.array([])
			while i<len(lineList):
				name_vec=np.append(name_vec,lineList[i])
				i+=1
		if 'k_vec' in dir() and 'p_vec_in' in dir():
			if 'conn_matrix' in dir():
				pass
			else:
				conn_matrix=np.zeros(shape=(len(k_vec),len(p_vec_in)))		
		if lineList[0] == '//connMatrix':
			j=1
			while j<len(lineList):
				conn_matrix[k][j-1]=lineList[j]
				j+=1
			k+=1	
				
	
	Na=6.02*10**23 #Avogadros constant [1/mol]
	print("repeats",repeats) # number of requested simulation runs to average
	print("t_interval",t_interval) #interval to store populations
	print("tmax",tmax) # maximum simulation length
	print("V",V) #volume [L]
	print("popVector",p_vec_in) #p_vec_in: initial population vector [unitless]
	print("rateVector",k_vec) #k_vec: rate constants vector [1/s] for type 1; [L/mol*s] for type 2 and 3
	print("nameVector",name_vec) #name_vec: molecule names vector - for plotting and exporting to file
	print("connMatrix",conn_matrix) #conn_matrix: reaction connectivity matrix - each row is a change of population vector if this reaction occurs

	
	
	
	
	FinaltVector,FinalPVector=runKMC(p_vec_in,k_vec,conn_matrix,V,repeats,t_interval,tmax)
	#write population vector output to text file
	outfile=infile.split(".")[0]+"_population.csv"
	OUT=open(outfile, 'w')
	OUT.write("t ")
	
	for name in name_vec:
		OUT.write(str(name)+" ")
	OUT.write("\n")
	for i in range(len(FinaltVector)):
		OUT.write(str(FinaltVector[i])+" ")
		for j in range(len(FinalPVector[i])):
			OUT.write(str(FinalPVector[i,j])+" ")
		OUT.write("\n")
	
	
	plot=plt.figure()
	ax=plot.add_subplot(111)
	i=0
	numberOfSpecies=len(FinalPVector[0,:])
	
	while i< numberOfSpecies:
		ax.scatter(FinaltVector,FinalPVector[:,i],label=name_vec[i])
		i+=1
	plt.legend(loc='best')
	plt.show()	
