'''
Abaqus py file to post process the odb file to generate txt file with J integral output
Aditya Jhunjhunwala - UC Davis
4/15/2022

To run this file independent of the toolbox - 
	1) uncomment the input type 1 
	2) run following from cmd - abaqus cae noGUI="postProcess.py"
	or) run following from cmd - abaqus cae script="postProcess.py"

'''

session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

from abaqus import *
from abaqusConstants import *
import regionToolset
from caeModules import *
# import odbAccess
import math 
import numpy as np
import re
import json
import os

def data_to_array(data):
	# function to convert abaqus data array to numpy array
	# it assumes that the abaqus array is n x 2
	n = len(data)
	arr = np.zeros((n,2))
	for i in range(n): 
		arr[i,0] = data[i][0]
		arr[i,1] = data[i][1]
	return arr

def stressAlongPath(myOdb, pth, stressType, stepNum):
	# this function extracts the stress along the specified path for all time steps.
	# myOdb - odb object (and not the name)
	# stressType - S11, S22, S33, S12 with COMPONENT; Mises, Pressure with INVARIANT
	# we will always have one step, i.e., step 0
	myStepKey = myOdb.steps.keys()[stepNum]
	myStep = myOdb.steps[myStepKey]
	numFrames = len(myStep.frames)

	for i in range(numFrames-1): 
		
		if (stressType =='Mises' or stressType =='Pressure'):
			session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
				variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, stressType), )
		else:
			session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
				variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, stressType))
		
		session.viewports['Viewport: 1'].odbDisplay.setFrame(step=stepNum, frame=i+1)
		dataName = 'stress_plot'
		xy_data = session.XYDataFromPath(name=dataName, path=pth, includeIntersections=False, 
			projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
			projectionTolerance=0, shape=UNDEFORMED, labelType=X_COORDINATE, 
			removeDuplicateXYPairs=True, includeAllElements=False)
		xy_arr = data_to_array(xy_data)
		if i == 0:
			outputData = xy_arr
		else:
			outputData = np.append(outputData, np.reshape(xy_arr[:,1],(-1,1)), axis=1)

	return outputData

def Jintegral2D(odbName, stepNum, maxStress, outputFileName):
	# max of all J integrals in 2D case
	# Usually the stepNum is 0
	if '.odb' not in odbName:
		odbName = odbName +'.odb'
	if '.txt' not in outputFileName:
		outputFileName = outputFileName + '.txt'

	myOdb = session.openOdb(name=odbName, readOnly=False)
	myStepKey = myOdb.steps.keys()[stepNum]
	hr = myOdb.steps[myStepKey].historyRegions.keys()
	J_keys = myOdb.steps[myStepKey].historyRegions[hr[1]].historyOutputs.keys()
	nJ = len(J_keys)		# no of contours that have been extracted

	# extract the maximum J values 
	for i in range(nJ):
		variableName = 'J-integral: ' + J_keys[i] + ' in ELSET  ALL ELEMENTS'
		xy_data = session.XYDataFromHistory(name='dataJ', odb=myOdb, outputVariableName=variableName, steps=(myStepKey, ))
		xy_arr = data_to_array(xy_data)
		if i == 0:
			stress = xy_arr[:,0]*maxStress
			Jint = xy_arr[:,1]
		else: 
			Jint = np.maximum(Jint,xy_arr[:,1])
	combinedData = np.vstack((stress, Jint))
	outputData = np.transpose(combinedData)
	np.savetxt(outputFileName, outputData, fmt='%.8f') 	
	myOdb.close()

def nodeListFromNodeSet(odbName, nodeSetName, outputFileName, sortByCoord):
	myOdb = session.openOdb(name=odbName, readOnly=False)
	myAssembly = myOdb.rootAssembly
	numNodes = len(myAssembly.nodeSets[nodeSetName].nodes[0])
	coord = np.zeros((numNodes,4))
	for i in range(numNodes):
		coord[i,0] = myAssembly.nodeSets[nodeSetName].nodes[0][i].label
		coord[i,1] = myAssembly.nodeSets[nodeSetName].nodes[0][i].coordinates[0]
		coord[i,2] = myAssembly.nodeSets[nodeSetName].nodes[0][i].coordinates[1]
		coord[i,3] = myAssembly.nodeSets[nodeSetName].nodes[0][i].coordinates[2]
	coord = coord[coord[:, sortByCoord].argsort()]
	# np.savetxt(outputFileName, coord, fmt='%.6f') 	
	np.savetxt(outputFileName, coord) 	
	myOdb.close()

def Jintegral3D(odbName, jOutputName, maxStress, outputFileName, nC = 0):
	jOutputName = jOutputName.upper()
	jStartKeyword = 'J at'
	myOdb = session.openOdb(name=odbName, readOnly=False)
	myStepKey = myOdb.steps.keys()[0]
	hr = myOdb.steps[myStepKey].historyRegions.keys()
	all_keys = myOdb.steps[myStepKey].historyRegions[hr[0]].historyOutputs.keys()
	J_keys = [s for s in all_keys if (jOutputName in s) & (jStartKeyword in s)]
	nOutput = len(J_keys)
	valTemp = []
	for i in range(nOutput): 
		valTemp.append(map(int, re.findall('\d+', J_keys[i])))
	intData = np.array(valTemp)
	nContour = max(intData[:,-1])
	nJ = nOutput/nContour 
	jData = []

	if nC == 0: 
		checkContour = nContour
	else:
		checkContour = nC

	for i in range(nJ):
		for ii in range(checkContour):
			key_val = J_keys[i*nContour+ii]
			xy_arr = data_to_array(myOdb.steps[myStepKey].historyRegions[hr[0]].historyOutputs[key_val].data)
			if ii == 0:
				stress = xy_arr[:,0]*maxStress
				jInt = xy_arr[:,1]
			else: 
				jInt = np.maximum(jInt,xy_arr[:,1])
		jData.append(jInt)

	jData = np.array(jData)
	outputData = np.vstack((stress,jData)).T
	np.savetxt(outputFileName, outputData, fmt='%.8f') 	
	myOdb.close()

def clickSnapShot(odbName, stressType, imgFileName=None, stepNum=0, fy=55): 
	'''
	give odbName without the .odb extension
	give imgFileName without png extension
	'''
	title = odbName
	odbName = odbName +'.odb'
	myOdb = session.openOdb(name=odbName, readOnly=False)
	myStepKey = myOdb.steps.keys()[stepNum]
	myStep = myOdb.steps[myStepKey]
	numFrames = len(myStep.frames)
	vps = session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=150, height=150)
	session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
	vps.makeCurrent()
	vps.setValues(displayedObject=myOdb)
	vps.odbDisplay.commonOptions.setValues(visibleEdges=FEATURE)
	vps.odbDisplay.contourOptions.setValues(maxAutoCompute=OFF, maxValue=fy, minAutoCompute=OFF, minValue=-fy)
	vps.odbDisplay.display.setValues(plotState=(CONTOURS_ON_UNDEF, ))
	if (stressType =='Mises' or stressType =='Pressure'):
		vps.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, stressType),)
	else:
		vps.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, stressType),)
	vps.viewportAnnotationOptions.setValues(compass=OFF, title=OFF, state=OFF)
	vps.viewportAnnotationOptions.setValues(legendDecimalPlaces=1, legendNumberFormat=FIXED)
	session.pngOptions.setValues(imageSize=(2048, 2048))
	if imgFileName is None: 
		for i in range(numFrames): 
			figureFolder = os.path.join(os.getcwd(), title)
			try:
			    os.mkdir(figureFolder)
			except OSError as error:
			    print('Note: Output figure folder already exists. Same folder will be used.')
			imgFileName = os.path.join(figureFolder, title + '-Frame'+ str(i+1) + '.png')
			vps.odbDisplay.setFrame(step=stepNum, frame=i )
			session.printToFile(fileName=imgFileName, format=PNG, canvasObjects=(vps, ))
	else:
		if '.png' not in imgFileName:
			imgFileName = imgFileName + '.png'
		session.printToFile(fileName=imgFileName, format=PNG, canvasObjects=(vps, ))
	vps.maximize()
	myOdb.close()

def stressAlongNodeSet(odbName, setName, maxStress, outputFileName, stepNum=0, frames = [-1]):
	'''
	function extracts the stress S22 along the nodeSet for all the frames in a step
	# odbName - name of odb file (with or without the extension .odb)
	# setName - name of set defined in the model before analysis
	# outputFileName - name of output file (with or without the extension .txt)
	# stepNum - step number for which the output is to be extracted (if only one load step then stepNum = 0 always)
	# frames - -1 : all frames (>0) or specify array of frame numbers for which the stress is to be extracted
	# output - outputFileName.txt - stress S22 values along the nodeset for all load factors
	# note that the stress values are stored for nodes from z = bf to z = 0 as this is how the J is stored by abaqus
	'''
	if '.odb' not in odbName:
		odbName = odbName + '.odb'
	myOdb = session.openOdb(name=odbName, readOnly=False)
	myAssembly = myOdb.rootAssembly
	myStepKey = myOdb.steps.keys()[stepNum]
	myStep = myOdb.steps[myStepKey]
	numFrames = len(myStep.frames)
	eleRegion = myAssembly.elementSets[setName]
	if -1 in frames:
		framesConsidered = np.arange(1, numFrames, 1)
	else:
		framesConsidered = frames

	# find all the node numbers along the flange width sorted in z direction - in descending - because J is stored in descending! 
	numNodes_alongCrack = len(myAssembly.nodeSets[setName].nodes[0])
	coord = np.zeros((numNodes_alongCrack,4))
	for i in range(numNodes_alongCrack):
		coord[i,0] = myAssembly.nodeSets[setName].nodes[0][i].label
		coord[i,1] = myAssembly.nodeSets[setName].nodes[0][i].coordinates[0]
		coord[i,2] = myAssembly.nodeSets[setName].nodes[0][i].coordinates[1]
		coord[i,3] = myAssembly.nodeSets[setName].nodes[0][i].coordinates[2]
	coord = coord[coord[:, 3].argsort()]

	outputData = np.zeros((numFrames-1, numNodes_alongCrack+1)) 

	# find the stress for the elementSet for all frames 
	# then extract stresses in each node in numNodes_alongCrack and take average and store it!

	for ii, frame_num in enumerate(framesConsidered): 
		outputData[ii,0] = myStep.frames[frame_num].frameValue * maxStress
		stressField = myStep.frames[frame_num].fieldOutputs['S']
		field = stressField.getSubset(region=eleRegion, position=ELEMENT_NODAL)
		fieldValues = field.values
		numNodes = len(fieldValues)
		stressData = np.zeros((numNodes,2))
		for i in range(numNodes):
			stressData[i,0] = fieldValues[i].nodeLabel
			stressData[i, 1] = fieldValues[i].data[1]
		for j, nodeLabel in enumerate(coord[::-1,0]):
			idx = np.where(stressData[:,0] == nodeLabel)
			outputData[ii,j+1] = np.average(stressData[idx, 1])

	np.savetxt(outputFileName, outputData, fmt='%.8f')
	myOdb.close()

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# %% POST PROCESS ODB FOR J INTEGRAL
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# %% INPUT FOR DIRECT RUN------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Sims = ['BNM-FM-PB']
# maxStress = [1.0]
# nC_sims = [1, 2, 3, 4, 5]
# fy = 55.0
# flag_generateImg = False

# %% INPUT FOR TOOLBOX---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
with open('geomData.json') as f:
    geomData = json.load(f)

modelName = str(geomData['modelName'])	

Sims = [modelName]
maxStress = [1.0]
nC_sims = [1, 2, 3, 4, 5]
fy = geomData['sigma_ys']
flag_generateImg = geomData['flag_generateImg']

# %% POST PROCESS--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
jOutputName = 'J-FLANGEOUTPUT'
nodeSetName = 'ALONGWIDTH'
crackNodeSetName = 'ALONGCRACK'

for i, sim in enumerate(Sims):
	odbName = sim + '.odb'
	outputFileNodeSet = sim + '-NodeCoord.txt'
	outputFileStress = sim + '-Stress.txt'
	stress = maxStress[i]
	nodeListFromNodeSet(odbName, nodeSetName, outputFileNodeSet, 3)
	stressAlongNodeSet(odbName, crackNodeSetName, stress, outputFileStress)
	for nC in nC_sims:
		outputFileJ = sim + '-J-' + str(nC) + '.txt'
		Jintegral3D(odbName, jOutputName, stress, outputFileJ, nC)
		if flag_generateImg == True:
			stressType = 'S22'
			imgFileName = None
			clickSnapShot(sim, stressType, imgFileName, stepNum=0, fy=fy)