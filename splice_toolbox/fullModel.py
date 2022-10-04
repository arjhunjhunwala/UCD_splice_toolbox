'''
Abaqus py file to generate splice full model
Aditya Jhunjhunwala - UC Davis
4/15/2022

To run this file independent of the toolbox - 
	1) uncomment the input type 1 
	2) run following from cmd - abaqus cae script="fullModel.py"
'''

from abaqus import *
from abaqusConstants import *
import regionToolset
from caeModules import *
import math
import numpy as np
import json

# %% Options required to store abaqus.rpy in readable format
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


# %% functions to create individiual parts and then create the full model

def createWSection(tf, bf, tw, bw, rtran, length, a_crack, modelName, partName, secName, neleAlongWidth=40): 
	'''
	Function creates the W section with meshing
	'''
	try: 
		myModel = mdb.models[modelName]
	except: 
		myModel = mdb.Model(name=modelName)

	# variables 
	dsec = tf*2 + rtran*2 + bw
	delta= 0.01

	#------------------------------------------------------------------------------------
	# a) create sketch of flange
	s = myModel.ConstrainedSketch(name=partName, sheetSize=20)
	s.ConstructionLine(point1=(0.0, 0.0), point2=(0.0, 1.0))
	s.ConstructionLine(point1=(0.0, 0.0), point2=(1.0, 0.0))
	s.Line((-(bw/2 + rtran + tf), 0.0), (-(bw/2 + rtran + tf), -bf/2))
	s.Line((-(bw/2 + rtran + tf), -bf/2), (-(bw/2 + rtran), -bf/2))
	s.Line((-(bw/2 + rtran), -bf/2), (-(bw/2 + rtran), -(tw/2 + rtran)))
	s.ArcByCenterEnds((-bw/2, -(tw/2 + rtran)), (-(bw/2 + rtran), -(tw/2 + rtran)), (-bw/2, -tw/2), direction=CLOCKWISE)
	s.Line((-bw/2, -tw/2), (0.0, -tw/2))
	allItems = s.geometry.items()
	geomItems = []
	for allItem in allItems:
		geomItems.append(allItem[1])
	s.copyMirror(mirrorLine=s.geometry.findAt((0.0, 1.0), ), objectList=geomItems)
	allItems = s.geometry.items()
	geomItems = []
	for allItem in allItems:
		geomItems.append(allItem[1])
	s.copyMirror(mirrorLine=s.geometry.findAt((1.0, 0.0), ), objectList=geomItems)

	myPart = myModel.Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
	myPart.BaseSolidExtrude(sketch=s, depth=length)
	c = myPart.cells

	#------------------------------------------------------------------------------------
	# Assign section
	cells = c.findAt(((0.0, 0.0, length/2), ))
	region = myPart.Set(cells=cells, name=partName)
	myPart.SectionAssignment(region=region, sectionName=secName, offset=0.0, offsetType=MIDDLE_SURFACE, 
		offsetField='', thicknessAssignment=FROM_SECTION)

	#------------------------------------------------------------------------------------
	# Assign mesh type 
	e, f, c, v = myPart.edges, myPart.faces, myPart.cells, myPart.vertices
	elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
	elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
	elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)

	myPart.setElementType(regions=(myPart.cells,), elemTypes=(elemType1, elemType2, elemType3))

	#------------------------------------------------------------------------------------
	# Create partition for the meshing
	# partition the web 
	pickedCells = c.findAt(((0.0, 0.0, length/2.0), ))
	e, v, d = myPart.edges, myPart.vertices, myPart.datums
	myPart.PartitionCellByPlanePointNormal(point=v.findAt(coordinates=(bw/2, tw/2, 0.0)), 
		normal=e.findAt(coordinates=(0.0, tw/2, 0.0)), cells=pickedCells)

	pickedCells = c.findAt(((0.0, 0.0, length/2.0), ))
	e, v, d = myPart.edges, myPart.vertices, myPart.datums
	myPart.PartitionCellByPlanePointNormal(point=v.findAt(coordinates=(-bw/2, tw/2, 0.0)), 
		normal=e.findAt(coordinates=(0.0, tw/2, 0.0)), cells=pickedCells)

	pickedCells = c.findAt(((0.0, 0.0, length/2.0), ))
	e, v, d = myPart.edges, myPart.vertices, myPart.datums
	myPart.PartitionCellByPlanePointNormal(normal=e.findAt(coordinates=(bw/2, tw/4, 0.0)), cells=pickedCells, 
		point=myPart.InterestingPoint(edge=e.findAt(coordinates=(bw/2, tw/4, 0.0)), rule=MIDDLE))

	# partititon the transition region
	pickedCells = c.findAt(((bw/2+rtran, 0.0, length/2.0), ))
	e, v, d = myPart.edges, myPart.vertices, myPart.datums
	myPart.PartitionCellByPlanePointNormal(point=v.findAt(coordinates=(bw/2+rtran, tw/2+rtran, 0.0)), 
		normal=e.findAt(coordinates=(0.0, tw/2, 0.0)), cells=pickedCells)

	pickedCells = c.findAt(((-bw/2-rtran, 0.0, length/2.0), ))
	e, v, d = myPart.edges, myPart.vertices, myPart.datums
	myPart.PartitionCellByPlanePointNormal(point=v.findAt(coordinates=(-bw/2-rtran, tw/2+rtran, 0.0)), 
		normal=e.findAt(coordinates=(0.0, tw/2, 0.0)), cells=pickedCells)

	# partition the compression flange for crack
	datum = myPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=(-bw/2-rtran-a_crack))
	datum_id = datum.id
	pickedCells = c.findAt(((-bw/2-rtran-tf/2, 0.0, length/2.0), ))
	myPart.PartitionCellByDatumPlane(datumPlane=myPart.datums[datum_id], cells=pickedCells)

	#------------------------------------------------------------------------------------
	# Meshing
	e, f, c, v = myPart.edges, myPart.faces, myPart.cells, myPart.vertices

	# global mesh size = thickness of flange 
	ele_size_global = tf
	myPart.seedPart(size=ele_size_global, deviationFactor=0.1, minSizeFactor=0.1)

	# Mesh the tension flange 
	neleTransition = int(math.ceil((tw+2*rtran)/bf*neleAlongWidth/2)*2)
	neleOther = neleAlongWidth - neleTransition
	ele_size_alongwidth = bf/neleAlongWidth

	neleThk = 2

	pickedEdges = e.findAt(((dsec/2.0, bf/4.0, 0.0), ), ((dsec/2.0, -bf/4.0, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleAlongWidth/2, constraint=FINER)	

	pickedEdges = e.findAt(((dsec/2.0-tf, 0.0, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleTransition, constraint=FINER)	

	pickedEdges = e.findAt(((dsec/2.0-tf, bf/2.0 - 0.1, 0.0), ), ((dsec/2.0-tf, -bf/2.0 + 0.1, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleOther/2, constraint=FINER)	

	pickedEdges = e.findAt(((dsec/2.0-tf/2.0, bf/2.0, 0.0), ), ((dsec/2.0-tf/2.0, -bf/2.0, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleThk, constraint=FINER)	

	# Mesh the tension transition region 
	ele_size_TranCurve = ele_size_alongwidth
	x_coord = bw/2.0 + rtran*math.cos(math.pi/4.0)
	y_coord = tw/2.0 + rtran - rtran*math.sin(math.pi/4.0)
	pickedEdges = e.findAt(((x_coord, y_coord, 0.0), ), ((x_coord, -y_coord, 0.0), ))
	myPart.seedEdgeBySize(edges=pickedEdges, size=ele_size_TranCurve, deviationFactor=0.1, constraint=FINER)

	# Mesh the web 
	if bw > 11: 
		neleWeb = 16
	else: 
		neleWeb = 8
	pickedEdges = e.findAt(((bw/4, tw/2, 0.0), ), ((bw/4, -tw/2, 0.0), ), ((-bw/4, tw/2, 0.0), ), ((-bw/4, -tw/2, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWeb/2, constraint=FINER)
	pickedEdges = e.findAt(((bw/4, tw/2, length), ), ((bw/4, -tw/2, length), ), ((-bw/4, tw/2, length), ), ((-bw/4, -tw/2, length), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWeb/2, constraint=FINER)
	pickedEdges = e.findAt(((0.0, 0.0, 0.0), ), ((0.0, 0.0, length), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWeb, constraint=FINER)

	neleWebThk = 2
	pickedEdges = e.findAt(((bw/2, tw/4, 0.0), ), ((bw/2, -tw/4, 0.0), ), ((-bw/2, tw/4, 0.0), ), ((-bw/2, -tw/4, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWebThk/2, constraint=FINER)

	# Mesh the compression transition region 
	pickedEdges = e.findAt(((-bw/2-rtran, 0.0, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWebThk, constraint=FINER)

	# Mesh the compression flange
	neleOther = 4
	pickedEdges = e.findAt(((-bw/2-rtran, bf/2-0.1, 0.0), ), ((-bw/2-rtran, -bf/2+0.1, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleOther, constraint=FINER)

	pickedEdges = e.findAt(((-dsec/2, bf/2-0.1, 0.0), ), ((-dsec/2, -bf/2+0.1, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleOther + neleWebThk/2, constraint=FINER)

	pickedEdges = e.findAt(((-dsec/2+tf-a_crack, 0.0, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=2*neleOther + neleWebThk, constraint=FINER)

	# Assign mesh control 
	coordMin = (-(dsec/2 + delta), -(bf/2 + delta), -delta)
	coordMax = ((dsec/2 + delta), (bf/2 + delta), length + delta)
	boundingCoords = coordMin + coordMax
	allCells = myPart.cells.getByBoundingBox(*boundingCoords)
	myPart.setMeshControls(regions=allCells, technique=SWEEP, algorithm=ADVANCING_FRONT)
	myPart.assignStackDirection(referenceRegion=f.findAt(coordinates=(delta, delta, 0.0)), cells=allCells)

	# Assign sweep path for each element - head ache 
	# from left flange to right flange 
	myPart.setSweepPath(region=c.findAt(coordinates=(-dsec/2+delta, 0.0, length/2)), edge=e.findAt(coordinates=(-dsec/2, bf/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(-bw/2-rtran-a_crack/2, 0.0, length/2)), edge=e.findAt(coordinates=(-dsec/2+tf, bf/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(-bw/2-rtran/2, 0.0, length/2)), edge=e.findAt(coordinates=(-bw/2, tw/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(0.0, delta, length/2)), edge=e.findAt(coordinates=(-bw/2, tw/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(0.0, -delta, length/2)), edge=e.findAt(coordinates=(-bw/2, -tw/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(bw/2+rtran/2, 0.0, length/2)), edge=e.findAt(coordinates=(bw/2, tw/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(dsec/2-tf/2, 0.0, length/2)), edge=e.findAt(coordinates=(dsec/2, bf/2, length/2)), sense=FORWARD)

	myPart.generateMesh()


def createWSectionCut(tf, bf, tw, bw, rtran, length, a_crack, modelName, partName, secName, ele_size_bottom, ele_size_top, flag): 
	'''
	Function creates W section with tension flange and tension transition zone removed
	flag : 'Top' or 'Bottom' - required for mesh transition direction 
	'''
	try: 
		myModel = mdb.models[modelName]
	except: 
		myModel = mdb.Model(name=modelName)

	# variables 
	dsec = tf*2 + rtran*2 + bw
	delta= 0.01

	#------------------------------------------------------------------------------------
	# a) create sketch of flange
	s = myModel.ConstrainedSketch(name=partName, sheetSize=20)
	s.ConstructionLine(point1=(0.0, 0.0), point2=(0.0, 1.0))
	s.ConstructionLine(point1=(0.0, 0.0), point2=(1.0, 0.0))
	s.Line((-(bw/2 + rtran + tf), 0.0), (-(bw/2 + rtran + tf), -bf/2))
	s.Line((-(bw/2 + rtran + tf), -bf/2), (-(bw/2 + rtran), -bf/2))
	s.Line((-(bw/2 + rtran), -bf/2), (-(bw/2 + rtran), -(tw/2 + rtran)))
	s.ArcByCenterEnds((-bw/2, -(tw/2 + rtran)), (-(bw/2 + rtran), -(tw/2 + rtran)), (-bw/2, -tw/2), direction=CLOCKWISE)
	s.Line((-bw/2, -tw/2), (bw/2, -tw/2))
	s.Line((bw/2, -tw/2), (bw/2, 0.0))
	allItems = s.geometry.items()
	geomItems = []
	for allItem in allItems:
		geomItems.append(allItem[1])
	s.copyMirror(mirrorLine=s.geometry.findAt((1.0, 0.0), ), objectList=geomItems)

	myPart = myModel.Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
	myPart.BaseSolidExtrude(sketch=s, depth=length)
	c = myPart.cells

	#------------------------------------------------------------------------------------
	# Assign section
	cells = c.findAt(((0.0, 0.0, length/2), ))
	region = myPart.Set(cells=cells, name=partName)
	myPart.SectionAssignment(region=region, sectionName=secName, offset=0.0, offsetType=MIDDLE_SURFACE, 
		offsetField='', thicknessAssignment=FROM_SECTION)

	#------------------------------------------------------------------------------------
	# Assign mesh type 
	e, f, c, v = myPart.edges, myPart.faces, myPart.cells, myPart.vertices
	elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
	elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
	elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)

	myPart.setElementType(regions=(myPart.cells,), elemTypes=(elemType1, elemType2, elemType3))

	#------------------------------------------------------------------------------------
	# Create partition for the meshing
	# partition the web 
	pickedCells = c.findAt(((0.0, 0.0, length/2.0), ))
	e, v, d = myPart.edges, myPart.vertices, myPart.datums
	myPart.PartitionCellByPlanePointNormal(point=v.findAt(coordinates=(-bw/2, tw/2, 0.0)), 
		normal=e.findAt(coordinates=(0.0, tw/2, 0.0)), cells=pickedCells)

	pickedCells = c.findAt(((0.0, 0.0, length/2.0), ))
	e, v, d = myPart.edges, myPart.vertices, myPart.datums
	myPart.PartitionCellByPlanePointNormal(normal=e.findAt(coordinates=(bw/2, tw/4, 0.0)), cells=pickedCells, 
		point=myPart.InterestingPoint(edge=e.findAt(coordinates=(-bw/2, tw/4, 0.0)), rule=MIDDLE))

	# partititon the transition region
	pickedCells = c.findAt(((-bw/2-rtran, 0.0, length/2.0), ))
	e, v, d = myPart.edges, myPart.vertices, myPart.datums
	myPart.PartitionCellByPlanePointNormal(point=v.findAt(coordinates=(-bw/2-rtran, tw/2+rtran, 0.0)), 
		normal=e.findAt(coordinates=(0.0, tw/2, 0.0)), cells=pickedCells)

	# partition the compression flange for crack
	datum = myPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=(-bw/2-rtran-a_crack))
	datum_id = datum.id
	pickedCells = c.findAt(((-bw/2-rtran-tf/2, 0.0, length/2.0), ))
	myPart.PartitionCellByDatumPlane(datumPlane=myPart.datums[datum_id], cells=pickedCells)

	#------------------------------------------------------------------------------------
	# Meshing
	e, f, c, v = myPart.edges, myPart.faces, myPart.cells, myPart.vertices

	# global mesh size = thickness of flange 
	ele_size_global = tf
	myPart.seedPart(size=ele_size_global, deviationFactor=0.1, minSizeFactor=0.1)

	# Mesh the web 
	if bw > 11: 
		neleWeb = 16
	else: 
		neleWeb = 8
	pickedEdges = e.findAt(((0.0, tw/2, 0.0), ), ((0.0, -tw/2, 0.0), ), ((0.0, tw/2, length), ), ((0.0, -tw/2, length), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWeb, constraint=FINER)
	pickedEdges = e.findAt(((0.0, 0.0, 0.0), ), ((0.0, 0.0, length), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWeb, constraint=FINER)

	neleWebThk = 2
	pickedEdges = e.findAt(((bw/2, tw/4, 0.0), ), ((bw/2, -tw/4, 0.0), ), ((-bw/2, tw/4, 0.0), ), ((-bw/2, -tw/4, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWebThk/2, constraint=FINER)

	# Mesh the compression transition region 
	pickedEdges = e.findAt(((-bw/2-rtran, 0.0, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleWebThk, constraint=FINER)

	# Mesh the compression flange
	neleOther = 4
	pickedEdges = e.findAt(((-bw/2-rtran, bf/2-0.1, 0.0), ), ((-bw/2-rtran, -bf/2+0.1, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleOther, constraint=FINER)

	pickedEdges = e.findAt(((-dsec/2, bf/2-0.1, 0.0), ), ((-dsec/2, -bf/2+0.1, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=neleOther + neleWebThk/2, constraint=FINER)

	pickedEdges = e.findAt(((-dsec/2+tf-a_crack, 0.0, 0.0), ))
	myPart.seedEdgeByNumber(edges=pickedEdges, number=2*neleOther + neleWebThk, constraint=FINER)

	# Assign mesh control 
	coordMin = (-(dsec/2 + delta), -(bf/2 + delta), -delta)
	coordMax = ((dsec/2 + delta), (bf/2 + delta), length + delta)
	boundingCoords = coordMin + coordMax
	allCells = myPart.cells.getByBoundingBox(*boundingCoords)
	myPart.setMeshControls(regions=allCells, technique=SWEEP, algorithm=ADVANCING_FRONT)
	myPart.assignStackDirection(referenceRegion=f.findAt(coordinates=(delta, delta, 0.0)), cells=allCells)

	# Assign sweep path for each element
	# from left flange to right flange 
	myPart.setSweepPath(region=c.findAt(coordinates=(-dsec/2+delta, 0.0, length/2)), edge=e.findAt(coordinates=(-dsec/2, bf/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(-bw/2-rtran-a_crack/2, 0.0, length/2)), edge=e.findAt(coordinates=(-dsec/2+tf, bf/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(-bw/2-rtran/2, 0.0, length/2)), edge=e.findAt(coordinates=(-bw/2, tw/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(0.0, delta, length/2)), edge=e.findAt(coordinates=(-bw/2, tw/2, length/2)), sense=FORWARD)
	myPart.setSweepPath(region=c.findAt(coordinates=(0.0, -delta, length/2)), edge=e.findAt(coordinates=(-bw/2, -tw/2, length/2)), sense=FORWARD)
	# myPart.setSweepPath(region=c.findAt(coordinates=(bw/2+rtran/2, 0.0, length/2)), edge=e.findAt(coordinates=(bw/2, tw/2, length/2)), sense=FORWARD)
	# myPart.setSweepPath(region=c.findAt(coordinates=(dsec/2-tf/2, 0.0, length/2)), edge=e.findAt(coordinates=(dsec/2, bf/2, length/2)), sense=FORWARD)

	# Assign biased meshing for the vertical height 
	pickedEdges = e.findAt(((-dsec/2, bf/2, length/2), ), ((-(bw/2 + rtran), bf/2, length/2), ), 
		((-(bw/2 + rtran), tw/2 + rtran, length/2), ), ((bw/2, tw/2, length/2), ), ((bw/2, -tw/2, length/2), ),)
	if flag == 'Top': 
		myPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=pickedEdges, minSize=ele_size_bottom, maxSize=ele_size_top, constraint=FINER)
	else: 
		myPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges, minSize=ele_size_bottom, maxSize=ele_size_top, constraint=FINER)

	myPart.generateMesh()


def createCrackedFlange(t_u, t_l, bf, tw_u, tw_l, bw, rtran, h_plate, a_crack, aw_crack, modelName, steelSectionName, weldSectionName, ele_size_bottom, ele_size_top, neleAlongWidth=40):
	'''
	Function creates the flange with crack and the transition regions
	'''
	try: 
		myModel = mdb.models[modelName]
	except: 
		myModel = mdb.Model(name=modelName)

	# ---------------------------------------------------------------------------
	# Inputs processing and constant parameters 
	# ---------------------------------------------------------------------------

	# transition details 
	r_transition = rtran                     # radius of transition for W14 & W24 typically = 1.25
	# crack details
	r_crack = 0.001 / 2                     # radius of crack
	a_del = r_crack * 0                     # crack increment (used in J verification)
	a_crack = a_crack + a_del

	# constants
	delta = 1e-5                            # will be useful when selecting faces
	crackWidthRatio = a_crack / t_u
	xi = t_u / t_l
	print("Model with eta : ", crackWidthRatio," and xi : ", xi)


	# ---------------------------------------------------------------------------
	# Initialization
	# ---------------------------------------------------------------------------

	# mesh input details
	nQuart = 10
	nRad = 15

	eleFac = 2.0
	ele_size_crack_min = ele_size_bottom/eleFac         # reduction from ele_size_bottom towards crack tip
	ele_size_alongwidth = bf/neleAlongWidth             # this is important - refer meshing technique

	# crack 1st contour radius details
	dmin = min(0.15, ele_size_bottom)                   # half of min distance between the inner edge and the circle
	r2max = 0.25                                        # max of the inner circle
	r2 = r2max
	if (r2 >= a_crack - 2*dmin) or r2 >= (t_u - a_crack - 2*dmin) : 
	    r2 = min(a_crack, t_u - a_crack)*0.6


	# ---------------------------------------------------------------------------
	# Create parts and assign section
	# ---------------------------------------------------------------------------

	# %% ------------------------ Flange portion ---------------------------------

	flangePartName = 'flange'

	# a) create sketch of flange
	s = myModel.ConstrainedSketch(name='flangeSketch', sheetSize=20)
	s.Line((0.0, r_crack), (a_crack, r_crack))
	s.Line((0.0, -r_crack), (a_crack, -r_crack))
	s.ArcByCenterEnds((a_crack, 0.0), (a_crack, r_crack), (a_crack, -r_crack), direction=CLOCKWISE)
	s.Line((0.0, r_crack), (0.0, h_plate))
	s.Line((0.0, h_plate), (t_u, h_plate))
	s.Line((t_u, h_plate), (t_u, t_u - a_crack))
	s.Line((t_u, t_u - a_crack), (t_u, 0.0))
	if t_u != t_l:
	    s.Line((t_u, 0.0), (t_l, 0.0))
	s.Line((t_l, 0.0), (t_l, -h_plate))
	s.Line((t_l, -h_plate), (0.0, -h_plate))
	s.Line((0.0, -h_plate), (0.0, -r_crack))

	# b) extrude it to part/ create a part from sketch
	flangePart = myModel.Part(name=flangePartName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
	flangePart.BaseSolidExtrude(sketch=s, depth=bf)
	e, f, c, v = flangePart.edges, flangePart.faces, flangePart.cells, flangePart.vertices

	# c1.1) create edges for circle around crack
	f1 = f.findAt((t_u / 2, t_u, bf), )
	e1 = e.findAt((t_u, t_u, bf), )
	t = flangePart.MakeSketchTransform(sketchPlane=f1, sketchPlaneSide=SIDE1, origin=(0.0, 0.0, bf), sketchUpEdge=e1)
	s1 = myModel.ConstrainedSketch(name='flangeCircle', sheetSize=20, transform=t)
	flangePart.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
	s1.CircleByCenterPerimeter((a_crack, 0.0), (a_crack + r2, 0.0))
	flangePart.PartitionFaceBySketch(faces=f1, sketch=s1, sketchUpEdge=e1)

	# c1.2) extrude circle to partition the cells
	pc1 = c.findAt((t_u / 2, t_u, bf / 2), )
	pe1 = (e.findAt((a_crack, r2, bf)), e.findAt((a_crack, -r2, bf)))
	flangePart.PartitionCellByExtrudeEdge(line=e.findAt(coordinates=(t_u, h_plate, bf / 2)),
	                                      cells=pc1, edges=pe1, sense=REVERSE)

	# c2.1) create edges in the circle around the crack
	f2 = f.findAt((a_crack + r2 / 2, r2 / 2, bf), )
	e2 = e.findAt((t_u, t_u, bf), )
	t2 = flangePart.MakeSketchTransform(sketchPlane=f2, sketchPlaneSide=SIDE1,
	                                    origin=(0.0, 0.0, bf), sketchUpEdge=e2)
	s2 = myModel.ConstrainedSketch(name='flangeCircleLines', sheetSize=20, transform=t2)
	flangePart.projectReferencesOntoSketch(sketch=s2, filter=COPLANAR_EDGES)
	s2.Line((a_crack, -r2), (a_crack, r2))
	s2.Line((a_crack, 0.0), (t_u, t_u - a_crack))
	flangePart.PartitionFaceBySketch(faces=f2, sketch=s2, sketchUpEdge=e2)

	# c2.2) extrude edges to partition the cells
	pc2 = c.findAt((a_crack + r2 / 2, delta, bf / 2), )
	pe2 = (e.findAt((a_crack, r2 / 2, bf)), e.findAt((a_crack, -r2 / 2, bf)), e.findAt((a_crack + r2 / 2, r2 / 2, bf)))
	flangePart.PartitionCellByExtrudeEdge(line=e.findAt(coordinates=(t_u, h_plate, bf / 2)),
	                                      cells=pc2, edges=pe2, sense=REVERSE)

	# c3.1) create edges outside the circle for weld part
	f3 = f.findAt((t_u / 2, t_u, bf), )
	e3 = e.findAt((t_u, t_u, bf), )
	t3 = flangePart.MakeSketchTransform(sketchPlane=f3, sketchPlaneSide=SIDE1,
	                                    origin=(0.0, 0.0, bf), sketchUpEdge=e3)
	s3 = myModel.ConstrainedSketch(name='flangeOuterLines', sheetSize=20, transform=t3)
	flangePart.projectReferencesOntoSketch(sketch=s3, filter=COPLANAR_EDGES)
	s3.Line((a_crack, 0.0), (t_u, t_u - a_crack))
	flangePart.PartitionFaceBySketch(faces=f3, sketch=s3, sketchUpEdge=e3)

	# c3.2) extrude edges to partition the cells
	pc3 = c.findAt((t_u / 2, t_u, bf / 2), )
	pe3 = (e.findAt((t_u - delta, t_u - a_crack - delta, bf)),)
	flangePart.PartitionCellByExtrudeEdge(line=e.findAt(coordinates=(t_u, h_plate, bf / 2)),
	                                      cells=pc3, edges=pe3, sense=REVERSE)

	# c4) create horizontal datum partition 
	pc4 = c.findAt(((t_l / 2, -t_l, bf / 2),), ((a_crack + r2 / 2, 0.0, bf / 2),))
	flangePart.PartitionCellByPlanePointNormal(point=v.findAt(coordinates=(t_l, 0.0, bf)),
	                                           normal=e.findAt(coordinates=(t_u, t_u, bf)), cells=pc4)

	# Assign section
	cells = c.findAt(((t_u / 2, t_u, bf / 2),), ((t_l / 2, -t_l, bf / 2),), ((a_crack - delta, r2 / 2, bf / 2),),
	                 ((a_crack - delta, -r2 / 2, bf / 2),), ((a_crack + delta, r2 / 2, bf / 2),),
	                 ((a_crack + delta, -r2 / 2, bf / 2),))
	region = flangePart.Set(cells=cells, name=flangePartName + 'Steel')
	flangePart.SectionAssignment(region=region, sectionName=steelSectionName, offset=0.0,
	                             offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

	cells = c.findAt(((a_crack + r2 / 2, delta, bf / 2),), ((a_crack + r2 + delta, delta, bf / 2),))
	region = flangePart.Set(cells=cells, name=flangePartName + 'Weld')
	flangePart.SectionAssignment(region=region, sectionName=weldSectionName, offset=0.0,
	                             offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

	# %% ------------------------ Transition portion --------------------------

	# top transition
	topTransitionPartName = 'TopTransition'
	len1 = tw_u + 2 * r_transition
	len2 = tw_u
	thk = r_transition
	ht = h_plate - r_crack
	s = myModel.ConstrainedSketch(name='TopTransitionSketch', sheetSize=20)
	s.Line((-len1 / 2, 0.0), (len1 / 2, 0.0))
	s.Line((-len2 / 2, thk), (len2 / 2, thk))
	s.ArcByCenterEnds((-len1 / 2, thk), (-len2 / 2, thk), (-len1 / 2, 0.0), direction=CLOCKWISE)
	s.ArcByCenterEnds((len1 / 2, thk), (len2 / 2, thk), (len1 / 2, 0.0), direction=COUNTERCLOCKWISE)
	topTransitionPart = myModel.Part(name=topTransitionPartName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
	topTransitionPart.BaseSolidExtrude(sketch=s, depth=ht)
	# e, f, c, v = topTransitionPart.edges, topTransitionPart.faces, topTransitionPart.cells, topTransitionPart.vertices
	region = (topTransitionPart.cells,)
	topTransitionPart.SectionAssignment(region=region, sectionName=steelSectionName, offset=0.0,
	                                    offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

	# bottom transition
	botTransitionPartName = 'BotTransition'
	len1 = tw_l + 2 * r_transition
	len2 = tw_l
	thk = r_transition
	ht = h_plate - r_crack
	s = myModel.ConstrainedSketch(name='BotTransitionSketch', sheetSize=20)
	s.Line((-len1 / 2, 0.0), (len1 / 2, 0.0))
	s.Line((-len2 / 2, thk), (len2 / 2, thk))
	s.ArcByCenterEnds((-len1 / 2, thk), (-len2 / 2, thk), (-len1 / 2, 0.0), direction=CLOCKWISE)
	s.ArcByCenterEnds((len1 / 2, thk), (len2 / 2, thk), (len1 / 2, 0.0), direction=COUNTERCLOCKWISE)
	botTransitionPart = myModel.Part(name=botTransitionPartName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
	botTransitionPart.BaseSolidExtrude(sketch=s, depth=ht)
	# e, f, c, v = botTransitionPart.edges, botTransitionPart.faces, botTransitionPart.cells, botTransitionPart.vertices
	region = (botTransitionPart.cells,)
	botTransitionPart.SectionAssignment(region=region, sectionName=steelSectionName, offset=0.0,
	                                    offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

	# %% ------------------------ Web portion ---------------------------------
	# note : this web part is not used in the full model 

	webPartName = 'Web'
	s = myModel.ConstrainedSketch(name='WebSketch', sheetSize=20)
	s.Line((-tw_u / 2, r_crack), (-tw_u / 2 + aw_crack, r_crack))
	s.Line((-tw_l / 2, -r_crack), (-tw_u / 2 + aw_crack, -r_crack))
	s.ArcByCenterEnds((-tw_u / 2 + aw_crack, 0.0), (-tw_u / 2 + aw_crack, r_crack), (-tw_u / 2 + aw_crack, -r_crack),
	                  direction=CLOCKWISE)
	s.Line((-tw_u / 2, r_crack), (-tw_u / 2, h_plate))
	s.Line((-tw_u / 2, h_plate), (tw_u / 2, h_plate))
	s.Line((tw_u / 2, h_plate), (tw_u / 2, tw_u - aw_crack))
	s.Line((tw_u / 2, tw_u - aw_crack), (tw_l / 2, 0.0))
	s.Line((tw_l / 2, 0.0), (tw_l / 2, -h_plate))
	s.Line((tw_l / 2, -h_plate), (-tw_l / 2, -h_plate))
	s.Line((-tw_l / 2, -h_plate), (-tw_l / 2, -r_crack))

	webPart = myModel.Part(name=webPartName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
	webPart.BaseSolidExtrude(sketch=s, depth=bw / 2)
	e, f, c, v = webPart.edges, webPart.faces, webPart.cells, webPart.vertices

	# create horizontal partition at the two webs interface
	pc = c.findAt((0.0, tw_u, bw / 4), )
	webPart.PartitionCellByPlanePointNormal(point=v.findAt((tw_l / 2, 0.0, 0.0)),
	                                        normal=e.findAt((tw_u / 2, tw_u, 0.0)), cells=pc)

	# create 45 degree partition 
	pc = c.findAt(((0.0, tw_u, bw / 4),))
	ip = webPart.InterestingPoint(edge=e.findAt(
	    (-tw_u / 2 + aw_crack + r_crack * math.cos(math.pi / 4), r_crack * math.sin(math.pi / 4), bw / 2)), rule=CENTER)
	webPart.PartitionCellByPlaneThreePoints(point1=v.findAt((tw_u / 2, tw_u - aw_crack, 0.0)),
	                                        point2=v.findAt((tw_u / 2, tw_u - aw_crack, bw / 2)), point3=ip, cells=pc)

	pc = c.findAt(((0.0, tw_u, bf / 4),), ((0.0, -tw_l, bf / 4),))
	region = webPart.Set(cells=pc, name=webPartName + 'Steel')
	webPart.SectionAssignment(region=region, sectionName=steelSectionName, offset=0.0,
	                          offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

	cells = c.findAt(((tw_u / 2 - delta, delta, bf / 4),))
	region = webPart.Set(cells=cells, name=webPartName + 'Weld')
	webPart.SectionAssignment(region=region, sectionName=weldSectionName, offset=0.0,
	                          offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

	# ---------------------------------------------------------------------------
	# Apply MESHING!!!
	# ---------------------------------------------------------------------------

	# %% assigning the element type
	elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
	elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
	elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)

	flangePart.setElementType(regions=(flangePart.cells,), elemTypes=(elemType1, elemType2, elemType3))
	topTransitionPart.setElementType(regions=(topTransitionPart.cells,), elemTypes=(elemType1, elemType2, elemType3))
	botTransitionPart.setElementType(regions=(botTransitionPart.cells,), elemTypes=(elemType1, elemType2, elemType3))
	webPart.setElementType(regions=(webPart.cells,), elemTypes=(elemType1, elemType2, elemType3))

	# %% meshing the flange!

	# due to inherent rules of abaqus, the along width seeds are provided everywhere
	# then the other seeds are assigned

	flangePart.seedPart(size=ele_size_alongwidth, deviationFactor=0.1, minSizeFactor=0.1)

	e, f, c, v = flangePart.edges, flangePart.faces, flangePart.cells, flangePart.vertices

	## assign the element type for the flange
	# for the top and bottom bulky portion
	pickedCells = c.findAt(((t_u / 2, 2*t_u, bf / 2),), ((t_l / 2, -2*t_l, bf / 2),))
	flangePart.assignStackDirection(referenceRegion=f.findAt(coordinates=(t_u / 2, 2*t_u, bf)), cells=pickedCells)
	flangePart.setMeshControls(regions=pickedCells, algorithm=ADVANCING_FRONT)

	# assign stack direction for the crack front
	coordMin = (-delta, -(t_l + delta), -delta )
	coordMax = (t_l + delta, t_u + delta, bf + delta )
	boundingCoords = coordMin + coordMax
	crackFrontCells = flangePart.cells.getByBoundingBox(*boundingCoords)
	flangePart.assignStackDirection(referenceRegion=f.findAt((a_crack + r2 / 2, delta, bf)), cells=crackFrontCells)

	# assign meshing type to the weld portion
	pickedCells = c.findAt(((t_u - delta, delta, bf / 2),), )
	flangePart.setMeshControls(regions=pickedCells, technique= SWEEP, algorithm=ADVANCING_FRONT)

	## local seeds
	c1 = (a_crack, 0.0, bf - delta)
	c2 = (a_crack, 0.0, bf + delta)
	crackFrontEdges = e.getByBoundingCylinder(c1, c2, r2 + delta)
	flangePart.seedEdgeByNumber(edges=crackFrontEdges, number=nQuart, constraint=FINER)

	# corrections to the half quarter edges
	coord1 = (a_crack - delta, -delta, bf - delta)
	coord2 = (a_crack + r2 + delta, r2 + delta, bf + delta)
	coord = coord1 + coord2
	crackFrontEdges = e.getByBoundingBox(*coord)
	flangePart.seedEdgeByNumber(edges=crackFrontEdges, number=nQuart / 2, constraint=FINER)

	# corrections to the st lines
	pickedEdges = e.findAt(((a_crack - r2 / 2, r_crack, bf),), ((a_crack - r2 / 2, -r_crack, bf),),
	                       ((a_crack, r2 / 2, bf),), ((a_crack, -r2 / 2, bf),), ((a_crack + r2 / 2, 0.0, bf),),
	                       ((a_crack + r2 / 2, r2 / 2, bf),), )
	flangePart.seedEdgeByNumber(edges=pickedEdges, number=nRad, constraint=FINER)

	# local seeds outside circular crackFront (biased mesh)
	pickedEdges1 = e.findAt(((a_crack + r2, r2, bf), ), ((a_crack - r2 - delta, r_crack, bf), ))
	pickedEdges2 = e.findAt(((a_crack + r2 + delta, 0.0, bf), ), ((a_crack - r2 - delta, -r_crack, bf), ))
	flangePart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
	    end2Edges=pickedEdges2, minSize=ele_size_crack_min, maxSize=ele_size_bottom, constraint=FINER)

	# # local seeds to additional rectangle and edge of weld - removed in R3
	# pickedEdges = e.findAt(((t_u, (t_u - a_crack) / 2, bf),), ((dmin, r_crack + delta, bf), ),
	#     ((dmin, -r_crack - delta, bf), ), ((t_l, -delta, bf), ), ((t_u, t_u - delta, bf), ),)
	# flangePart.seedEdgeBySize(edges=pickedEdges, size=ele_size_global / eleFac, constraint=FINER)

	# local seeds to the no transition part of the bottom flange
	if t_u != t_l: 
	    pickedEdges = e.findAt(((t_l - delta, 0.0, bf), ),)
	    flangePart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges, 
	        minSize=ele_size_crack_min, maxSize=ele_size_bottom, constraint=FINER)

	# local seeds along height
	pickedEdges1 = e.findAt(((0.0, t_u, bf),), ((t_l, -2*t_l, bf),))
	pickedEdges2 = e.findAt(((t_u, 2*t_u, bf),), ((0, -t_l, bf),))
	flangePart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
	    end2Edges=pickedEdges2, minSize=ele_size_bottom, maxSize=ele_size_top, constraint=FINER)
	pickedEdges = e.findAt(((t_u, delta, bf), ))
	flangePart.seedEdgeBySize(edges=pickedEdges, size=ele_size_bottom, constraint=FINER)

	# local seeds along the width 
	# pickedEdges = e.findAt(((t_u, h_plate, bf/2),), ((t_l, -h_plate, bf/2),),
	#                        ((t_l, 0.0, bf/2),), )
	# flangePart.seedEdgeByNumber(edges=pickedEdges, number=neleAlongWidth, constraint=FINER)

	# seeds as per ele_size_global
	ntop = max(math.ceil(t_u/ele_size_top),2)
	nbot = max(math.ceil(t_u/ele_size_top),2)
	pickedEdges = e.findAt(((t_u/ 2, h_plate, bf), ), )
	flangePart.seedEdgeByNumber(edges=pickedEdges, number=ntop, constraint=FINER)
	pickedEdges = e.findAt(((t_l/2, -h_plate, bf), ),)
	flangePart.seedEdgeByNumber(edges=pickedEdges, number=nbot, constraint=FINER)

	flangePart.generateMesh()

	# %% meshing the top transition!

	topTransitionPart.seedPart(size=ele_size_alongwidth, deviationFactor=0.1, minSizeFactor=0.1)
	e, f, c, v = topTransitionPart.edges, topTransitionPart.faces, topTransitionPart.cells, topTransitionPart.vertices

	region = c.findAt(((delta, r_transition/2, h_plate/2),))
	topTransitionPart.setMeshControls(regions=region, technique=SWEEP, algorithm=ADVANCING_FRONT)
	topTransitionPart.assignStackDirection(referenceRegion=f.findAt((0.0, r_transition / 2, 0.0)),
	                                       cells=c.findAt(((0.0, r_transition / 2, delta),), ))

	nft = int(math.ceil(((tw_u + 2 * r_transition) / ele_size_alongwidth) / 2) * 2)
	pickedEdges = e.findAt(((delta, 0.0, 0.0),), )
	topTransitionPart.seedEdgeByNumber(edges=pickedEdges, number=nft, constraint=FINER)

	nwt = int(math.ceil(tw_u / ele_size_alongwidth))
	pickedEdges = e.findAt(((delta, r_transition, 0.0),), )
	topTransitionPart.seedEdgeByNumber(edges=pickedEdges, number=nwt, constraint=FINER)

	pickedEdges1 = e.findAt(((tw_u / 2, r_transition, delta),), )
	topTransitionPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1,
	                                 minSize=ele_size_bottom, maxSize=ele_size_top, constraint=FINER)

	topTransitionPart.generateMesh()

	# %% meshing the bottom transition!

	botTransitionPart.seedPart(size=ele_size_alongwidth, deviationFactor=0.1, minSizeFactor=0.1)

	e, f, c, v = botTransitionPart.edges, botTransitionPart.faces, botTransitionPart.cells, botTransitionPart.vertices

	region = c.findAt(((delta, r_transition/2, h_plate/2),))
	botTransitionPart.setMeshControls(regions=region, technique=SWEEP, algorithm=ADVANCING_FRONT)
	botTransitionPart.assignStackDirection(referenceRegion=f.findAt((0.0, r_transition / 2, 0.0)),
	                                       cells=c.findAt(((0.0, r_transition / 2, delta),), ))

	nfl = int(math.ceil(((tw_l + 2 * r_transition) / ele_size_alongwidth) / 2) * 2)
	pickedEdges = e.findAt(((delta, 0.0, 0.0),), )
	botTransitionPart.seedEdgeByNumber(edges=pickedEdges, number=nfl, constraint=FINER)

	nwl = int(math.ceil(tw_l / ele_size_alongwidth))
	pickedEdges = e.findAt(((delta, r_transition, 0.0),), )
	botTransitionPart.seedEdgeByNumber(edges=pickedEdges, number=nwl, constraint=FINER)

	pickedEdges1 = e.findAt(((tw_l / 2, r_transition, bw / 4),), )
	botTransitionPart.seedEdgeByBias(biasMethod=SINGLE, end2Edges=pickedEdges1,
	                                 minSize=ele_size_bottom, maxSize=ele_size_top, constraint=FINER)

	botTransitionPart.generateMesh()

	# %% meshing the web portion

	webPart.seedPart(size=ele_size_alongwidth, deviationFactor=0.1, minSizeFactor=0.1)
	e, f, c, v = webPart.edges, webPart.faces, webPart.cells, webPart.vertices

	pickedCells = c.findAt(((0.0, tw_u, bw / 4),), ((0.0, -tw_l, bw / 4),))
	webPart.assignStackDirection(referenceRegion=f.findAt(coordinates=(0.0, tw_u, bw / 2)), cells=pickedCells)
	webPart.setMeshControls(regions=pickedCells, algorithm=MEDIAL_AXIS)

	pickedCells = webPart.cells.findAt(((tw_u / 2 - delta, delta, bw / 4),), )
	webPart.assignStackDirection(referenceRegion=f.findAt((tw_u / 2 - delta, delta, bw / 2)), cells=pickedCells)

	pickedEdges1 = e.findAt(((tw_l / 2, -tw_l, 0.0),), ((-tw_u / 2, tw_u, 0.0),))
	pickedEdges2 = e.findAt(((tw_u / 2, tw_u, 0.0),), ((-tw_l / 2, -tw_l, 0.0),))
	webPart.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1,
	                       end2Edges=pickedEdges2, minSize=ele_size_bottom, maxSize=ele_size_top,
	                       constraint=FINER)

	# adding seeds to the weld part of the web
	# neleWeldWeb = int(math.ceil((tw_u - aw_crack)/(ele_size_bottom)))
	pickedEdges = e.findAt((((tw_u+tw_l)/4, (tw_u - aw_crack)/2.0, 0.0),), )
	webPart.seedEdgeBySize(edges=pickedEdges, size=ele_size_bottom, constraint=FINER)

	webPart.generateMesh()


	return r2, r_crack

    
def createFullModel(tu, tl, bf, twu, twl, bw, rtran, a_crack, aw_crack, length_crack, length_bottom, length_top, modelName, sigma_uniform, sigma_primary, sigma_secondary, tau_x, tau_z):
	'''
	Gx and Gy are gradient in flange wrt to flange average stress - in decimel - not percentage
	'''

	# section input details 
	r_transition = rtran

	# analysis step details 
	max_inc = 0.1
	init_inc = 0.02
	# analysis and contour evaluation details
	flag_NLGEOM = 1                         # flag for large deformations - 1: NLGEOM on, 0: NLGEOM off (ON is by default)
	flag_True = 1                           # flag for true stress strain; 0 for using engg stress strain (which is incorrect)
	flag_NLFEM = 1
	nContour = 5                            # No of contours to be evaluated from disc edge

	# mesh details
	ele_size_bottom = round(0.075*tu,3)                             # this was figured out after lots of trial and error
	ele_size_top = tu
	neleAlongWidth=40

	# some variables
	depthBotSec = 2*tl + bw + 2*rtran
	depthTopSec = 2*tu + bw + 2*rtran

	# ---------------------------------------------------------------------------
	# Create model
	# ---------------------------------------------------------------------------

	myModel = mdb.Model(name=modelName)
	# del mdb.models['Model-1']


	# # ---------------------------------------------------------------------------
	# # Create materials
	# # ---------------------------------------------------------------------------

	Es = 29000
	nu = 0.3

	steel_plasticity_data = []
	if flag_True == 1:
	    fhand = open('steelPlasticTrue.txt', 'r')
	else:
	    fhand = open('steelPlastic.txt', 'r')
	for line in fhand:
	    value_list = (line.strip()).split()
	    value_tuple = (float(value_list[0]), float(value_list[1]))
	    steel_plasticity_data.append(value_tuple)
	steel_plasticity_table = tuple(steel_plasticity_data)
	steelHardeningMaterial = myModel.Material(name='Steel_Hardening')
	steelHardeningMaterial.Plastic(table=steel_plasticity_table)
	steelHardeningMaterial.Elastic(table=((Es, nu),))
	fhand.close()
	# steelHardeningMaterial 'Steel_Hardening'

	weld_plasticity_data = []
	if flag_True == 1:
	    fhand = open('weldPlasticTrue.txt', 'r')
	else:
	    fhand = open('weldPlastic.txt', 'r')
	for line in fhand:
	    value_list = (line.strip()).split()
	    value_tuple = (float(value_list[0]), float(value_list[1]))
	    weld_plasticity_data.append(value_tuple)
	weld_plasticity_table = tuple(weld_plasticity_data)
	weldHardeningMaterial = myModel.Material(name='Weld_Hardening')
	weldHardeningMaterial.Plastic(table=weld_plasticity_table)
	weldHardeningMaterial.Elastic(table=((Es, nu),))
	fhand.close()
	# weldHardeningMaterial 'Weld_Hardening'

	matName = 'Steel_Elastic'
	steelElasticMaterial = myModel.Material(name=matName)
	steelElasticMaterial.Elastic(table=((Es, nu),))
	# steelElasticMaterial 'Steel_Elastic'

	matName = 'Weld_Elastic'
	weldElasticMaterial = myModel.Material(name=matName)
	weldElasticMaterial.Elastic(table=((Es, nu),))
	# weldElasticMaterial 'Weld_Elastic'


	# # ---------------------------------------------------------------------------
	# # Create sections
	# # ---------------------------------------------------------------------------

	if flag_NLFEM == 1:
	    steelMatName = 'Steel_Hardening'
	    weldMatName = 'Weld_Hardening'
	else:
	    steelMatName = 'Steel_Elastic'
	    weldMatName = 'Weld_Elastic'

	steelSectionName = 'Steel_Section'
	steelSection = myModel.HomogeneousSolidSection(name=steelSectionName, material=steelMatName, thickness=None)
	weldSectionName = 'Weld_Section'
	weldSection = myModel.HomogeneousSolidSection(name=weldSectionName, material=weldMatName, thickness=None)


	# # ---------------------------------------------------------------------------
	# # Create parts
	# # ---------------------------------------------------------------------------


	flangePartName = 'flange'
	topTransitionPartName = 'TopTransition'
	botTransitionPartName = 'BotTransition'
	topSectionPartName = 'TopSection'
	topSectionCutPartName = 'TopSectionCut'
	botSectionPartName = 'BotSection'
	botSectionCutPartName = 'BotSectionCut'

	createWSection(tl, bf, twl, bw, rtran, length_bottom, a_crack, modelName, botSectionPartName, steelSectionName)

	flag = 'Bot'
	createWSectionCut(tl, bf, twl, bw, rtran, length_crack, a_crack, modelName, botSectionCutPartName, steelSectionName, ele_size_bottom, ele_size_top, flag)

	createWSection(tu, bf, twu, bw, rtran, length_top, a_crack, modelName, topSectionPartName, steelSectionName)

	flag = 'Top'
	createWSectionCut(tu, bf, twu, bw, rtran, length_crack, a_crack, modelName, topSectionCutPartName, steelSectionName, ele_size_bottom, ele_size_top, flag)

	r2, r_crack = createCrackedFlange(tu, tl, bf, twu, twl, bw, rtran, length_crack, a_crack, aw_crack, modelName, steelSectionName, weldSectionName, ele_size_bottom, ele_size_top, neleAlongWidth=40)

	flangePart = myModel.parts[flangePartName]
	topTransitionPart = myModel.parts[topTransitionPartName]
	botTransitionPart = myModel.parts[botTransitionPartName]
	topSectionPart = myModel.parts[topSectionPartName]
	topSectionCutPart = myModel.parts[topSectionCutPartName]
	botSectionPart = myModel.parts[botSectionPartName]
	botSectionCutPart = myModel.parts[botSectionCutPartName]


	# ---------------------------------------------------------------------------
	# Create assembly
	# ---------------------------------------------------------------------------

	myAssembly = myModel.rootAssembly
	flangeInstName = 'flangeInst'
	topTransitionInstName = 'topTransitionInst'
	botTransitionInstName = 'botTransitionInst'
	topSectionInstName = 'topSectionInst'
	topSectionCutInstName = 'topSectionCutInst'
	botSectionInstName = 'botSectionInst'
	botSectionCutInstName = 'botSectionCutInst'
	webInstName = 'webInst'

	flangeInst = myAssembly.Instance(name=flangeInstName, part=flangePart, dependent=ON)
	topTransitionInst = myAssembly.Instance(name=topTransitionInstName, part=topTransitionPart, dependent=ON)
	botTransitionInst = myAssembly.Instance(name=botTransitionInstName, part=botTransitionPart, dependent=ON)
	# webInst = myAssembly.Instance(name=webInstName, part=webPart, dependent=ON)
	topSectionInst = myAssembly.Instance(name=topSectionInstName, part=topSectionPart, dependent=ON)
	topSectionCutInst = myAssembly.Instance(name=topSectionCutInstName, part=topSectionCutPart, dependent=ON)
	botSectionInst = myAssembly.Instance(name=botSectionInstName, part=botSectionPart, dependent=ON)
	botSectionCutInst = myAssembly.Instance(name=botSectionCutInstName, part=botSectionCutPart, dependent=ON)


	# setting the transition parts
	h_plate = length_crack
	ht = h_plate - r_crack
	x0 = np.array([0.0, 0.0, ht])
	x1 = np.array([0.0, h_plate, bf / 2])
	dx = tuple(x1 - x0)
	myAssembly.translate(instanceList=(topTransitionInstName,), vector=dx)
	myAssembly.rotate(instanceList=(topTransitionInstName,), axisPoint=(0.0, h_plate, 0.0),
	                  axisDirection=(0.0, 0.0, 1.0), angle=90.0)
	myAssembly.rotate(instanceList=(topTransitionInstName,), axisPoint=(0.0, h_plate, bf / 2),
	                  axisDirection=(-1.0, 0.0, 0.0), angle=90.0)

	x0 = np.array([0.0, 0.0, ht])
	x1 = np.array([0.0, h_plate, bf / 2])
	dx = tuple(x1 - x0)
	myAssembly.translate(instanceList=(botTransitionInstName,), vector=dx)
	myAssembly.rotate(instanceList=(botTransitionInstName,), axisPoint=(0.0, h_plate, 0.0),
	                  axisDirection=(0.0, 0.0, 1.0), angle=90.0)
	myAssembly.rotate(instanceList=(botTransitionInstName,), axisPoint=(0.0, h_plate, bf / 2),
	                  axisDirection=(-1.0, 0.0, 0.0), angle=90.0)
	dx = (0.0, -(h_plate + r_crack), 0.0)
	myAssembly.translate(instanceList=(botTransitionInstName,), vector=dx)

	# # setting the web parts
	# x0 = np.array([0.0, h_plate, bw / 2])
	# x1 = np.array([-r_transition, h_plate, bf / 2])
	# dx = tuple(x1 - x0)
	# myAssembly.translate(instanceList=(webInstName,), vector=dx)
	# myAssembly.rotate(instanceList=(webInstName,), axisPoint=(-r_transition, h_plate, bf / 2),
	#                   axisDirection=(0.0, 1.0, 0.0), angle=90.0)

	# # add the four sections here --------------
	myAssembly.rotate(instanceList=(topSectionInstName, topSectionCutInstName, botSectionInstName, botSectionCutInstName), axisPoint=(0.0, 0.0, 0.0),
	                  axisDirection=(1.0, 0.0, 0.0), angle=90.0)

	x0 = np.array([bw/2+rtran, 0.0, 0.0])
	x1 = np.array([0.0, 0.0, bf/2])
	dx = tuple(x1-x0)
	myAssembly.translate(instanceList=(topSectionInstName, topSectionCutInstName, botSectionInstName, botSectionCutInstName), vector=dx)

	x0 = np.array([0.0, 0.0, 0.0])
	x1 = np.array([0.0, length_crack, 0.0])
	dx = tuple(x1-x0)
	myAssembly.translate(instanceList=(topSectionCutInstName,), vector=dx)

	x0 = np.array([0.0, 0.0, 0.0])
	x1 = np.array([0.0, length_crack+length_top, 0.0])
	dx = tuple(x1-x0)
	myAssembly.translate(instanceList=(topSectionInstName,), vector=dx)

	x0 = np.array([0.0, 0.0, 0.0])
	x1 = np.array([0.0, -length_crack, 0.0])
	dx = tuple(x1-x0)
	myAssembly.translate(instanceList=(botSectionInstName,), vector=dx)


	# ---------------------------------------------------------------------------
	# Create Constraints
	# ---------------------------------------------------------------------------

	# Tie constraint 1 and 2 - vertical between flange and transition and web

	surfFace1 = flangeInst.faces.findAt(((0.0, tu, bf / 2),), ((0.0, -tl, bf / 2),))
	mSurfFlange = myAssembly.Surface(side1Faces=surfFace1, name='m_surf_flange')

	surfFace1 = topTransitionInst.faces.findAt(((0.0, tu, bf / 2),))
	surfFace2 = botTransitionInst.faces.findAt(((0.0, -tl, bf / 2),))
	sSurfTransition = myAssembly.Surface(side1Faces=surfFace1 + surfFace2, name='s_surf_transition')

	surfFace1 = topTransitionInst.faces.findAt(((-r_transition, tu, bf / 2),))
	surfFace2 = botTransitionInst.faces.findAt(((-r_transition, -tl, bf / 2),))
	sSurfTransition2 = myAssembly.Surface(side1Faces=surfFace1 + surfFace2, name='s_surf_transition2')

	surfFace1 = topSectionCutInst.faces.findAt(((-r_transition, tu, bf / 2 + twu/4),), ((-r_transition, tu, bf / 2 - twu/4),))
	surfFace2 = botSectionCutInst.faces.findAt(((-r_transition, -tl, bf / 2 + twl/4),), ((-r_transition, -tl, bf / 2 - twl/4),))
	mSurfWeb = myAssembly.Surface(side1Faces=surfFace1 + surfFace2, name='m_surf_web')

	myModel.Tie(name='Constraint-1', master=mSurfFlange, slave=sSurfTransition,
	            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
	myModel.Tie(name='Constraint-2', master=mSurfWeb, slave=sSurfTransition2,
	            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

	# Tie constraint 3
	coordMin = (-(depthBotSec-tl)-0.001, -length_crack, -0.001 )
	coordMax = (tl + 0.001, -length_crack, bf + 0.001 )
	boundingCoords = coordMin + coordMax	
	surfFace1 = botSectionInst.faces.getByBoundingBox(*boundingCoords)
	mSurf1 = myAssembly.Surface(side1Faces=surfFace1, name='m_surf_1')

	coordMin = (-(depthBotSec-tl)-0.001, -length_crack, -0.001 )
	coordMax = (tl + 0.001, -length_crack, bf + 0.001 )
	boundingCoords = coordMin + coordMax	
	surfFace1 = botSectionCutInst.faces.getByBoundingBox(*boundingCoords)
	surfFace2 = flangeInst.faces.getByBoundingBox(*boundingCoords)
	surfFace3 = botTransitionInst.faces.getByBoundingBox(*boundingCoords)
	sSurf1 = myAssembly.Surface(side1Faces=surfFace1 + surfFace2 + surfFace3, name='s_surf_1')

	myModel.Tie(name='Constraint-3', master=mSurf1, slave=sSurf1,
	            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

	# Tie constraint 4 
	surfFace1 = botSectionCutInst.faces.findAt(((-r_transition-bw/2, 0.0, bf/2 - twl/2 + 0.001),), ((-(depthBotSec-tl)+0.001, 0.0, bf/2),))
	mSurf2 = myAssembly.Surface(side1Faces=surfFace1, name='m_surf_2')

	surfFace1 = topSectionCutInst.faces.findAt(((-r_transition-bw/2, 0.0, bf/2 - twu/2 + 0.001),), ((-(depthTopSec-tu)+0.001, 0.0, bf/2),))
	sSurf2 = myAssembly.Surface(side1Faces=surfFace1, name='s_surf_2')

	myModel.Tie(name='Constraint-4', master=mSurf2, slave=sSurf2,
	            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

	# Tie constraint 5
	coordMin = (-(depthBotSec-tl)-0.001, length_crack, -0.001 )
	coordMax = (tl + 0.001, length_crack, bf + 0.001 )
	boundingCoords = coordMin + coordMax	
	surfFace1 = topSectionCutInst.faces.getByBoundingBox(*boundingCoords)
	surfFace2 = flangeInst.faces.getByBoundingBox(*boundingCoords)
	surfFace3 = topTransitionInst.faces.getByBoundingBox(*boundingCoords)
	mSurf3 = myAssembly.Surface(side1Faces=surfFace1 + surfFace2 + surfFace3, name='m_surf_3')

	coordMin = (-(depthBotSec-tl)-0.001, length_crack, -0.001 )
	coordMax = (tl + 0.001, length_crack, bf + 0.001 )
	boundingCoords = coordMin + coordMax	
	surfFace1 = topSectionInst.faces.getByBoundingBox(*boundingCoords)
	sSurf3 = myAssembly.Surface(side1Faces=surfFace1, name='s_surf_3')

	myModel.Tie(name='Constraint-5', master=mSurf3, slave=sSurf3,
	            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

	# Interaction 1 (gap element to simulate compression flange behavior)
	myModel.ContactProperty('IntProp-1')
	myModel.interactionProperties['IntProp-1'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, 
		constraintEnforcementMethod=DEFAULT)

	surfFace1 = botSectionCutInst.faces.findAt(((-r_transition-bw/2, 0.0, bf/2 + twl/2 - 0.001),), ((-(depthBotSec-2*tl)-0.001, 0.0, bf/2),), 
		((-(depthBotSec-2*tl)+0.001, 0.0, bf/2),))
	surface2 = botTransitionInst.faces.findAt(((-r_transition/2, -r_crack, bf/2),),)
	surface3 = flangeInst.faces.findAt(((a_crack - r2 - 0.001, -r_crack, bf/2),), ((a_crack - 0.001, -r_crack, bf/2),))
	mSurfInt = myAssembly.Surface(side1Faces=surfFace1+surface2+surface3, name='m_surf_int')

	surfFace1 = topSectionCutInst.faces.findAt(((-r_transition-bw/2, 0.0, bf/2 + twu/2 - 0.001),), ((-(depthTopSec-2*tu)-0.001, 0.0, bf/2),),
		((-(depthTopSec-2*tu)+0.001, 0.0, bf/2),))
	surface2 = topTransitionInst.faces.findAt(((-r_transition/2, r_crack, bf/2),),)
	surface3 = flangeInst.faces.findAt(((a_crack - r2 - 0.001, r_crack, bf/2),), ((a_crack - 0.001, r_crack, bf/2),))
	sSurfInt = myAssembly.Surface(side1Faces=surfFace1+surface2+surface3, name='s_surf_int')

	myModel.SurfaceToSurfaceContactStd(name='Int-1', createStepName='Initial', master=mSurfInt, slave=sSurfInt, sliding=FINITE, 
		thickness=ON, interactionProperty='IntProp-1', adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

	# ---------------------------------------------------------------------------
	# Create step
	# ---------------------------------------------------------------------------

	loadStepName = 'Loading'
	if flag_NLGEOM == 1:
	    myModel.StaticStep(name=loadStepName, previous='Initial', initialInc=init_inc, maxInc=max_inc, nlgeom=ON)
	else:
	    myModel.StaticStep(name=loadStepName, previous='Initial', initialInc=init_inc, maxInc=max_inc, nlgeom=OFF)


	# ---------------------------------------------------------------------------
	# Create J integral region
	# ---------------------------------------------------------------------------

	J_integral_name = 'J-integralFlange'

	c1 = (a_crack, 0.0, 0.0)
	c2 = (a_crack, 0.0, bf)
	crackFrontCells = flangeInst.cells.getByBoundingCylinder(c1, c2, r2 + 0.001)

	# coordMin = (-delta, -(t_l + delta), -delta )
	# coordMax = (t_l + delta, t_u + delta, bf + delta )
	# boundingCoords = coordMin + coordMax
	# crackFrontCells = flangeInst.cells.getByBoundingBox(*boundingCoords)

	crackFront = regionToolset.Region(cells=crackFrontCells)

	crackTipEdge = flangeInst.edges.findAt(((a_crack + r_crack, 0.0, bf / 2),))
	crackTip = regionToolset.Region(edges=crackTipEdge)

	myAssembly.engineeringFeatures.ContourIntegral(name=J_integral_name, symmetric=OFF,
	                                               crackFront=crackFront, crackTip=crackTip,
	                                               extensionDirectionMethod=Q_VECTORS,
	                                               qVectors=(((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),), midNodePosition=0.5,
	                                               collapsedElementAtTip=NONE)


	# ---------------------------------------------------------------------------
	# Create field output request & history output request
	# ---------------------------------------------------------------------------

	J_integral_out = 'J-flangeOutput'
	myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'U'))
	myModel.historyOutputRequests['H-Output-1'].setValues(variables=('ALLIE', 'ETOTAL'))
	myModel.HistoryOutputRequest(name=J_integral_out, createStepName=loadStepName, contourIntegral=J_integral_name,
	                             sectionPoints=DEFAULT, rebar=EXCLUDE, numberOfContours=nContour)
	del myModel.historyOutputRequests['H-Output-1']


	# # ---------------------------------------------------------------------------
	# # Apply loads
	# # ---------------------------------------------------------------------------

	# top surface
	coordMin = (-(depthTopSec-tu)-0.001, length_crack+length_top, -0.001 )
	coordMax = (tu + 0.001, length_crack+length_top, bf + 0.001 )
	boundingCoords = coordMin + coordMax	
	surfFace1 = topSectionInst.faces.getByBoundingBox(*boundingCoords)
	topSurface = myAssembly.Surface(side1Faces=surfFace1, name='topSurface')

	expression = ''
	# all stresses are extreme fiber stress for corresponding loading scenario
	# input sign convention: tension +ve; compression -ve
	# input converted to abaqus convention: tension -ve pressure; compression +ve pressure
	# the direction of sigma_primary and tau_x has to be pre determined such that the flange on the right side has tension at splice location
	# sigma_total = -(sigma_uniform + sigma_primary + sigma_secondary)
	if sigma_uniform != 0.0 or sigma_primary != 0.0 or sigma_secondary != 0.0: 
		print('inside the loop - check is okay!')
		if sigma_primary == 0.0 and sigma_secondary == 0.0: 
			myModel.Pressure(name='StressField', createStepName=loadStepName,
		        region=topSurface, distributionType=UNIFORM, field='', magnitude=-sigma_uniform, amplitude=UNSET)
		else: 
			uniform_factor = str(float(-sigma_uniform))
			primary_factor = str(float(-sigma_primary))
			secondary_factor = str(float(sigma_secondary))
			expr_primary = primary_factor + ' * ( 1.0 - 2 * ' + str(tu) + ' / ' + str(depthTopSec) + ' + 2 * X / ' + str(depthTopSec) + ' )'
			expr_secondary = secondary_factor + ' * ( 2 * Z / ' + str(bf) + ' - 1.0)'
			expression += uniform_factor + ' + ' + expr_primary + ' + ' + expr_secondary
			fieldName = 'non-uniform'
			myModel.ExpressionField(name=fieldName, localCsys=None, description='', expression=expression)
			myModel.Pressure(name='StressField', createStepName=loadStepName, region=topSurface, distributionType=FIELD, 
				field=fieldName, magnitude=1.0, amplitude=UNSET)

	# top surface shear - x direction 
	if tau_x != 0.0: 
		myModel.SurfaceTraction(name='ShearX', createStepName=loadStepName, region=topSurface, magnitude=tau_x, 
			directionVector=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), distributionType=UNIFORM, 
			field='', localCsys=None, resultant=ON)

	# top surface shear - y direction
	if tau_z != 0.0:
		myModel.SurfaceTraction(name='ShearZ', createStepName=loadStepName, region=topSurface, magnitude=tau_z, 
			directionVector=((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), distributionType=UNIFORM, 
			field='', localCsys=None, resultant=ON)

	# ---------------------------------------------------------------------------
	# Apply boundary conditions
	# ---------------------------------------------------------------------------

	# bottom surface
	coordMin = (-(depthBotSec-tl)-0.001, -(length_crack + length_bottom), -0.001 )
	coordMax = (tl + 0.001, -(length_crack + length_bottom), bf + 0.001 )
	boundingCoords = coordMin + coordMax	
	surfFace1 = botSectionInst.faces.getByBoundingBox(*boundingCoords)
	botSurface = myAssembly.Set(faces=surfFace1, name='botSurface')
	myModel.DisplacementBC(name='flangeBottom', createStepName='Initial',
	                       region=botSurface, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
	                       amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

	edges1 = flangeInst.edges.findAt(((tl, 0.0, bf / 2),))
	myAssembly.Set(edges=edges1, name='alongWidth')

	edges1 = flangeInst.edges.findAt(((a_crack+r_crack, 0.0, bf / 2),))
	myAssembly.Set(edges=edges1, name='alongCrack')


	# ---------------------------------------------------------------------------
	# Create JOB
	# ---------------------------------------------------------------------------

	jobName = modelName

	if flag_NLGEOM == 0:
	    jobName += 'NLGEOM_OFF'

	mdb.Job(name=jobName, model=modelName, description='', type=ANALYSIS,
	        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
	        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
	        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
	        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
	        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
	        numGPUs=0)
	mdb.jobs[jobName].writeInput(consistencyChecking=OFF)


	session.viewports['Viewport: 1'].setValues(displayedObject=myAssembly)
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, 
		interactions=OFF, constraints=OFF, connectors=OFF, engineeringFeatures=OFF)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# %% CREATE MODEL
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# %% INPUT FOR DIRECT RUN------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# modelName = 'BNM-Loading-Case2'					# Model name
# tu = 3.5
# tl = 4.52
# twu = 2.1875
# twl = 2.8125
# bf = 17.000
# bw = 10.000
# rtran = 1.25 
# a_crack = 0.5*tu
# aw_crack = 0.5*twu

# story_ht = 10.0*12.0
# splice_ht = 5.0*12.0

# # sign convention - tension positive 
# # directions as per abaqus 

# Mz1 = 45000.0
# Mz2 = 45000.0
# Py = 0.0
# Mx1 = 0.0
# Mx2 = 0.0

# %% INPUT FOR TOOLBOX ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

with open('geomData.json') as f:
    geomData = json.load(f)
with open('loadData.json') as f:
    loadData = json.load(f)

modelName = str(geomData['modelName'])	

tu = geomData['tu']
tl = geomData['tl']
twu = geomData['twu']
twl = geomData['twl']
bf = geomData['bf']
bw = geomData['bw']
rtran = geomData['rtran']
a_crack = geomData['a_crack']
aw_crack = geomData['aw_crack']

story_ht = geomData['story_ht']
splice_ht = geomData['splice_ht']

# sign convention - tension positive 
# directions as per abaqus 

Mz1 = loadData['Mz1']
Mz2 = loadData['Mz2']
Py = loadData['Py']
Mx1 = loadData['Mx1']
Mx2 = loadData['Mx2']

# ###############
# %% CALCULATIONS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ###############

depth_top_sec = 2*tu + 2*rtran + bw
Iz = 2*((bf*tu**3)/12 + (bf*tu)*(bw/2+rtran+tu/2)**2) + (twu*(bw+rtran*2)**3)/12
Sz = Iz/(depth_top_sec/2)
Ix = 2*((tu*bf**3)/12) + ((bw+rtran*2)*twu**3)/12
Sx = Ix/(bf/2)
Area = 2*(bf*tu) + (bw+2*rtran)*twu + (4*rtran**2 - math.pi*rtran**2)

splice_ht_top = story_ht - splice_ht

length_crack = 10.0*tu
length_bottom = splice_ht - length_crack
length_top = (story_ht - splice_ht) - length_crack

# load calculation
Vx1 = (Mz1-Mz2)/story_ht
Vz1 = (Mx2-Mx1)/story_ht

Mz_splice = Mz1 - Vx1*splice_ht_top
Mx_splice = Mx1 + Vz1*splice_ht_top

"""
Comments R1 - 
Earlier the following check was performed to change sign of loads such that the right flange is always critical. 
User has to now be sure of the loading directions. 
"""
# if Mz_splice < 0.0: 
# 	Mz1 = -Mz1
# 	Mz2 = -Mz2
# 	Vx1 = -Vx1
# 	Mz_splice = -Mz_splice

sigma_uniform = round(Py/Area,3)
sigma_primary = round(Mz1/Sz,3)
sigma_secondary = round(Mx1/Sx,3)
tau_x = round(Vx1/Area,3)
tau_z = round(Vz1/Area,3)

# print('sigma unifrom is ',sigma_uniform)
# print('sigma primary is ',sigma_primary)
# print('sigma secondary is ',sigma_secondary)
# print('tau_x is ',tau_x)

# %% Splice flange in complete tension check and Create Model -------------------------------------------------------------------------------------------------------------------------------
x1 = (depth_top_sec/2)
x2 = x1 - tu
s1 = Py/Area + Mz_splice*x1/Iz + Mx_splice/Sx
s2 = Py/Area + Mz_splice*x1/Iz - Mx_splice/Sx
s3 = Py/Area + Mz_splice*x2/Iz + Mx_splice/Sx
s4 = Py/Area + Mz_splice*x2/Iz - Mx_splice/Sx

"""
Comments R1 - 
Earlier the following checks gave a warning to the user if any part of the right flange was under compression.
Now this is taken care of in the post processing by checking the stresses at crack tip to remove portion of flange in compression in Pf calculations.
"""
# if s1 < 0.0 or s2 < 0.0 or s3 < 0.0 or s4 < 0.0: 
# 	print('[WARNING] The tension flange is in partial or full compression - J results will be useless!')
# 	createFullModel(tu, tl, bf, twu, twl, bw, rtran, a_crack, aw_crack, length_crack, length_bottom, length_top, modelName, sigma_uniform, sigma_primary, sigma_secondary, tau_x, tau_z)
# else: 
# 	createFullModel(tu, tl, bf, twu, twl, bw, rtran, a_crack, aw_crack, length_crack, length_bottom, length_top, modelName, sigma_uniform, sigma_primary, sigma_secondary, tau_x, tau_z)

createFullModel(tu, tl, bf, twu, twl, bw, rtran, a_crack, aw_crack, length_crack, length_bottom, length_top, modelName, sigma_uniform, sigma_primary, sigma_secondary, tau_x, tau_z)

caeFileName = modelName
mdb.saveAs(pathName=caeFileName)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------