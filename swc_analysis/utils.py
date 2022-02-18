import numpy as np
import copy
import params
import xml.etree.ElementTree as ET

from shutil import copyfile
from os.path import join, exists
from os import system, mkdir
from Node import Node
from sys import exit
from scipy.stats import ranksums, kruskal


def angleBetweenVectors(vecA, vecB):

	uvecA = vecA / np.linalg.norm(vecA)
	uvecB = vecB / np.linalg.norm(vecB)
	dotProduct = np.dot(uvecA, uvecB)

	if dotProduct < -1:
		dotProduct = -1

	if dotProduct > 1:
		dotProduct = 1	

	angle = np.arccos(dotProduct)

	return angle



## statistical test functions

def p_value_to_string(p):

	"""
	Converts p value to asterisks based on significance

	"""

	if p >= .1:
		return ''
	elif .01 <= p < .1:
		return '*'
	elif .001 <= p < .01:
		return '**'
	elif .0001 <= p < .001:
		return '***'
	elif .00001 <= p < .0001:
		return '****'
	else:
		return '*****'

def wilcoxon(dataA, dataB):

	"""
	Wrapper for wilcoxon rank sum test
	Returns p valye and significance asterisks
	
	"""

	_, p = ranksums(dataA, dataB)

	pString = pValueToString(p)

	return p, pString

def kruskalWallis(data):
	
	"""
	Wrapper for kruskal wallis test

	"""

	if len(data) <= 1:
		pString = ''
	else:

		if len(data) == 2:
			_, p = kruskal(data[0], data[1])
		elif len(data) == 3:
			_, p = kruskal(data[0], data[1], data[2])

		pString = pValueToString(p)

	return pString
		

## Persistent Homology Functions

def arcCosDissimilarity(vectorA, vectorB):
	
	"""
	Calculates arccosin dissimilarity metric between two vectors

	"""

	cosTheta = vectorA.dot(vectorB) / (np.linalg.norm(vectorA)*np.linalg.norm(vectorB))
	if cosTheta > 1:
		cosTheta = 1

	return np.arccos(cosTheta) / (np.pi/2)





###### Utility functions for Oblique -> Coronal Transformation #######

# These functions assume transformations are occuring on images
# therefore, the size of the image must be provided and the
# input coordinates are assumed to be in pixel coordinates of the oblique image

def oblique_to_coronal(coords, fusedDims, shearFactor):

	if not isinstance(fusedDims,dict):
		exit('Error: fuseDims must be of type dict')

	# convert to numpy arrays
	coords = np.array(coords)
	imageLengths = np.array([fusedDims['x']['length'],fusedDims['y']['length'],fusedDims['z']['length']])

	# shift coordinates 
	coords -= np.array([fusedDims['x']['min'],fusedDims['y']['min'],fusedDims['z']['min']])		

	# perform transformations
	coords, imageLengths = reslice(coords,imageLengths)
	coords = verticalFlip(coords, imageLengths)
	coords, imageLengths = shear(coords,imageLengths,shearFactor)
	coords, imageLengths = reslice(coords,imageLengths)
	coords, imageLengths = rotate(coords,imageLengths)
	
	return list(coords)
		

def reslice(coords, imageLengths):

	return coords[[1,2,0]], imageLengths[[1,2,0]]


def verticalFlip(coords, imageLengths):
	
	coords[1] = imageLengths[1] - coords[1] - 1

	return coords


def shear(coords, imageLengths, shearFactor):

	imageLengthsShear = copy.deepcopy(imageLengths)
	imageLengthsShear[1] = imageLengths[1] + abs(shearFactor)*imageLengthsShear[0]
	imageLengthsShear = np.round(imageLengthsShear).astype(int)

	coords[1] = coords[1] - imageLengths[1]/2 + shearFactor*(coords[0]-imageLengths[0]/2) + imageLengthsShear[1]/2

	return coords, imageLengthsShear

def rotate(coords,imageLengths):

	coords = coords[[1,0,2]]
	coords[0] = imageLengths[1] - coords[0] - 1
	
	return coords, imageLengths[[1,0,2]]

#################################################################

######### Utility Functions for Branch Markup #####################

def parseBranchMarkupFile(path):

	coords, ids = [], [1]

	with open(path,'r') as fp:
		lines = [x for x in fp.readlines() if x != '\n']
		for i,line in enumerate(lines):
			if i == 0:
				coord = [int(x) for x in line.split(': ')[1].split(',')]
			else:
				coord = [int(x) for x in line.split(': ')[-1].split(' ')]
			coords.append(coord)
			if i > 0:
				ids.append(line.split(' :')[:2])

	return np.array(coords), ids

def writeBranchMarkupFile(path, coords, ids):

	with open(path, 'w') as fp:
		for i,coord in enumerate(coords):
			if i == 0:
				fp.write('Soma: ' + str(coord[0]) + ',' + str(coord[1]) + ',' + str(coord[2]) + '\n\n')
			else:
				fp.write(ids[i][0] + ' : ' + ids[i][1] +' : ' + str(coord[0]) + ',' + str(coord[1]) + ',' + str(coord[2]) + '\n')

def read_branch_markups(path):

	label_info = []

	if not exists(path):
		return False

	with open(path,'r') as fp:
		lines = fp.readlines()
		for line in lines:
			if '\t' in line:
				id_ = int(line.split(':')[1])
				label = int([l for l in line.split('\t') if l != '' and l != '\n'][-1][0])
				label_info.append((id_,label))
	
	return label_info

#############################################################

def zyToVolumeID(z,y):
	return 'Z' + '{:02}'.format(z) + '_Y' + '{:02}'.format(y)


def loadCropInfo():

	cropInfo = np.genfromtxt(params.CROP_INFO_PATH,dtype=object)
	cropInfo[:,0] = cropInfo[:,0].astype(str)
	cropInfo[:,1] = cropInfo[:,1].astype(int)

	cropDict = {}

	for cropping in cropInfo:
		cropDict[cropping[0]] = cropping[1]

	return cropDict


def writeCoordsForRegistration(path, coords):

	with open(path,'w') as fp:
		fp.write('index\n')
		fp.write(str(len(coords))+'\n')
		for i,coord in enumerate(coords):
			fp.write(str(coord[0])+' ')
			fp.write(str(coord[1])+' ')
			fp.write(str(coord[2]))
			if i != len(coords)-1:
				fp.write('\n')

def parseTransformixOutput(path):

	registeredCoords = np.zeros(shape=(0,3))
	with open(path,'r') as fp:
		lines = fp.readlines()
		for line in lines:
			coords = np.array([x for x in line.split(';') if 'OutputPoint' in x][0].split(' ')[4:7],dtype=np.float)
			registeredCoords = np.vstack((registeredCoords,coords))

	return registeredCoords

def runTransformix(outPath, tpPath, inPath):
	
	system('transformix -out ' + outPath  + ' -tp ' + tpPath + ' -def ' + inPath)



def Bresenham3D(coordA, coordB): 

	"""
	Bresenhams algorithm to get pixel coordinates of line between 2 pixels in 3D

	"""

	x1, y1, z1 = coordA
	x2, y2, z2 = coordB
	

	ListOfPoints = [] 
	ListOfPoints.append((x1, y1, z1)) 
	dx = abs(x2 - x1) 
	dy = abs(y2 - y1) 
	dz = abs(z2 - z1) 
	if (x2 > x1): 
		xs = 1
	else: 
		xs = -1
	if (y2 > y1): 
		ys = 1
	else: 
		ys = -1
	if (z2 > z1): 
		zs = 1
	else: 
		zs = -1
  
	# Driving axis is X-axis" 
	if (dx >= dy and dx >= dz):         
		p1 = 2 * dy - dx 
		p2 = 2 * dz - dx 
		while (x1 != x2):
			x1 += xs 
			if (p1 >= 0): 
				y1 += ys 
				p1 -= 2 * dx 
			if (p2 >= 0): 
				z1 += zs 
				p2 -= 2 * dx 
			p1 += 2 * dy 
			p2 += 2 * dz 
			ListOfPoints.append((x1, y1, z1)) 

	# Driving axis is Y-axis" 
	elif (dy >= dx and dy >= dz):        
		p1 = 2 * dx - dy 
		p2 = 2 * dz - dy 
		while (y1 != y2): 
			y1 += ys 
			if (p1 >= 0): 
				x1 += xs 
				p1 -= 2 * dy 
			if (p2 >= 0): 
				z1 += zs 
				p2 -= 2 * dy 
			p1 += 2 * dx 
			p2 += 2 * dz 
			ListOfPoints.append((x1, y1, z1)) 

	# Driving axis is Z-axis" 
	else:         
		p1 = 2 * dy - dz 
		p2 = 2 * dx - dz 
		while (z1 != z2): 
			z1 += zs 
			if (p1 >= 0): 
				y1 += ys 
				p1 -= 2 * dz 
			if (p2 >= 0): 
				x1 += xs 
				p2 -= 2 * dz 
			p1 += 2 * dy 
			p2 += 2 * dx 
			ListOfPoints.append((x1, y1, z1)) 
	return ListOfPoints 


def get_swc_info(swc_path):


	"""
	Extracts info from name of swc, return dict

	"""

	# remove path, only keeping file name
	swc_name = swc_path.split('/')[-1]

	# remove extension
	if swc_name.endswith('.swc') or swc_name.endswith('.txt'):
		swc_name = swc_name.split('.')[0]

	# split swc name
	swc_name_parts  = swc_name.split('_')

	# check if layer is included
	if len(swc_name_parts[0]) == 1:
		layer = int(swc_name_parts.pop(0))
	else:
		layer = None

	swc_name = '_'.join(swc_name_parts)


	# special case for 170329_500
	if swc_name_parts[0] == '170329' and swc_name_parts[1] == '500':
		brain = '_'.join(swc_name_parts[:2])
		swc_name_parts = swc_name_parts[2:]
	else:
		brain = swc_name_parts.pop(0)

	swc_id = swc_name_parts[0]

	return {'layer':layer, 'brain':brain, 'id':swc_id, 'name':swc_name}





