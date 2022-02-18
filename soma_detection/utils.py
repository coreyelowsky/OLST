import numpy as np
import params
from skimage import io
from os import system
import xml.etree.ElementTree as ET
import copy

RES_FACTORS = np.array([0.406,0.406,2.5])

def loadCropInfo():

	cropInfo = np.genfromtxt(params.CROP_INFO_PATH,dtype=object)
	cropInfo[:,0] = cropInfo[:,0].astype(str)
	cropInfo[:,1] = cropInfo[:,1].astype(int)

	cropDict = {}

	for cropping in cropInfo:
		cropDict[cropping[0]] = cropping[1]

	return cropDict

###### Utility functions for Oblique -> Coronal Transformation #######

# These functions assume transformations are occuring on images
# therefore, the size of the image must be provided and the
# input coordinates are assumed to be in pixel coordinates of the oblique image

def obliqueToCoronal(coords, fusedDims, shearFactor):

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

def runTransformix(outPath, tpPath, inPath):
	
	system('transformix -out ' + outPath  + ' -tp ' + tpPath + ' -def ' + inPath)

def parseTransformixOutput(path):

	registeredCoords = np.zeros(shape=(0,3))
	with open(path,'r') as fp:
		lines = fp.readlines()
		for line in lines:
			coords = np.array([x for x in line.split(';') if 'OutputPoint' in x][0].split(' ')[4:7],dtype=np.float)
			registeredCoords = np.vstack((registeredCoords,coords))

	return registeredCoords


def pad(image,seg_or_raw,cropped_size):

	x,y,z = croppped_size

	if seg_or_raw == 'seg': 
		if image.shape != (z,y,x):
			image = np.pad(image,((0,z-image.shape[0]),(0,y-image.shape[1]),(0,x-image.shape[2])),'constant')
	elif seg_or_raw == 'raw': 
		if image.shape != (z,y,x):
			image = np.pad(image,((0,z-image.shape[0]),(0,y-image.shape[1]),(0,x-image.shape[2])),'constant',constant_values=np.amin(image))
	return image

# Reads in tif to numpy array
def tif_to_numpyArray(tiff_file_name):
	return io.imread(tiff_file_name)

def distance(a,b,input_units,output_units):
	if input_units == 'voxels' and output_units == 'microns':
		return np.sqrt(np.sum(((np.array(a) - np.array(b))*RES_FACTORS)**2))
	else:
		return 'ERROR: BAD INPUT TYPES'

# Extracts stitching information from xml file generated from BigStitcher
def extract_stitching_parameters(xml_file):

	# read in xml
	tree = ET.parse(xml_file)
	root = tree.getroot() 

	# output lists
	files = []
	registrations = []
	stitchings = []

	# find nodes for files, registrations, stitchings
	for child in root.iter():
		if child.tag == 'files':
			files_node = child
		elif child.tag == 'ViewRegistrations':
			registrations_node = child
		elif child.tag == 'StitchingResults':
			stitchings_node = child

	for child in files_node:
		setup_number = child.attrib['view_setup']
		file_name = child[0].text
		dict_data = {"setup number":setup_number,"file name":file_name}
		files.append(dict_data)

	for child in registrations_node:
		setup_number = child.attrib['setup']
		stitching_transform = np.fromstring(child[0][1].text, sep=' ')
		translation_regular_grid = np.fromstring(child[1][1].text, sep=' ')
		calibration = np.fromstring(child[2][1].text, sep=' ')
		dict_data = {"setup number":setup_number,"stitching transform":stitching_transform,"translation to regular grid":translation_regular_grid,"calibration":calibration}
		registrations.append(dict_data)

	for child in stitchings_node:
		setup_a_number = child.attrib['view_setup_a']
		setup_b_number = child.attrib['view_setup_b']	
		shift = np.fromstring(child[0].text, sep=' ')
		bounding_box = np.fromstring(child[3].text, sep=' ')
		dict_data = {"setup number a":setup_a_number,"setup number b":setup_b_number,"shift":shift,"bounding box":bounding_box}
		stitchings.append(dict_data)

	return files, registrations, stitchings

def get_stitching_matrices(files,registrations,volume):

	# Get associated setup id
	setup_id = -1
	for f in files:
		if volume in f['file name']:
			setup_id = f['setup number']
			break

	# Make sure volume was found in XML
	if setup_id == -1:
		print("Error: Volume name not found in XML")
		return -1

	# Get associated stitching parameters
	stitch_params = -1
	for r in registrations:
		if r['setup number'] == setup_id:
			stitch_params = r
			break

	# Make sure volume was found in XML
	if stitch_params == -1:
		print("Error: Volume stitching parameters name not found in XML")
		return -1

	# Extract affine transformation matrices and add fourth row (matrices are 3x4)
	translation_to_grid_matrix = np.vstack([np.reshape(stitch_params['translation to regular grid'], (3,4)),np.array([0,0,0,1])])
	stitching_matrix = np.vstack([np.reshape(stitch_params['stitching transform'], (3,4)),np.array([0,0,0,1])])
	calibration_matrix = np.vstack([np.reshape(stitch_params['calibration'], (3,4)),np.array([0,0,0,1])])

	return translation_to_grid_matrix, stitching_matrix, calibration_matrix
