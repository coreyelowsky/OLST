from sys import exit
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import utils
from os.path import dirname, join


"""
Stitching XML 
	- provides functions to parse Stitching XML file from Big Stitcher
"""


class StitchingXML:

	def __init__(self, xmlPath, * , numZVolumes=None, numYVolumes=None, verbose=True):
		
		self.xmlPath = xmlPath
		self.tree = ET.parse(xmlPath)
		self.root =  self.tree.getroot()

		# parse XML for setups, registrations, and stitchings
		self.setupsAndRegistrations = self.parseSetups()
		self.parseRegistrations()
		self.pairwiseStitchings = self.parsePairwiseStitchings()

		# calculate volume dimension info
		if numYVolumes == None or numZVolumes == None:
			self.numYVolumes, self.numZVolumes = self.inferYZVolumes()
		else:
			self.numYVolumes = numYVolumes
			self.numZVolumes = numZVolumes
		self.numVolumes = self.numYVolumes*self.numZVolumes

		# get volume sizes
		self.volumeSize = self.getVolumeSize()

		# get resolution
		self.resolution = self.getResolution()

		# get anisotropy factor
		self.anisotropyFactor = self.getAnisotropyFactor()
		
		if verbose:
			print(self)



	def getAnisotropyFactor(self):

		"""
		gets anisotropy factor (Z resolution / XY Resolution) 
		"""

		voxelSize = self.setupsAndRegistrations[0]['voxel size']
		anisotropyFactor =  voxelSize[2]/voxelSize[1]

		# check if all anisotropy factors are the same
		for _, setupDict in self.setupsAndRegistrations.items():
			
			voxelSizeLocal = setupDict['voxel size']
			anisotropyFactorLocal = voxelSizeLocal[2]/voxelSizeLocal[1]
	
			if anisotropyFactor != anisotropyFactorLocal:
				print()
				print('###############################################')
				print('WARNING: All anisotropy factors are not the same!!!')
				print('###############################################')
				print()
			
		return anisotropyFactor

	def getResolution(self):

		"""
		gets resolution
		"""

		resolution = self.setupsAndRegistrations[0]['voxel size']

		# check if all anisotropy factors are the same
		for _, setupDict in self.setupsAndRegistrations.items():
			
			resolutionLocal = setupDict['voxel size']
	
			if resolution != resolutionLocal:
				print()
				print('###############################################')
				print('WARNING: Resolutions are not the same!!!')
				print('###############################################')
				print()
			
		return resolution

	def getVolumeSize(self):

		"""
		gets volumes size

		"""

		volumeSize = self.setupsAndRegistrations[0]['size']

		# check if all volume sizes are the same
		for _, setupDict in self.setupsAndRegistrations.items():
			
			volumeSizeLocal = setupDict['size']

			if volumeSize != volumeSizeLocal:
				print()
				print('###############################################')
				print('WARNING: All Volume Sizes are not the same!!!')
				print('###############################################')
				print()
			
		return volumeSize


	def getAdjacencyType(self, setupIDa, setupIDb):

		"""
		Returns adjacency type (y, z, yz, other) of two volumes
		"""

		zA = self.setupIDtoZ(setupIDa)
		yA = self.setupIDtoY(setupIDa)
	
		zB = self.setupIDtoZ(setupIDb)
		yB = self.setupIDtoY(setupIDb)
		
		if zA == zB and abs(yA - yB) == 1:
			return 'y'
		elif yA == yB and abs(zA - zB) == 1:
			return 'z'
		elif abs(yA - yB) == 1 and abs(zA - zB) == 1:
			return 'yz'
		else:
			return 'other'


	def setupIDtoY(self, setupID):

		"""
		Converts setup ID to Y
		"""

		return (setupID % self.numYVolumes) + 1

	def setupIDtoZ(self, setupID):

		"""
		Converts setup ID to Z
		"""

		return (setupID // self.numYVolumes) + 1

	def setupIDtoVolume(self, setupID):

		"""
		Converts setup ID to Volume (Z12_Y08)
		"""

		y = '{:02}'.format(self.setupIDtoY(setupID))
		z = '{:02}'.format(self.setupIDtoZ(setupID))

		return 'Z' + z + '_Y' + y

	def volumeToSetupID(self, volume):

		"""
		Converts Volume (Z12_Y08) to setup ID
		"""

		z = int(volume.split('_')[0][1:])
		y = int(volume.split('_')[1][1:])

		return (z-1)*self.numYVolumes + (y-1)


	def parseSetups(self):

		"""
		Parses setup information into dictionary where setup ids are keys
		"""

		setups = {}

		for setup in self.root.iter('ViewSetup'):

			setupID = int(setup.find('id').text)
			size = setup.find('size').text.split(' ')
			voxelSize = setup.find('voxelSize').find('size').text.split(' ')
			
			setups[setupID] = {'size':[int(x) for x in size], 'voxel size':[float(x) for x in voxelSize]}

		return setups

	def parseRegistrations(self):

		"""	
		Parses registration information which are three matrices
		and stores them in setupsAndRegistrations dictionary 
		"""
		
		for registration in self.root.iter('ViewRegistration'):

			setupID = int(registration.attrib['setup'])

			for transform in registration.iter('ViewTransform'):

				name = transform.find('Name').text
				matrix = np.array([float(x) for x in transform.find('affine').text.split()]).reshape(3,4)
				self.setupsAndRegistrations[setupID][name] = np.vstack((matrix,np.array([0,0,0,1])))

	def parsePairwiseStitchings(self):
		
		"""
		Parses pairwise stitchings and stores in list of dictionaries
		"""


		pairwiseStitchings = []

		for pairwise in self.root.iter('PairwiseResult'):

			setupID_a = int(pairwise.attrib['view_setup_a'])
			setupID_b = int(pairwise.attrib['view_setup_b'])
			shift = np.array([float(x) for x in pairwise.find('shift').text.split()]).reshape(3,4)
			shift = np.vstack((shift,np.array([0,0,0,1])))
			correlation = np.round(float(pairwise.find('correlation').text),2)
			bbox = [float(x) for x in pairwise.find('overlap_boundingbox').text.split()]

			pairwiseStitchings.append({'setupID_a':setupID_a,'setupID_b':setupID_b,'shift':shift,'correlation':correlation,'bbox':bbox})
			

		return pairwiseStitchings

	
	def inferYZVolumes(self):
		
		"""
		Calculates number of Y and Z volumes based on translsation to grid matrices
		"""

		x1 = self.getTranslationToGridMatrix(0, 'x')

		for setupID in sorted(self.setupsAndRegistrations.keys()):
	
			if setupID > 0:			
				x2 = self.getTranslationToGridMatrix(setupID, 'x')

				# handle 1 dimension case in Y
				if setupID == len(self.setupsAndRegistrations)-1:
					numYVolumes = setupID + 1

				# reached next Z
				if x1 == x2:
					numYVolumes = setupID
					break

		numZVolumes = int(len(self.setupsAndRegistrations)/numYVolumes)

		if numZVolumes*numYVolumes != len(self.setupsAndRegistrations):
			print('Warning: Number of y and z volumes might not be correct infered correctly!')	

		return numYVolumes, numZVolumes


	def viewAllCorrelations(self, * , zPair=None, yPair=None, verbose=True):
		
	
		"""
		Provides breakdown of correlations in sorted order

		Parmams: zPair, set that contains pair of Zs to get correlations across, (e.g. {4,5})
			 yPair, set that contains pair of Ys to get correlations across, (e.g. {4,5})

		"""

		s = '\n\n'
		s += '#############\n'
		s += 'Correlations\n'
		s += '#############\n'

		yCorrs, zCorrs, yzCorrs, otherCorrs = [], [], [], []

		for pairwise in self.pairwiseStitchings:

			setupID_a = pairwise['setupID_a']
			setupID_b = pairwise['setupID_b']

			if zPair:
				if type(zPair) != set:
					exit('Error: zPair must be of type set!')
				if {self.setupIDtoZ(setupID_a), self.setupIDtoZ(setupID_b)} != zPair:
					continue

			if yPair:
				if type(yPair) != set:
					exit('Error: zPair must be of type set!')
				if {self.setupIDtoY(setupID_a), self.setupIDtoY(setupID_b)} != zPair:
					continue
			
			correlation = pairwise['correlation']

			if self.getAdjacencyType(setupID_a,setupID_b) == 'y':
				yCorrs.append(correlation)
			elif self.getAdjacencyType(setupID_a,setupID_b) == 'z':
				zCorrs.append(correlation)
			elif self.getAdjacencyType(setupID_a,setupID_b) == 'yz':
				yzCorrs.append(correlation)
			elif self.getAdjacencyType(setupID_a,setupID_b) == 'other':
				otherCorrs.append((correlation,self.setupIDtoVolume(setupID_a),self.setupIDtoVolume(setupID_b)))

		yCorrs.sort()
		zCorrs.sort()
		yzCorrs.sort()
		otherCorrs.sort()

		yCorrs = yCorrs[::-1]
		zCorrs = zCorrs[::-1]
		yzCorrs = yzCorrs[::-1]
		otherCorrs = otherCorrs[::-1]
		

		s += 'Y Corrs: ' + str(yCorrs) + '\n\n'
		s += 'Z Corrs: ' + str(zCorrs) + '\n\n'
		s += 'YZ Corrs: ' + str(yzCorrs) + '\n\n'
		s += 'Other Corrs: ' + str(otherCorrs) + '\n\n'
		s += '# Y: ' + str(len(yCorrs)) + '\n'
		s += '# Z: ' + str(len(zCorrs)) + '\n'
		s += '# YZ: ' + str(len(yzCorrs)) + '\n'
		s += '# Other: ' + str(len(otherCorrs)) + '\n\n'

		if verbose:
			print(s)

		return s


	def getStitchingMatrices(self, volume):

		"""
		Gets transformation matrices for given volume

		Params: volume, e.g. Z12_Y06
		Return: translation to grid matrix, stitching matrix, calibration matrix stored in dictionary 

		"""

		stitchingDict = self.setupsAndRegistrations[self.volumeToSetupID(volume)]

		stitchingMatrices = {}
		stitchingMatrices['Translation to Regular Grid'] = stitchingDict['Translation to Regular Grid']
		stitchingMatrices['Stitching Transform'] = stitchingDict['Stitching Transform']
		stitchingMatrices['calibration'] = stitchingDict['calibration']

		return stitchingMatrices

	def getFusedDimensions(self, * , downsampling=1):
		
		"""
		Calculate dimensions of big stitcher image after fusing

		Params: downamspling, int that determines how much image will be downsampled
		Return: dictionary containing dimension information

		"""

		xCoords, yCoords, zCoords = [], [], []

		# iterate through volumes
		for setupID in self.setupsAndRegistrations:

			volumeSize = self.setupsAndRegistrations[setupID]['size']

			volume = self.setupIDtoVolume(setupID)
			stitchingMatrices = self.getStitchingMatrices(volume)
			translationToGridMatrix = stitchingMatrices['Translation to Regular Grid']
			stitchingMatrix = stitchingMatrices['Stitching Transform']

			# append min and max bounds
			xCoords.append(translationToGridMatrix[0,3] + stitchingMatrix[0,3])
			xCoords.append(translationToGridMatrix[0,3] + stitchingMatrix[0,3] + volumeSize[0])
			yCoords.append(translationToGridMatrix[1,3] + stitchingMatrix[1,3])
			yCoords.append(translationToGridMatrix[1,3] + stitchingMatrix[1,3] + volumeSize[1])
			zCoords.append((translationToGridMatrix[2,3] + stitchingMatrix[2,3])/self.anisotropyFactor)
			zCoords.append((translationToGridMatrix[2,3] + stitchingMatrix[2,3])/self.anisotropyFactor + volumeSize[2])

		minX, maxX = min(xCoords)/downsampling, max(xCoords)/downsampling
		minY, maxY = min(yCoords)/downsampling, max(yCoords)/downsampling
		minZ, maxZ = min(zCoords)/downsampling, max(zCoords)/downsampling
		lengthX = int(np.round(maxX - minX))
		lengthY = int(np.round(maxY - minY))
		lengthZ = int(np.round(maxZ - minZ))

		dimensions = {'anisotropy factor':self.anisotropyFactor ,'x':{'min':minX,'max':maxX,'length':lengthX},'y':{'min':minY,'max':maxY,'length':lengthY},'z':{'min':minZ,'max':maxZ,'length':lengthZ}}
		
		return dimensions

	def printCorrelationGrid_Z(self, verbose=True):

		"""
		Prints grid of Z correlations

		"""

		s = '\n\n'
		s += '######################\n'
		s += ' Z Correlations\n'
		s += '######################\n\n'


		for z in range(1,self.numZVolumes):
			corrList = []
			for y in range(1,self.numYVolumes+1):

				volumeID = utils.zyToVolumeID(z,y)
				setupID = self.volumeToSetupID(volumeID)
				setupIDBelow = setupID + self.numYVolumes

				pairwiseStitching = self.getPairwiseStitching(setupID,setupIDBelow)

				if pairwiseStitching == -1:
					corrList.append(' . ')
				else:
					if pairwiseStitching['correlation'] == 1.0:
						corrList.append('1.0')
					else:
						corrList.append('{:.2f}'.format(abs(pairwiseStitching['correlation']))[1:])
			corrList.append('('+str(z)+') ')	
			corrList = corrList[::-1]
			s += ' '.join(corrList) + '\n'

		if verbose:
			print(s)

		return s

	def printCorrelationGrid_ZY(self, verbose=True):

		"""
		Prints grid of ZY correlations

		"""


		s = '\n'
		s += '######################\n'
		s += ' ZY Correlations\n'
		s += '######################\n\n'


		for z in range(1,self.numZVolumes):
			corrList = []
			for y in range(1,self.numYVolumes+1):

				volumeID = utils.zyToVolumeID(z,y)
				setupID = self.volumeToSetupID(volumeID)
				setupIDBelowRight = setupID + self.numYVolumes - 1
				setupIDBelowLeft = setupID + self.numYVolumes + 1

		
				if self.setupIDtoY(setupID) == 1:
					pairwiseStitchingBelowLeft = self.getPairwiseStitching(setupID,setupIDBelowLeft)
					pairwiseStitchingBelowRight = -1

				elif self.setupIDtoY(setupID) == self.numYVolumes:
					pairwiseStitchingBelowLeft = -1
					pairwiseStitchingBelowRight = self.getPairwiseStitching(setupID,setupIDBelowRight)
				else:
					pairwiseStitchingBelowLeft = self.getPairwiseStitching(setupID,setupIDBelowLeft)
					pairwiseStitchingBelowRight = self.getPairwiseStitching(setupID,setupIDBelowRight)
			
				if pairwiseStitchingBelowRight == -1:
					corrList.append(' . ')
				else:
					if pairwiseStitchingBelowRight['correlation'] == 1.0:
						corrList.append('1.0')
					else:
						corrList.append('{:.2f}'.format(abs(pairwiseStitchingBelowRight['correlation']))[1:])

				if pairwiseStitchingBelowLeft == -1:
					corrList.append(' . ')
				else:
					if pairwiseStitchingBelowLeft['correlation'] == 1.0:
						corrList.append('1.0')
					else:
						corrList.append('{:.2f}'.format(abs(pairwiseStitchingBelowLeft['correlation']))[1:])
		
			corrList.append('('+str(z)+') ')
			corrList = corrList[::-1]
			s += ' '.join(corrList) + '\n'

		if verbose:
			print(s)

		return s

	def printCorrelationGrid_Y(self, verbose=True):
		
		"""
		Prints grid of Y correlations

		"""

		s = '\n'
		s += '######################\n'
		s += ' Y Correlations\n'
		s += '######################\n\n'

		for z in range(1,self.numZVolumes+1):
			corrList = []
			for y in range(1,self.numYVolumes):

				volumeID = utils.zyToVolumeID(z,y)
				setupID = self.volumeToSetupID(volumeID)
				setupIDLeft = setupID + 1

				pairwiseStitching = self.getPairwiseStitching(setupID,setupIDLeft)

				if pairwiseStitching == -1:
					corrList.append(' . ')
				else:

					if pairwiseStitching['correlation'] == 1.0:
						corrList.append('1.0')
					else:
						corrList.append('{:.2f}'.format(abs(pairwiseStitching['correlation']))[1:])
			corrList.append('('+str(z)+') ')	
			corrList = corrList[::-1]
			s += ' '.join(corrList) + '\n'
		
		if verbose:
			print(s)
		
		return s
	
	def getPairwiseStitching(self,setupIDa, setupIDb):

		"""
		Returns pairwise stitching dictionary for pair of volumes
		If pairwise stitching doesnt exist, reutrns -1

		"""		

		for pairwiseStitching in self.pairwiseStitchings:
			
			sA = pairwiseStitching['setupID_a']
			sB = pairwiseStitching['setupID_b']
			
		
			if {sA, sB} == {setupIDa, setupIDb}:
				return pairwiseStitching


		return -1

	def printOverlapGrid_Z(self, verbose=True):

		"""
		Prints grid of  Z overlap percentages 

		"""

		s = '\n\n'
		s += '######################\n'
		s += 'Z Overlap Percentages\n'
		s += '######################\n\n'

		for z in range(1,self.numZVolumes):
			overlapList = []
			for y in range(1,self.numYVolumes+1):


				volumeID = utils.zyToVolumeID(z,y)
				setupID = self.volumeToSetupID(volumeID)
				setupIDBelow = setupID + self.numYVolumes

				volAbove = self.setupsAndRegistrations[setupID]
				volBelow = self.setupsAndRegistrations[setupIDBelow]

				if 'Stitching Transform' not in volAbove:
					topZ = self.getTranslationToGridMatrix(setupIDBelow, 'y')
					bottomZ =  self.getTranslationToGridMatrix(setupID, 'y') + volAbove['size'][1]
				else:	
					topZ = self.getTranslationToGridMatrix(setupIDBelow, 'y') + self.getStitchingMatrix(setupIDBelow, 'y')
					bottomZ = self.getTranslationToGridMatrix(setupID, 'y') + self.getStitchingMatrix(setupID, 'y') + volAbove['size'][1]

			
				zOverlap = np.round(100*(bottomZ-topZ) / volAbove['size'][1],1)
				overlapList.append('{:02}'.format(zOverlap))

			overlapList.append('('+str(z)+') ')
			overlapList = overlapList[::-1]
			s += ' '.join(overlapList) + '\n'

		if verbose:
			print(s)

		return s	

	def printOverlapGrid_Y(self, verbose=True):

		"""
		Prints grid of Y overlap percentages 

		"""

		s = '\n\n'
		s += '######################\n'
		s += 'Y Overlap Percentages\n'
		s += '######################\n\n'
	
		for z in range(1,self.numZVolumes+1):
			overlapList = []
			
			for y in range(1,self.numYVolumes):

				volumeID = utils.zyToVolumeID(z,y)
				setupID = self.volumeToSetupID(volumeID)
				setupIDLeft = setupID + 1

				volRight = self.setupsAndRegistrations[setupID]
				volLeft = self.setupsAndRegistrations[setupIDLeft]

				if 'Stitching Transform' not in volRight:
					rightY = volLeft['Translation to Regular Grid'][0,-1] +  volLeft['size'][0]
					leftY = volRight['Translation to Regular Grid'][0,-1] 
				else:
					rightY = volLeft['Translation to Regular Grid'][0,-1] + volLeft['Stitching Transform'][0,-1] + volLeft['size'][0]
					leftY = volRight['Translation to Regular Grid'][0,-1] + volRight['Stitching Transform'][0,-1] 

				yOverlap = np.round(100*(rightY-leftY) / volLeft['size'][0],1)

				overlapList.append('{:02}'.format(yOverlap))
					
			overlapList.append('('+str(z)+') ')
			overlapList = overlapList[::-1]
			s += ' '.join(overlapList) + '\n'

		if verbose:
			print(s)
	
		return s

	def printOverlapGrid_X(self, verbose=True):

		"""
		Prints grid of  X overlap percentages

		"""

		s = '\n\n'
		s += '######################\n'
		s += 'X Overlap Percentages\n'
		s += '######################\n\n'
	
		for z in range(1,self.numZVolumes):
			overlapList = []
			for y in range(1,self.numYVolumes+1):

				volumeID = utils.zyToVolumeID(z,y)
				setupID = self.volumeToSetupID(volumeID)
				setupIDOut = setupID + self.numYVolumes

				volIn = self.setupsAndRegistrations[setupID]
				volOut = self.setupsAndRegistrations[setupIDOut]

				if 'Stitching Transform' not in volIn:
					inX = volIn['Translation to Regular Grid'][2,-1]
					outX = volOut['Translation to Regular Grid'][2,-1] + volOut['size'][2]*self.anisotropyFactor
				else:	
					inX = volIn['Translation to Regular Grid'][2,-1] + volIn['Stitching Transform'][2,-1]
					outX = volOut['Translation to Regular Grid'][2,-1] + volOut['Stitching Transform'][2,-1] + volOut['size'][2]*self.anisotropyFactor

			
				xOverlap = np.round(100*(outX-inX) / (volIn['size'][2]*self.anisotropyFactor),1)
				overlapList.append('{:02}'.format(xOverlap))

			overlapList.append('('+str(z)+') ')
			overlapList = overlapList[::-1]
			s += ' '.join(overlapList) + '\n'

		if verbose:
			print(s)

		return s	

	def removeNonAdjacentCorrelationsAndSaveXML(self):
		
		# iterate through pairwise stitchings and collect all non adjacent pairwise stitchings
		stitchingsToRemove = []
		for pairwise in self.root.iter('PairwiseResult'):

			setupID_a = int(pairwise.attrib['view_setup_a'])
			setupID_b = int(pairwise.attrib['view_setup_b'])

			if self.getAdjacencyType(setupID_a, setupID_b) == 'other':
				stitchingsToRemove.append(pairwise)

		# remove all non adjacent pairwise stitchings from xml
		for pairwise in stitchingsToRemove:
			self.root.find('StitchingResults').remove(pairwise)

		print('# Non Adjacent Pairwise Stitchings Removed: ', len(stitchingsToRemove))

		# save xml
		self.writeXML('removed_non_adjacent.xml')

	def modifyImageLoaderForSavingAsN5Parallelized(self):
		
		# get image loader element
		imageLoader = self.root.find('SequenceDescription').find('ImageLoader')

		# set attributes
		imageLoader.set('format','bdv.n5')
		imageLoader.set('version','1.0')

		# remove all children elements
		children = []
		for child in imageLoader:
			children.append(child)
		for child in children:
			imageLoader.remove(child)

		# add n5 child element
		n5Child = ET.SubElement(imageLoader,'n5')
		n5Child.set('type','relative')
		n5Child.text = 'dataset.n5'

		# save xml
		self.writeXML('translate_to_grid.xml')
	

	def analyze(self):
	
		self.viewAllCorrelations()
		self.printCorrelationGrid_Y()
		self.printCorrelationGrid_Z()
		self.printCorrelationGrid_ZY()
		self.printOverlapGrid_Y()
		self.printOverlapGrid_Z()
		self.printOverlapGrid_X()
		
		print()

	def generateReport(self):

		outPath = self.xmlPath[:-4] + '.txt'
		with open(outPath, 'w') as fp:
			fp.write(self.__str__())
			fp.write(self.viewAllCorrelations(verbose=False))
			fp.write(self.printCorrelationGrid_Y(verbose=False))
			fp.write(self.printCorrelationGrid_Z(verbose=False))
			fp.write(self.printCorrelationGrid_ZY(verbose=False))
			fp.write(self.printOverlapGrid_Y(verbose=False))
			fp.write(self.printOverlapGrid_Z(verbose=False))
			fp.write(self.printOverlapGrid_X(verbose=False))


	def setStitchingMatricesToIdentity(self):

		print()
		print('Setting Stitching Matrices to Identity...')
		print()

		setupID = 0
		while setupID < len(self.setupsAndRegistrations):

			setupsDict = self.setupsAndRegistrations[setupID]

			setupsDict['Stitching Transform'] = np.identity(4)

			setupID += 1

		


	def setTranslationToGrid(self, *,  yOverlapPercentage, zOverlapPercentage, zOverlapPercentageSectioning, zSectioningInterval, xOverlapPercentage, verbose=True):
		
		# y is left and right
		# z is up and down
		# x is into and out of screen

		if verbose:
			print()
			print('Set Translation to Grid...')
			print()
			print('New Y Overlap:',yOverlapPercentage,'%')
			print('New Z Overlap:',zOverlapPercentage,'%')
			print('New Z Overlap (Sectioning):',zOverlapPercentageSectioning,'%')
			print('Z Sectioning Interval:',zSectioningInterval)
			print('New X Overlap:',xOverlapPercentage,'%')
			print()

		setupID = 0
		while setupID < len(self.setupsAndRegistrations):

			# get y and z
			y = self.setupIDtoY(setupID)
			z = self.setupIDtoZ(setupID)
			volume = self.setupIDtoVolume(setupID)
	
			# calculate new y coordinate
			# fix all Y01's
			if y != 1:

				rightDict = self.setupsAndRegistrations[setupID-1]
				leftDict = self.setupsAndRegistrations[setupID]

				yRightSize = rightDict['size'][0]
				yOverlapPixels = yRightSize*yOverlapPercentage/100
				yRightLeftBorder = rightDict['Translation to Regular Grid'][0,-1]
				newYLeftBorder = yRightLeftBorder + yOverlapPixels - yRightSize
	
				leftDict['Translation to Regular Grid'][0,-1] = newYLeftBorder

			# calculate new x and z coordinate
			# fix all Z01
			if z!= 1:

				# z calculations
				aboveDict = self.setupsAndRegistrations[setupID-self.numYVolumes]
				belowDict = self.setupsAndRegistrations[setupID]

				zAboveSize = rightDict['size'][1]

				if (z-1) % zSectioningInterval == 0:
					zOverlapPixels = zAboveSize*zOverlapPercentageSectioning/100
				else:
					zOverlapPixels = zAboveSize*zOverlapPercentage/100
				zAboveBottomBorder = aboveDict['Translation to Regular Grid'][1,-1] + zAboveSize
				newZ = zAboveBottomBorder - zOverlapPixels
	
				belowDict['Translation to Regular Grid'][1,-1] = newZ

				# x calulations
				outDict = self.setupsAndRegistrations[setupID-self.numYVolumes]
				inDict = self.setupsAndRegistrations[setupID]

				xOutSize = outDict['size'][2]*self.anisotropyFactor
				xOverlapPixels = xOutSize*xOverlapPercentage/100
				xOut = outDict['Translation to Regular Grid'][2,-1]
				newX = xOut + xOverlapPixels - xOutSize
	
				inDict['Translation to Regular Grid'][2,-1] = newX

			setupID += 1

			
	def writeXML(self, xmlName, writeStitchingTransform=False):

		# writes xml file with updated registrations matrices

		# set registrations

		for registration in self.root.iter('ViewRegistration'):
			
			setupID = int(registration.attrib['setup'])
			
			# add Stitching Transform if not there
			if writeStitchingTransform:
				if not self.isStitchingTransformInViewRegistration(registration):
					viewTransformStitching = ET.SubElement(registration,'ViewTransform')
					viewTransformStitching.set('type','affine')
					nameChild = ET.SubElement(viewTransformStitching,'Name')
					nameChild.text = 'Stitching Transform'
					matrixChild = ET.SubElement(viewTransformStitching,'affine')
					matrixChild.text = self.matrixToString(np.identity(4))

			for transform in registration.iter('ViewTransform'):
				
				name = transform.find('Name').text

				if name == 'Translation to Regular Grid':
					translationMatrixNode = transform.find('affine')
					translationMatrix = self.setupsAndRegistrations[setupID]['Translation to Regular Grid']
					translationMatrixNode.text = self.matrixToString(translationMatrix)

				if writeStitchingTransform:
					if name == 'Stitching Transform':
						stitchingMatrixNode = transform.find('affine')
						stitchingMatrix = self.setupsAndRegistrations[setupID]['Stitching Transform']
						stitchingMatrixNode.text = self.matrixToString(stitchingMatrix)


		outPath = join(dirname(self.xmlPath),xmlName)

		print('Writing XML:', outPath)
		print()
		self.tree.write(outPath)


	def isStitchingTransformInViewRegistration(self,viewRegistrationNode):
		
		for transform in viewRegistrationNode.iter('ViewTransform'):
			
			name = transform.find('Name').text
			
			if name == 'Stitching Transform':
				return True
		
		return False

	def matrixToString(self, matrix):
		
		# converts 2d numpy array into string that can be placed in xml

		return ' '.join(matrix[:-1].flatten().astype(str))

	def getCoordsAdjacentVolume(self, sourceVolume, imageCoordsSource, targetVolume):
		
		print()
		print('Source Volume:', sourceVolume)
		print('Target Volume:', targetVolume)
		print()
		print('Image Coords (Source):', imageCoordsSource)

		# get setup IDs
		setupIDSource = self.volumeToSetupID(sourceVolume)
		setupIDTarget = self.volumeToSetupID(targetVolume)

		# get pairwise shift
		pairwiseStitching = self.getPairwiseStitching(setupIDSource, setupIDTarget)

		shiftMatrix = pairwiseStitching['shift']
		bbox = pairwiseStitching['bbox']
		correlation = pairwiseStitching['correlation']

		# shift describes how setup b should move
		viewSetupA = pairwiseStitching['setupID_a']
		viewSetupB = pairwiseStitching['setupID_b']

		# invert shift if volume are swapped	
		# since shift is how to go from view_setup b to view setup a
		# so source volume must be view setup b
		if setupIDSource != viewSetupB:
			shiftMatrix = np.linalg.inv(shiftMatrix)

		calibrationMatrix = self.getCalibrationMatrix(setupIDSource)
		translationMatrixSource = self.getTranslationToGridMatrix(setupIDSource)
		translationMatrixTarget = self.getTranslationToGridMatrix(setupIDTarget)
		
		# transform coords
		imageCoordsSource = np.append(np.array(list(imageCoordsSource)),1).reshape(-1,1)
		gridCoordsSource = translationMatrixSource @ calibrationMatrix @ imageCoordsSource
		gridCoordsShifted = shiftMatrix @ gridCoordsSource
		imageCoordsTarget = np.linalg.inv(calibrationMatrix) @ np.linalg.inv(translationMatrixTarget) @ gridCoordsShifted

		imageCoordsTarget = np.round(imageCoordsTarget.flatten()[:-1])
		print('Image Coords (Target) :',imageCoordsTarget)	
		print()
		print('Correlation:', correlation)

		return imageCoordsTarget
		
	
	def getStitchingMatrix(self, setupID, dim=None):

		if dim == 'x':
			return self.setupsAndRegistrations[setupID]['Stitching Transform'][0,-1]
		elif dim == 'y':
			return self.setupsAndRegistrations[setupID]['Stitching Transform'][1,-1]
		elif dim == 'z':
			return self.setupsAndRegistrations[setupID]['Stitching Transform'][2,-1]
		else:	
			return self.setupsAndRegistrations[setupID]['Stitching Transform']

		
	def getTranslationToGridMatrix(self, setupID, dim=None):

		if dim == 'x':
			return self.setupsAndRegistrations[setupID]['Translation to Regular Grid'][0,-1]
		elif dim == 'y':
			return self.setupsAndRegistrations[setupID]['Translation to Regular Grid'][1,-1]
		elif dim == 'z':
			return self.setupsAndRegistrations[setupID]['Translation to Regular Grid'][2,-1]
		else:
			return self.setupsAndRegistrations[setupID]['Translation to Regular Grid']
	
	def getCalibrationMatrix(self, setupID):
	
		return self.setupsAndRegistrations[setupID]['calibration']		
		

	def __str__(self):

		s = '\n###############\n'
		s += 'Stitching Info\n'
		s += '###############\n\n'
		s += 'XML Path: ' + self.xmlPath + '\n'
		s += '# Volumes: ' + str(self.numVolumes) + '\n'
		s += '# Y: ' + str(self.numYVolumes) + '\n'
		s += '# Z: ' + str(self.numZVolumes) + '\n'
		s += 'Resolution: ' + str(self.resolution) + '\n'
		s += 'Volume Size: ' + str(self.volumeSize) + '\n'
		s += 'Anisotropy Factor: ' + str(self.anisotropyFactor) + '\n'
		s += '# Correlations: ' + str(len(self.pairwiseStitchings)) + '\n'
		
		return s
		 


if __name__ == '__main__':

	xml = StitchingXML('/data/elowsky/OLSTv2/data/VIP-GFP-M2/pairwise_shift_tests/pairwise_shifts.xml')
	#xml.analyze()

	volumeA = 'Z06_Y11'
	imageCoordsA = [415, 1559, 1710]
	volumeB = 'Z07_Y11'
	xml.getCoordsAdjacentVolume(volumeA, imageCoordsA, volumeB)
	



	#xml = StitchingXML('/mnt/nfs/grids/hpc_norepl/elowsky/PV-GFP-M2_stitching_output/no_stitching/translate_to_grid.xml')
	#xml.setTranslationToGrid(yOverlapPercentage=yOverlapPercentage, zOverlapPercentage=zOverlapPercentage, zOverlapPercentageSectioning=zOverlapPercentageSectioning, zSectioningInterval=zSectioningInterval, xOverlapPercentage=xOverlapPercentage)
	#xml.setStitchingMatricesToIdentity()
	#xml.writeXML('modified_translation_to_grid_for_no_stitching.xml', writeStitchingTransform=True)
	#xml.analyze()


	#xml = StitchingXML('/mnt/nfs/grids/hpc_norepl/elowsky/SST-GFP_stitching_output/xmls_50%_overlap/modified_transaltion_to_grid_for_no_stitching.xml')
	#xml.analyze()



	



	


