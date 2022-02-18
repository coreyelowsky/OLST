import numpy as np
import params
import utils
import tifffile as tif

from StitchingXML import StitchingXML
from sys import exit
from os.path import join

class Centroids:

	def __init__(self, inpath, brain):
		
		self.inpath = inpath
		self.brain = brain

		# make sure brain is in inpath
		if self.brain not in self.inpath:
			exit('Error: brain must be in inpath')

		# load centroids
		self.load_centroids()



	def load_centroids(self):

		"""
		load centroids from csv
	
		"""
		
		if self.inpath.endswith('.csv'):
			centroids_array = np.genfromtxt(self.inpath, delimiter=',', dtype=object)[:,3:7]
			centroids_array[:,0] = centroids_array[:,0].astype(str)
			centroids_array[:,1:] = centroids_array[:,1:].astype(int)
			self.centroids_array = centroids_array
		else:
			exit('Centroids file must be of type csv')

	def register(self):

		"""
		oblique to coronal
		
		"""

		# load stitching xml
		xml_path = join(params.BASE_PATH,self.brain,'stitching_parameters.xml')
		stitching_xml = StitchingXML(xml_path, verbose=False)
		fused_dims = stitching_xml.getFusedDimensions()

		# load cropping infor
		crop_dict = utils.loadCropInfo()

		coords_list = []
	
		for centroid in self.centroids_array:

			# extract volume and centroid
			volume = centroid[0]
			centroid_coords = centroid[1:]

			# get stitching matrices 
			stitching_matrices = stitching_xml.getStitchingMatrices(volume)
			
	
			# apply stitching matrices
			centroid_coords = np.concatenate((centroid_coords, [1]))
			transformed_coords = np.linalg.inv(stitching_matrices['calibration']) @ stitching_matrices['Stitching Transform'] @ stitching_matrices['Translation to Regular Grid'] @ stitching_matrices['calibration'] @ centroid_coords.reshape(-1,1)
			transformed_coords = transformed_coords.flatten()[:-1]
			
			# transform from oblique to coronal
			coronal_coords = utils.obliqueToCoronal(transformed_coords, fused_dims, params.SHEAR_FACTOR_ANISOTROPIC)
		

			# apply cropping shift
			shift = crop_dict[self.brain]
			coronal_coords[1] -= shift*params.CROP_FACTOR

			coords_list.append(coronal_coords)

			
		coords_list = np.array(coords_list)
		
		
		brain_path = join(params.REGISTRATION_BASE_PATH, self.brain)

		# write coords to file
		outPath = join(brain_path,'swc_in.txt')
		utils.writeCoordsForRegistration(outPath, coords_list)

		# run transformix
		tpPath = join(brain_path,'TransformParameters_labels.1_full_scale.txt')
		inPath = join(brain_path,'swc_in.txt')
		utils.runTransformix(brain_path, tpPath, inPath)

		# read in registered coordinates
		transformedCoordsPath = join(brain_path,'outputpoints.txt')
		registeredCoords = utils.parseTransformixOutput(transformedCoordsPath)

		registeredCoords = registeredCoords.round().astype(int)

		# create tif image
		ccf = np.zeros(shape=(528,320,456), dtype=np.uint8)
		for coord in registeredCoords:
			ccf[coord[2], coord[1],coord[0]] = 1

		tif.imsave(join(params.BASE_PATH, self.brain, 'somas_ccf_25_um.tif'), ccf)
		







if __name__ == '__main__':
	

	for brain in params.BRAINS:

		if brain == '190306':
			continue
		inpath = join(params.BASE_PATH, brain, 'inside_brain_somas.csv')

		c = Centroids(inpath, brain)
		c.register()



