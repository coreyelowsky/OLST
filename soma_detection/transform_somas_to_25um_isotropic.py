import numpy as np
import xml.etree.ElementTree as ET
import os
import sys
from utils import *
from BRAINS import BRAINS,BRAINS_4000,BRAINS_3600

##################################################
BASE_PATH = '/data/elowsky/OLST/soma_detection/'
##################################################

#####  Initialize Parameters #########
DOWNSIZE_FACTOR = 10
SHEAR_FACTOR = -.7
XY_RES = .406
Z_RES = 2.5
ANISOTROPY_FACTOR = Z_RES/XY_RES
X_SIZE = 2048
Y_SIZE = 1024
######################################

for brain in BRAINS:

	brain_path = os.path.join(BASE_PATH,brain)
	xml_file = os.path.join(brain_path,'stitching_parameters.xml')                      
	soma_path = os.path.join(brain_path,'triaged_somas_close_200.npy')              

	
	if not os.path.exists(soma_path):
		sys.exit('Soma Path Does Not Exist!!')

	print('Brain:',brain)


	# Lookup Z Size
	if brain in BRAINS_4000:
		Z_SIZE = 4000
	elif brain in BRAINS_3600:
		Z_SIZE = 3600
	else:
		sys.exit('Brain ID not in  4000 or 3600 list')


	# load somas and stitching parameters
	somas = np.load(soma_path, allow_pickle=True)
	files, registrations, stitchings = extract_stitching_parameters(xml_file)

	# get image length
	x_list, y_list, z_list = [],[],[]

	for r in registrations:
		x = r['translation to regular grid'][3] + r['stitching transform'][3]
		y = r['translation to regular grid'][7] + r['stitching transform'][7]
		z = r['translation to regular grid'][11] + r['stitching transform'][11]

		x_list.append(x)
		y_list.append(y)
		z_list.append(z)
		
	minX, maxX = min(x_list), max(x_list) + X_SIZE
	x_len_stitched = maxX - minX + 1
	x_len_downsized = int(round(x_len_stitched / DOWNSIZE_FACTOR ))
	x_len_downsized_isotropic = int(round(x_len_downsized / ANISOTROPY_FACTOR))
	
	minY, maxY = min(y_list), max(y_list) + Y_SIZE
	y_len_stitched = maxY - minY + 1
	y_len_downsized = int(round(y_len_stitched / DOWNSIZE_FACTOR ))
	y_len_downsized_isotropic = int(round(y_len_downsized / ANISOTROPY_FACTOR))

	minZ, maxZ = min(z_list)/ANISOTROPY_FACTOR, max(z_list)/ANISOTROPY_FACTOR + Z_SIZE
	z_len_stitched = maxZ - minZ + 1
	z_len_downsized = int(round(z_len_stitched / DOWNSIZE_FACTOR))
	z_len_downsized_isotropic = int(round(z_len_downsized))

	# loop through somas
	for i,soma in enumerate(somas):
		num_duplicates = int(sum(soma!=None)/7)

		for j in range(num_duplicates):

			volume = soma[3 + 7*j]

			translation_to_grid_matrix, stitching_matrix, calibration_matrix = get_stitching_matrices(files,registrations,volume)	

			# Perform Stitching Transforms
			soma_coord = np.array([soma[4+7*j],soma[5+7*j],soma[6+7*j],1])
			transformed_soma_coord = np.matmul(calibration_matrix,soma_coord) 
			transformed_soma_coord = np.matmul(translation_to_grid_matrix,transformed_soma_coord)
			transformed_soma_coord = np.matmul(stitching_matrix,transformed_soma_coord)
			transformed_soma_coord = np.matmul(np.linalg.inv(calibration_matrix),transformed_soma_coord)

			# Shift coordinates
			x_shift = transformed_soma_coord[0] - minX
			y_shift = transformed_soma_coord[1] - minY
			z_shift = transformed_soma_coord[2] - minZ

			# Downsize
			x_downsized = int(round(x_shift/DOWNSIZE_FACTOR)) 
			y_downsized = int(round(y_shift/DOWNSIZE_FACTOR))
			z_downsized = int(round(z_shift/DOWNSIZE_FACTOR)) 

			# Make isotropic
			x_isotropic = int(round(x_downsized/ANISOTROPY_FACTOR))
			y_isotropic = int(round(y_downsized/ANISOTROPY_FACTOR))
			z_isotropic = z_downsized

			# Reslice Dimensions
			x_reslice = y_isotropic
			y_reslice = z_isotropic
			z_reslice = x_isotropic
			x_len_reslice = y_len_downsized_isotropic
			y_len_reslice = z_len_downsized_isotropic
			z_len_reslice = x_len_downsized_isotropic

			# Vertical Flip
			x_vert_flip = x_reslice
			y_vert_flip = y_len_reslice - y_reslice - 1
			z_vert_flip = z_reslice
	
			# shear_tranform
			x_len_shear = int(round(x_len_reslice))
			y_len_shear = int(round(y_len_reslice + abs(SHEAR_FACTOR)*x_len_reslice))
			z_len_shear = int(round(z_len_reslice))
			x_shear = x_vert_flip
			y_shear = int(round(y_vert_flip-y_len_reslice/2 + SHEAR_FACTOR*(x_vert_flip-x_len_reslice/2) + y_len_shear/2))  
			z_shear = z_vert_flip

			# reslicing
			x_reslice_2 = y_shear
			y_reslice_2 = z_shear
			z_reslice_2 = x_shear
			x_len_reslice_2 = y_len_shear
			y_len_reslice_2 = z_len_shear
			z_len_reslice_2 = x_len_shear

			# rotate
			x_rotate = y_len_reslice_2 - y_reslice_2 - 1
			y_rotate = x_reslice_2
			z_rotate = z_reslice_2
			x_len_rotate = y_len_reslice_2
			y_len_rotate = x_len_reslice_2
			z_len_rotate = z_len_reslice_2
	
			somas[i,7*j:3+7*j] = np.array([int(round(x_rotate)),int(round(y_rotate)),int(round(z_rotate))])
	print(somas[0:10])
	input()
	#np.save(os.path.join(brain_path,'transformed_somas_close.npy'),somas)
	#np.savetxt(os.path.join(brain_path,'transformed_somas_close.csv'),somas,delimiter=",",fmt="%s")
	
