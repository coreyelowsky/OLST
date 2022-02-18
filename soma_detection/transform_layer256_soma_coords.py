import numpy as np
import params
import sys

from os import listdir
from os.path import join, exists
from utils import extract_stitching_parameters, get_stitching_matrices

##### PARAMETERS ###############
SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
XML_BASE_PATH = '/data/anarasim/data/mop/'
BRAINS = ['180606']
LAYERS = [2,5,6]
SAVE = False
################################

# iterate through brains
for brain in BRAINS:

	# load somas file
	brain_somas_path = join(SOMA_BASE_PATH, brain)
	somas = np.load(join(brain_somas_path,'motor_cortex_somas.npy'),allow_pickle=True)

	# iterate through somas
	for soma_number,soma in enumerate(somas,1):

		layer = soma[0]

		print('###########################')
		print('Soma #:',soma_number)
		print('Layer:',layer)

		# only transform somas in specified layers
		if layer in LAYERS:

			####### Logic to find XML file ##########

			# get raw and prediction path for soma
			# if neither raw noor prediction path exists then go to next soma
			predictions_path = join(XML_BASE_PATH, 'Layer' + str(layer) + '_predictions',brain,str(soma_number))
			print('Predictions:',predictions_path)
			if not exists(predictions_path):
				print('Predictions Path Does Not Exist: '+predictions_path)
				somas[soma_number-1,1:4] = np.array([-1,-1,-1])
				continue

			raw_path = join(XML_BASE_PATH, 'Layer' + str(layer),brain,str(soma_number))
			print('Raw:',raw_path)
			if not exists(raw_path):
				print('Raw Path Does Not Exist: '+raw_path)
				somas[soma_number-1,1:4] = np.array([-1,-1,-1])
				continue

			# if nothing in directory or dataset.xml doesnt exist in path
			# look for dataset.xml in other layers (could be "same as")
			if not exists(join(predictions_path,'dataset.xml')):
	
				# look for 'Same' for soma number with identical xml
				file_list = [f for f in listdir(raw_path) if f[0:4] == 'Same']
				if len(file_list) == 0:
					print('XML isnt in predictions directory and no Same As in raw directory')
					somas[soma_number-1,1:4] = np.array([-1,-1,-1])
					continue
				else:
					sameAsFile = file_list[0].replace('.txt','')
					soma_num_same_xml = int(sameAsFile.split('_')[-1])
					layer_same_xml = somas[soma_num_same_xml-1][0]
					predictions_path = join(XML_BASE_PATH, 'Layer' + str(layer_same_xml)+'_predictions',brain,str(soma_num_same_xml))
					print('Predictions path (Same As):', predictions_path)
		
			###########################################
	
			# xml path
			xml_path = join(predictions_path,'dataset.xml')
			print('XML Path:',xml_path)
			volume = soma[4]

			# parse stitching parameters 
			files, registrations, stitchings = extract_stitching_parameters(xml_path)
			translation_to_grid_matrix, stitching_matrix, calibration_matrix = get_stitching_matrices(files,registrations,volume)

			# get min bounds of fused image
			minX = min([x['translation to regular grid'][3] + x['stitching transform'][3] for x in registrations])
			minY = min([x['translation to regular grid'][7] + x['stitching transform'][7] for x in registrations])
			minZ = min([x['translation to regular grid'][11] + x['stitching transform'][11] for x in registrations])

			# extract soma
			soma_coord = np.array([soma[5],soma[6],soma[7],1])

			# transform to stitching coordinates
			transformed_soma_coord = np.matmul(calibration_matrix,soma_coord) 
			transformed_soma_coord = np.matmul(translation_to_grid_matrix,transformed_soma_coord)
			transformed_soma_coord = np.matmul(stitching_matrix,transformed_soma_coord)

			# get coordinates in fused image and downsize Z to original resolution
			x_new = int((transformed_soma_coord[0] - minX))
			y_new = int((transformed_soma_coord[1] - minY))
			z_new = int((transformed_soma_coord[2] - minZ)/params.ANISOTROPY_FACTOR)

			somas[soma_number-1,1:4] = np.array([x_new,y_new,z_new])
		else:
			print('DONT TRANSFORM SOMA - Not in specified layer')
			somas[soma_number-1,1:4] = np.array([-1,-1,-1])


	if SAVE:		
		somas = somas[:,0:8]
		np.save(join(brain_somas_path,'transformed_motor_cortex_somas.npy'),somas)
		np.savetxt(join(brain_somas_path,'transformed_motor_cortex_somas.csv'),somas,delimiter=",",fmt="%s")

