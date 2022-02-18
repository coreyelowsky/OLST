import numpy as np
from shutil import copyfile
import os

###############
err_type = 'fp'
###############

##################
BRAIN_ID = 180206
THRESHOLD = 1000
##################

##############################################################################################
brain_path = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '_Emxcre_Reconstruction/'
model_path = brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/'
err_path = model_path + err_type + '/'
truth_path = brain_path + 'soma_detection/truth/'
##############################################################################################

#--------#
truth_somas_path = truth_path + 'detected_somas_truth.npy'
truth_somas = np.load(truth_somas_path,allow_pickle=True)

model_somas_path = model_path + 'detected_somas_model.npy'
model_somas = np.load(model_somas_path,allow_pickle=True)
#--------#

#--------#
err_somas_path = err_path + err_type + '.npy'
err_somas = np.load(err_somas_path,allow_pickle=True)



#--------#

if err_type == 'fn':
	source = truth_somas
	destination = err_somas
elif err_type == 'tp':
	source = model_somas
	destination = err_somas
elif err_type == 'fp':
	source = model_somas
	destination = err_somas

for i in range(1,len(destination)+1):

	soma = destination[i-1,0:7]


	index = np.where((source[:,0:7] == soma).all(axis=1))[0][0] + 1

	if err_type == 'fn':
		source_path_seg = truth_path + 'cropped_segmented/' + soma[3] + '/' + str(index) + '.tif'
		destination_path_seg = err_path + 'cropped_segmented/'  + soma[3] + '/'  + str(i) + '.tif'

		source_path_raw = truth_path + 'cropped_raw/' + soma[3] + '/' + str(index) + '.tif'
		destination_path_raw = err_path + 'cropped_raw/' + soma[3] + '/'  + str(i) + '.tif'

		if not os.path.isdir(err_path + 'cropped_segmented/'  + soma[3] + '/'):
			os.mkdir(err_path + 'cropped_segmented/'  + soma[3] + '/')	

		if not os.path.isdir(err_path + 'cropped_raw/'  + soma[3] + '/'):
			os.mkdir(err_path + 'cropped_raw/'  + soma[3] + '/')
	
	elif err_type == 'tp':
		source_path_seg = model_path + 'cropped_segmented/' + soma[3] + '/' + str(index) + '.tif'
		destination_path_seg = err_path + 'cropped_segmented/'  + soma[3] + '/'  + str(i) + '.tif'

		source_path_raw = model_path + 'cropped_raw/' + soma[3] + '/' + str(index) + '.tif'
		destination_path_raw = err_path + 'cropped_raw/' + soma[3] + '/'  + str(i) + '.tif'

		if not os.path.isdir(err_path + 'cropped_segmented/'  + soma[3] + '/'):
			os.mkdir(err_path + 'cropped_segmented/'  + soma[3] + '/')	

		if not os.path.isdir(err_path + 'cropped_raw/'  + soma[3] + '/'):
			os.mkdir(err_path + 'cropped_raw/'  + soma[3] + '/')
	
	elif err_type == 'fp':
		source_path_seg = model_path + 'cropped_segmented/' + soma[3] + '/' + str(index) + '.tif'
		destination_path_seg = err_path + 'cropped_segmented/'  + soma[3] + '/'  + str(i) + '.tif'

		source_path_raw = model_path + 'cropped_raw/' + soma[3] + '/' + str(index) + '.tif'
		destination_path_raw = err_path + 'cropped_raw/' + soma[3] + '/'  + str(i) + '.tif'

		if not os.path.isdir(err_path + 'cropped_segmented/'  + soma[3] + '/'):
			os.mkdir(err_path + 'cropped_segmented/'  + soma[3] + '/')	

		if not os.path.isdir(err_path + 'cropped_raw/'  + soma[3] + '/'):
			os.mkdir(err_path + 'cropped_raw/'  + soma[3] + '/')
	

	copyfile(source_path_seg,destination_path_seg)
	copyfile(source_path_raw,destination_path_raw)

	
