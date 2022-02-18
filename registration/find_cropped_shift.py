from skimage import io
import numpy as np
import os
import tifffile as tif
from BRAINS import BRAINS

######################################################################
REGISTRATION_BASE_PATH = '/data/elowsky/OLST/registration/'
ARUN_REGISTRATION_PATH = '/data/anarasim/data/registration_images/'
######################################################################

transformed_labels_basepath = os.path.join(ARUN_REGISTRATION_PATH,'MOp_registeredFiles')
out_path = os.path.join(REGISTRATION_BASE_PATH,'crop_info.txt')

crop_info = np.genfromtxt(out_path,delimiter=' ',dtype=str)[:,0]

print()
for brain in BRAINS:

	if brain in crop_info:
		continue

	brain_path = os.path.join(REGISTRATION_BASE_PATH,brain)

	print('Brain:',brain)
	
	if brain == '170329_500':
		brain = '170329'

	# read in transformed labels cropped
	transformed_cropped_labels_brain_path = os.path.join(brain_path,'mop_labels_registered.nrrd')
	if not os.path.exists(transformed_cropped_labels_brain_path):
		continue
	print('Read in Cropped Labels...')
	transformed_cropped_labels = io.imread(transformed_cropped_labels_brain_path)

	# read in transformed labels
	print('Read in Transformed Labels...')
	transformed_labels_path = os.path.join(transformed_labels_basepath,'MOp_'+brain+'_Emxcre_25um_isotropic_removedStripes.tif')
	transformed_labels = io.imread(transformed_labels_path)

	# get highest y value to start search
	y_top_label = np.min(np.where(transformed_labels > 0)[1])

	# iterate through possible shifts
	for i in range(y_top_label-75,y_top_label):
		print(i)

		# place labels in image
		out_image = np.zeros(shape=transformed_labels.shape)
		start_y = i
		end_y = start_y + transformed_cropped_labels.shape[1]
		out_image[:,start_y:end_y,:] = transformed_cropped_labels

		if np.all(out_image == transformed_labels):
			print('y shift:',i)
			with open(out_path,'a') as fp:
				fp.write(brain + ' ' + str(i)+'\n')
			break













