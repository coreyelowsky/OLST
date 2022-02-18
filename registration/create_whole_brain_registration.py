import os
import numpy as np
import BRAINS
from skimage.io import imread
import tifffile as tif

###########################################################
REGISTRATION_BASE_PATH = '/data/elowsky/OLST/registration/'
ARUN_REGISTRATION_PATH = '/data/anarasim/data/registration_images/MOp_registeredFiles'
###########################################################

# load crop info
crop_info_path = os.path.join(REGISTRATION_BASE_PATH,'crop_info.txt')
crop_info = np.genfromtxt(crop_info_path,delimiter=' ',dtype=str)

for brain in BRAINS.BRAINS:

	if brain == '190306':
		continue

	print()
	print('Brain:',brain)

	crop_shift = int(crop_info[crop_info[:,0] == brain][0][1])

	# load mop registered to get size
	if brain == '170329_500':
		brain = '170329'
	print('Load Mop Registered...')
	mop_registered_path = os.path.join(ARUN_REGISTRATION_PATH,'MOp_'+brain+'_Emxcre_25um_isotropic_removedStripes.tif')
	mop_registered = imread(mop_registered_path)
	mop_shape = mop_registered.shape

	# load whole brain registered labels
	if brain == '170329':
		brain = '170329_500'
	print('Load Whole Brain Labels...')
	whole_brain_path = os.path.join(REGISTRATION_BASE_PATH,brain,'whole_brain_labels_registered.nrrd')
	whole_brain_labels = imread(whole_brain_path)

	# create output label image
	out_image = np.zeros(shape=mop_shape,dtype=np.uint32)
	out_image[:,crop_shift:crop_shift+whole_brain_labels.shape[1],:] = whole_brain_labels

	# save image
	print('Save Image...')
	out_path = os.path.join(REGISTRATION_BASE_PATH,brain,'whole_brain_labels_registered_padded.nrrd')
	tif.imsave(out_path,out_image.astype(np.uint32))
	









	
