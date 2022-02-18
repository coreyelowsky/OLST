import numpy as np
import tifffile as tif
from timeit import default_timer as timer
from utils import *
import os

print()
print('###########')
print('Crop Somas')
print('###########')
print()

################################################################################
BRAIN = 180206
SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
RAW_PATH = '/mnt/brainstore8/anarasim/data/180206_Emxcre_Reconstruction/100pc/'
SEG_PATH = '/mnt/brainstore8/anarasim/data/180206_Emxcre_Reconstruction/predictions/'
IM_TYPE = 'raw' #'seg'
CROP_RADIUS = 25 # microns
RES_FACTORS = np.array([0.406,0.406,2.5])
################################################################################

soma_brain_path = os.path.join(SOMA_BASE_PATH, str(BRAIN))
soma_path = os.path.join(soma_brain_path,'detected_somas_clustering.npy')

print('Brain:',BRAIN)
print('Raw Path:',RAW_PATH)
print('Seg Path:',SEG_PATH)
print('Crop Radius:',CROP_RADIUS,'(um)')

# load somas
somas = np.load(soma_path, allow_pickle=True)

crop_radii = (np.array([CROP_RADIUS,CROP_RADIUS,CROP_RADIUS])/RES_FACTORS).astype(np.uint16)
print('Cropping Radii (voxels):', crop_radii)
print()


volume_prior = None
for soma_number,soma in enumerate(somas,1):

	print('Soma #:',soma_number)

	volume = soma[3]
	soma = soma[4:7]

	volume_full_path_raw = os.path.join(RAW_PATH,volume + '.tif')
	volume_full_path_seg = os.path.join(SEG_PATH,volume + '.tif')

	print(volume_full_path_raw)
	input()

	if volume != volume_prior:
		out_vol_path_raw = output_path + 'cropped_raw/'+ volume + '/'
		out_vol_path_seg = output_path + 'cropped_segmented/' + volume + '/'

		if im_type == 'raw':
			if not os.path.isdir(out_vol_path_raw):
				os.mkdir(out_vol_path_raw)
			# load volume
			print("Loading Volume: " + volume_full_path_raw)
			load_start = timer()
			full_volume_raw = tif_to_numpyArray(volume_full_path_raw)
			load_end = timer()
			print(str(round(load_end - load_start)) + " seconds")
			volume_prior = volume
		elif im_type == 'seg':
			if not os.path.isdir(out_vol_path_seg):
				os.mkdir(out_vol_path_seg)				

			# load volume
			print("Loading Volume: " + volume_full_path_seg)
			load_start = timer()
			full_volume_seg = tif_to_numpyArray(volume_full_path_seg)
			load_end = timer()
			print(str(round(load_end - load_start)) + " seconds")
			volume_prior = volume
		


	(x,y,z) = (int(local_coords[0]), int(local_coords[1]),int(local_coords[2]))

	# SMALLER
	if z-z_pix_radius < 0:
		min_z = 0
	else:
		min_z = z-z_pix_radius

	if y-y_pix_radius < 0:
		min_y = 0
	else:
		min_y = y-y_pix_radius

	if x-x_pix_radius < 0:
		min_x = 0
	else:
		min_x = x-x_pix_radius


	# crop volumes
	print('Crop Volumes...')
	if im_type == 'raw':
		cropped_volume_raw = full_volume_raw[min_z:z+z_pix_radius,min_y:y+y_pix_radius,min_x:x+x_pix_radius]
	elif im_type == 'seg':
		cropped_volume_seg = full_volume_seg[min_z:z+z_pix_radius,min_y:y+y_pix_radius,min_x:x+x_pix_radius]

	
	# save cropped volume
	print('Save Volumes...')
	if im_type == 'raw':
		tif.imsave(out_vol_path_raw + str(soma_number) + '.tif',cropped_volume_raw)
	elif im_type == 'seg':
		tif.imsave(out_vol_path_seg + str(soma_number) + '.tif',cropped_volume_seg)




	
