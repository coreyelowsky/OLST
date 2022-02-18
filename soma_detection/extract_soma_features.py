import numpy as np
import math
import tifffile as tif
from timeit import default_timer as timer
from utils import *
import os

BRAIN_ID = 180206
THRESHOLD = 1000

brain_path = '/data/elowsky/OLST/reconstruction/'+ str(BRAIN_ID) +'_Emxcre_Reconstruction/'
model_path = brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/'

#soma_array = np.load(brain_path + 'soma_detection/truth/detected_somas_truth.npy',allow_pickle=True)
#soma_array = np.load(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/detected_somas_model.npy',allow_pickle=True)
#soma_array = np.load(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/fp/fp.npy',allow_pickle=True)

soma_array_1 = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/tp/tp.npy',allow_pickle=True)
soma_array_2 = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/fp/fp.npy',allow_pickle=True)

soma_array = np.vstack((soma_array_1,soma_array_2))
tp_somas = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/tp/tp.npy',allow_pickle=True)
fp_somas = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/fp/fp.npy',allow_pickle=True)

num_cols = soma_array.shape[1]

if num_cols < 18:
	soma_array = np.hstack((soma_array,np.full((len(soma_array),18-num_cols),None)))


for soma_number in range(1, len(soma_array)+1):
	if soma_number % 200 == 0:
		print('Soma #: ' + str(soma_number))
	soma = soma_array[soma_number-1,:]

	volume = soma[3]
	
	
	if any((soma[0:7] == tp_somas[:,0:7]).all(1)):
		soma_id = np.where((soma[0:7] == tp_somas[:,0:7]).all(axis=1))[0][0] + 1
		source_seg = model_path + 'tp/cropped_segmented/' + volume + '/' + str(soma_id) + '.tif'
		source_raw = model_path + 'tp/cropped_raw/' + volume + '/' + str(soma_id) + '.tif'
	elif any((soma[0:7] == fp_somas[:,0:7]).all(1)):
		soma_id = np.where((soma[0:7] == fp_somas[:,0:7]).all(axis=1))[0][0] + 1
		source_seg = model_path + 'fp/cropped_segmented/' + volume + '/' + str(soma_id) + '.tif'
		source_raw = model_path + 'fp/cropped_raw/' + volume + '/' + str(soma_id) + '.tif'
	else:
		continue

	cropped_volume_raw = tif_to_numpyArray(source_raw)
	cropped_volume_seg = tif_to_numpyArray(source_seg)

	# MODEL
	#cropped_volume_raw = tif_to_numpyArray(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/cropped_raw/' + volume + '/' + str(soma_number) + '.tif')
	#cropped_volume_seg = tif_to_numpyArray(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/cropped_segmented/' + volume + '/' + str(soma_number) + '.tif')

	# ERROR
	#cropped_volume_raw = tif_to_numpyArray(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/fp/cropped_raw/' + volume + '/' + str(soma_number) + '.tif')
	#cropped_volume_seg = tif_to_numpyArray(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/fp/cropped_segmented/' + volume + '/' + str(soma_number) + '.tif')

	raw_mean = np.mean(cropped_volume_raw)
	raw_std = np.std(cropped_volume_raw)
	raw_max_intensity = np.amax(cropped_volume_raw)
	
	# raw info
	soma_array[soma_number-1,7] = int(round(raw_mean))
	soma_array[soma_number-1,8] = int(round(raw_std))
	soma_array[soma_number-1,9] = int(raw_max_intensity)

	seg_mean_1 = np.mean(cropped_volume_seg[cropped_volume_seg > 0])
	seg_std_1 = np.std(cropped_volume_seg[cropped_volume_seg > 0])
	seg_mean_2 = np.mean(cropped_volume_seg)
	seg_std_2 = np.std(cropped_volume_seg)
	seg_max_intensity = np.amax(cropped_volume_seg)
	seg_count = (cropped_volume_seg > 0).sum()

	if math.isnan(seg_mean_1):
		seg_mean_1 = 0
	if math.isnan(seg_std_1):
		seg_std_1 = 0
	if math.isnan(seg_mean_2):
		seg_mean_2 = 0
	if math.isnan(seg_std_2):
		seg_std_2 = 0

	# seg info
	soma_array[soma_number-1,10] = int(round(seg_mean_1))
	soma_array[soma_number-1,11] = int(round(seg_std_1))
	soma_array[soma_number-1,12] = int(round(seg_mean_2))
	soma_array[soma_number-1,13] = int(round(seg_std_2))
	soma_array[soma_number-1,14] = int(seg_max_intensity)
	soma_array[soma_number-1,15] = int(seg_count)

	# z_y
	soma_array[soma_number-1,16] = int(volume[1:3])
	soma_array[soma_number-1,17] = int(volume[5:7])

np.save('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/X_features',soma_array)
np.savetxt('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/X_features.csv',soma_array,delimiter=",",fmt="%s")
	
# Model
#np.save(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/detected_somas_model.npy',soma_array)
#np.savetxt(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/detected_somas_model.csv',soma_array,delimiter=",",fmt="%s")

# Error
#np.save(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/fp/fp_features.npy',soma_array)
#np.savetxt(brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/fp/fp_features.csv',soma_array,delimiter=",",fmt="%s")



