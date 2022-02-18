import numpy as np
from utils import *
from sklearn.model_selection import train_test_split
from shutil import copyfile
import tifffile as tif



###################
BRAIN_ID = 180206
THRESHOLD = 1000
###################

#########################################################################################################
brain_path = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '_Emxcre_Reconstruction/'
model_path = brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/'
out_path = model_path + 'train_test_set/'
#########################################################################################################

######## Pads image with zeros #################
def pad(image,seg_or_raw):
	if seg_or_raw == 'seg': 
		if image.shape != (20,122,122):
			image = np.pad(image,((0,20-image.shape[0]),(0,122-image.shape[1]),(0,122-image.shape[2])),'constant')
	if seg_or_raw == 'raw': 
		if image.shape != (20,122,122):
			image = np.pad(image,((0,20-image.shape[0]),(0,122-image.shape[1]),(0,122-image.shape[2])),'constant',constant_values=np.amin(image))
	return image
################################################


########### Load data ##########################
tp_somas_path = model_path + '/tp/tp.npy'
tp_somas = np.load(tp_somas_path,allow_pickle=True)

fp_somas_path = model_path + '/fp/fp.npy'
fp_somas = np.load(fp_somas_path,allow_pickle=True)

all_fp_path = model_path + '/fp/fp.npy'
all_fp_somas = np.load(all_fp_path,allow_pickle=True)
##################################################

n_tp_somas = len(tp_somas)
n_fp_somas = len(fp_somas)

if False:
	print('Splitting into train and test sets....')
	all_somas = np.vstack((tp_somas,fp_somas))
	all_labels = np.concatenate([np.ones(n_tp_somas),np.zeros(n_fp_somas)])

	# split into train and test
	X_train_features, X_test_features, y_train_features, y_test_features = train_test_split(all_somas,all_labels,test_size=.00000001,shuffle=True)

	X = np.vstack((X_train_features,X_test_features))
	y = np.concatenate((y_train_features,y_test_features))
	np.save(out_path + 'X_features',X)
	np.savetxt(out_path + 'X_features.csv',X,delimiter=",",fmt="%s")
	np.save(out_path + 'y',y)
else:
	# Load train/test sets
	print('Loading train and test sets....')
	X_train_features = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/X_features.npy',allow_pickle=True)

# go through all train and test and see whether its in tp or fp, and then create image datasets

if True:
	X_train_images = np.empty((len(X_train_features),20,122,122,2),dtype=np.int32)

	for i in range(len(X_train_features)):
		if i % 100 == 0:		
			print(i)
		soma_features = X_train_features[i]
		volume = soma_features[3]

		if any((soma_features[0:7] == tp_somas[:,0:7]).all(1)):
			soma_number = np.where((soma_features[0:7] == tp_somas[:,0:7]).all(axis=1))[0][0] + 1
			source_seg = model_path + 'tp/cropped_segmented/' + volume + '/' + str(soma_number) + '.tif'
			source_raw = model_path + 'tp/cropped_raw/' + volume + '/' + str(soma_number) + '.tif'
		if any((soma_features[0:7] == all_fp_somas[:,0:7]).all(1)):
			soma_number = np.where((soma_features[0:7] == all_fp_somas[:,0:7]).all(axis=1))[0][0] + 1
			source_seg = model_path + 'fp/cropped_segmented/' + volume + '/' + str(soma_number) + '.tif'
			source_raw = model_path + 'fp/cropped_raw/' + volume + '/' + str(soma_number) + '.tif'

		image_seg = tif_to_numpyArray(source_seg)
		image_seg = pad(image_seg,'seg')

		image_raw = tif_to_numpyArray(source_raw)
		image_raw = pad(image_raw,'raw')

		X_train_images[i,:,:,:,0] = image_seg
		X_train_images[i,:,:,:,1] = image_raw

	np.save(out_path + 'X_images',X_train_images)
	
if False:
	X_test_images = np.empty((len(X_test_features),20,122,122,2), dtype=np.int32)
	for i in range(len(X_test_features)):
		if i % 100 == 0:		
			print(i)
		soma_features = X_test_features[i]
		volume = soma_features[3]
	
		if any((soma_features[0:7] == tp_somas[:,0:7]).all(1)):
			soma_number = np.where((soma_features[0:7] == tp_somas[:,0:7]).all(axis=1))[0][0] + 1
			source_seg = model_path + 'tp/cropped_segmented/' + volume + '/' + str(soma_number) + '.tif'
			source_raw = model_path + 'tp/cropped_raw/' + volume + '/' + str(soma_number) + '.tif'
		if any((soma_features[0:7] == fp_somas[:,0:7]).all(1)):
			soma_number = np.where((soma_features[0:7] == fp_somas[:,0:7]).all(axis=1))[0][0] + 1
			source_seg = model_path + 'fp/cropped_segmented/' + volume + '/' + str(soma_number) + '.tif'
			source_raw = model_path + 'fp/cropped_raw/' + volume + '/' + str(soma_number) + '.tif'

		image_seg = tif_to_numpyArray(source_seg)
		image_seg = pad(image_seg,'seg')

		image_raw = tif_to_numpyArray(source_raw)
		image_raw = pad(image_raw,'raw')

		X_test_images[i,:,:,:,0] = image_seg
		X_test_images[i,:,:,:,1] = image_raw


	np.save(out_path + 'X_images',X_test_images)





