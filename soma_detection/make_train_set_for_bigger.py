import numpy as np
from utils import *
import math

##################
BRAIN_ID = 180926
THRESHOLD = 1000
##################

##############################################################################################
brain_path = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '_Emxcre_Reconstruction/'
############################################################################################## 

path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/training/X_train_features.npy'
data = np.load(path, allow_pickle=True)

path_labels = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/training/y_train.npy'
labels = np.load(path_labels, allow_pickle=True)


###############################################################################################################
tp_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/tp/tp.npy'
tp_somas = np.load(tp_path, allow_pickle=True)

fp_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/fp/fp.npy'
fp_somas = np.load(fp_path, allow_pickle=True)
######################################################################################################
print(data.shape)

out_data = np.empty(shape=data.shape,dtype=object)


for soma_number in range(1, len(data)+1):
	print('Soma #: ' + str(soma_number))
	soma = data[soma_number-1,:]

	out_data[soma_number-1,0:7] = soma[0:7]

	volume = soma[3]
	
	label = labels[soma_number-1]
	
	
	if label == 1:
		index = np.where((tp_somas[:,0:7] == soma[0:7]).all(axis=1))[0][0]
		out_data[soma_number-1,:] = tp_somas[index]



	elif label == 0:
		index = np.where((fp_somas[:,0:7] == soma[0:7]).all(axis=1))[0][0]
		out_data[soma_number-1,:] = fp_somas[index]




np.save('/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/X_features_bigger.npy',out_data)	
	
