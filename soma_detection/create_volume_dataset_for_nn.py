import numpy as np
from utils import *

#########################
BRAIN = 170329
CROPPED_SIZE = (122,122,20)
#########################

# Load Soma location File
soma_path = '/data/anarasim/data/soma_detections/' + str(brain_id) + '_Emxcre/soma_detection/model/1000/detected_somas_model.npy'
somas = np.load(soma_path,allow_pickle=True)

image_path = '/data/anarasim/data/croppedraw_croppedsegmented/' + str(brain_id) + '_Emxcre/'

output_path = '/data/elowsky/OLST/reconstruction/' + str(brain_id) + '/'

# Create volume dataset
X = np.empty((len(somas),CROPPED_SIZE[2],CROPPED_SIZE[1],CROPPED_SIZE[0],2),dtype=np.int32)

for soma_number in range(1,len(somas)+1):
	if soma_number % 500 == 0:		
		print(soma_number)
	soma = somas[soma_number-1]

	volume = soma[3]

	path_seg = image_path + 'cropped_segmented/' + volume + '/' + str(soma_number) + '.tif'
	path_raw = image_path + 'cropped_raw/' + volume + '/' + str(soma_number) + '.tif'

	# load and pad images
	image_seg = tif_to_numpyArray(path_seg)
	image_seg = pad(image_seg,'seg',CROPPED_SIZE)

	image_raw = tif_to_numpyArray(path_raw)
	image_raw = pad(image_raw,'raw',CROPPED_SIZE)

	X[soma_number-1,:,:,:,0] = image_seg
	X[soma_number-1,:,:,:,1] = image_raw

np.save(output_path + 'input_somas_cnn',X)

# normalize

print('Normalizing...')

#nomalized between 0 and 255
mins = X.min(axis=(1,2,3), keepdims=True)
ptps = np.ptp(X, axis=(1,2,3), keepdims=True)
ptps[ptps == 0] = 1 # dont want to divide by zero
X = (255*(X-mins)/ptps).astype(np.uint8)

np.save(output_path + 'input_somas_cnn_normalized',X)




