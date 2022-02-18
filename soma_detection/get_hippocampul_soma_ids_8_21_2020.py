import numpy as np
from skimage.io import imread
from os.path import join

###########################################################
BRAINS = ['180606', '190522']
REGISTRATION_BASE_PATH = '/data/elowsky/OLST/registration/'
SOMAS_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
SAVE = True
############################################################

# load shift infomration
shift_path = join(REGISTRATION_BASE_PATH,'crop_info.txt')
shift_info = np.genfromtxt(shift_path,dtype=str)

for brain in BRAINS:

	print()
	print('Brain:', brain)
	print()

	somas_brain_path = join(SOMAS_BASE_PATH, brain)
	registration_brain_path = join(REGISTRATION_BASE_PATH, brain)

	labels_path = join(registration_brain_path,'registered_hippocampul_formation_25x25x25.nrrd')
	somas_path = join(somas_brain_path,'transformed_somas.npy')
	fullsize_somas_path = join(somas_brain_path,'triaged_somas_close_200.npy')
	
	out_path = join(somas_brain_path, 'hippocampul_formation_soma_ids.csv')
	out_path_center = join(somas_brain_path, 'hippocampul_formation_soma_ids_center.csv')	


	# get shift
	shift = int(shift_info[shift_info[:,0] == brain,1][0])

	# load somas
	somas = np.load(somas_path, allow_pickle=True)
	full_size_somas = np.load(fullsize_somas_path, allow_pickle=True)

	# get 25um coords
	# shift y cooridnate
	coords = somas[:,:3]
	coords[:,1] -= shift
	coords = coords.astype(np.int64)

	# load labels
	labels = imread(labels_path)	
	coords[:,1] = np.clip(coords[:,1],0,labels.shape[1]-1)

	# get hippocampul somas
	hippocampul_soma_ids = np.where(labels[coords[:,2],coords[:,1],coords[:,0]])[0] + 1
	out_data = np.hstack((hippocampul_soma_ids.reshape(-1,1),full_size_somas[hippocampul_soma_ids-1][:,3:7]))
	out_data_center = out_data[(out_data[:,2] > 900) & (out_data[:,2] < 1100) & (out_data[:,3] > 400) & (out_data[:,3] < 600)]

	if SAVE:
		np.savetxt(out_path,out_data,delimiter=',',fmt=['%d','%s','%d','%d','%d'])
		np.savetxt(out_path_center,out_data_center,delimiter=',',fmt=['%d','%s','%d','%d','%d'])



