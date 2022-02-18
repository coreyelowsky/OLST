from BRAINS import BRAINS
from os.path import join
import numpy as np
import tifffile as tif
from skimage.io import imread
import csv

#####################################################
SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
ANNOTATION_IMAGE_BASE_PATH = '/data/elowsky/OLST/registration/'
ANNOTATION_CSV_PATH = '/data/anarasim/data/registration_images/20170803_Annotation/ARA2_annotation_structure_info.csv'
SAVE = False
#####################################################

# remove brain without registration
BRAINS.remove('190306')

# load and parse annotation csv
# ignore first row of headers
# save id and names in dictionary
# add in zero for outside brain ('OUTSIDE BRAIN')
annotation_dict = {}
with open(ANNOTATION_CSV_PATH, newline='') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for i,row in enumerate(reader):
		if i != 0:
			annotation_dict[int(row[0])] = row[1]
annotation_dict[0] = 'OUTSIDE BRAIN'

not_found_set = set()

# iterate through brains
for brain in BRAINS:

	print()
	print('Brain:', brain)
	print()

	# set up paths
	somas_brain_path = join(SOMA_BASE_PATH, brain)
	somas_path = join(somas_brain_path,'transformed_somas.npy')

	annotation_brain_path = join(ANNOTATION_IMAGE_BASE_PATH,brain)
	annotation_image_path = join(annotation_brain_path,'whole_brain_labels_registered_padded.nrrd')

	# load annotation image
	annotation_image = tif.imread(annotation_image_path)

	# load somas
	somas = np.load(somas_path,allow_pickle=True)

	# get labels
	annotation_labels = annotation_image[somas[:,2].astype(int),somas[:,1].astype(int),somas[:,0].astype(int)]

	# get annotation names
	annotation_names = []
	for label in annotation_labels:
		if label not in annotation_dict:
			not_found_set.add(label)
			annotation_names.append('NOT IN ANNOTAITON CSV')
		else:
			annotation_names.append(annotation_dict[label])

	# Oustide Brain: 0
	# Not Layer 1-6: -1
	# Not in Annotation CSV: -2
	annotation_layers = []
	for name in annotation_names:
		if 'layer 1' in name or 'Layer 1' in name:
			annotation_layers.append(1)
		elif 'layer 2/3' in name or 'Layer 2/3' in name:
			annotation_layers.append(2)
		elif 'layer 4' in name or 'Layer 4' in name:
			annotation_layers.append(4)
		elif 'layer 5' in name or 'Layer 5' in name:
			annotation_layers.append(5)
		elif 'layer 6' in name or 'Layer 6' in name or '6a' in name or '6b' in name:
			annotation_layers.append(6)	
		elif name == 'OUTSIDE BRAIN':
			annotation_layers.append(0)
		elif name == 'NOT IN ANNOTAITON CSV':
			annotation_layers.append(-2)
		else:
			annotation_layers.append(-1)

	# create output somas
	somas = somas[:,[3,4,5,6,10,11,12,13,17,18,19,20,24,25,26,27]]
	somas = np.hstack((np.array(annotation_layers).reshape(-1,1),somas))
	
	# save somas
	if SAVE:
		np.savetxt(join(somas_brain_path,'whole_brain_isocortex_somas.csv'),somas,delimiter=',',fmt='%s')
		np.save(join(somas_brain_path,'whole_brain_isocortex_somas.npy'),somas)

	






	

