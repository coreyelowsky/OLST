import numpy as np
from utils import *
import os
import sys
import tifffile as tif
import params

##########################################################
SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
REGISTRATION_BASE_PATH = '/data/elowsky/OLST/registration/'
ARUN_REGISTRATION_PATH = '/data/anarasim/data/registration_images/forArun/'
REGION = 'piriform'
REGION_IDS = params.REGION_IDS[REGION]
SAVE = True
##########################################################

total_somas = 0    
total_layers = [0,0,0,0,0,0]  

print()
print('###################')
print('Region:', REGION)
print('###################')
print()

# loop through brains
for brain in params.BRAINS:

	# set up paths
	brain_path = os.path.join(SOMA_BASE_PATH,brain) 
	soma_path = os.path.join(brain_path,'transformed_somas.npy')
	labels_path = os.path.join(ARUN_REGISTRATION_PATH, brain+'_Emxcre_25um_isotropic_removedStripes_crop_ARA.tif')                     
	
	# make sure somas files exist
	if not os.path.exists(soma_path):
		print(soma_path)
		sys.exit('Soma File Does Not Exist')	

	if not os.path.exists(labels_path):
		print('Label File Does Not Exist')
		continue			


	# load somas
	somas = np.load(soma_path,allow_pickle=True)
	
	# load labels
	labels = tif.imread(labels_path)

	# set up output
	# if value are dicts then there are subregion
	if isinstance(list(REGION_IDS.values())[0],dict):
		output_somas = np.empty((0,30),dtype=object)
	else:
		output_somas = np.empty((0,29),dtype=object)

	
	local_layers = [0,0,0,0,0,0]

	for i,soma in enumerate(somas):

		(x,y,z) = (int(soma[0]),int(soma[1]),int(soma[2]))

		soma_label = labels[z,y,x]
		
		# if value are dicts then there are subregion
		if isinstance(list(REGION_IDS.values())[0],dict):

			for region, region_dict in REGION_IDS.items():
				if soma_label in region_dict:
					layer = region_dict[soma_label]
					soma = np.concatenate(([layer],soma, [region]))
					output_somas = np.vstack((output_somas,soma))
					local_layers[layer-1] += 1

		else:
			if soma_label in REGION_IDS:
				layer = REGION_IDS[soma_label]
				soma = np.concatenate(([layer], soma))
				output_somas = np.vstack((output_somas,soma))
				local_layers[layer-1] += 1


	print(brain,':',len(output_somas), local_layers)
	total_somas += len(output_somas)
	total_layers  = [x + y for x,y in zip(total_layers, local_layers)]
	

	if SAVE:	
		np.savetxt(os.path.join(brain_path,'regions', REGION + '.csv'),output_somas,delimiter=",",fmt="%s")

print()
print('Total Somas:', total_somas)
print('Total Layer:', total_layers)


	
			
