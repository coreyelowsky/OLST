import os
import numpy as np
from BRAINS import BRAINS

################################################
BASE_PATH = '/data/elowsky/OLST/soma_detection/'
################################################

for brain in BRAINS:
	
	print()
	print('Brain:',brain)	

	# set up paths
	brain_path = os.path.join(BASE_PATH,brain)
	somas_path = os.path.join(brain_path,'triaged_somas_duplicates.npy')
	triaged_somas_close_path = os.path.join(brain_path,'triaged_somas_close_200.npy') 
	
	# load somas
	somas = np.load(somas_path,allow_pickle=True)
	triaged_somas = np.load(triaged_somas_close_path,allow_pickle=True)

	# get non common somas from two somas lists
	non_common = np.array(list(set(map(tuple,somas.tolist()))^set(map(tuple,triaged_somas.tolist()))))

	# sort by volume	
	non_common = non_common[non_common[:,3].argsort()]

	np.save(os.path.join(brain_path,'close_somas.npy'),non_common)
	np.savetxt(os.path.join(brain_path,'close_somas.csv'),non_common,delimiter=",",fmt="%s")
	
