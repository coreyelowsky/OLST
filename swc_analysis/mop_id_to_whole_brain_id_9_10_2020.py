import os
import numpy as np


###########################################################
BRAINS = ['190522', '190327']
MOP_IDS = {'190522':[4, 32, 16, 23, 27, 33], '190327':[42, 21, 18, 31, 45, 46]}
SOMAS_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
############################################################

for brain in BRAINS:

	print('Brain:',brain)

	somas_brain_path = os.path.join(SOMAS_BASE_PATH, brain)
	mop_somas_path = os.path.join(somas_brain_path, 'motor_cortex_somas.npy')
	whole_brain_somas_path = os.path.join(somas_brain_path, 'transformed_somas.npy')

	mop_somas = np.load(mop_somas_path, allow_pickle=True)
	whole_brain_somas = np.load(whole_brain_somas_path, allow_pickle=True)

	brain_mop_ids = MOP_IDS[brain]

	for mop_id in brain_mop_ids:

		mop_soma = mop_somas[mop_id-1][4:8]	
		whole_brain_id = np.where(np.all(whole_brain_somas[:,3:7] == mop_soma,axis=1))[0][0] + 1

		print('MOP:', mop_id, 'Whole Brain:', whole_brain_id)




