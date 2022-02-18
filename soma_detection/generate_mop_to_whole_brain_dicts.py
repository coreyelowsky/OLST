from BRAINS import BRAINS
from os.path import join, exists
import numpy as np


SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'

for brain in BRAINS:

	print()
	print('Brain:', brain)
	print()

	brainPath = join(SOMA_BASE_PATH, brain)

	wholeBrainPath = join(brainPath,'whole_brain_isocortex_somas.npy')
	mopPath = join(brainPath, 'motor_cortex_somas.npy')

	if not exists(wholeBrainPath) or not exists(mopPath):
		print('Soma file doesnt not exist!')
		continue

	mopToWholeBrainDict = {}

	wholeBrainSomas = np.load(wholeBrainPath, allow_pickle=True)
	mopSomas = np.load(mopPath,allow_pickle=True)

	# modify mop somas to remove stitching coord columns
	mopSomas = mopSomas[:,[4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28]]

	for somaID, soma in enumerate(mopSomas, start=1):
		
		wholeBrainID = np.where(np.all(wholeBrainSomas[:,1:] == soma,axis=1))[0][0] + 1
	
		mopToWholeBrainDict[somaID] = wholeBrainID

	outPath = join(brainPath, 'mopToWholeBrainDict.npy')
	np.save(outPath, mopToWholeBrainDict)
