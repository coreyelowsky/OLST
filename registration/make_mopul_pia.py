import numpy as np
from skimage.io import imread
import tifffile as tif

mopul_path = '/data/elowsky/OLST/registration/MOpul_Layers_528.tif'

labels = imread(mopul_path)

layer_1 = np.where(labels == 1)

for z,y,x in zip(layer_1[0],layer_1[1],layer_1[2]):
	
	#above is zero
	if labels[z,y-1,x] == 0:
		# below is 1
		if labels[z,y+1,x] == 1:
			labels[z,y-1,x] = 100

print('Write tif')
tif.imwrite('/data/elowsky/OLST/registration/MOpul_Layers_528_with_pia.tif',labels)
		
