import numpy as np
from skimage.io import imread, imsave
from os import listdir
from os.path import join
from sklearn.metrics import pairwise_distances

##################################################################
CCF_PATH = '/data/elowsky/OLST/registration/MOpul_Layers_528_with_pia.tif'
SWC_BASE_PATH = '/data/elowsky/OLST/swc_analysis/automatically_traced/flagship/'
LAYERS = [2,5]
RES_SCALE_FACTOR = 25
PIA_INTENSITY = 100
##################################################################

# load ccf
mop_labels = imread(CCF_PATH)

# get pia
pia_voxels = np.argwhere(mop_labels == 100)

out_distances = []

olga_id = 1

for layer in LAYERS:

	print()
	print('Layer:',layer)
	print()

	layer_path = join(SWC_BASE_PATH,'layer_'+str(layer))
	swc_path = join(layer_path,'registered_smoothed_microns')
	
	swc_names = listdir(swc_path)

	for swc_name in sorted(swc_names):

		swc_id = swc_name[:-4]
		print(swc_id)

		swc = np.genfromtxt(join(swc_path,swc_name))

		soma = (swc[0,2:5]/RES_SCALE_FACTOR)[::-1].round().astype(int)

		# get pairwise distances between pia and soma and find minimum
		distances = pairwise_distances(soma.reshape(1,-1),pia_voxels)
		min_dist = np.round(np.min(distances)*RES_SCALE_FACTOR,1)


		out_distances.append([olga_id,swc_id,min_dist])

		olga_id += 1

out_distances = np.array(out_distances,dtype=object)
print(out_distances)
np.savetxt(join(SWC_BASE_PATH,'soma_distances_to_pia_layers_2_5.csv'),out_distances,delimiter=',',fmt=['%d','%s','%g'])
	

	

