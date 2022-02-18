from os import listdir, remove
from os.path import join
from skimage.io import imread
from scipy.ndimage import gaussian_filter
import tifffile as tif
import numpy as np
import shutil
from SWC import SWC

def a():
	path = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/registered/'
	outpath = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/test.tif'
	downsampling = 5
	referenceImage= np.zeros(shape=(2640,1600,2280), dtype=np.uint8)

	swcNames = listdir(path)


	for id_, swcName in enumerate(swcNames,start=1):

		swc = SWC(join(path, swcName))

	
		swcArray = (swc.generateSWCArray(onlyCoords=True)/downsampling).round().astype(int)

	
		referenceImage[swcArray[:,2],swcArray[:,1],swcArray[:,0]] = id_

	# save
	tif.imsave(outpath, referenceImage)

def move_swc():

	path =  '/data/elowsky/OLST/swc_analysis/automatically_traced/motor_cortex/registered/'
	out_path = '/data/elowsky/OLST/swc_analysis/automatically_traced/motor_cortex/for_uygar_6_9_2021/'

	swc_names = listdir(path)

	for swc_name in swc_names:
	
		layer = swc_name[0]
		
		swc_name_no_layer = '_'.join(swc_name.split('_')[1:])
		
		source = join(path,swc_name)
		target = join(out_path,'layer_' + layer,swc_name_no_layer)

		shutil.copyfile(source,target)

def create_markup_csv():

	path = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/'
	out_data = []

	for i in range(5):
		
		norm_path = join(path, 'normalized_oblique_' + str(i))
		filenames = listdir(norm_path)
		for filename in filenames:
			out_data.append([filename[:-4], i])

	out_data = np.array(out_data)
	print(out_data)
	np.savetxt(join(path, 'initial_markup.csv'), out_data, delimiter=',', fmt=['%s','%s'])

def create_markup_csv_gcut():

	path = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/'
	out_data = []

	for i in [0,1,2,4]:
		
		norm_path = join(path, 'normalized_oblique_gcut_' + str(i))
		filenames = listdir(norm_path)
		for filename in filenames:
			out_data.append([filename[:-4], i])

	out_data = np.array(out_data)
	print(out_data)
	np.savetxt(join(path, 'initial_markup_gcut.csv'), out_data, delimiter=',',fmt=['%s','%s'])


if __name__ == '__main__':
	
	create_markup_csv_gcut()
	




	

