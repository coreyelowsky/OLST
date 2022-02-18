from SWC import SWC

import numpy as np
import tifffile as tif
import matplotlib.pyplot as plt

from os import listdir, rename, mkdir
from os.path import join, exists
from shutil import copyfile

import re

def plot_swc_clusters():

	# plot swcs from cluster

	FONT_SIZE=7
	BASAL_COLOR='red'
	APICAL_COLOR='blue'
	DOT_SIZE=.1
	X_SPACING = .1
	REGION = 'barrel_cortex'
	SWC_PATH = join('/data/elowsky/OLST/swc_analysis/analysis/', REGION, 'pia_white_matter_normalized/swcs/removed_nan_nodes/base/')
	OUT_PATH = join('/data/elowsky/OLST/swc_analysis/analysis/', REGION, 'pia_white_matter_normalized/analysis/removed_nan_nodes/scaled_1000/soma_centered/')
	#CLUSTER_FILENAME = 'swcs_by_cluster_all_5-15_mc4.csv'
	#CLUSTER_FILENAME = 'swcs_by_cluster_all_5-25_mc4.csv'
	CLUSTER_FILENAME = 'swcs_by_cluster_all_2-10_mc4.csv'
	#CLUSTER_FILENAME = 'swcs_by_cluster_all_2-7_mc4.csv'
	CLUSTER_PATH = join('/data/palmer/consensus_clustering/output/', REGION)
	DATA_TYPES = ['morphometrics']
	#LAYERS = ['all']
	LAYERS = ['4']
	NORM_TYPES = ['pca']
	STRUCTURE_TYPES = ['basal_apical', 'basal_apical_concat']

	for data_type in DATA_TYPES:
		
		print('Data Type:', data_type)
		
		cluster_path_data_type = join(CLUSTER_PATH, data_type)
		out_path_data_type = join(OUT_PATH, data_type)

		for layer in LAYERS:

			print('Layer:', layer)

			cluster_path_layer = join(cluster_path_data_type, 'layer_' + layer)
			out_path_layer = join(out_path_data_type, 'layer_' + layer, 'consensus_clustering')
			if not exists(out_path_layer):
				mkdir(out_path_layer)

			for norm_type in NORM_TYPES:

				print('Normalization Type:', norm_type)		

				cluster_path_norm = join(cluster_path_layer, norm_type)
				out_path_norm = join(out_path_layer, norm_type)
				if not exists(out_path_norm):
					mkdir(out_path_norm)

				for structure_type in STRUCTURE_TYPES:

					print('Structure Type:', structure_type)

					out_path_structure = join(out_path_norm, structure_type)
					if not exists(out_path_structure):
						mkdir(out_path_structure)
			
					cluster_path_structure = join(cluster_path_norm, structure_type, CLUSTER_FILENAME)

					cluster_info = np.genfromtxt(cluster_path_structure, delimiter=' ', dtype=object)
					cluster_info[:,0] = cluster_info[:,0].astype(str)
					cluster_info[:,1] = cluster_info[:,1].astype(str)
					cluster_info[:,2] = cluster_info[:,2].astype(int)


					# iterate through cluster
					for cluster_id in np.unique(cluster_info[:, 2]):
	
						print('Cluster ID:', cluster_id)
		
						# get all swcs in cluster
						swc_ids_in_cluster = cluster_info[cluster_info[:,2] == cluster_id][:,:2]

						# instaniate figure
						fig = plt.figure(figsize=(18,10))
						ax = plt.gca()
						ax.set_title('Layer: '+ layer + '\nNormalization Type: ' + norm_type + '\nStructure Type: ' + structure_type + '\nCluster: ' + str(cluster_id))
						ax.axes.xaxis.set_visible(False)
						ax.axes.yaxis.set_visible(False)
						ax.axis('off')

						# pia and white matter lines
						plt.axhline(y=.5, color='black', linestyle='-')
						plt.axhline(y=-.5, color='black', linestyle='-')
		
						# iterate through swcs to plot
						x_pos = 0
						for i, swc_id_list in enumerate(swc_ids_in_cluster):
			
							swc_id = '_'.join(swc_id_list[::-1])
							swc_full_path = join(SWC_PATH, swc_id + '.swc')

							# load swc
							swc = SWC(swc_full_path)

							# extract only basal and apical
							swc.extract_type_tree(['basal dendrite', 'apical dendrite'])
			
							# generate swc array
							swc_array = swc.generate_swc_array()
							basal_array = swc_array[swc_array[:,1] == 3, 2:5]
							apical_array = swc_array[swc_array[:,1] == 4, 2:5]
							plt.scatter(basal_array[:,0] + x_pos - swc_array[:,2].min() , -basal_array[:,1], c=BASAL_COLOR, s=DOT_SIZE)
							plt.scatter(apical_array[:,0] + x_pos - swc_array[:,2].min() , -apical_array[:,1], c=APICAL_COLOR, s=DOT_SIZE)

							# place swc name
							if i % 2 == 0:
								plt.text(x_pos, .51, swc_id, fontsize=FONT_SIZE, fontweight='bold')
							else:
								plt.text(x_pos, -.53, swc_id, fontsize=FONT_SIZE, fontweight='bold')
			
							# update x position
							x_pos += swc_array[:,2].max() - swc_array[:,2].min() + X_SPACING

			
						out_path = join(out_path_structure, 'cluster_' + str(cluster_id) + '.png')
						plt.savefig(out_path)
						#plt.show()




def give_soma_radius_batch(in_path, out_path, radius):

	swc_names = listdir(in_path)

	for swc_name in swc_names:
		
		swc_full_path = join(in_path, swc_name)
		
		# load swc but dont build tree
		swc = SWC(swc_full_path, build_tree=False)

		# set soma radius
		swc.swc_array[0,-2] = radius

		swc.save_swc(out_path, from_array=True)

def scale_swcs(in_path, out_path, scale_factor):


	swc_names = listdir(in_path)

	for swc_name in swc_names:
		
		swc_full_path = join(in_path, swc_name)
		
		# load swc but dont build tree
		swc = SWC(swc_full_path, build_tree=False)

		# set soma radius
		swc.swc_array[:,2:5] *= scale_factor 

		swc.save_swc(out_path, from_array=True)


def remove_nan_nodes(in_path, out_path, nan_swc_names):

	swc_names = ['2_171012_34','2_180926_2','2_190123_14', '5_190123_11']
	
	for nan_swc_name in nan_swc_names:
	
		swc_path = join(in_path, nan_swc_name + '.swc')
		swc = SWC(swc_path)
		swc.remove_nan_nodes()
		swc.save_swc(out_path)

def soma_center_swcs(in_path, out_path):

	swc_names = listdir(in_path)

	for swc_name in swc_names:
		
		# load swc
		swc_full_path = join(in_path, swc_name)
		swc = SWC(swc_full_path)

		# center around soma
		swc.center_around_soma()

		# save
		swc.save_swc(out_path)

def rename_persistent_vectors():

	path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/removed_nan_nodes/persistent_homology/vectors/'
	dirs = ['basal', 'apical']

	for dir_ in dirs:
		
		full_dir = join(path, dir_)
		
		file_names = listdir(full_dir)
		
		for file_name in file_names:

			new_name = '_'.join(file_name.split('_')[1:])
			
			source = join(full_dir, file_name)
			dest = join(full_dir, new_name)

			rename(source, dest)

def concat_persistent_vectors():

	path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/removed_nan_nodes/persistent_homology/vectors/'

	apical_path = join(path, 'apical')
	basal_path = join(path, 'basal')
	out_path = join(path, 'basal_apical_concat')

	file_names = listdir(basal_path)

	for file_name in sorted(file_names):
		
		# load apical and basal vector
		basal_vector = np.genfromtxt(join(basal_path, file_name))
		apical_vector = np.genfromtxt(join(apical_path, file_name))
		
		# concatenate
		basal_apical_concat_vector = np.concatenate((basal_vector, apical_vector)).reshape(-1,1)

		np.savetxt(join(out_path,file_name), basal_apical_concat_vector)

def find_regions_of_old_barrel_cortex():
	
	old_barrel_path = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex_old_annotation/registered/'
	soma_base_path = '/data/elowsky/OLST/soma_detection/'
	labels_base_path = '/data/anarasim/data/registration_images/forArun/'

	soma_labels = []

	for old_barrel_file_name in sorted(listdir(old_barrel_path)):
		break		

		old_barrel_file_name = '_'.join(old_barrel_file_name.split('_')[1:])
		
		brain = old_barrel_file_name.split('_')[0]
		soma_id = old_barrel_file_name[:-4].split('_')[1]

		soma_brain_path = join(soma_base_path, brain, 'regions', 'barrel_cortex_old_annotation.csv')
		labels_path = join(labels_base_path, brain+'_Emxcre_25um_isotropic_removedStripes_crop_ARA.tif') 

		# load somas
		somas = np.genfromtxt(soma_brain_path, delimiter=',', dtype=object)

		soma_row = somas[int(soma_id)-1]

		x,y,z = soma_row[1:4].astype(int)

		# load labels
		labels = tif.imread(labels_path)

		soma_label = labels[z,y,x]
		
		soma_labels.append(soma_label)

	#soma_labels = np.unique(soma_labels)
	#soma_labels.sort()
	
	soma_labels = [0,9,201,233,241,308,346,532,625,635,670,683,686,793,854,862,865,921,981,1006,1038,1047,1070,1086,1111]

	regions = load_annotation()


	for soma_label in soma_labels:

		if soma_label == 0:
			print('Out of Brain')
		else:
			print(regions[np.where(regions[:,0] == soma_label)[0][0],1])





def load_annotation():

	path = '/data/elowsky/OLST/registration/new_annotation_from_rmc.csv'
	regions = []
		
	with open(path, 'r') as file:
		lines = file.readlines()
		for i,line in enumerate(lines):
			if i > 0:
	
				line = re.split(r',(?=")', line)
				id_ = line[0]
				name = line[1][1:-1]
			

				regions.append([int(id_), name])	

	regions = np.array(regions, dtype=object)
	regions[:,0] = regions[:,0].astype(np.uint32)
	regions[:,1] = regions[:,1].astype(str)
	
	return regions



		
def pool_pfc():
	
	pfc_path = '/data/elowsky/OLST/swc_analysis/automatically_traced/pfc/'
	out_path = join(pfc_path, 'all_registered')

	dir_names = ['anterior_cingulate', 'limbic_cortex', 'orbital_cortex']

	for dir_name in dir_names:

		swcs_path = join(pfc_path, dir_name,'registered')
		swc_names = listdir(swcs_path)

		for swc_name in swc_names:
			
			pooled_name = dir_name[0] + '_' + swc_name
			
			source = join(swcs_path, swc_name)
			dest = join(out_path, pooled_name)

			copyfile(source, dest)
	

	



	
		


	

if __name__ == '__main__':
	
	#in_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/with_soma_radius/removed_nan_swcs/all_structures/'
	#out_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/with_soma_radius/removed_nan_swcs/all_structures_/'
	#radius = 1
	#give_soma_radius_batch(in_path, out_path, radius)

	#in_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/all_structures_not_modified/'
	#out_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/all_structures_scaled_by_1000/'
	#scale_factor = 1000
	##scale_swcs(in_path, out_path, scale_factor)

	#in_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/all_structures_scaled_by_1000/'
	#out_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/removed_nan_nodes/all_structures/'
	#nan_swc_names = ['2_171012_34','2_180926_2','2_190123_14', '5_190123_11']
	#remove_nan_nodes(in_path, out_path, nan_swc_names)


	#in_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/swcs/removed_nan_nodes/apical/'
	#out_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/swcs/removed_nan_nodes/soma_centered/apical/'
	#soma_center_swcs(in_path, out_path)

	#rename_persistent_vectors()

	#concat_persistent_vectors()

	#find_regions_of_old_barrel_cortex()

	#load_annotation()

	#pool_pfc()

	plot_swc_clusters()



