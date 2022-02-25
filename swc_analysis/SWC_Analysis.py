from SWC import SWC
from os.path import join, exists, isfile
from os import listdir, mkdir, system, rename, remove
from sys import exit
from shutil import copyfile
from scipy.stats import zscore
from sklearn.decomposition import PCA
import scipy.cluster.hierarchy as shc
from dynamicTreeCut import cutreeHybrid
from scipy.spatial.distance import pdist
import numpy as np
from scipy.cluster.hierarchy import linkage
from skimage.io import imread, imsave
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import kruskal, ranksums, mannwhitneyu, wilcoxon
import matplotlib as mpl
import seaborn as sns
from scipy.spatial import distance_matrix
	

import pprint
import pandas as pd
import json
import params
import utils
import numpy as np
import filecmp
import pickle
from pylab import *
import Clusters

class SWC_Analysis():

	def __init__(self, region, remove_nan_type, scale_type, soma_type, swc_type):

		# base path for analysis
		self.analysis_base_path = '/data/elowsky/OLST/swc_analysis/analysis/'

		# info for type of dataset
		self.region = region
		self.remove_nan_type = remove_nan_type
		self.scale_type = scale_type
		self.soma_type = soma_type
		self.swc_type = swc_type

		# set paths
		self.region_path = join(self.analysis_base_path, self.region)
		self.swc_type_path = join(self.region_path, self.swc_type)
		self.swc_path = join(self.swc_type_path, 'swcs')
		self.analysis_path = join(self.swc_type_path, 'analysis', remove_nan_type, scale_type, soma_type)
		self.morphometrics_json_path = join(self.analysis_base_path, 'lmeasure_morphometrics.json')

		# layers
		self.layers_map = {
			'mop': ['2', '5', '6', 'all'],
			'barrel_cortex': ['2', '4', '5', '6', 'all'],
			'pfc': ['2','5','6']
		}

		self.layers = self.layers_map[self.region]

		# structures
		self.structures = ['basal', 'apical', 'basal_apical']
		self.analysis_structures = self.structures + ['basal_apical_concat']


		# cluster file names
		self.cluster_file_name_dict = {

			'mop': {

				'basal_apical':{},
				'basal_apical_concat':{
					'2': 'swcs_by_cluster_all_2-10_mc4.csv',
					'5': 'swcs_by_cluster_all_2-10_mc4.csv',
					'6': 'swcs_by_cluster_all_2-10_mc4.csv',
					'all': 'swcs_by_cluster_all_5-30_mc4.csv'
				}
			},

			'barrel_cortex':{

				'basal_apical':{},
				'basal_apical_concat':{
					'2': 'swcs_by_cluster_all_2-10_mc4.csv',
					'4': 'swcs_by_cluster_all_2-10_mc4.csv',
					'5': 'swcs_by_cluster_all_2-10_mc4.csv',
					'6': 'swcs_by_cluster_all_2-10_mc4.csv',
					'all': 'swcs_by_cluster_all_5-25_mc4.csv'
				}
			
			}
		}


		# other helpful variables
		self.remove_nan_types = ['removed_nan_nodes', 'removed_nan_swcs']

	@staticmethod
	def check_if_two_folders_are_same(path_a, path_b):

		"""
		checks if 2 paths have the exact same files in them
		
			- number of files must be same
			- file contents must be same
			- does not check directories inside paths
		"""

		path_a_files = [f for f in listdir(path_a) if isfile(join(path_a, f))]
		path_b_files = [f for f in listdir(path_b) if isfile(join(path_b, f))]

		# check if folder contain same number of files
		if len(path_a_files) != len(path_b_files):
			exit('Folders do not have the same number of files....')

		# iterate through every file and check if same
		path_a_files = [f for f in listdir(path_a) if isfile(join(path_a, f))]

		for path_a_file_name in path_a_files:
			
			# check if exists in b
			if not exists(join(path_b, path_a_file_name)):
				exit(path_a_file_name + ' does not exist in path b')

			# check content
			if not filecmp.cmp(join(path_a, path_a_file_name), join(path_b, path_a_file_name)):
				exit(path_a_file_name + ' is different in both paths')

		print('Paths contain the same files!')
				

			
	def consolidate_swcs_from_olga(self):
		
		"""
		consolidate swcs from olga into one directory
			
			- looks to see if swcs were sent in batches
			- assumes swcs were sent in separate layers

		"""

		
		from_olga_path = join(self.swc_type_path, 'swcs', 'from_olga')
		consolidated_path = join(self.swc_type_path, 'swcs', 'from_olga_consolidated')

		# make out path if doesnt exist
		if not exists(consolidated_path):
			mkdir(consolidated_path)

		# check if came in batches
		if any(['batch' in x for x in listdir(from_olga_path)]):
			
			# iterate through batches
			for batch_name in sorted(listdir(from_olga_path)):
				
				batch_path = join(from_olga_path, batch_name)	
				
				# iterate through layers
				for layer_name in sorted(listdir(batch_path)):
		
					layer_path = join(batch_path, layer_name)
			
					# extract_layer
					layer = layer_name.split('_')[-1]
					
					for swc_name in sorted(listdir(layer_path)):
						
						swc_name_with_layer = '_'.join([layer, swc_name])
						
						# copy file
						source = join(layer_path, swc_name)
						dest = join(consolidated_path, swc_name_with_layer)

						copyfile(source, dest)
						
		else:
			print('No Batches')

			# check if came in layers
			if any(['layer' in x for x in listdir(from_olga_path)]):
	
				# iterate through layers
				for layer_name in sorted(listdir(from_olga_path)):
	
					layer_path = join(from_olga_path, layer_name)
		
					# extract_layer
					#layer = layer_name.split('_')[-1]
				
					for swc_name in sorted(listdir(layer_path)):
					
						#swc_name_with_layer = '_'.join([layer, swc_name])
					
						# copy file
						source = join(layer_path, swc_name)
						dest = join(consolidated_path, swc_name)

						copyfile(source, dest)
			else:
				print('No Layers')
				
				# just modify name and copy


	def handle_nan_nodes(self):
		
		"""
		deals with nan nodes due to registration errors
			
			- saves 2 separate data sets
				- removed nan nodes
				- removed swc nodes
		"""

		swcs_path = join(self.swc_type_path, 'swcs')
		consolidated_path = join(swcs_path, 'from_olga_consolidated')
		removed_nan_nodes_path = join(swcs_path, 'removed_nan_nodes')
		removed_nan_swcs_path = join(swcs_path, 'removed_nan_swcs')

		# make directories if they dont exist
		if not exists(removed_nan_nodes_path):
			mkdir(removed_nan_nodes_path)
		if not exists(removed_nan_swcs_path):
			mkdir(removed_nan_swcs_path)
		if not exists(join(removed_nan_nodes_path, 'base')):
			mkdir(join(removed_nan_nodes_path, 'base'))
		if not exists(join(removed_nan_swcs_path, 'base')):
			mkdir(join(removed_nan_swcs_path, 'base'))

		
		# figure out which swcs contain nans
		swc_names = listdir(consolidated_path)
		
		# keep track of removed nan nodes
		removed_nan_nodes_list = []

		for swc_name in sorted(swc_names):
			
			# check if contains nan
			# build tree and resave so duplicates will be removed
			swc_path = join(consolidated_path, swc_name)
			swc = SWC(swc_path)
		
			# if contains nan nodes remove nan nodes
			# it not then save to remove nan swcs
			if swc.contains_nan_coords()[0]:
				swc.remove_nan_nodes()
				removed_nan_nodes_list.append(swc_name)
			else:
				swc.save_swc(join(removed_nan_swcs_path,'base'))

			# all swcs are save to removed nan nodes
			swc.save_swc(join(removed_nan_nodes_path, 'base'))

		print()
		print('Removed nan nodes:', removed_nan_nodes_list)


	def scale_swcs(self, scale_factor=1000):

		"""
		scale swcs coordinates

		"""

		for remove_nan_type in self.remove_nan_types:
	
			in_path = join(self.swc_path, remove_nan_type)
			out_path = join(in_path, 'scaled_' + str(scale_factor))
	
			# make out directories
			if not exists(out_path):
				mkdir(out_path)
			if not exists(join(out_path, 'base')):
				mkdir(join(out_path, 'base'))


			# scale removed nan nodes
			for swc_name in sorted(listdir(join(in_path, 'base'))):
			
				swc = SWC(join(in_path, 'base', swc_name))
				swc.scale_coords([scale_factor, scale_factor, scale_factor])
				swc.save_swc(join(out_path, 'base'))

	def soma_center_swcs(self):

		"""
		soma center swcs

		"""
	
		for remove_nan_type in self.remove_nan_types:
	
			in_path = join(self.swc_path, remove_nan_type, 'scaled_1000')
			out_path = join(in_path, 'soma_centered')

	
			# make out directories
			if not exists(out_path):
				mkdir(out_path)
			if not exists(join(out_path, 'base')):
				mkdir(join(out_path, 'base'))


			# soma center
			for swc_name in sorted(listdir(join(in_path, 'base'))):
			
				swc = SWC(join(in_path, 'base', swc_name))
				swc.center_around_soma()
				swc.save_swc(join(out_path, 'base'))

	def extract_branch_types(self, endpoint_path):

		"""
		extract basal, apical, and basal_apical

		"""

		structure_types = ['basal', 'apical', 'basal_apical']
		structure_dict = {
				'basal':['basal dendrite'],
				'apical':['apical dendrite'],
				'basal_apical':['basal dendrite', 'apical dendrite'],
				}

		for remove_nan_type in self.remove_nan_types:
		
			in_path = join(self.swc_path, remove_nan_type, endpoint_path)
			swcs_path = join(in_path, 'base')

			for structure in structure_types:
			
				out_path = join(in_path, structure)

				if not exists(out_path):
					mkdir(out_path)

				# extract from swcs
				for swc_name in sorted(listdir(swcs_path)):
	
					swc = SWC(join(swcs_path, swc_name))
					swc.extract_type_tree(structure_dict[structure])
					swc.save_swc(out_path)
			
	def run_lmeasure_morphometrics(self):

		"""
		runs lmeasure to calculate Morphometrics	

		"""
	
		# set up paths
		swc_in_path = join(self.swc_path, self.remove_nan_type, self.scale_type, self.soma_type)
		morphometrics_out_path = join(self.swc_type_path, 'morphometrics', self.remove_nan_type, self.scale_type, self.soma_type)

		# structure types
		structure_types = ['basal', 'apical', 'basal_apical']

		# load morphometrics json
		with open(self.morphometrics_json_path,'r') as fp:
			morphometrics_dict = json.load(fp)

		# sort morphometric names
		morphometric_names = sorted(morphometrics_dict)

		# create function specifier for lmeasure
		function_specifier = ' '.join(['-f' + str(params.LMEASURE_FUNCTION_IDS[m])+',0,0,10.0' for m in morphometric_names])

		# run lmeasure to obtain morphometrics for all structure types
		for structure_type in structure_types:

			# set up paths
			swc_path_full = join(swc_in_path, structure_type)
			swc_paths = ' '.join([join(swc_path_full, f) for f in sorted(listdir(swc_path_full))])

			# form full lmeasure command
			out_path_full = join(morphometrics_out_path, 'lmout_morphometrics_' + structure_type + '.txt')
			lmeasure_command = ' '.join([params.LMEASURE_EXE, function_specifier, '-s' + out_path_full, swc_paths])

			# run lmeasure from command line
			system(lmeasure_command)


	def parse_lmeasure_morphometrics(self):
	
		"""
		parse lmeasure morphometrics		

		"""

		morphometrics_path = join(self.swc_type_path, 'morphometrics', self.remove_nan_type, self.scale_type, self.soma_type)
		out_path_base = join(self.analysis_path, 'morphometrics')

		# structure types
		structure_types = ['basal', 'apical', 'basal_apical']
		
		# get layers for region
		layers = self.layers

		# load morphometrics json to decide which morphometrics to calculate
		with open(self.morphometrics_json_path,'r') as fp:
			morpho_dict = json.load(fp)

		# dict to store morphometrics
		morphometrics_dict = {}
		for layer in layers:
			morphometrics_dict[layer] = {}


		# load basal and apical and concatenate
		for structure_type in structure_types:

			print(structure_type)
		
			# save concatrnated basal apical morphometrics
			if structure_type == 'basal_apical':				
		
				# save to separate layers
				for layer in layers:

					out_layer_path = join(out_path_base, 'layer_'+ str(layer))
					if not exists(out_layer_path):
						mkdir(out_layer_path)

					# load and concatenate
					basal_morpho = np.loadtxt(join(out_path_base, 'layer_'+ str(layer), 'basal.txt'),dtype=object)
					apical_morpho = np.loadtxt(join(out_path_base, 'layer_'+ str(layer), 'apical.txt'),dtype=object)

					# append structure type to headers
					basal_headers = ['basal:' + x for x in basal_morpho[0,2:]]
					basal_morpho[0,2:] = basal_headers
					apical_headers = ['apical:' + x for x in apical_morpho[0,2:]]
					apical_morpho[0,2:] = apical_headers
			
					# horizontal stack data
					morphometrics = np.hstack((basal_morpho, apical_morpho[:,2:]))

					out_path =  join(out_layer_path, 'basal_apical_concat.txt')
				

					np.savetxt(out_path, morphometrics, fmt='%s')
			

			# get lmeasure output path
			lmout_morpho_path = join(morphometrics_path, 'lmout_morphometrics_' + structure_type + '.txt')

			# read in data into dictionary
			with open(lmout_morpho_path,'r') as fp:
				lines =  fp.readlines()
				for line in lines:

					line = line.split()

					if len(line):
						swc_name = line[0].split('/')[-1]

						# extract layer and swcID
						swc_info = utils.get_swc_info(swc_name)
						layer = str(swc_info['layer'])
						swc_name = swc_info['name']

						# extract metric
						metric_name = line[1]

						# extract all values to consider associates with metric
						data = []
						for v in morpho_dict[metric_name]:
							index = params.LMEASURE_INDICES[v]
							data.append(float(line[index]))

						# create new dict for swc if needed
						if not swc_name in morphometrics_dict[layer]:
							morphometrics_dict[layer][swc_name] = {}

						# append data to dict
						morphometrics_dict[layer][swc_name][metric_name] = data


			# store data in array
			morphometrics = []

			# create labels row
			labels = ['SWC_id','layer']
			for metric in sorted(morpho_dict):
				for value in morpho_dict[metric]:
					labels.append( metric + ':'+ value)
			morphometrics.append(labels)
		
			# append data
			for layer in sorted(morphometrics_dict):
				for swc_name in sorted(morphometrics_dict[layer]):
					row_data = [swc_name,layer]
					for metric in sorted(morphometrics_dict[layer][swc_name]):
						for item in morphometrics_dict[layer][swc_name][metric]:
							row_data.append(item)
					morphometrics.append(row_data)

			# convert to np array
			morphometrics = np.array(morphometrics,dtype=object)

			# save to separate layers
			for layer in layers:
				out_layer_path = join(out_path_base, 'layer_'+ str(layer))
				if not exists(out_layer_path):
					mkdir(out_layer_path)
				out_path =  join(out_layer_path, structure_type + '.txt')
				
				if layer == 'all':
					morphometrics_out = morphometrics
				else:
					morphometrics_out = morphometrics[morphometrics[:,1] == layer]
					morphometrics_out = np.vstack((morphometrics[0],morphometrics_out))
		
				np.savetxt(out_path, morphometrics_out, fmt='%s')

	def consolidate_arbor_densities(self):

		arbor_density_in_path = join(self.swc_type_path, 'arbor_density', 'from_olga')
		arbor_density_out_path = join(self.swc_type_path, 'arbor_density', 'from_olga_consolidated')
		layers = self.layers[self.region]

		if 'all' in layers:
			layers.remove('all')

		# iterate through layers
		for layer in layers:

			layer_path = join(arbor_density_in_path, 'layer_' + str(layer) + '_arbor_densities', self.soma_type)
	
			for arbor_density_name in sorted(listdir(layer_path)):

				if 'apical' in arbor_density_name:
					structure = 'apical'
				elif 'basal' in arbor_density_name:
					structure = 'basal'
				else:
					structure = 'basal_apical'
			
				arbor_density_name_extracted = arbor_density_name.split('_')[-3:]
				arbor_density_name_extracted = '_'.join(arbor_density_name_extracted)

				source_path = join(layer_path, arbor_density_name)
				dest_path = join(arbor_density_out_path, self.soma_type, structure, arbor_density_name_extracted)

				copyfile(source_path, dest_path)

	def concat_arbor_densities(self):

		ARBOR_DENSITY_SIZE = 4800

		# inpaths
		arbor_density_path = join(self.swc_type_path, 'arbor_density', 'from_olga_consolidated', self.soma_type)

		# layers
		layers = self.layers

		# dict to hold arbor densities
		arbor_densities = {}

		# iterate through structures
		for structure in self.structures:
			
			print(structure)

			structure_in_path = join(arbor_density_path, structure)

			# get sorted list of all neurons
			arbor_density_names = listdir(structure_in_path)

			# first row should be headers
			headers = ['SWC_id', 'layer'] +  [str(i+1) for i in range(ARBOR_DENSITY_SIZE)]
			arbor_densities[structure] = np.array([headers])
						
			for arbor_density_name in sorted(arbor_density_names):
				
				vector_path = join(structure_in_path, arbor_density_name)
				vector = np.loadtxt(vector_path, delimiter=',', dtype=np.float64).flatten()

				# prepend swc id and layers
				swc_info = [arbor_density_name[2:-4], arbor_density_name[0]]
				vector = np.hstack((swc_info, vector))
	
				# concatenate to all vectors
				arbor_densities[structure] = np.vstack([arbor_densities[structure], vector])

		# concat for basal_apical_concat
		arbor_densities['basal_apical_concat'] = np.hstack((arbor_densities['basal'],arbor_densities['apical'][:,2:]))
		arbor_densities['basal_apical_concat'][0] = ['SWC_id', 'layer'] +  [str(i+1) for i in range(ARBOR_DENSITY_SIZE*2)]

		# save

		# iterate through stuctures
		for structure in arbor_densities:
	
			arbor_densities_structure = arbor_densities[structure]
		
			# iteraye through layers
			for layer in layers:
				
				if layer != 'all':
					vectors_out = arbor_densities_structure[arbor_densities_structure[:,1] == layer]
					vectors_out = np.vstack((arbor_densities_structure[0],vectors_out))
				else:
					vectors_out = arbor_densities_structure

				
				out_path = join(self.analysis_path, 'arbor_density', 'layer_' + str(layer), structure + '.txt')
				np.savetxt(out_path, vectors_out, fmt='%s')


	def rename_persistent_vectors(self):


		vectors_path = join(self.swc_type_path, 'persistent_homology', self.remove_nan_type, self.scale_type, self.soma_type, 'persistent_vectors')

		for structure in self.structures:
			
			structure_path = join(vectors_path, structure)

			for file_name in listdir(structure_path):
		
				file_name_fixed = '_'.join(file_name.split('_')[1:])
			
				source = join(structure_path, file_name)
				dest = join(structure_path, file_name_fixed)

				rename(source, dest)

	def concat_persistent_homology_vectors(self):

		# inpaths
		persistent_vectors_path = join(self.swc_type_path, 'persistent_homology', self.remove_nan_type, self.scale_type, self.soma_type, 'persistent_vectors')

		# layers
		layers = self.layers

		# dict to hold arbor densities
		persistent_vectors = {}

		# iterate through structures
		for structure in self.structures:
			
			print(structure)
			
			structure_in_path = join(persistent_vectors_path, structure)

			# get sorted list of all neurons
			vector_names = listdir(structure_in_path)

			# get size by loading first vector
			p_v_size = len(np.loadtxt(join(structure_in_path, vector_names[0]), delimiter=',', dtype=np.float64))			

			# first row should be headers
			headers = ['SWC_id', 'layer'] +  [str(i+1) for i in range(p_v_size)]
			persistent_vectors[structure] = np.array([headers])
						
			for vector_name in sorted(vector_names):
				
				vector_path = join(structure_in_path, vector_name)
				vector = np.loadtxt(vector_path, delimiter=',', dtype=np.float64).flatten()

				# prepend swc id and layers
				swc_info = [vector_name[2:-4], vector_name[0]]
				vector = np.hstack((swc_info, vector))
	
				# concatenate to all vectors
				persistent_vectors[structure] = np.vstack([persistent_vectors[structure], vector])

		# concat for basal_apical_concat
		persistent_vectors['basal_apical_concat'] = np.hstack((persistent_vectors['basal'],persistent_vectors['apical'][:,2:]))
		persistent_vectors['basal_apical_concat'][0] = ['SWC_id', 'layer'] +  [str(i+1) for i in range((persistent_vectors['basal'].shape[1]-2 + persistent_vectors['apical'].shape[1]-2))]

		# save

		# iterate through stuctures
		for structure in persistent_vectors:
	
			persistent_vectors_structure = persistent_vectors[structure]
		
			# iteraye through layers
			for layer in layers:
				
				if layer != 'all':
					vectors_out = persistent_vectors_structure[persistent_vectors_structure[:,1] == layer]
					vectors_out = np.vstack((persistent_vectors_structure[0],vectors_out))
				else:
					vectors_out = persistent_vectors_structure

				
				out_path = join(self.analysis_path, 'persistent_homology', 'layer_' + str(layer), structure + '.txt')
				np.savetxt(out_path, vectors_out, fmt='%s')
			




	def pca(self, data_type, explained_variance_thresh = .95):


		data_path = join(self.analysis_path, data_type)

		layers = self.layers

		if data_type == 'arbor_density' and self.region == 'mop':
			structures = ['basal', 'apical', 'basal_apical_concat']
		else:
			structures = ['basal', 'apical', 'basal_apical', 'basal_apical_concat']	
	
		for layer in layers:
			print('Layer:', layer)
			layer_path = join(data_path, 'layer_' + str(layer))

			for structure in structures:

				print('Structure:', structure)

				# load data
				normalized_data_path = join(layer_path,'normalized', structure + '.txt')
				data = pd.read_csv(normalized_data_path, sep=' ', header=0)

				# only run pca on data
				X = data[data.columns[2:]]
	
				# pca object
				if explained_variance_thresh >=1:
					exit('Error: variance threshold must be 0 < x < 1')
				pca = PCA(explained_variance_thresh)

				# fit pca
				pca.fit(X)

				num_components = len(pca.explained_variance_ratio_)
				print('# Features:', X.shape[1])
				print('# Features (after PCA):', num_components)

				# transform data
				y = pca.transform(X)

				print('Output Shape:', y.shape)
		
				# append columns and headers
				y = np.hstack((np.array(data['SWC_id']).reshape(-1,1),np.array(data['layer']).reshape(-1,1),y))
				headers = ['SWC_id', 'layers']
				pc_headers = ['pc_' + str(i) for i in range(1,num_components+1)]
				y = np.vstack((headers+pc_headers,y))
		
				# save vectors
				out_path = join(layer_path, 'pca', structure + '.txt')
				if not exists(join(layer_path, 'pca')):
					mkdir(join(layer_path, 'pca'))
				np.savetxt(out_path, y, delimiter=' ', fmt='%s')


	def pool_data(pca=False, explained_variance_thresh=.95):

		base_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/removed_nan_nodes/soma_centered/'
	
		# load data
		arbor_density_data = load_data('arbor_density')
		morphometrics_data = load_data('morphometrics')
		persistent_homology_data = load_data('persistent_homology')

		# concat data
		pooled_data = np.hstack((arbor_density_data, morphometrics_data, persistent_homology_data))
		print(pooled_data.shape)

		# save pooled data
		np.save(join(base_path, 'pooled_data_raw'), pooled_data)

		# remove features with 0 variance and normalize data
		pooled_data = pooled_data[:, ~(pooled_data == pooled_data[0,:]).all(0)]
		pooled_data_normalized = zscore(pooled_data)
		print(pooled_data_normalized.shape)
	
		# save normazlied data
		np.save(join(base_path, 'pooled_data_normalized'), pooled_data_normalized)

		# run pca
		if pca:
			# pca object
			if explained_variance_thresh >=1:
				exit('Error: variance threshold must be 0 < x < 1')
			pca = PCA(explained_variance_thresh)

			# fit pca
			pca.fit(pooled_data_normalized)

			num_components = len(pca.explained_variance_ratio_)
			print('# Features:', pooled_data_normalized.shape[1])
			print('# Features (after PCA):', num_components)

			# transform data
			y = pca.transform(pooled_data_normalized)

			print('Output Shape:', y.shape)

			# save vectors
			np.save(join(base_path, 'pooled_data_pca_reduced'), y)
	

	def hierarchical_clustering(self, plot_dendro=False, save_dendro=False, cut_method=None):


		base_path = self.analysis_path
		data_types = ['morphometrics', 'arbor_density', 'persistent_homology']
		norm_types = ['pca', 'normalized']

		for data_type in data_types:
			print(data_type)
			
			data_type_path = join(base_path, data_type)

			for layer in self.layers:
				print(layer)

				layer_path = join(data_type_path, 'layer_' + str(layer))
				out_layer_path = join(layer_path, 'hierarchical_clustering')
				if not exists(out_layer_path):
					mkdir(out_layer_path)

				for norm_type in norm_types:
					
					print(norm_type)
					
					norm_path = join(layer_path, norm_type)

					out_norm_path = join(out_layer_path, norm_type)
					if not exists(out_norm_path):
						mkdir(out_norm_path)

					for structure in self.analysis_structures:
	
						print(structure)

						data_path = join(norm_path, structure + '.txt')

						# make sure data exists		
						if exists(data_path):

	
							# load data
							data = pd.read_csv(data_path, sep=' ', header=0)
		
							# only on data
							X = data[data.columns[2:]]
							num_data = len(X)

							# calculate pairwise distances
							distances = pdist(X, 'euclidean')

							# use wards algorithm
							link = linkage(distances, method='ward', metric='euclidean')	

							# cut tree to form clusters
							if cut_method ==' cut_tree_hybrid':
								clusters = cutreeHybrid(link, distances, minClusterSize=min_cluster_size)
								cluster_labels= clusters['labels']

								colors = ['k','b', 'g', 'r', 'c', 'm', 'y']
								if len(np.unique(cluster_labels)) > len(colors):
									exit('Error: More clusters than colors')


								# for all links calculate colors
								# key will be link id
								# value will be color

								link_colors = {}
								cluster_id = num_data
								for link_row in link:
			
									cluster_id_a = int(link_row[0])
									cluster_id_b = int(link_row[1])

									# if either of cluster are a specific data point
									# then just color based on cluster of that point
									if cluster_id_a < num_data:
										link_colors[cluster_id] = colors[cluster_labels[cluster_id_a]]
									elif cluster_id_b < num_data:
										link_colors[cluster_id] = colors[cluster_labels[cluster_id_b]]

									else:
										# neither of the clusters are original data points
										# so have to color from dictionary
										# but make sure left and right are same
										# if not just color black
										if link_colors[cluster_id_a] == link_colors[cluster_id_b]:
											link_colors[cluster_id] = link_colors[cluster_id_a]
										else:
											link_colors[cluster_id] = 'k'

									cluster_id += 1

							elif cut_method == 'fcluster':
								
								out_path_struct = join(out_norm_path, structure)
								if not exists(out_path_struct):
									mkdir(out_path_struct)

								# iterate through many max clusters
								for max_clusters in range(5,30+1):
									clusters = shc.fcluster(link, t=max_clusters, criterion='maxclust')-1
									out_path = join(out_path_struct, 'fcluster_max_clust_' + str(max_clusters) + '.csv')
									np.savetxt(out_path, clusters, fmt='%d', delimiter='\n')


	

							# instantiate plot
							plt.figure(figsize=(10, 7))
							labels = np.arange(1, len(data)+1)
				
							if cut_method ==' cut_tree_hybrid':
								dend = shc.dendrogram(link, labels=labels, link_color_func=lambda k:link_colors[k])
							else:
								dend = shc.dendrogram(link, labels=labels, color_threshold=-1, link_color_func=lambda k:'black')

							# plot ltitle
							plt.title(data_path, fontsize=6)

							if plot_dendro:			
								plt.show()

							if save_dendro:
								out_path = join(out_norm_path, structure + '.png')
								plt.savefig(out_path)

							plt.close()
						else:
							print('DATA DOES NOT EXIST')
	


	def normalize_data(self, data_type):


		data_base_path = join(self.analysis_path, data_type)
		
		layers = self.layers
		
		if data_type == 'arbor_density' and self.region == 'mop':
			structure_types = ['basal', 'apical', 'basal_apical_concat']
		else:
			structure_types = ['basal', 'apical', 'basal_apical', 'basal_apical_concat']

		for layer in layers:

			print('Layer:', layer)			

			layer_path = join(data_base_path, 'layer_' + str(layer))
		
			for structure in structure_types:

				print('Structure:', structure)

				data_path = join(layer_path, structure + '.txt')

				# load data
				data = pd.read_csv(data_path, sep=' ', header=0)
				

		
				# remove columns where all data is same	
				# ignoring first 2 columns
				nunique = data.nunique()[2:]
				cols_to_drop = nunique[nunique == 1].index
				data = data.drop(cols_to_drop, axis=1)
				print('Removed Columns:', cols_to_drop)
		
				# z score data columns
				morpho_cols = data.columns[2:]

				data[morpho_cols] = data[morpho_cols].apply(zscore)
		
				# save
				out_path = join(layer_path, 'normalized', structure + '.txt')
				if not exists(join(layer_path, 'normalized')):
					mkdir(join(layer_path, 'normalized'))
				data.to_csv(out_path, sep=' ', index=False)


	def create_2d_swc_plots(self):

		# swc path
		swc_in_path = join(self.swc_type_path, 'swcs', self.remove_nan_type, self.scale_type, self.soma_type, 'base')
		out_path = join(self.analysis_path, '2d_swc_plots')


		# get all swc names
		swc_names = listdir(swc_in_path)


		# find x lims
		min_x_list, max_x_list = [], []
		for swc_name in swc_names:
			swc_full_path = join(swc_in_path, swc_name)
			swc = SWC(swc_full_path, build_tree=False)
			min_x_list.append(swc.swc_array[:,2].min())
			max_x_list.append(swc.swc_array[:,2].max())
		x_lim = [np.min(min_x_list) - 10, np.max(max_x_list) + 10]
		print(x_lim)

		# find y limes
		min_y_list, max_y_list = [], []
		for swc_name in swc_names:
			swc_full_path = join(swc_in_path, swc_name)
			swc = SWC(swc_full_path, build_tree=False)

			min_y_list.append(swc.swc_array[:,3].min())
			max_y_list.append(swc.swc_array[:,3].max())
		y_lim = [-np.max(max_y_list) - 10, -np.min(min_y_list) + 10]
		print(y_lim)


		for swc_name in swc_names:

			swc_full_path = join(swc_in_path, swc_name)
			out_path_full = join(out_path, swc_name[:-4] + '.png')
	
			swc = SWC(swc_full_path, build_tree=False)

			swc.plot(save=True, out_path=out_path_full, x_lim=x_lim, y_lim=y_lim, from_array=True)

	def arbor_densities_to_binary_images(self):

		in_path = join(self.swc_type_path,'arbor_density', 'from_olga_consolidated')
		soma_types = ['soma_centered']
		out_path = join(self.swc_type_path,'arbor_density', 'binary_images')
		structures = ['basal', 'apical', 'basal_apical']
	

		for soma_type in soma_types:
		
			print(soma_type)
			in_path_soma_type = join(in_path, soma_type)
			out_path_soma_type = join(out_path, soma_type)
		
			for structure in structures:
				
				in_path_structure = join(in_path_soma_type, structure)
				file_names = listdir(in_path_structure)
	
				for file_name in file_names:

					file_path = join(in_path_structure, file_name)

					# load arbor density
					arbor_density = np.loadtxt(file_path, delimiter=',', dtype=np.float64)

					# make binary for 8 bit image
					arbor_density[arbor_density > 0] = 255
					arbor_density = arbor_density.astype(np.uint8)

					# save as image
					out_path_full = join(out_path_soma_type, structure, file_name[:-4] + '.tif')
					imsave(out_path_full, arbor_density)

	def plot_swc_clusters(self, cluster_path, data_type):

		# plot swcs from cluster

		FONT_SIZE=5
		BASAL_COLOR='red'
		APICAL_COLOR='blue'
		DOT_SIZE=.1
		X_SPACING = .1
		swc_path = join('/data/elowsky/OLST/swc_analysis/analysis/', self.region, 'pia_white_matter_normalized/swcs/removed_nan_nodes/base/')
		out_path = join('/data/elowsky/OLST/swc_analysis/analysis/', self.region, 'pia_white_matter_normalized/analysis/removed_nan_nodes/scaled_1000/soma_centered/', data_type)
		#out_path = join('/data/elowsky/OLST/swc_analysis/analysis/', self.region, 'pia_white_matter_normalized/analysis/removed_nan_nodes/scaled_1000/soma_centered/all_metrics/')	
		layers = self.layers
		norm_types = ['pca']
		structure_types = ['basal_apical_concat']

		layer_borders = {
			'mop':{ 
				'2': .1,
				'5': .36,
				'6': .60,},
			'barrel_cortex':{
				'2': .1,
				'4': .37,
				'5': .43,
				'6': .64,},
			'pfc':{}
		}


		if not exists(out_path):
			mkdir(out_path)

		for layer in layers:

			print('Layer:', layer)
	
			cluster_path_layer = join(cluster_path, f'layer_{layer}')
			out_path_layer = join(out_path, f'layer_{layer}')
			if not exists(out_path_layer):
				mkdir(out_path_layer)
			out_path_layer_consensus = join(out_path, f'layer_{layer}', 'consensus_clustering_5_25')
			if not exists(out_path_layer_consensus):
				mkdir(out_path_layer_consensus)
	
			for norm_type in norm_types:

				print('Normalization Type:', norm_type)		

				cluster_path_norm = join(cluster_path_layer, norm_type)
				out_path_norm = join(out_path_layer_consensus, norm_type)
				if not exists(out_path_norm):
					mkdir(out_path_norm)

				for structure_type in structure_types:

					print('Structure Type:', structure_type)

					out_path_structure = join(out_path_norm, structure_type)
					if not exists(out_path_structure):
						mkdir(out_path_structure)
		
					# load clusters
					cluster_path_structure = join(cluster_path_norm, structure_type, self.cluster_file_name_dict[self.region][structure_type][layer])
					clusters = self.load_clusters(cluster_path_structure)


					# iterate through cluster
					for cluster_id in np.unique(clusters[:, 2]):

						print('Cluster ID:', cluster_id)
	
						# get all swcs in cluster
						swc_ids_in_cluster = clusters[clusters[:,2] == cluster_id][:,0]

						# instaniate figure
						fig = plt.figure(figsize=(18,10))
						ax = plt.gca()
						#ax.set_title('Layer: '+ layer + '\nNormalization Type: ' + norm_type + '\nStructure Type: ' + structure_type + '\nCluster: ' + str(cluster_id))
						ax.axes.xaxis.set_visible(False)
						ax.axes.yaxis.set_visible(False)
						ax.axis('off')

						ax.set_ylim([-.6,.7])

						# pia and white matter lines
						#plt.axhline(y=.5, color='black', linestyle='-')
						#plt.axhline(y=-.5, color='black', linestyle='-')

						# create layer boundaries
						#for border_val in layer_borders[self.region].values():
						#	plt.axhline(y=.5-border_val, color='black', linestyle='--', linewidth=1)
	
						# iterate through swcs to plot
						x_pos = 0
						for i, swc_id in enumerate(swc_ids_in_cluster):
		
							swc_full_path = join(swc_path, swc_id + '.swc')

							# load swc
							swc = SWC(swc_full_path)

							# extract only basal and apical
							swc.extract_type_tree(['basal dendrite', 'apical dendrite'])
		
							# generate swc array
							swc_array = swc.generate_swc_array()
							basal_array = swc_array[swc_array[:,1] == 3, 2:5]
							apical_array = swc_array[swc_array[:,1] == 4, 2:5]
							plt.scatter(basal_array[:,0] + x_pos - swc_array[:,2].min() , -basal_array[:,1] + basal_array[0,1], c=BASAL_COLOR, s=DOT_SIZE)
							plt.scatter(apical_array[:,0] + x_pos - swc_array[:,2].min() , -apical_array[:,1] + apical_array[0,1], c=APICAL_COLOR, s=DOT_SIZE)

							# place swc name
							#if i % 2 == 0:
							#	plt.text(x_pos, .51, swc_id, fontsize=FONT_SIZE, fontweight='bold')
							#else:
							plt.text(basal_array[0,0] + x_pos - swc_array[:,2].min(), -.6, swc_id, fontsize=FONT_SIZE, fontweight='bold', ha='center')
		
							# update x position
							x_pos += swc_array[:,2].max() - swc_array[:,2].min() + X_SPACING

		
						out_path_full = join(out_path_structure, 'cluster_' + str(cluster_id) + '.png')
						plt.savefig(out_path_full)
						plt.close('all')



	def logistic_regression(
			self, 
			bootstrap=False, 
			bootstrap_samples=1000, 
			generate_bootstrap=False, 
			bootstrap_method='individual_features',
			bootstrap_layer=None,
			max_iter=1000, 
			c_precision=2, 
			remove_zero_features=True, 
			plot=False,
			save_plot=False,
			cluster_base_path=None,
			save_model=False,
			only_sig_features=False,
			layers_to_process=None):

		"logistic regression"

		base_path = join(self.analysis_path, 'morphometrics')
		structures = ['basal_apical_concat']

		if layers_to_process is None:
			layers = self.layers


		swc_path = join('/data/elowsky/OLST/swc_analysis/analysis/', self.region, 'pia_white_matter_normalized/swcs/removed_nan_nodes/base/')
		# load full list of features
		feature_path = '/data/elowsky/OLST/swc_analysis/analysis/barrel_cortex/pia_white_matter_normalized/analysis/removed_nan_nodes/scaled_1000/soma_centered/morphometrics/layer_2/basal_apical_concat.txt'

		significant_features = {
			'barrel_cortex':[
				'apical:EucDistance:Maximum',
				'basal:Depth:Total_Sum',
				'apical:Terminal_degree:Maximum',
				'basal:Length:Total_Sum',
				'apical:Branch_Order:Maximum',
				'apical:Length:Total_Sum',
				'apical:Height:Total_Sum',
				'apical:N_stems:Total_Sum',
				'basal:PathDistance:Average',
				'apical:PathDistance:Average',
				'apical:TerminalSegment:Total_Sum']
		}
	



		sig=set()

		print()
		print('Logistic Regression...')
		print()
		
		fig_size = (12,9)
		point_size = 15
		sig_size = 7
		x_tick_font_size = 8

		# iterate through layers
		for layer in layers:

			print('Layer:', layer)
	
			layer_path = join(base_path, f'layer_{layer}')
			log_reg_path = join(layer_path, 'logistic_regression')
			
			# make output path if it doesnt exist
			if not exists(log_reg_path):
				mkdir(log_reg_path)
	
			# iterate through structures
			for structure in structures:

				print(f'Structure: {structure}')
	
				# load clusters
				cluster_path = join(
					cluster_base_path, 
					f'layer_{layer}', 
					'pca',
					structure, 
					self.cluster_file_name_dict[self.region][structure][layer])
				print(f'Cluster Path: {cluster_path}')
				
				clusters = self.load_clusters(cluster_path)
		
				# make sure data exists	
				data_path = join(layer_path, 'normalized', f'{structure}.txt')	
				
				if exists(data_path):

					# load data
					data_matrix = np.genfromtxt(data_path, delimiter=' ', dtype=object)
					feature_names = data_matrix[0, 2:].astype(str)
					X = data_matrix[1:,2:].astype(float)

					# get clusters	
					y = clusters[:,2].astype(int)

					# bootstrap
					if bootstrap:

						# generate bootstrap coefficients
						if generate_bootstrap:
						
							if layer == bootstrap_layer:											
								bootstrap_coeffs = self.bootstrap_logistic_regression(
												X, 
												y, 
												method=bootstrap_method,
												c_precision=c_precision, 
												max_iter=max_iter, 
												bootstrap_samples=bootstrap_samples,
												save=True,
												out_path=log_reg_path)
							continue

						# load bootstrap coeffs
						bootstrap_coeffs_path = join(log_reg_path, f'bootstrap_coeffs_{bootstrap_samples}.npy')
						bootstrap_coeffs = np.load(bootstrap_coeffs_path)

			
					# load model if exists
					model_path = join(log_reg_path, 'log_reg_model.sav')
					if exists(model_path):
						model = pickle.load(open(model_path, 'rb'))
					else:

						# find highest c value that still gives 100% accuraacy
						# this will ensure strongest regulariszation that still has 100% accuracy
						optimal_c = self.logistic_regression_find_optimal_c(X, y, max_iter=max_iter, c_precision=c_precision)
						print(f'\t\t\t\tC: {optimal_c}')

						# train logistic regression
						model = LogisticRegression(
								penalty='l1',
								solver='liblinear',
								multi_class='auto',
								C=optimal_c,
								max_iter=max_iter)

						# fit model
						model.fit(X,y)

						# get accuracy
						accuracy = model.score(X,y)
						print(f'\t\t\t\tAccuracy: {accuracy}')

						if save_model:
							pickle.dump(model, open(model_path,'wb'))

			
					# remove all features where all models have zero coefficients
					if remove_zero_features:
						coefs_bool_map = ~np.all(model.coef_ == 0, axis=0)
						coeffs = model.coef_[:,coefs_bool_map]
						feature_names = feature_names[coefs_bool_map]
						if bootstrap:
							bootstrap_coeffs = bootstrap_coeffs[coefs_bool_map,:,:]
					else:
						coeffs = model.coef_



					full_features = np.genfromtxt(feature_path, dtype=str)[0,2:]

					# if feature size create new 
					if len(full_features) != len(feature_names):

						for f_idx, feature in enumerate(full_features):
							if feature not in feature_names:
								coeffs = np.insert(coeffs, f_idx, np.zeros(coeffs.shape[0]), axis=1)
								bootstrap_coeffs = np.insert(bootstrap_coeffs, f_idx, np.zeros((bootstrap_coeffs.shape[1],bootstrap_coeffs.shape[2])), axis=0)

						feature_names = full_features

					if only_sig_features:
						sig_feature_mask = np.in1d(feature_names,significant_features[self.region])
						coeffs = coeffs[:, sig_feature_mask]
						feature_names = feature_names[sig_feature_mask]
						bootstrap_coeffs = bootstrap_coeffs[sig_feature_mask,:,:]


				
					
					# create plot
					fig, axs = plt.subplots(
							len(model.classes_),
							1, 
							sharex=True, 
							figsize=fig_size)


					# get p values
					if bootstrap:
						p_values = np.zeros_like(coeffs)
						# for every coefficient 
						for i, cluster_coeffs in enumerate(coeffs):
							for j, coeff in enumerate(cluster_coeffs):
								
								greater_bool_matrix = bootstrap_coeffs[j,:,i] >= abs(coeff)
								less_bool_matrix = bootstrap_coeffs[j,:,i] <= -abs(coeff)
							
								bool_matrix = np.logical_or(greater_bool_matrix, less_bool_matrix)
							
								p = bool_matrix.sum() / bool_matrix.size
								p_values[i,j] = p
							
		
						# multiple comparisons corrections
						p_values = multipletests(p_values.flatten(), method='bonferroni')[1]
						p_values = p_values.reshape(coeffs.shape)



					# order clusters in same order as plots
					# swcs in cluster lists come in pairs
					# individual layers have order preserved 
					# layer all is reordered based on height
					cluster_plot_list = Clusters.clusters[self.region][layer]	
					unique_clusters = model.classes_.copy()

					cluster_order = []
					
					if layer == 'all':

						# sort clusters by height
						heights = []
						for swc_id in cluster_plot_list:

							swc_full_path = join(swc_path, swc_id + '.swc')
							swc = SWC(swc_full_path, build_tree=False)
							heights.append(-(swc.swc_array[:,3] - swc.swc_array[0,3]).min() )
						heights = np.array(heights)
						max_heights = [max(a,b) for a,b in zip(heights[::2],heights[1::2])]
						indices = np.argsort(max_heights)
					
						# skip every other swc index by indices
						cluster_plot_list = np.array(cluster_plot_list[::2])[indices]

						#for each swc get whigh cluster its in
						for swc_id in cluster_plot_list:
							cluster_order.append(clusters[clusters[:,0] == swc_id][0,2])

						# order clusters
						cluster_order = np.array(cluster_order)
						unique_clusters = unique_clusters[cluster_order]
				
	
					else:		

						# skip every other swc
						cluster_plot_list = cluster_plot_list[::2]
			
						#for each swc get whigh cluster its in
						for swc_id in cluster_plot_list:
							cluster_order.append(clusters[clusters[:,0] == swc_id][0,2])

						# order clusters
						cluster_order = np.array(cluster_order)
						unique_clusters = unique_clusters[cluster_order]


				

					p_strings = []

					# iterate through clusters
					for i, cluster in enumerate(unique_clusters):
	

						# get coefficients for cluster
						cluster_coeffs = coeffs[cluster]


						axs[i].axhline(y=0, c="grey",linewidth=0.5, alpha=.4)

						if layer == '6' and i < len(model.classes_)-1:
							axs[i].xaxis.set_visible(False)
						
						p_strings_row = []
						for j, coeff in enumerate(cluster_coeffs):

							if bootstrap:
								p = p_values[cluster,j]
								p_string = utils.p_value_to_string(p)
								p_strings_row.append(p_string)
								if p_string != '':
									print(f'Cluster: {cluster}, Feature: {feature_names[j]}, p: {p}') 
									sig.add(feature_names[j])
							if coeff == 0:	
								axs[i].scatter(j, coeff,s=point_size,color='black',facecolors='none')
							elif coeff > 0:
								axs[i].scatter(j, coeff,s=point_size,color='blue')
								if bootstrap:
									axs[i].text(j, coeff + 1, p_string,size=sig_size, ha='center')
							else:
								axs[i].scatter(j, coeff,s=point_size,color='red')
								if bootstrap:
									axs[i].text(j, coeff + 1, p_string, size=sig_size, ha='center')

						p_strings.append(p_strings_row)

						# find coeff with maximum magnitude
						max_coeff = max(abs(cluster_coeffs))


						# get accuracies (should all be 100%)
						cluster_acc = model.score(X[y==cluster],y[y==cluster])
						if cluster_acc != 1:
							print(f'Error: Cluster Accuracy is not 100% is actually {np.round(cluster_acc*100,2)}')

						# axis options
						axs[i].set_ylim([-10, 10])
						axs[i].set_xticks([])
						axs[i].set_yticklabels([])
						axs[i].set_yticks([])


						axs[i].yaxis.grid()
						axs[i].set_facecolor('lightgray')
						axs[i].tick_params(axis='y', which='major', labelsize=5)


					f_type_map = {'maximum': 'max', 'total_sum':'sum', 'average':'avg'}
					s_type_map = {'apical':'api','basal':'bas'}
					# modify feature names
					for i, feature_name in enumerate(feature_names):
						parts = feature_name.lower().split(':')
						
						feature_name_modified = '_'.join([s_type_map[parts[0]],parts[1],f_type_map[parts[2]]])
						
						feature_names[i] = feature_name_modified

					# plot options
					if layer in ['6', 'all']:
						axs[-1].set_xticks(np.arange(len(feature_names)))
						axs[-1].set_xticklabels(feature_names, rotation=90,fontsize=x_tick_font_size)
						

						# add space below 
						fig.subplots_adjust(bottom=0.2)

					# show plot
					if plot:
						plt.show()
		
					if save_plot:
						if only_sig_features:
							plt.savefig(
								join(log_reg_path, f'{structure}_only_significant_features_{fig_size}.png'), 
								dpi=300)
						else:
							plt.savefig(
								join(log_reg_path, f'{structure}.png'), 
								dpi=300)


					# heatmap

					# create plot
					fig, ax = plt.subplots(figsize=fig_size)
				
					coeffs = coeffs[unique_clusters]
					cmap = sns.diverging_palette(240, 10, s=100, as_cmap=True)
					ax = sns.heatmap(coeffs.T, cbar=False, square=True, cmap=cmap, annot=np.array(p_strings).T, fmt="", annot_kws={'fontsize':5})

					ax.set_xticklabels(np.arange(1,len(unique_clusters)+1))
					ax.set_yticklabels(feature_names, rotation=0)

					if save_plot:
						plt.savefig(
							join(log_reg_path, f'{structure}_heatmap.png'), 
							dpi=300)
	
				

	def logistic_regression_find_optimal_c(self, 
			X, 
			y, 
			c_init=1000, 
			c_precision=2, 
			max_iter=1000):


		# set c to inital
		c = c_init
		initial_precision = c_precision

		# set accuracy to 1
		acc = 1
		
		# set initial step
		step_size = c_init/10 

		# run logistic regression until accuracy is 100%
		while True:

			if c == 0:
				c = step_size
				step_size = c/10
				continue
			
			# instantiate model
			model = LogisticRegression(
				penalty='l1',
				solver='liblinear',
				multi_class='auto',
				C=c,
				max_iter=max_iter)
	
			# fit model
			model.fit(X,y)

			# get accuracy
			acc = model.score(X,y)

			# if initial c cant get full accuracy then just choose that c
			if c == c_init and acc != 1:
				return np.round(c, initial_precision)

			# if accuracy is still 100% then decrement c and try again
			# otherwise 
			if acc == 1:
				c -= step_size
			else:
				if c_precision == 1:
					if c != c_init:
						c += step_size
						break
					else:
						exit('INITIAL C VALUE DID NOT ALLOW FOR 100% ACCURACY')
				else:
					c += step_size
					if step_size <= .1:
						c_precision -= 1
					step_size /= 10
			
		return np.round(c, initial_precision)


	def bootstrap_logistic_regression(
			self,
			X, 
			y, 
			method,			
			bootstrap_samples=1000, 
			c_precision=2, 
			max_iter=1000, 
			save=False, 
			out_path=None):
	
		"""
		Performs bootstrapping for logistic regression to find p value
			
			- assigns each sample to a random cluster

		"""

		if method == 'clusters_ids':
			coeffs = []

			for i in range(bootstrap_samples):
		
				print('i:',i)

				# assign to random clusters
				y = np.random.randint(0, max(y)+1, size=len(X))

				# find optimal c value
				c = self.logistic_regression_find_optimal_c(X, y, c_precision=c_precision, max_iter=max_iter)

				# train logistic regression
				clf = LogisticRegression(penalty='l1',solver='liblinear',multi_class='auto',C=c, max_iter=max_iter).fit(X,y)

				coeffs.append(clf.coef_)

			coeffs = np.array(coeffs)
			
		elif method == 'individual_features':
	
			# output array		
			coeffs = np.zeros((X.shape[1], bootstrap_samples, len(np.unique(y))))
	
			# iterate through features
			print(X.shape, bootstrap_samples)
			
			for feature_idx in range(X.shape[1]):
				print(f'Feature : {feature_idx}')				

				# iterate through bootstrap runs
				for i in range(bootstrap_samples):

					print(f'i: {i}')

					# get features data
					feature_data = X[:,feature_idx]
			
					# randomly sample with replacement
					sampled_data = np.random.choice(feature_data, size=len(feature_data), replace=True)
	
					# place sampled data in new array
					X_sampled = np.copy(X)
					X_sampled[:,feature_idx] = sampled_data

					# find optimal c value
					c = self.logistic_regression_find_optimal_c(
							X_sampled, 
							y, 
							c_precision=c_precision, 
							max_iter=max_iter)

					print(f'\tc: {c}')

					# run logistic regression
					model = LogisticRegression(
							penalty='l1',
							solver='liblinear',
							multi_class='auto',
							C=c, 
							max_iter=max_iter).fit(X_sampled,y)

					# get coefficients for specific feature
					feature_coeffs = model.coef_[:,feature_idx]
					
					# place in output array
					coeffs[feature_idx,i,:] = feature_coeffs

	
		else:
			exit(f'Error: invalid bootstrap method - {bootstrap_method}')
		

		if save:
			np.save(join(out_path,f'bootstrap_coeffs_{bootstrap_samples}_{method}'), coeffs)		

		return coeffs

	def cluster_morphometrics_statistics(self):

		base_path = join(self.analysis_path, 'morphometrics')
		structures = ['basal_apical_concat']
		layers = self.layers
		
		# iterate through layers
		for layer in layers:
			print('\tLayer:', layer)

			layer_path = join(base_path, 'layer_' + str(layer))
	
			for structure in structures:
				
				# load morphometrics
				morpho_path = join(layer_path, structure + '.txt')
				morphometrics = np.genfromtxt(morpho_path, dtype=object)
				feature_names = morphometrics[0,2:].astype(str)
				morphometrics = morphometrics[1:,2:].astype(np.float64)

				# load clusters
				cluster_path_structure = join(self.cluster_path, 'layer_' + layer, 'pca', structure, self.cluster_file_name_dict[self.region][structure][layer])

				cluster_info = np.genfromtxt(cluster_path_structure, delimiter=' ', dtype=object)
				cluster_info[:,0] = cluster_info[:,0].astype(str)
				cluster_info[:,1] = cluster_info[:,1].astype(str)
				cluster_info[:,2] = cluster_info[:,2].astype(int)
				cluster_info[:,1] = cluster_info[:,1] + '_' + cluster_info[:,0]
				cluster_info = cluster_info[:,1:]
				cluster_info = cluster_info[cluster_info[:,0].argsort()]
				unique_clusters = np.unique(cluster_info[:,1])



				means = [] 
				stds = []
				for cluster_id in unique_clusters:
					cluster_bool_map = cluster_info[:,1] == cluster_id
					cluster_morpho = morphometrics[cluster_bool_map,:]
					means.append(list(np.mean(cluster_morpho, axis=0).round(2)))
					stds.append(list(np.std(cluster_morpho, axis=0).round(2)))
				
				out_path_full = join(layer_path, f'morpho_stats_means_{structure}.txt')
				with open(out_path_full, 'w') as fp:
					fp.write(','.join([str(x) for x in feature_names]) + '\n')
					for cluster_id, mean_list in enumerate(means):
						fp.write(','.join([str(x) for x in mean_list]) + '\n')

				out_path_full = join(layer_path, f'morpho_stats_stds_{structure}.txt')
				with open(out_path_full, 'w') as fp:
					fp.write(','.join([str(x) for x in feature_names]) + '\n')
					for cluster_id, std_list in enumerate(stds):
						fp.write(','.join([str(x) for x in std_list]) + '\n')


	
				

	def load_clusters(self, path):
		
		"load clusters"

		# make sure cluster path exists
		if not exists(path):
			exit(f'CLUSTER PATH DOES NOT EXIST: {path}' )

		# load file
		clusters = np.genfromtxt(path, delimiter=' ', dtype=object)

		# prepend cluster id
		clusters[:,0] = clusters[:,1] + b'_' + clusters[:,0]

		# change data types
		clusters[:,0] = clusters[:,0].astype(str)
		clusters[:,1:] = clusters[:,1:].astype(int)


		# sort by id
		clusters = clusters[clusters[:,0].argsort()]

		return clusters


	def box_plots(self, cluster_path, plot=False, save=False):


		base_path = join(self.analysis_path, 'morphometrics')
		structures = ['basal_apical_concat']
		layers = self.layers

		print()
		print('Box Plots...')
		print()
		
		# iterate through layers
		for layer in layers:

			print('Layer:', layer)
	
			layer_path = join(base_path, f'layer_{layer}')
			out_path = join(layer_path, 'box_plots')
			
			# make output path if it doesnt exist
			if not exists(out_path):
				mkdir(out_path)
	
			# iterate through structures
			for structure in structures:

				print(f'Structure: {structure}')
	
				# load clusters
				cluster_path_layer = join(
					cluster_path, 
					f'layer_{layer}', 
					'pca',
					structure, 
					self.cluster_file_name_dict[self.region][structure][layer])

				
				clusters = self.load_clusters(cluster_path_layer)
				unique_clusters = np.unique(clusters[:,2])

				out_path_structure = join(out_path, structure)
				if not exists(out_path_structure):
					mkdir(out_path_structure)
		
				# make sure data exists	
				data_path = join(layer_path, 'normalized', f'{structure}.txt')	
				
				if exists(data_path):

					morphometrics = np.genfromtxt(data_path, dtype=object)
					feature_names = morphometrics[0,2:].astype(str)
					morphometrics = morphometrics[1:,2:].astype(np.float64)


					for feature_idx, feature_name in enumerate(feature_names):
					
						feature_data = morphometrics[:,feature_idx]
		
						fig, ax  = plt.subplots()

						b_plot_data = []
						for cluster_id in unique_clusters:
							cluster_bool_map = clusters[:,2] == cluster_id
							cluster_morpho = feature_data[cluster_bool_map]
							b_plot_data.append(cluster_morpho)
				

						#ax.set_xticklabels(unique_clusters)
						ax.set_xlabel('Cluster')
						ax.set_title(f'{structure} : {feature_name}')
		
						bplot = ax.boxplot(b_plot_data, whis=float('infinity'), positions=unique_clusters)

						# wilcoxon
						for i in range(len(unique_clusters)):
							
							group_a = b_plot_data[i]
							group_b = np.concatenate([x for idx,x in enumerate(b_plot_data) if idx != i])
		
							_, p = ranksums(group_a, group_b)
							p_string = utils.p_value_to_string(p)
						
							
							ax.text(i, bplot['caps'][2*i+1].get_ydata()[0], p_string, ha='center')


						if plot:
							plt.show()

						if save:
							plt.savefig(join(out_path_structure, feature_name + '.png'))
		
						plt.close('all')

	def find_zero_morphometrics(self):
	
		regions = ['mop', 'barrel_cortex']
		base_path = '/data/elowsky/OLST/swc_analysis/analysis/'
		feature_dict = {}

		for region in regions:
			region_path = join(base_path, region, 'pia_white_matter_normalized/analysis/removed_nan_nodes/scaled_1000/soma_centered/morphometrics/layer_all/basal_apical_concat.txt')

			data, feature_names = self.load_feature_data(region_path)
			
			zero_features = feature_names[(data == 0).any(axis=0)]
			swc_counts = (data == 0).sum(axis=0)[(data == 0).any(axis=0)]
			print(swc_counts)
			for zero_feature, swc_count in zip(zero_features, swc_counts):
				if zero_feature in feature_dict:
					feature_dict[zero_feature] += swc_count
				else:
					feature_dict[zero_feature]  = swc_count
		pprint.pprint(feature_dict)
			
	def load_feature_data(self, path):

		morphometrics = np.genfromtxt(path, delimiter=' ', dtype=object)
		swc_names = morphometrics[1:,0].astype(str)
		feature_names = morphometrics[0,2:].astype(str)
		data = morphometrics[1:, 2:].astype(np.float64)
		
		return data, feature_names, swc_names


	def cross_layer_cluster_analysis(self, cluster_path):
		
		base_path = join(self.analysis_path, 'morphometrics')
		structures = ['basal_apical_concat']
		layers = self.layers
		out_base_path = join(base_path, 'cross_layer_analysis')

		# make out path
		if not exists(out_base_path):
			mkdir(out_base_path)

		global_mapping_data = []		

		# iterate through all 2-permutations of layers
		for structure in structures:
			for layer_input in layers:
				for layer_model in layers:
			
					if layer_input != layer_model:

						print('Input Layer: ', layer_input)
						print('Model Layer: ', layer_model)
					
						# load input data
						data_path = join(base_path, f'layer_{layer_input}', 'normalized', f'{structure}.txt')
						data_matrix = np.genfromtxt(data_path, delimiter=' ', dtype=object)
						feature_names = data_matrix[0, 2:].astype(str)
						X = data_matrix[1:,2:].astype(float)

						# load model data
						data_path_model = join(base_path, f'layer_{layer_model}', 'normalized', f'{structure}.txt')
						data_matrix_model = np.genfromtxt(data_path_model, delimiter=' ', dtype=object)
						feature_names_model = data_matrix_model[0, 2:].astype(str)


						# handle case when different features in model input data
						diff_features = set(feature_names).symmetric_difference(set(feature_names_model))
						if len(diff_features):
							
							X_new = np.zeros(shape=(X.shape[0],0))

							# iterate through feature names in model
							# if feature in dataset then concat
							# if not concat zeros
							for feature_name_model in feature_names_model:
								if feature_name_model in feature_names:
									index = np.argwhere(feature_names == feature_name_model)[0][0]
									X_new = np.hstack((X_new, X[:,index].reshape(-1,1)))
								else:
									X_new = np.hstack((X_new, np.zeros(shape=(X_new.shape[0],1))))
	
							X = X_new
							feature_names = feature_names_model


						# load model
						model_path = join(base_path, f'layer_{layer_model}', 'logistic_regression', 'log_reg_model.sav')
						model = pickle.load(open(model_path, 'rb'))

						# predict on input data
						predictions = model.predict(X)

						# sort by clusters in input layer
						cluster_path = join(cluster_base_path, f'layer_{layer_input}', 'pca', structure, self.cluster_file_name_dict[self.region][structure][layer_input])
						cluster_data = self.load_clusters(cluster_path)
						sorted_order = np.argsort(cluster_data[:,-1], kind='mergesort')
						cluster_data_sorted = cluster_data[sorted_order]

						# sort predictions
						predictions_sorted = predictions[sorted_order]

						# create output array
						output = np.hstack((cluster_data_sorted, predictions_sorted.reshape(-1,1)))
						headers = ['swc_id', 'layer', 'input_cluster', 'output_cluster']
						output = np.vstack((headers, output))

						out_path = join(out_base_path,f'{layer_input}_{layer_model}.csv')
						np.savetxt(out_path, output, fmt='%s')


						# iterate through cluster and store percentage mappings
						data = output[1:]
						for in_cluster in np.unique(data[:,2]):
							c_data = data[data[:,2] == in_cluster]
							for out_cluster in np.unique(c_data[:,3]):
								perc = np.round((c_data[:,3] == out_cluster).sum() / len(c_data),2)
								global_mapping_data.append((layer_input, in_cluster, layer_model, out_cluster,  len(c_data), perc))
						

		# find clusters that map most to other clusters
		global_mapping_data.sort(key=lambda y: y[-1], reverse=True)
		global_mapping_data = np.array(global_mapping_data)

		out_path = join(out_base_path,f'global_cluster_maping.csv')
		np.savetxt(out_path, global_mapping_data, fmt='%s')


	def plot_clusters_flexible(self, clusters):
		

		clusters = clusters[self.region]


		# plot parameters
		DOT_SIZE= 1.5
		X_SPACING = .1

		cmap_a = mpl.cm.get_cmap(name='tab10')
		cmap_b = mpl.cm.get_cmap(name='Dark2')
		cmap_c = mpl.cm.get_cmap(name='Set1')
		colors = cmap_a.colors + cmap_b.colors + cmap_c.colors


		color_idx = 0


		# paths
		swc_path = join('/data/elowsky/OLST/swc_analysis/analysis/', self.region, 'pia_white_matter_normalized/swcs/removed_nan_nodes/base/')
		out_path = join('/data/elowsky/OLST/swc_analysis/analysis/', self.region, 'pia_white_matter_normalized/analysis/removed_nan_nodes/scaled_1000/soma_centered/cluster_images/')

		layer_borders = {
			'mop':{ 
				'2': .1,
				'5': .36,
				'6': .60,},
			'barrel_cortex':{
				'2': .1,
				'4': .37,
				'5': .43,
				'6': .64,},
			'pfc':{}
		}



	
		for layer, swc_ids in clusters.items():
			

			if layer == 'all':
				color_idx = 0

				# sort clusters by height
				heights = []
				for swc_id in swc_ids:

					swc_full_path = join(swc_path, swc_id + '.swc')
					swc = SWC(swc_full_path, build_tree=False)
					heights.append(-(swc.swc_array[:,3] - swc.swc_array[0,3]).min() )
				heights = np.array(heights)
				max_heights = [max(a,b) for a,b in zip(heights[::2],heights[1::2])]
	
				indices = np.argsort(max_heights)*2
				sorted_inds = []
				for x in indices:
					sorted_inds.append(x)
					sorted_inds.append(x+1)

				# sort ids
				swc_ids = np.array(swc_ids)[sorted_inds]

				# split into 4 groups
				swc_id_groups = np.split(swc_ids, 4)
			else:
				swc_id_groups = [swc_ids]	
	
			for group_idx, swc_ids in enumerate(swc_id_groups):
				
				# instaniate figure
				fig = plt.figure(figsize=(18,10))
				ax = plt.gca()

				ax.axes.xaxis.set_visible(False)
				ax.axes.yaxis.set_visible(False)
				ax.axis('off')	
				ax.set_ylim([-.6,.6])	


				# iterate through swcs to plot
				x_pos = 0
				for i, swc_id in enumerate(swc_ids):

					swc_full_path = join(swc_path, swc_id + '.swc')

					# load swc
					swc = SWC(swc_full_path)

					# extract only basal and apical
					swc.extract_type_tree(['basal dendrite', 'apical dendrite'])

					# generate swc array
					swc_array = swc.generate_swc_array()
					basal_array = swc_array[swc_array[:,1] == 3, 2:5]
					apical_array = swc_array[swc_array[:,1] == 4, 2:5]

					basal_color = [*colors[color_idx],1]
					apical_color = [*colors[color_idx],.1]
					basal_color[0] *= .6
					basal_color[1] *= .6
					basal_color[2] *= .6


					plt.scatter(basal_array[:,0] + x_pos - swc_array[:,2].min() , -basal_array[:,1]+basal_array[0,1], color=basal_color, s=DOT_SIZE)
					plt.scatter(apical_array[:,0] + x_pos - swc_array[:,2].min() , -apical_array[:,1]+apical_array[0,1], color=apical_color, s=DOT_SIZE)

					if i % 2 == 1:
						color_idx += 1
					
					# update x position
					x_pos += swc_array[:,2].max() - swc_array[:,2].min() + X_SPACING

				out_path_full = join(out_path, f'layer_{layer}_group_{group_idx}.png')
				plt.savefig(out_path_full, dpi=300)
				plt.close('all')


	def soma_density_plots(self, cluster_path):

		# parameters
		pia_pos = 0
		white_matter_pos = 1
		num_bins = 20
		bins = np.linspace(pia_pos, white_matter_pos, num_bins)

		# output path
		out_path = join(self.analysis_path, 'soma_density_plots')
	
		# path were swcs are
		swcs_path = join(self.swc_path, self.remove_nan_type, 'base')

		# iterate though layers
		for layer in self.layers:

			# load clusters
			cluster_layer_path = join(
				cluster_path, 
				f'layer_{layer}',
				'pca', 
				'basal_apical_concat',
				self.cluster_file_name_dict[self.region]['basal_apical_concat'][layer])

			clusters = self.load_clusters(cluster_layer_path)
			
			layer_somas = []

			# iterate through cluster
			for cluster_id in np.unique(clusters[:, 2]):
	
				# get all swcs names from cluster
				cluster_swc_names = clusters[clusters[:,2] == cluster_id][:,0]
				
				# iterate through swcs
				cluster_somas = []
				for swc_name in cluster_swc_names:
				
					# get soma coordinate from swc
					# negate soma cooridnate and subtract from .5
					swc_full_path = join(swcs_path, f'{swc_name}.swc')
					swc = SWC(swc_full_path, build_tree=False)
					soma_y_pos = swc.swc_array[0,3]
					cluster_somas.append(.5-soma_y_pos)
				
				layer_somas.append(cluster_somas)

			# create density plot
			
			fig, ax = plt.subplots()
			ax.set_ylim([pia_pos, white_matter_pos])
			ax.set_xticks([])
			ax.set_yticks([])
			ax.axis('off')	

			# pia and white matter lines
			plt.axhline(y=0, color='black', linestyle='-')
			plt.axhline(y=1, color='black', linestyle='-')		

			#for cluster_somas in layer_somas:
								
				#ax.hist(
				#	cluster_somas, 
				#	bins=bins, 
				#	orientation='horizontal',
				#	histtype='step')
				#
				# plot point at middle of every bin
				#hist, bin_edges = np.histogram(cluster_somas, bins=bins)
				#bin_centers = [(a+b)/2 for a,b in zip(bins[:-1],bins[1:])]
				#hist =[0] + list(hist) + [0]
				#bin_centers = [0] + bin_centers + [1]
				#ax.plot(hist, bin_centers)

				#ax.hist(cluster_somas, bins=bins, orientation='horizontal')

			ax.hist(layer_somas, bins=bins, orientation='horizontal')
			out_path_full = join(out_path, f'layer_{layer}_bins={num_bins}.png')
			plt.savefig(out_path_full)
			plt.close('all')
			

			
		
				
	def feature_type_pairwise_comparison(self, cluster_path):
		
		"""
		- comparison of clusters obtained from morphometrics
		with other feature types

			- arbor density
			- persistent homology
	
		- performs intra/inter cluster comparisons

		"""

		swc_path = join('/data/elowsky/OLST/swc_analysis/analysis/', self.region, self.swc_type, 'swcs', self.remove_nan_type, 'base')
		feature_types = ['arbor_density', 'persistent_homology']
		fig_size = (20,6)
		x_tick_font_size = 12
		sig_font_size = 15
		inter_intra_x_diff = .7

		# output path
		out_path = join(self.analysis_path, 'cross_feature_comparison')
		if not exists(out_path):
			mkdir(out_path)
	
		# iterate through feature types
		for feature_type in feature_types:

			print('Feature Type:', feature_type)

			out_path_feature = join(out_path, feature_type)
			if not exists(out_path_feature):
				mkdir(out_path_feature)


			# iterate though layers
			for layer in self.layers:

				print('Layer:', layer)

				out_path_layer = join(out_path_feature, f'layer_{layer}')
				if not exists(out_path_layer):
					mkdir(out_path_layer)

				# load data from feature type
				feature_data_path = join(
						self.analysis_path,
						feature_type,
						f'layer_{layer}',
						'pca',
						'basal_apical_concat.txt')
			
				feature_data, feature_names, feature_swc_names = self.load_feature_data(feature_data_path)


				# load clusters
				cluster_layer_path = join(
					cluster_path, 
					f'layer_{layer}',
					'pca', 
					'basal_apical_concat',
					self.cluster_file_name_dict[self.region]['basal_apical_concat'][layer])

				clusters = self.load_clusters(cluster_layer_path)

			
				# array to hold all feature data
				unique_clusters = np.unique(clusters[:, 2]) 
				
				# order clusters in same order as plots
				# swcs in cluster lists come in pairs
				# individual layers have order preserved 
				# layer all is reordered based on height
				cluster_plot_list = Clusters.clusters[self.region][layer]

				cluster_order = []

				if layer == 'all':

					# sort clusters by height
					heights = []
					for swc_id in cluster_plot_list:

						swc_full_path = join(swc_path, swc_id + '.swc')
						swc = SWC(swc_full_path, build_tree=False)
						heights.append(-(swc.swc_array[:,3] - swc.swc_array[0,3]).min() )
					heights = np.array(heights)
					max_heights = [max(a,b) for a,b in zip(heights[::2],heights[1::2])]
					indices = np.argsort(max_heights)
					
					# skip every other swc index by indices
					cluster_plot_list = np.array(cluster_plot_list[::2])[indices]

					#for each swc get whigh cluster its in
					for swc_id in cluster_plot_list:
						cluster_order.append(clusters[clusters[:,0] == swc_id][0,2])

					# order clusters
					cluster_order = np.array(cluster_order)
					unique_clusters = unique_clusters[cluster_order]
				
	
				else:		

					# skip every other swc
					cluster_plot_list = cluster_plot_list[::2]
			
					#for each swc get whigh cluster its in
					for swc_id in cluster_plot_list:
						cluster_order.append(clusters[clusters[:,0] == swc_id][0,2])

					# order clusters
					cluster_order = np.array(cluster_order)
					unique_clusters = unique_clusters[cluster_order]


				cluster_features_data_global = []
				global_layer_data = []	
				global_p_values = []
				global_ax_labels = []
				x_pos = 1
				x_positions = []	


				# iterate through clusters
				for cluster_id in unique_clusters:

					# get all swcs names from cluster
					cluster_swc_names = clusters[clusters[:,2] == cluster_id][:,0]
					cluster_swc_names = ['_'.join(x.split('_')[1:]) for x in cluster_swc_names]
			
					# get all feature data from cluster
					swc_name_mask = np.in1d(feature_swc_names, cluster_swc_names)
					cluster_feature_data = feature_data[swc_name_mask]

					cluster_features_data_global.append(cluster_feature_data)
				
					x_positions.append(x_pos)
					x_positions.append(x_pos+inter_intra_x_diff)	
					x_pos += 2
					

				# for each cluster computer intra and inter pairwise distances
				for i, data_i in enumerate(cluster_features_data_global):

					data_dict = {
						'intra':None,
						'inter':None,
					}

					for j, data_j in enumerate(cluster_features_data_global):

						pairwise_distances = distance_matrix(data_i, data_j).flatten()

						# if intra remove duplicate distances
						if i == j:
	
							# if intra then matrix will be square
							# remove zeros on diagonal
							# each distance is duplicated so just take every other
					
							pairwise_distances.sort()
							pairwise_distances = pairwise_distances[len(data_i):][::2]
							data_dict['intra'] = pairwise_distances
						else:
							if data_dict['inter'] is None:
								data_dict['inter'] = pairwise_distances
							else:
								data_dict['inter'] = np.concatenate((data_dict['inter'],pairwise_distances))

					global_layer_data.append(data_dict['intra'])
					global_layer_data.append(data_dict['inter'])
					global_ax_labels.append(f'intra')
					global_ax_labels.append(f'inter')

					# run wilxocon test
					p = mannwhitneyu(data_dict['intra'], data_dict['inter']).pvalue
					p_string = utils.p_value_to_string(p)
					global_p_values.append(p_string)

			
				# create box and whisker plots
				fig, ax = plt.subplots(figsize=fig_size)

				ax.boxplot(global_layer_data, whis=float('infinity'), positions=x_positions)

				x_tick_positions = [(a + b) / 2 for a, b in zip(x_positions[::2], x_positions[1::2])]

				x_labels = np.arange(1,len(x_tick_positions)+1)
				ax.set_xticks(x_tick_positions)

				ax.set_xticklabels([])

				ax.yaxis.set_visible(False)
	

				y_max = np.max([np.max(x) for x in global_layer_data])
				ax.set_ylim([0,y_max+10])

				# place p value strings
				for i in range(0, len(global_p_values)*2, 2):
					if global_p_values[i//2] == '':
						global_p_values[i//2] = 'ns'
					ax.text((x_positions[i] + x_positions[i+1])/2, y_max+3, global_p_values[i//2], size=sig_font_size, ha='center')
	

				# save figure
				plt.savefig(join(out_path_layer, f'all_clusters_{fig_size}'))


		
						

				


					
					


		



if __name__ == '__main__':

	region = 'barrel_cortex'
	remove_nan_type = 'removed_nan_nodes'
	scale_type = 'scaled_1000'
	soma_type = 'soma_centered'
	swc_type = 'pia_white_matter_normalized'

	swc_analysis = SWC_Analysis(region, remove_nan_type, scale_type, soma_type, swc_type)

	#swc_analysis.consolidate_swcs_from_olga()

	#swc_analysis.handle_nan_nodes()

	#swc_analysis.scale_swcs()

	#swc_analysis.soma_center_swcs()

	#swc_analysis.extract_branch_types('scaled_1000/soma_centered/')

	#swc_analysis.run_lmeasure_morphometrics()

	#swc_analysis.parse_lmeasure_morphometrics()

	#swc_analysis.consolidate_arbor_densities()

	#swc_analysis.concat_arbor_densities()

	#swc_analysis.rename_persistent_vectors()

	#data_type = 'persistent_homology'
	#swc_analysis.normalize_data(data_type)

	#data_type = 'persistent_homology'
	#swc_analysis.pca(data_type)

	#swc_analysis.concat_persistent_homology_vectors()


	###########################################################

	#swc_analysis.hierarchical_clustering(cut_method='fcluster')

	#swc_analysis.create_2d_swc_plots()

	#swc_analysis.arbor_densities_to_binary_images()

	#path_a = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/swcs/removed_nan_nodes/scaled/basal_apical/'
	#path_b = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/swcs/removed_nan_nodes_old/soma_positioned/basal_apical/'
	#SWC_Analysis.check_if_two_folders_are_same(path_a, path_b)

	"""
	path = '/data/elowsky/OLST/color_test/'
	a =240
	b = 10

	fig, ax = plt.subplots()
	sns.palplot(sns.diverging_palette(a, b,s=100, n=9))
	plt.savefig(join(path, f'{a}_{b}.png'))
	"""
	
	swc_analysis.logistic_regression(
			bootstrap=True,
			bootstrap_samples=10000,
			generate_bootstrap=False,
			bootstrap_method='individual_features',
			bootstrap_layer='all', 
			max_iter=1000, 
			c_precision=2, 
			remove_zero_features=False,
			plot=False,
			save_plot=True,
			cluster_base_path = join('/data/palmer/consensus_clustering/output_swcs_removed/', region , 'morphometrics'),
			save_model=False,
			only_sig_features=False,
			layers_to_process=None)
	
	
	#swc_analysis.cluster_morphometrics_statistics()
	
	#cluster_path = join('/data/palmer/consensus_clustering/output_swcs_removed/', region , 'morphometrics')
	#swc_analysis.box_plots(cluster_path, plot=False, save=True)

	#data_type = 'morphometrics'
	#cluster_path = join('/data/palmer/consensus_clustering/output_swcs_removed/', region , 'morphometrics')
	#data_type = 'morphometrics'
	#swc_analysis.plot_swc_clusters(cluster_path, data_type)

	#swc_analysis.find_zero_morphometrics()

	#swc_analysis.cross_layer_cluster_analysis()


	#swc_analysis.plot_clusters_flexible(Clusters.clusters)


	#cluster_path = join('/data/palmer/consensus_clustering/output_swcs_removed/', region , 'morphometrics')	
	#swc_analysis.soma_density_plots(cluster_path)
	
	#cluster_path = join('/data/palmer/consensus_clustering/output_swcs_removed/', region , 'morphometrics')	
	#swc_analysis.feature_type_pairwise_comparison(cluster_path)
	




	





	
 







