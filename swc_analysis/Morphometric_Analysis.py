from SWC import SWC
from sklearn.decomposition import PCA
from os import listdir, system
from os.path import join, exists
from os import mkdir
from sys import exit
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D
from skimage.io import imread, imsave
from scipy.stats import zscore
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import pairwise_distances
from statsmodels.sandbox.stats.multicomp import multipletests
from shutil import copyfile

import scipy.cluster.hierarchy as shc
import params
import json
import utils
import numpy as np
import matplotlib.pyplot as plt
import math
import tifffile as tif

np.seterr(divide='ignore', invalid='ignore')


class Morphometric_Analysis:

	def __init__(self):

		print()
		print('#####################')
		print('Morphometric Analysis')		
		print('#####################')	
		print()

	@staticmethod
	def consolidate_swcs(in_path, in_dir_list, out_path, layers=(2,5,6)):

		
		# iterate through input directories
		for in_dir in in_dir_list:
			
			in_path_full = join(in_path, in_dir)
			
			# iterate through layers
			for layer in layers:
			
				layer_path = join(in_path_full, 'layer_' + str(layer))
				swc_names = listdir(layer_path)

				# iterate through swcs
				for swc_name in swc_names:
					
					source = join(layer_path, swc_name)
					destination = join(out_path, str(layer) + '_' + swc_name)
					copyfile(source, destination)
				
			
		
		
	
	def load_swcs(self, swc_path, *, registered=False, orientation=False):

		"""
		Loads all SWCs from layers specified in constructor

		"""

		# check if any SWC in swcPath, if not looks for layer folders
		swc_names = [f for f in listdir(swc_path) if f.endswith('.swc')]

		if len(swc_names) > 0:

			print()
			print('Loading all SWCs from swc path...')
			print()
			print('# SWCs:', len(swc_names))
			print()

			# instantiate swc d
			swcs = []

			for swc_name in sorted(swc_names):
				
				# get swc info 
				swc_info = utils.get_swc_info(swc_name)

				# load swc	
				swc_path_full = join(swc_path, swc_name)
				swc = SWC(swc_path_full)			
				swcs.append(swc)
		
		else:
			print()
			print('Loading SWCs from layer folders...')
			print()
	
			layerFolders = [f for f in listdir(swcpath) if f.startswith('layer_')]
			
			if len(layerFolders) == 0:
				exit('Error: No SWCs or Layer Folders exist in directory!')

			# instantiate swc dict
			swcs = {}
			for layer in params.LAYERS:
				swcs[layer] = []

			for layerFolder in sorted(layerFolders):

				layer = int(layerFolder.split('_')[-1])
				print('Layer:', layer)

				layerString = '_registered' if registered else '_not_registered' 
				layerPath = join(join(swcpath, layerFolder, 'all_swcs_' + orientation + layerString))
				print(layerPath)
	
				swcFilenames = listdir(layerPath)
				swcFilenames.sort()
	
				for swcFilename in swcFilenames:
					swc = SWC(join(layerPath,swcFilename))
					swcs[layer].append(swc)	

		self.swcs = swcs



	def extract_structures(self, structures, normalize=False):

		for swc in self.swcs:

			# extract structure trees
			swc.extract_type_tree(structures)

			if normalize:
				swc.center_around_soma()
	
	def save_swcs(self, out_path, *,  layer_folders=False, layer_prefix=True):

		if layer_folders:
			# to do
			pass
		else:

			for swc in self.swcs:
				swc_info = utils.get_swc_info(swc.swc_path)

				if layer_prefix:
					swc_out_path = join(out_path, str(swc_info['layer']) + '_' + swc_info['name'] + '.swc')
				else:
					swc_out_path = join(out_path, swc_info['name'] + '.swc')
				
				swc.save_swc(swc_out_path)
			


	def run_lmeasure_morphometrics(self, swc_path, out_path):

		"""
		Runs lmeasure to calculate Morphometrics	
	
		"""

		print()
		print('############################')
		print('Run Lmeasure (Morphometrics)')
		print('############################')
		print()

		# load morphometrics json to decide which morphometrics to calculate
		morphometrics_json_path = join(self.in_path, 'lmeasure_morphometrics.json')

		with open(morphometrics_json_path,'r') as fp:
			morphometrics_dict = json.load(fp)

		morphometric_names = sorted(morphometrics_dict)

		# create function specifier for lmeasure
		function_specifier = ' '.join(['-f' + str(params.LMEASURE_FUNCTION_IDS[m])+',0,0,10.0' for m in morphometric_names])

		for structure_type in ['basal', 'apical', 'basal_apical']:

			# swc paths
			swc_path_full = join(swc_path, structure_type)
			swc_paths = ' '.join([join(swc_path_full, f) for f in sorted(listdir(swc_path_full))])

			# form full lmeasure command
			out_path_full = join(out_path, 'lmout_morphometrics_' + structure_type + '.txt')
			lmeasure_command = ' '.join([params.LMEASURE_EXE, function_specifier, '-s' + out_path_full, swc_paths])

			# run lmeasure from command line
			system(lmeasure_command)

	def parse_lmeasure_morphometrics(self, lmout_path, json_path):
		
		"""
		Parse lmeasure morphometrics		

		"""
		
		print()
		print('#####################################')
		print('Parse Lmeasure Output (Morphometrics)')
		print('#####################################')
		print()
	
		# load morphometrics json to decide which morphometrics to calculate
		with open(json_path,'r') as fp:
			morpho_dict = json.load(fp)

		# dict to store morphometrics
		morphometrics_dict = {}
		for layer in params.LAYERS:
			morphometrics_dict[layer] = {}
		

		# iterate through structures
		for structure_type in ['basal', 'apical', 'basal_apical']:
			
			# save concatrnated basal apical morphometrics
			if structure_type == 'basal_apical':				
			
				# load and concatenate
				basal_morpho = np.loadtxt(join(lmout_path, 'morphometrics_basal.txt'),dtype=object)
				apical_morpho = np.loadtxt(join(lmout_path, 'morphometrics_apical.txt'),dtype=object)

				morphometrics = np.hstack((basal_morpho, apical_morpho[:,2:]))

				# save
				out_path =  join(lmout_path, 'morphometrics_'+structure_type+'_concat.txt')
				np.savetxt(out_path, morphometrics, fmt='%s')

	
			# get lmeasure output path
			lmout_morpho_path = join(lmout_path, 'lmout_morphometrics_' + structure_type + '.txt')

			# read in data into dictionary
			with open(lmout_morpho_path,'r') as fp:
				lines =  fp.readlines()
				for line in lines:
	
					line = line.split()

					if len(line):
						swc_name = line[0].split('/')[-1]
	
						# extract layer and swcID
						swc_info = utils.get_swc_info(swc_name)
						layer = swc_info['layer']
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
			for layer in morphometrics_dict:
				for swc_name in sorted(morphometrics_dict[layer]):
					row_data = [swc_name,layer]
					for metric in sorted(morphometrics_dict[layer][swc_name]):
						for item in morphometrics_dict[layer][swc_name][metric]:
							row_data.append(item)
					morphometrics.append(row_data)
	
			morphometrics = np.array(morphometrics,dtype=object)

			# save
			outpath =  join(lmout_path, 'morphometrics_'+ structure_type+ '.txt')
			np.savetxt(outpath, morphometrics, fmt='%s')

	def load_morphometrics(self, in_path, layers=None):
	
		morphometrics = np.loadtxt(in_path, dtype=object)

		# get morphometric labels
		labels = morphometrics[0,2:]

		# remove first row of labels
		morphometrics = morphometrics[1:,:]

		# only keep swcs in layer
		if layers is not None:
			morphometrics = morphometrics[np.isin(morphometrics[:,1].astype(int),layers)]

		swc_ids = morphometrics[:,0]
		layer_info = morphometrics[:,1].astype(int)
		data = morphometrics[:,2:].astype(float)

		return data, swc_ids, layer_info, labels
		


	def append_distance_to_pia_to_morphometrics(self, in_path, swc_path, pia_coord=-500):

		for structure in ['basal', 'apical', 'basal_apical', 'basal_apical_concat']:

			morphometrics_path = join(in_path, 'morphometrics_' + structure + '.txt')
			morphometrics = np.loadtxt(morphometrics_path, dtype=object)

			distances = ['Distance_to_Pia']
			for idx, swc_morphos in enumerate(morphometrics):
	
				if idx > 0:
					swc_name = '_'.join([swc_morphos[1],swc_morphos[0]])
				
					# soma coordinates are same regardless of structure
					swc_full_path = join(swc_path, 'all_structures', swc_name + '.swc')
					swc = SWC(swc_full_path, build_tree=False)
					soma_y = swc.swc_array[0,3]
					distances.append(abs(pia_coord-soma_y).round(2))
			
			distances =np.array(distances).reshape(-1,1)
			morphometrics = np.hstack((morphometrics, distances))

			
			np.savetxt(morphometrics_path, morphometrics, fmt='%s')


	def runLmeasureSholl(self, independentVar, dependentVar, structure, *, absolute, bins=20):

		"""
		Runs lmeasure to run sholl analysis	
	
		"""
		print()
		print('##################')
		print('Run Lmeasure Sholl')
		print('##################')
		print()

		if structure == 'all':
			structures = params.SWC_STRUCTURE_TYPES
		else:
			structures = [structure]

		# check parameters
		if independentVar not in params.SHOLL_INDEPENDENT_VARIABLES:
			exit('ERROR: Independent Variable must be one of the following ' + str(params.SHOLL_INDEPENDENT_VARIABLES))

		if dependentVar not in params.SHOLL_DEPENDENT_VARIABLES:
			exit('ERROR: Dependent Variable must be one of the following ' + str(params.SHOLL_DEPENDENT_VARIABLES))

		print('Independent Variable:', independentVar)
		print('Dependent Variable:', dependentVar)
		print('Absolute:', absolute)
		print('# bins / bin width:', bins)

		# function specifier
		functionSpecifier = '-' + ','.join(['f' +str(params.LMEASURE_FUNCTION_IDS[dependentVar]),'f' + str(params.LMEASURE_FUNCTION_IDS[independentVar]),'0', str(int(absolute)), str(bins)])

		# output path
		# create directories if non existent
		outpath = join(self.inpath, 'sholl', independentVar + '_' + dependentVar)
		if not exists(outpath):
			mkdir(outpath)
	
		if absolute:
			absString = 'absolute'
		else:
			absString = 'relative'
		outpath = join(outpath, absString)
		if not exists(outpath):
			mkdir(outpath)
		
		# iterate through structures
		for structureType in structures:

			outPathFinal = join(outpath, 'lmout_' + structureType + '.txt')
	
		
			# lmeasure command
			lmeasureCommand = ' '.join([params.LMEASURE_EXE, functionSpecifier, '-s' + outPathFinal])
	
			# swc path
			swcPath = join(self.inpath, 'swcs', structureType, 'normalized')
			swcPaths = ' '.join([join(swcPath,f) for f in listdir(swcPath)])

			lmeasureCommand = ' '.join([lmeasureCommand, swcPaths])

			# run lmeasure from command line
			system(lmeasureCommand)
		
	def PCA_persistent_vectors(self, in_path, structures, layers, *,  normalize=True, n_components=2, annotate=False, save=True, dpi=300, scree=False):


		print()
		print('###')
		print('PCA - Persistent Homology')
		print('###')
		print()

		for structure_type in structures:
			
			structure_path = join(in_path, 'vectors', structure_type)

			for layer in layers:

				vector_files = [x for x in listdir(structure_path) if x.startswith(str(layer))]
				vectors = []

				for file_name in sorted(vector_files):
			
					full_path = join(structure_path, file_name)

					vector = np.genfromtxt(full_path)
					vectors.append(vector)

				data = np.array(vectors)
		
				# center data and divide by sd			
				if normalize:
					data -= data.mean(0)
					data /= data.std(0)


				# scree plot
				if scree:
					num_components_scree = np.min(data.shape)
					pca_scree = PCA(n_components=num_components_scree)
					pca_scree.fit(data)
					x_labels = np.arange(1, num_components_scree+1)
					explained_variances = pca_scree.explained_variance_ratio_
					fig, ax = plt.subplots()
					ax.plot(x_labels,explained_variances,'ro-')
					ax.set_xlabel('Principal Component')
					ax.set_ylabel('Proportion of Variance Explained')
					ax.set_facecolor('lightgray')
					ax.set_title('Scree Plot (' + params.SWC_STRUCTURE_DISPLAY[structure_type] + ')')

					if save:
						out_path = join(in_path, 'pca')
						if not exists(out_path):
							mkdir(out_path)
						out_path_full = join(out_path, structure_type + '_layer_' + str(layer)  + '_scree.png')
						plt.savefig(out_path_full, dpi=dpi)

				# instantiate PCA object
				pca = PCA(n_components=n_components)

				# fit pca
				pca.fit(data)

				# transform data
				transformed_data = pca.transform(data)

				# instantiate plot
				if n_components == 3:
					fig = plt.figure()
					ax = fig.add_subplot(111, projection='3d')
				else:
					fig, ax = plt.subplots()


				if n_components == 3:
					ax.scatter(transformed_data[:,0],transformed_data[:,1],transformed_data[:,2],)
				else:
					ax.scatter(transformed_data[:,0],transformed_data[:,1])

				# labels
				ax.set_xlabel('PC 1')
				ax.set_ylabel('PC 2')

				if n_components == 3:
					ax.set_zlabel('PC 3')

				if n_components == 2:
					ax.set_facecolor('lightgray')
					ax.set_title('Structure: ' + params.SWC_STRUCTURE_DISPLAY[structure_type] + '\nLayer: ' + str(layer) )

			
				ax.grid()
				ax.set_aspect('equal', adjustable='box')
	


				if save:
					out_path = join(in_path, 'pca')
					if not exists(out_path):
						mkdir(out_path)

					out_path_full = join(out_path, structure_type + '_layer_' + str(layer)  + '.png')
					plt.savefig(out_path_full, dpi=dpi)
				else:
					plt.show()

		


		

	def PCA(self, in_path, structures, layer, *,  normalize=True, n_components=2, annotate=False, save=True, dpi=300, scree=False):

		"""
		Principal Component Analysis

		"""

		print()
		print('###')
		print('PCA')
		print('###')
		print()


		# set struture list
		if structures == 'all':
			structures = params.SWC_STRUCTURE_TYPES
			structures.append('basal_apical_concat')
		else:
			structures = [structures]

	
		# iterate through structures
		for structure_type in structures:
	
			# path to morphometrics
			morphometrics_path = join(in_path,'morphometrics_' + structure_type + '.txt')
		
			# load morphometrics
			data, swc_ids, layer_info, labels = self.load_morphometrics(morphometrics_path, layers = [layer])

			# remove any column where all values are identical
			indices_to_remove = np.where(np.all(data == data[0], axis=0))[0]
			data = np.delete(data, indices_to_remove, axis=1)

			# center data and divide by sd			
			if normalize:
				data -= data.mean(0)
				data /= data.std(0)



			# scree plot
			if scree:
				num_components_scree = np.min(data.shape)
				pca_scree = PCA(n_components=num_components_scree)
				pca_scree.fit(data)
				x_labels = np.arange(1, num_components_scree+1)
				explained_variances = pca_scree.explained_variance_ratio_
				fig, ax = plt.subplots()
				ax.plot(x_labels,explained_variances,'ro-')
				ax.set_xlabel('Principal Component')
				ax.set_ylabel('Proportion of Variance Explained')
				ax.set_facecolor('lightgray')
				ax.set_title('Scree Plot (' + params.SWC_STRUCTURE_DISPLAY[structure_type] + ')')

				if save:
					out_path = join(in_path, 'pca')
					if not exists(out_path):
						mkdir(out_path)
					out_path_full = join(out_path, structure_type + '_layer_' + str(layer)  + '_scree.png')
					plt.savefig(out_path_full, dpi=dpi)
	
			# instantiate PCA object
			pca = PCA(n_components=n_components)

			# fit pca
			pca.fit(data)

			# transform data
			transformed_data = pca.transform(data)


			# instantiate plot
			if n_components == 3:
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
			else:
				fig, ax = plt.subplots()

			# iterate through layer to plot
			#for layer in layers:

				#layerData = transformedData[morphometrics[:,1].astype(int) == layer]
				#if n_components == 3:
				#	ax.scatter(layerData[:,0],layerData[:,1],layerData[:,2], color=params.LAYERS_COLORS_DICT[layer],label=params.LAYERS_DISPLAY[layer])
				#else:
				#	ax.scatter(layerData[:,0],layerData[:,1], color=params.LAYERS_COLORS_DICT[layer],label=params.LAYERS_DISPLAY[layer])

				#if annotate:
				#	for idx, coord in enumerate(layerData,start=1):
		    		#		ax.annotate(idx, (coord[0], coord[1]))

			if n_components == 3:
				ax.scatter(transformed_data[:,0],transformed_data[:,1],transformed_data[:,2],)
			else:
				ax.scatter(transformed_data[:,0],transformed_data[:,1])

			#if annotate:
			#	for idx, coord in enumerate(layerData,start=1):
	    		#		ax.annotate(idx, (coord[0], coord[1]))

			# labels
			ax.set_xlabel('PC 1')
			ax.set_ylabel('PC 2')

			if n_components == 3:
				ax.set_zlabel('PC 3')

			if n_components == 2:
				ax.set_facecolor('lightgray')
				ax.set_title('Structure: ' + params.SWC_STRUCTURE_DISPLAY[structure_type] + '\nLayer: ' + str(layer) )

			# legend
			#if layer == [2,5,6]:
			#	ax.legend(loc='lower center',bbox_to_anchor=(0.5, -.4))
			#else:
			#	ax.legend(loc='lower center',bbox_to_anchor=(0.5, -.3))
			#fig.subplots_adjust(bottom=0.25)
			
			ax.grid()
			ax.set_aspect('equal', adjustable='box')
	
			#if n_components == 3:
			#	azimuth = 156
			#	elevation = -96
			#	ax.view_init(elev=elevation, azim=azimuth)
			#	outpath = join(self.inpath, 'pca', structureType + '_' + ''.join([str(l) for l in layers]) + '_' + str(azimuth) + '_'+str(elevation) + '.png')
			#	plt.savefig(outpath,dpi=dpi)



			if save:
				
				out_path = join(in_path, 'pca')
				if not exists(out_path):
					mkdir(out_path)

				out_path_full = join(out_path, structure_type + '_layer_' + str(layer)  + '.png')
				plt.savefig(out_path_full, dpi=dpi)
			else:
				plt.show()
			


	def boxPlots(self, structure, *, normalizeByLayer=None, save=True, dpi=300, legend=True, tickSize=15):

		"""
		box plots of morphometrics
		
		"""
		print()
		print('#########')
		print('Box Plots')
		print('#########')
		print()

		# set struture list
		if structure == 'all':
			structures = params.SWC_STRUCTURE_TYPES
			structures.append('basal_apical_concat')
		else:
			structures = [structure]

		# iterate through structures
		for structureType in structures:

			# load morphometrics for structure
			morphometricsPath = join(self.inpath,'morphometrics_'+structureType + '.txt')

			morphometrics = np.loadtxt(morphometricsPath,dtype=object)
			morphometricLabels = morphometrics[0,2:]
			morphometrics = morphometrics[1:,:]

			# remove first 2 columns (swc and layer)
			data = np.array(morphometrics[:,2:], dtype=np.float64)

			# normalize by layer
			if not normalizeByLayer is None:
				data /= data[morphometrics[:,1].astype(int) == normalizeByLayer].mean(axis=0)


			# iterate through metrics
			for i,metric in enumerate(morphometricLabels):

				bplotData = {}

				# box plot for each layer
				for layer in params.LAYERS:
					bplotData[layer] = data[morphometrics[:,1].astype(int) == layer,i]

				# wilcoxon rank sum test (pairwise between layers 2-5 and 2-6)
				_, pString5 = utils.wilcoxon(bplotData[2], bplotData[5])
				_, pString6 = utils.wilcoxon(bplotData[2], bplotData[6])

				# create plot
				fig, ax  = plt.subplots()

				# set title
				metricName = metric.split(':')[0]
				value = metric.split(':')[1]
				ax.set_title(params.BOX_PLOT_METRIC_DISPLAY_DICT[metricName][value])

				# make box plot
				bplot = ax.boxplot( [v for (k,v) in sorted(bplotData.items())],whis=float('inf'),patch_artist=True)

				for j,layer in enumerate(sorted(bplotData)):
					bplot['boxes'][j].set_facecolor(params.LAYERS_COLORS_DICT[layer])

		
				# place siginifiance 
				ax.text(2-params.BOXPLOT_PSHIFT[pString5],max(bplotData[5]),pString5)
				ax.text(3-params.BOXPLOT_PSHIFT[pString6],max(bplotData[6]),pString6)

				# plot options
				ax.set_xticks([])
				ax.set_facecolor('lightgray')
				ax.grid()
				for tick in ax.yaxis.get_major_ticks():
					tick.label1.set_fontsize(tickSize)
					tick.label1.set_fontweight('bold')

				# legend
				if legend:
					legendElements = []
					for layer in params.LAYERS:
						legendElements.append(Patch(facecolor=params.LAYERS_COLORS_DICT[layer], edgecolor='black',label='Color Patch')) 
		
					layerStringArray = []
					for layer in params.LAYERS:
						layerStringArray.append(params.LAYERS_DISPLAY[layer])
					ax.legend(legendElements,layerStringArray, loc='lower center',bbox_to_anchor=(0.5, -.3))
					fig.subplots_adjust(bottom=0.25)
	
				if save:
					if normalizeByLayer is None:
						outpath = join(self.inpath,'box_plots','absolute', structureType, metric + '.png')
					else:
						outpath = join(self.inpath,'box_plots','normalized_by_layer_' + str(normalizeByLayer),structureType,metric + '.png' )
					plt.savefig(outpath,dpi=dpi)
				else:
					plt.show()



	def parseLmoutSholl(self, dataPath):

		"""
		parses output from lmeasure and store in dictionary

		"""

		# set up dictionary to store data
		data = {}
		for layer in params.LAYERS:
			data[layer] = {}

		# parse lmeasure output
		with open(dataPath, 'r') as fp:
			lines = fp.readlines()

			# iterate through lines in fine
			for i,line in enumerate(lines):

				lineList = line.split()
			
				# ignore empty lines
				if len(lineList) == 0:
					continue
			
				# lines alternate independent - dependent
				if i % 2 == 0:
					swcPath = lineList.pop(0)
					layer, _, _, brainSWCID  = utils.getSWCBrainAndID(swcPath, layerIncluded=True)
					metric = lineList.pop(0).split(':')[-1]
					data[layer][brainSWCID] = {}
				else:
					metric = lineList.pop(0)

				lineData = np.array(lineList, dtype=float)
				data[layer][brainSWCID][metric] = lineData

		return data
				

	def sholl(self, independentVar, dependentVar, structure, absolute, bins=20, * , save=False, barplot=False, dpi=300):

		"""
		Sholl Analysis		

		"""

		print()
		print('##############')
		print('Sholl Analysis')
		print('##############')
		print()

		# set struture list
		if structure == 'all':
			structures = params.SWC_STRUCTURE_TYPES
		else:
			structures = [structure]


		# check parameters
		if independentVar not in params.SHOLL_INDEPENDENT_VARIABLES:
			exit('ERROR: Independent Variable must be one of the following ' + str(params.SHOLL_INDEPENDENT_VARIABLES))

		if dependentVar not in params.SHOLL_DEPENDENT_VARIABLES:
			exit('ERROR: Dependent Variable must be one of the following ' + str(params.SHOLL_DEPENDENT_VARIABLES))
	
		
		# get lmeasure outpath
		if absolute:
			absString = 'absolute'
		else:
			absString = 'relative'

		print('Independent Variable:', independentVar)
		print('Dependent Variable:', dependentVar)
		print('Absolute/Relative:', absString)
		print('Structure:', structure)
		print('Bar Plot:', barplot)

		# iterate through structures
		for structureType in structures:

			dataPath = join(self.inpath, 'sholl', independentVar + '_' + dependentVar,absString,'lmout_' + structureType + '.txt' )

			# parse lmeasure output
			parsedData = self.parseLmoutSholl(dataPath)

			# organize data int list of arrays
			data = {}
			for layer in parsedData:
				data[layer] = {'indVar': [], 'depVar': []} 

			for layer in parsedData:
				for swc in parsedData[layer]:
					data[layer]['indVar'].append(np.array(parsedData[layer][swc][independentVar]))
					data[layer]['depVar'].append(np.array(parsedData[layer][swc][dependentVar]))
		
			# if relative, then normalize by sum
			if not absolute:
				for layer in data:
					for i in range(len(data[layer]['indVar'])):
						data[layer]['indVar'][i] = np.linspace(0,100,bins+1)
					for rowData in data[layer]['depVar']:
						rowData /= rowData.sum()
						rowData *= 100

			# reorganize data so 'timepoints' grouped
			dataGrouped = {}
			for layer in data:
				dataGrouped[layer] = {} 

			for layer in data:
				layerData = data[layer]
			
				# find max number of data points
				maxLength = max([len(l) for l in layerData['depVar']])
				dataGrouped[layer]['depVar'] = [[] for i in range(maxLength)]
				dataGrouped[layer]['indVar'] = layerData['indVar'][np.argmax(np.array([len(l) for l in layerData['indVar']]))]
			
				for swcRow in layerData['depVar']:
					for i, d in enumerate(swcRow):
						dataGrouped[layer]['depVar'][i].append(d)
				

			# plot
			fig, ax = plt.subplots()

			if barplot:

				# kruskal wallace
				pStrings = []
				maxNum= max([len(v['depVar']) for (k,v) in dataGrouped.items()])
				for i in range(maxNum):
					pLists = []
					for layer in dataGrouped:
						if len(dataGrouped[layer]['depVar']) > i:
							if np.any(dataGrouped[layer]['depVar'][i]):
								pLists.append(dataGrouped[layer]['depVar'][i])
					pStrings.append(utils.kruskalWallis(pLists))			

				# calculate medians/means/stds

				stds = {}
				for layer in dataGrouped:
					stds[layer] = ([0]*len(dataGrouped[layer]['depVar']), np.array([np.std(np.array(l)) for l in dataGrouped[layer]['depVar']]))
					dataGrouped[layer]['depVar'] = np.array([np.median(np.array(l)) for l in dataGrouped[layer]['depVar']])


				if absolute:
					width = 5
				else:
					width = 1
		
				bars = {}
				for i,layer in enumerate(dataGrouped):
					layerData = dataGrouped[layer]
					bars[layer] =ax.bar(layerData['indVar'] + i*width,layerData['depVar'],width, yerr=stds[layer],error_kw={'elinewidth':.5,'capsize':0,'capthick':.5,'lolims':True } , color=params.LAYERS_COLORS_DICT[layer], label=params.LAYERS_DISPLAY[layer])


					bars[layer].errorbar.lines[1][0].set_marker('_')
					bars[layer].errorbar.lines[1][0].set_markersize(2)



				# display significance above bars
				for i in range(max([len(v) for (k,v) in bars.items()])):
				
					heights = []
					xs = []
					for layer in dataGrouped:
						if len(dataGrouped[layer]['depVar']) > i:
							heights.append(bars[layer][i].get_height() + stds[layer][1][i])
							xs.append(bars[layer][i].get_x() + bars[layer][i].get_width() / 2)
	
					maxH = np.array(heights).max()
					xCoord = np.array(xs).mean()

					ax.annotate('{}'.format(pStrings[i]), xy=(xCoord, maxH), xytext=(0, 0),textcoords="offset points",ha='center', va='bottom')


			else:

				# calculate medians/means
				for layer in dataGrouped:
					dataGrouped[layer]['depVar'] = np.array([np.array(l).mean() for l in dataGrouped[layer]['depVar']])
			

				for layer in dataGrouped:
					layerData = dataGrouped[layer]
					plt.plot(layerData['indVar'], layerData['depVar'], linestyle='-', marker='o', color=params.LAYERS_COLORS_DICT[layer], label=params.LAYERS_DISPLAY[layer])
			
			# plot settings
			ax = plt.gca()
			ax.set_facecolor('lightgray')

			if not barplot:
				ax.grid()
		
			# set title and labels
			if absolute:
				#ax.set_title('Sholl Diagram : ' + params.SWC_STRUCTURE_DISPLAY[structure])
				ax.set_xlabel(params.METRIC_DISPLAY_DICT[independentVar] + ' (um)')	
				ax.set_ylabel(params.METRIC_DISPLAY_DICT[dependentVar] + ' (um)')
			else:
				#ax.set_title('Sholl-like Diagram : ' + params.SWC_STRUCTURE_DISPLAY[structure])
				ax.set_xlabel('Relative ' + params.METRIC_DISPLAY_DICT[independentVar] + ' (%)')		
				ax.set_ylabel('Relative ' + params.METRIC_DISPLAY_DICT[dependentVar] + ' (%)')

			plt.legend()

			if save:
				if barplot:
					outpath = join(self.inpath, 'sholl', independentVar + '_' + dependentVar,absString, structureType + '_barplot.png')
				else:
					outpath = join(self.inpath, 'sholl', independentVar + '_' + dependentVar,absString, structureType + '.png')
				plt.savefig(outpath, dpi=dpi)
			else:		
				plt.show()





	def calculateDistances(self, layerA,layerB=None,*,vectorDict):


		"""
		Calculate pairwise distances for persistent homology analysis

		"""

		distances = []
	
		if layerB is None:

			vectorList = list(vectorDict[layerA].values())

			for i,vectorA in enumerate(vectorList):
				for j in range(i+1,len(vectorList)):
					vectorB = vectorList[j]
					distances.append(utils.arcCosDissimilarity(vectorA, vectorB))
		else:
			vectorListA = list(vectorDict[layerA].values())
			vectorListB = list(vectorDict[layerB].values())

			for vectorA in vectorListA:
				for vectorB in vectorListB:
					distances.append(utils.arcCosDissimilarity(vectorA, vectorB))

		return distances

	def persistentHomology(self, inpath, outpath=None):

		# vector files
		vectorFiles = [f for f in listdir(inpath) if int(f.split('_')[1]) in self.layers]

		# create vector dictionary
		vectorDict = {}
		for layer in self.layers:
			vectorDict[layer] = {}	

		# load vectors
		for vectorFile in vectorFiles:
		
			vector = np.genfromtxt(join(inpath,vectorFile))
			layer = int(vectorFile.split('_')[1])
			brainSWCID = '_'.join(vectorFile.split('_')[2:])[:-4]
			vectorDict[layer][brainSWCID] = vector

		# intra layer pairwise distances
		intraLayerDistances = {}
		for layer in self.layers:
			intraLayerDistances[layer] = self.calculateDistances(layer,vectorDict=vectorDict)
		

		# inter layer distances
		interLayerDistances = {}

		# need to get all possible pairs of layers
		for i in range(len(self.layers)-1):
			for j in range(i+1, len(self.layers)):
				if i != j:
					layersKey = (self.layers[i],self.layers[j])
					layersKey = frozenset(layersKey)			
					interLayerDistances[layersKey] = self.calculateDistances(self.layers[i],self.layers[j],vectorDict=vectorDict)
		
		# create plot
		fig,axs = plt.subplots(1,len(self.layers)+1,sharex=True, figsize=(25,10))
		plt.xlim(0,1)


	
		# All layers

		# append all inter and intra layers
		allInterLayer = []

		for layer, distanceList in interLayerDistances.items():
			allInterLayer += distanceList

		allIntraLayer = []
		for layer, distanceList in intraLayerDistances.items():
			allIntraLayer += distanceList

			
		# All Layers
		fig.axes[0].set_title('All Layers')
		bplotAll = fig.axes[0].boxplot([allInterLayer,allIntraLayer],whis=float('inf'),patch_artist=True,vert=False,medianprops=dict(color='black'))
		bplotAll['boxes'][0].set_facecolor('white')
		bplotAll['boxes'][1].set_facecolor('lightcoral')
		p, p_string = utils.wilcoxon(allInterLayer,allIntraLayer)
		fig.axes[0].text(.5,1.5,p_string)

		# Individual Layers
		for i, layer in enumerate(self.layers, start=1):

			# get inter layer distances
			localInterLayerDistances = []
			for localLayer in self.layers:
				if localLayer != layer:
					localInterLayerDistances += interLayerDistances[frozenset((layer,localLayer))]

			fig.axes[i].set_title('Layer ' + str(layer) )
			bplot = fig.axes[i].boxplot([localInterLayerDistances,intraLayerDistances[layer]],whis=float('inf'),patch_artist=True,vert=False,medianprops=dict(color='black'))
			bplot['boxes'][0].set_facecolor('white')
			bplot['boxes'][1].set_facecolor('lightcoral')
			_, pString = utils.wilcoxon(localInterLayerDistances,intraLayerDistances[layer])
			fig.axes[i].text(.3,1.5,pString)



		for i in range(len(axs)):
			fig.axes[i].set_yticks([])
			fig.axes[i].set_facecolor('lightgray')
			fig.axes[i].grid()


		fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.2)
		legend_elements = [Patch(facecolor='lightcoral', edgecolor='black',label='Color Patch'),Patch(facecolor='white', edgecolor='black',label='Color Patch')]	
		fig.legend(legend_elements,['within','between'],loc=(.41,.04))
		
		if outpath is None:
			plt.show()
		else:
			plt.savefig(outpath)

	def hierarchicalClustering(self, outpath=None):
		
		data = np.array(self.morphometrics[:,2:], dtype=np.float64)
		swcLayers = self.morphometrics[:,1]		

		plt.figure(figsize=(10, 7))
		plt.title("Dendograms")
		dend = shc.dendrogram(shc.linkage(data, method='ward'),labels=swcLayers, )	
		
		if outpath is None:
			plt.show()
		else:
			plt.savefig(outpath)


	def getMopulSWCs(self, inpath, outpath):
		
		swcPaths = [join(inpath, f) for f in listdir(inpath)]
		
		for swcPath in sorted(swcPaths):
	

			swc  = SWC(swcPath)
			swc.extractTypeTree(['basal dendrite'])

			layer = swc.getMOPulLayer()

			if layer != 0:
				swc.saveSWC(join(outpath,str(layer) + '_' + swcPath.split('/')[-1]))

	def collapseSoma(self, inpath, outpath):

		swcPaths = [join(inpath, f) for f in listdir(inpath)]
		
		for swcPath in sorted(swcPaths):

			swcName = swcPath.split('/')[-1]
			swc = SWC(swcPath)
			swc.collapseSoma()
			swc.saveSWC(join(outpath,swcName))
			print()

	def saveInfoForLabelMarkup(self, inpath, outpath, downsampling, pixelConnected):

		swcPaths = [join(inpath, f) for f in listdir(inpath)]
		
		imageShape = [int(x/downsampling) for x in params.CCF_TIF_SHAPE]

		for swcPath in sorted(swcPaths):
			
			swc = SWC(swcPath)
			swc.saveBranchInfoForStructureLabelMarkup(outpath, downsampling=downsampling, normalize=True)
			swc.saveSkeleton(outpath,downsampling=downsampling,compressIds=True, pixelConnected=pixelConnected)
		
	def UMAP(self, structures='all', layers='all',*, normalize=True, save=False, dpi=300, n_neighbors=3, min_dist=.2, n_components=2, metric='euclidean', annotate=False, random_state=42):

		"""
		UMAP

		"""

		import umap

		print()
		print('####')
		print('UMAP')
		print('####')
		print()


		# set struture list
		if structures == 'all':
			structures = params.SWC_STRUCTURE_TYPES
			structures.append('basal_apical_concat')

		# set layers list
		if layers == 'all':
			layers = [[2],[5],[6],[2,5,6]]


		for structure in structures:

			for layer in layers:

				data, swcIDs, layerInfo, labels = self.loadMorphometrics(structure, layer)

				# remove any column where all values are identical
				indicesToRemove = np.where(np.all(data == data[0], axis=0))[0]
				data = np.delete(data, indicesToRemove, axis=1)

				# zscore
				if normalize:
					layerData = zscore(data)

				# run umap
				reducer = umap.UMAP(n_neighbors=n_neighbors,min_dist=min_dist,n_components=n_components,metric=metric,random_state=random_state)
				embedding = reducer.fit_transform(layerData)
				
				fig, ax = plt.subplots()
				
				idx = 1
				for local_layer in layer:
	

					localLayerData = embedding[layerInfo == local_layer]
					ax.scatter(localLayerData[:,0],localLayerData[:,1],color=params.LAYERS_COLORS_DICT[local_layer],label=params.LAYERS_DISPLAY[local_layer])

					if annotate:
						for coord in localLayerData:
							ax.annotate(idx, (coord[0], coord[1]), fontsize=5)
							idx += 1


				# labels
				ax.set_xlabel('Component 1')
				ax.set_ylabel('Component 2')
				ax.set_facecolor('lightgray')
				ax.set_title('UMAP - Morphometrics (' + params.SWC_STRUCTURE_DISPLAY[structure] + ')')

				# legend
				if layer == [2,5,6]:
					ax.legend(loc='lower center',bbox_to_anchor=(0.5, -.4))
				else:
					ax.legend(loc='lower center',bbox_to_anchor=(0.5, -.3))

				fig.subplots_adjust(bottom=0.25)
				ax.grid()
				ax.set_aspect('equal', adjustable='box')

				if save:
					outpath = join(self.inpath, 'umap','morphometrics', structure + '_' + ''.join([str(l) for l in layer]) + '_annot.png')
					plt.savefig(outpath, dpi=dpi)
				else:
					plt.show()


	def logRegFindOptimalC(self, X,y):
		c = params.MAX_C_L1

		acc = 1
		while True:
			clf = LogisticRegression(penalty='l1',solver='liblinear',multi_class='auto',C=c,max_iter=1000).fit(X,y)
			acc = clf.score(X,y)

			if acc == 1:
				c -= .1
			else:
				if c != params.MAX_C_L1:
					c += .1
				break
			
		return c
	
	def bootstrapLogisticRegression(self,X,y,n=10000):
	
		"""
		Performs bootstrapping for logistic regression to find p value

		"""

		print()
		print('Bootstrapping...')
		print()

		coeffs = []

		for i in range(n):		

			# sample randomly with replacement
			randIndices = np.random.randint(len(X),size=len(X))
			X = X[randIndices]

			# z score data
			X = zscore(X)

			# Replace any column that contains NA with 0
			colsWithNan = np.where(np.any(np.isnan(X),axis=0))
			if len(colsWithNan) > 0:
				X[:,colsWithNan[0]] = 0

			# find optimal c value
			c = self.logRegFindOptimalC(X,y)

			# train logistic regression
			clf = LogisticRegression(penalty='l1',solver='liblinear',multi_class='auto',C=c, max_iter=1000).fit(X,y)

			coeffs.append(clf.coef_)

		coeffs = np.array(coeffs)

		return coeffs



	def logisticRegression(self, structure, layer, *, normalize=True, clusters=None):

		"""
		Logistic Regression

		"""
		
		print()
		print('####################')
		print('Logistic Regression')
		print('####################')
		print()
		print('Structure:', structure)
		print('Layer:', layer)
		print()

		# load morphometrics for structure
		data, swcIDs, layerInfo, labels = self.loadMorphometrics(structure, layer)
		
		# z score data
		if normalize:
			X = zscore(data)
	
		# Replace any column that contains NA with 0

		colsWithNan = np.where(np.any(np.isnan(X),axis=0))
		if len(colsWithNan) > 0:
			X[:,colsWithNan[0]] = 0


		# If running all three layers, user layer info as clusters
		# Otherwise need to get cluster info

		if layer == [2,5,6]:
			#y = layerInfo
			y = [0]*len(layerInfo)
			for clusterID, cluster in enumerate(clusters,1):
				for swcID in cluster:
					y[swcID-1] = clusterID
		else:
			y = [0]*len(layerInfo)
			for clusterID, cluster in enumerate(params.L1_CLUSTERING[layer],1):
				for swcID in cluster:
					y[swcID-1] = clusterID
		y = np.array(y)

		# bootstrap
		bootstrapCoeffs = self.bootstrapLogisticRegression(data,y)

		# find highest c value that still gives 100% accuraacy
		c = self.logRegFindOptimalC(X,y)

		# train logistic regression
		clf = LogisticRegression(penalty='l1',solver='liblinear',multi_class='auto',C=c,max_iter=1000).fit(X,y)	
	
		# create plot
		fig, axs = plt.subplots(len(clf.classes_),1)
		
		# get p values
		pValues = []
		for i,cluster in enumerate(clf.classes_):
			coeffs = clf.coef_[i]
			for j, coeff in enumerate(clf.coef_[i]):
				bCoeffs = bootstrapCoeffs[:,i,j]
				p = sum(abs(bCoeffs) >= abs(coeff))/len(bCoeffs)
				pValues.append(p)

		# multiple comparisons corrections
		pValues = multipletests(pValues, method='bonferroni')[1]
		pValues = pValues.tolist()
	

		# iterate through clusters
		for i,cluster in enumerate(clf.classes_):
		
			# get coefficients for cluster
			coeffs = clf.coef_[i]
			

			for j, coeff in enumerate(clf.coef_[i]):


				p = pValues.pop(0)
				pString = utils.pValueToString(p)
		

				if coeff == 0:	
					axs[i].scatter(j,coeff,s=10,color='black',facecolors='none')
				elif coeff > 0:
					axs[i].scatter(j,coeff,s=10,color='blue')
					axs[i].text(j,coeff,pString,size=8)
				else:
					axs[i].scatter(j,coeff,s=10,color='red')
					axs[i].text(j,coeff-2.2,pString,size=8)

			# find coeff with maximum magnitude
			maxCoeff = max(abs(coeffs))


			# get accuracies (should all be 100%)
			clusterAcc = clf.score(X[y==cluster],y[y==cluster])
			axs[i].text(-1.6, -maxCoeff-3, 'Accuracy: ' + str(clusterAcc.round(2)) + ' (' +str(sum(y==cluster))+ ')', fontsize=7)
		

			# axis options
			axs[i].set_ylim([-maxCoeff-3, maxCoeff+3])
			axs[i].set_xticks([])
			axs[i].yaxis.grid()
			axs[i].set_facecolor('lightgray')


		# plot options
		axs[0].set_title(params.SWC_STRUCTURE_DISPLAY[structure] + ' C=' + str(np.round(c,1)))
		axs[-1].set_xticks(np.arange(len(labels)))
		axs[-1].set_xticklabels(labels,rotation=90)

		# morphometric labels
		for i,xtick in enumerate(axs[-1].get_xticklabels()):
			if i < 18:
				xtick.set_color('green')
			else:
				xtick.set_color('orange')
		#fig.text(0.25, 0.05, 'Basal', fontsize=14,color='green')
		#fig.text(0.75, 0.05, 'Apical', fontsize=14, color='orange')
		fig.subplots_adjust(bottom=0.5)

		# show plot
		plt.show()



	def boxPlotsClustering(self, inpath, layer, structure, save=False, dpi=300):

		# load morphometrics
		data, swcIDs, layerInfo, morphoLabels = self.loadMorphometrics(structure, layer)
		
		# split into clusters
		clusterData = {}
		for clusterID, clusterIds in enumerate(params.L1_CLUSTERING[layer], start=1):
			clusterIDXs = [x-1 for x in clusterIds]
			clusterData[clusterID] = data[clusterIDXs]

		# iterate through metrics
		for idx, morphoLabel in enumerate(morphoLabels):
			
			# create box plot
			bplotData = {}

			# create plot
			fig, ax  = plt.subplots()

			# set title
			metricName = morphoLabel.split(':')[0]
			value = morphoLabel.split(':')[1]
			ax.set_title(params.BOX_PLOT_METRIC_DISPLAY_DICT[metricName][value])

			# box plot for each layer
			for clusterID in sorted(clusterData):
					
				localClusterData = clusterData[clusterID]

				bplotData[clusterID] = localClusterData[:,idx]


			# pairwise wilcoxon rank sum
			sigString = ''
			for i in range(1,len(clusterData)+1):
				for j in range(i+1,len(clusterData)+1):
					_, pString = utils.wilcoxon(clusterData[i][:,idx], clusterData[j][:,idx])
					if pString == '':
						pString = 'ns'
					sigString += 'C' + str(i) + '-C' + str(j) + ' : ' + pString
					if i != len(clusterData)-1 or j != len(clusterData):
						sigString += '\n'
			
			# text box for significance levels
			props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
			plt.subplots_adjust(right=.8)
			ax.text(.82, .6, sigString, transform=plt.gcf().transFigure,fontsize=14, verticalalignment='top', bbox=props)
			
			# make box plot
			bplot = ax.boxplot( [v for (k,v) in sorted(bplotData.items())],whis=float('inf'),patch_artist=True)

			plt.xticks(np.arange(len(clusterData))+1,['C'+str(x) for x in list(np.arange(len(clusterData))+1)])
			ax.set_facecolor('lightgray')
			ax.yaxis.grid()

			if save:
				structString = '(Basal)' if idx < 18 else '(Apical)'	
				outpath = join(inpath,'box_plots_layer_' + str(layer), morphoLabel + '_' + structString)
				plt.savefig(outpath, dpi=dpi)
			else:
				plt.show()


		
	def somaDistances(self, swcPath, layer, * , dpi=300):

		# get swc filenames
		swcFilenames = sorted([f for f in listdir(swcPath) if f.startswith(str(layer))])

		# get cluster info
		clusterInfo = params.L1_CLUSTERING[layer]

		# separate into clusters
		somas = {}

		for i, cluster in enumerate(clusterInfo, start = 1):

			somas[i] = []			

			for id_ in cluster:

				swcArray = np.genfromtxt(join(swcPath,swcFilenames[id_-1]))
				soma = swcArray[0,2:5]
				if soma[0] > params.CCF_WIDTH/2:
					soma[0] = params.CCF_WIDTH - soma[0]
				somas[i].append(swcArray[0,2:5])

		for clusterID, clusterSomas in somas.items():

			fig, ax = plt.subplots()

			bplotData = []

			xlabels = []

			# intra cluster (remove diagonal)
			pairwiseDists = pairwise_distances(clusterSomas)
			pDists = pairwiseDists[~np.eye(pairwiseDists.shape[0],dtype=bool)] 
			bplotData.append(pDists)

			xlabels.append(str(clusterID)+'-'+str(clusterID)) 
			
			for clusterID_other, clusterSomas_other in somas.items():

				if clusterID != clusterID_other:
					pairwiseDists = pairwise_distances(clusterSomas, clusterSomas_other)
					pDists = pairwiseDists.flatten()
				
					bplotData.append(pDists)

					xlabels.append(str(clusterID)+'-'+str(clusterID_other))

			ax.boxplot(bplotData, whis=float('inf'), patch_artist=True)

			ax.set_facecolor('lightgray')
			ax.yaxis.grid()

			plt.xticks(np.arange(1,len(bplotData)+1),xlabels)
			
			plt.xlabel('Cluster')
			plt.ylabel('Distance (um)')
		
			outpath = join(self.outpath,'soma_distances','layer_' + str(layer) + '_cluster_' + str(clusterID))
			plt.savefig(outpath, dpi=dpi)


	def arbor_densities_to_binary_images(self):

		arbor_densities_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/arbor_density/from_olga/'
		arbor_density_types = ['soma_centered', 'soma_positioned']
		out_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/arbor_density/binary_images/'
		layers = [2, 5, 6]

		for arbor_density_type in arbor_density_types:
		
			for layer in layers:
				
				arbor_density_path = join(arbor_densities_path, 'layer_' + str(layer) + '_arbor_densities', arbor_density_type,)
				arbor_density_file_names = listdir(arbor_density_path)
	
				for arbor_density_file_name in arbor_density_file_names:
		
					# extract swc name
					structure_type = arbor_density_file_name[:-4].split('_')[2]
					swc_name = str(layer) + '_' + '_'.join(arbor_density_file_name[:-4].split('_')[3:])
		
					arbor_density_full_path = join(arbor_density_path, arbor_density_file_name)

					# load arbor density
					arbor_density = np.loadtxt(arbor_density_full_path, delimiter=',', dtype=np.float64)

					# make binary for 8 bit image
					arbor_density[arbor_density > 0] = 255
					arbor_density = arbor_density.astype(np.uint8)

					# save as image
					out_path_full = join(out_path, arbor_density_type, structure_type, swc_name + '.tif')
					imsave(out_path_full, arbor_density)

	def create_2d_swc_plots(self, *, center_soma_y):

		swc_base_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/swcs/'
		out_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/'
		dirs = ['removed_nan_swcs', 'removed_nan_nodes']
		if center_soma_y:
			out_dir = 'soma_centered'
		else:
			out_dir = 'soma_positioned'
			y_lim = [-500, 500]
		

		for dir_name in dirs:

			swc_dir = join(swc_base_path, dir_name, 'all_structures')

			swc_names = listdir(swc_dir)

			# find x lims
			min_x_list, max_x_list = [], []
			for swc_name in swc_names:
				swc_full_path = join(swc_dir, swc_name)
				swc = SWC(swc_full_path, build_tree=False)
				min_x_list.append(swc.swc_array[:,2].min())
				max_x_list.append(swc.swc_array[:,2].max())
			x_lim = [np.min(min_x_list), np.max(max_x_list)]

			if center_soma_y:
				min_y_list, max_y_list = [], []
				for swc_name in swc_names:
					swc_full_path = join(swc_dir, swc_name)
					swc = SWC(swc_full_path, build_tree=False)
					min_y_list.append(swc.swc_array[:,3].min() - swc.swc_array[0,3])
					max_y_list.append(swc.swc_array[:,3].max() - swc.swc_array[0,3])
				y_lim = [np.min(min_y_list), np.max(max_y_list)]

			for swc_name in swc_names:
		
				swc_full_path = join(swc_dir, swc_name)
				out_path_full = join(out_path, dir_name, '2d_images', out_dir, swc_name[:-4] + '.png')
			
				swc = SWC(swc_full_path)

				if center_soma_y:
					swc.center_around_soma()				

				swc.plot(save=True, out_path=out_path_full, x_lim=x_lim, y_lim=y_lim)

		
				 

					

				
if __name__ == '__main__':

	# Consolidate SWCs sent from Olga
	#in_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/'
	#in_dir_list = ['from_olga_batch_1', 'from_olga_batch_2']
	#out_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/'
	#Morphometric_Analysis.consolidate_swcs(in_path, in_dir_list, out_path)

	# create analaysis object
	morpho = Morphometric_Analysis()

	# load SWCs
	#swc_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/removed_nan_swcs/all_structures/'
	#morpho.load_swcs(swc_path, registered=True, orientation='coronal')

	# extract structures
	#structures=['apical dendrite']
	#morpho.extract_structures(structures=structures, normalize=False)

	# save
	#out_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/removed_nan_swcs/apical/'
	#morpho.save_swcs(out_path=out_path)

	# get morphometrics from lmeasure
	#swc_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/removed_nan_swcs/'
	#out_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/removed_nan_swcs/'
	#morpho.run_lmeasure_morphometrics(swc_path, out_path)

	# parse morphometrics output from lmeasure
	lmout_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/removed_nan_nodes/'
	json_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/lmeasure_morphometrics.json'
	morpho.parse_lmeasure_morphometrics(lmout_path, json_path)

	# append distance to pia morphometric
	#in_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/removed_nan_swcs/'
	#swc_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/removed_nan_swcs/'
	#morpho.append_distance_to_pia_to_morphometrics(in_path, swc_path)

	# PCA
	#in_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/removed_nan_nodes/without_soma_depth/'
	#structures = 'basal_apical_concat'
	#layer = 6
	##n_components = 3
	#save = False
	#scree = False
	#morpho.PCA(in_path, structures=structures, layer=layer, n_components=n_components, save=save, scree=scree)


	in_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/removed_nan_nodes/persistent_homology/'
	structures = ['basal', 'apical', 'basal_apical_concat']
	layers = [2,5,6]
	n_components = 2
	save = True
	scree = True
	morpho.PCA_persistent_vectors(in_path, structures=structures, layers=layers, n_components=n_components, save=save, scree=scree)

	# Box Plots
	#swcAnalysis.boxPlots(structure='all', normalizeByLayer=2, legend=False)

	# Sholl 
	#swcAnalysis.runLmeasureSholl('EucDistance','Length',structure='all', absolute=True) 
	#swcAnalysis.sholl('EucDistance','Length',structure='all', absolute=True, save=True, barplot=True)

	#swcAnalysis.getMopulSWCs(inpath='/data/elowsky/OLST/swc_analysis/paper/analysis/janelia/original_swcs/',outpath='/data/elowsky/OLST/swc_analysis/paper/analysis/janelia/mopul_swcs/')

	#swcAnalysis.saveInfoForLabelMarkup(inpath='/data/elowsky/OLST/swc_analysis/paper/analysis/janelia/mopul_swcs_dendrites_collapsed_soma/', outpath='/data/elowsky/OLST/swc_analysis/paper/analysis/janelia/mopul_swcs_dendrites_collapsed_soma_structure_label_markup/', downsampling=1, pixelConnected=True)

	# UMAP
	#morpho.UMAP(structures='all', layers=[[2,5,6]], save=True, annotate=True)


	# Logistic Regression
	#morpho.logisticRegression(structure='basal', layer=[2,5,6], clusters = params.L1_CLUSTERING_ALL_LAYERS_BASAL)
	
	#inpath = '/data/elowsky/OLST/swc_analysis/paper/analysis/coronal_registered/logistic_regression/arbor_density/soma_not_centered/'
	#swcAnalysis.boxPlotsClustering(inpath=inpath, layer=6, structure='basal_apical_concat', save=True)

	#swcPath = '/data/elowsky/OLST/swc_analysis/paper/analysis/coronal_registered/swcs/basal_apical/not_normalized/'
	#morpho.somaDistances(swcPath, layer=6)

	#morpho.arbor_densities_to_binary_images()

	#morpho.create_2d_swc_plots(center_soma_y=True)


	
	






