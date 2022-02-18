from SWC import SWC
from os import listdir, system
from os.path import join, exists
from os import mkdir, remove, stat
from shutil import copyfile, move
from skimage.io import imread
from StitchingXML import StitchingXML

import tifffile as tif
import numpy as np
import params
import sys
import utils


"""
	SWC_Processing contains functions used to batch process SWCs

"""

class SWC_Processing:

	def __init__(self, region, *, brain=None):
		
		print()
		print('##############')
		print('SWC Processing')		
		print('##############')	
		print()
		print('Region:', region)

		self.ORIGINAL_SWC_PATH = '/data/palmer/data/swc_files/'
		self.SWC_CROPPED_IMAGE_PATH = '/data/palmer/data/reconstructed_brains/'
		self.MARKUP_BASE_PATH = '/data/anarasim/data/'	
		self.SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'

		self.region = region
		self.brain = brain
		self.region_path = join(params.SWC_BASE_PATH, region)

		# allow for sub regions
		# by default no sub regions
		self.sub_region_list = ['']
		if self.region == 'pfc':
			self.sub_region_list = ['anterior_cingulate', 'limbic_cortex', 'orbital_cortex']

		if brain is not None:
			print('Brain: ', brain)

		print()


	def load_swcs(self, in_path, *, in_dir='', out_dir=None, sub_region=None):

		"""
		Loads all SWCs from input path and stores in list

		"""

		print()
		print('###############')
		print('Loading SWCs...')
		print('###############')
		print()


		# dictionary to store swcs
		# keys will be subregions
		# if no subregions key will be '/'
		swc_dict = {}

		if sub_region is not None:
			self.sub_region_list = [sub_region]

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			print('Sub Region:',sub_region_name)

			# create outpath if doesnt exist
			out_path = join(self.region_path, sub_region_name, out_dir)
			if not exists(out_path):
				mkdir(out_path)

			# in path
			in_path_sub_region = join(in_path, sub_region_name,in_dir)

			# get all swc file names
			swc_filenames = [f for f in listdir(in_path_sub_region) if f.endswith('.swc')]
	
			# get only for specific brain
			if self.brain is not None:
				swc_filenames = [f for f in swc_filenames if self.brain in f]				

			# remove all swcs that have been processed
			if out_path is not None:
				already_processed = [utils.get_swc_info(f)['name'] for f in listdir(out_path) if (f.endswith('.swc') or f.endswith('.txt'))]

				if self.brain is not None:
					already_processed = [f for f in already_processed if self.brain in f]
				
				swc_filenames = [swc for swc in swc_filenames if swc[:-4] not in already_processed]
			
		
			print('# SWCs:', len(swc_filenames))

			swcs = []

			for swc_filename in sorted(swc_filenames):

				swc_path = join(in_path_sub_region, swc_filename)
				print(swc_path)
				swc = SWC(swc_path)
				swcs.append(swc)
				

			swc_dict[sub_region_name] = swcs

			print()
		
		return swc_dict


	def reorder_and_invert_y_from_vaa3D(self, gcut):

		"""
		- Reorders SWCs automatically when loaded and resaves
		- This is needed since Vaa3D has incorrect ordering
		- Doesnt preserve radius and structure id from Vaa3d
		- Also inverts Y coordinate based on image shape

		"""

		print()
		print('###############################')
		print('Reorder And Invert Y From Vaa3D')
		print('###############################')
		print()

		# set up paths
		if gcut:
			in_path = join(self.region_path, 'swc_files')
			in_dir = 'gcut_neurons'
			out_dir = 'oblique_reordered_y_inverted_gcut'
		else:
			in_path = join(self.ORIGINAL_SWC_PATH, self.region)
			in_dir = ''
			out_dir = 'oblique_reordered_y_inverted'

		# load swcs
		swc_dict = self.load_swcs(in_path, in_dir=in_dir, out_dir=out_dir)

		# path for cropped images
		cropped_image_path = join(self.SWC_CROPPED_IMAGE_PATH, self.region)

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			print('Sub Region:', sub_region_name)

			# set up output path
			out_path = join(self.region_path, sub_region_name, out_dir)

			for swc in swc_dict[sub_region_name]:

				# get cropped image path
				swc_info = swc.get_swc_info()
				brain_id = swc_info['brain']
				soma_id = swc_info['id']
				cropped_image_soma_path = join(cropped_image_path, sub_region_name, brain_id + '_reconstructed', soma_id)

				raw_cropped_files = [f for f in listdir(cropped_image_soma_path) if 'raw_cropped' in f]
				if len(raw_cropped_files) != 1:
					exit('Error: Amount of raw_cropped files is ' + str(len(raw_cropped_file)) + ' but needs to be 1')
				cropped_path = join(cropped_image_soma_path, raw_cropped_files[0])

				# get cropped image size 
				image_y_size = imread(cropped_path).shape[1]
				if image_y_size <= 0:
					exit('Error: Image y size must be greater than 0')

				# invert Y
				swc.invert_y_pixel_coords(image_y_size)

				# save
				swc.save_swc(out_path, fmt=params.SWC_FORMAT_INT, preserve_sid=False, preserve_radius=False)

			print()

	def normalize(self, gcut):
		
		"""
		Normalize SWCs to microns, centers around zero

		"""

		print()
		print('#######################')
		print('Normalize and Invert Y')
		print('#######################')
		print()


		# setup up paths
		in_path = self.region_path
		if gcut:
			in_dir = 'oblique_reordered_y_inverted_gcut'
			out_dir = 'normalized_oblique_gcut'
		else:
			in_dir = 'oblique_reordered_y_inverted'
			out_dir = 'normalized_oblique'

		# load swcs
		swc_dict = self.load_swcs(in_path, in_dir=in_dir, out_dir=out_dir)

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			print('Sub Region:', sub_region_name)
			out_path = join(self.region_path, sub_region_name, out_dir)

			for swc in swc_dict[sub_region_name]:
	
				# center around soma and convet to microns
				swc.pixels_to_microns(orientation='oblique', center_around_soma=True)

				# save
				swc.save_swc(out_path)

			print()
		
	def copy_swcs_after_initial_markup(self, gcut):

		"""
		Copys SWCs to respective markup folders after initial markup
			- 0: Not usable
			- 1: Useable
			- 2: Potentially Useable (needs branches removed)
			- 3: Need to be gcut
			- 4: Potentially useable but needs more work (if need more in future)

		"""

		print()
		print('#########################')
		print('Copy After Initial Markup')
		print('#########################')
		print()

		
		# keep track of swcs needed for gcut
		swcs_to_gcut = {}
		
		
		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			print('Sub Region:', sub_region_name)

			# set up swcs to gcut dict
			swcs_to_gcut[sub_region_name] = []

			# create output folders if they dont already exist
			if gcut:
				for i in [0,1,2,4]:
					marked_up_path = join(self.region_path, sub_region_name,'normalized_oblique_gcut_' + str(i))
					if not exists(marked_up_path):
						mkdir(marked_up_path)
			else:
				for i in range(5):
					marked_up_path = join(self.region_path, sub_region_name, 'normalized_oblique_' + str(i))
					if not exists(marked_up_path):
						mkdir(marked_up_path)


			# copy Arun's markup csv
			if gcut:
				source = join(self.MARKUP_BASE_PATH, self.region, sub_region_name, 'markups','swc_markups_gcut.csv')
				destination = join(self.region_path, sub_region_name, 'initial_markup_gcut.csv')
			else:
				source = join(self.MARKUP_BASE_PATH, self.region, sub_region_name, 'markups','swc_markups.csv')
				destination = join(self.region_path, sub_region_name, 'initial_markup.csv')

			copyfile(source, destination)

			# load markups csv
			try:
				if gcut:
					markup_path = join(self.region_path, sub_region_name, 'initial_markup_gcut.csv')
				else:
					markup_path = join(self.region_path, sub_region_name, 'initial_markup.csv')

				# check if file is empty
				if stat(markup_path).st_size == 0:
					print()
					continue
					

				markups = np.genfromtxt(markup_path, delimiter=',', dtype=object)
				if markups.ndim == 1:
					markups = markups.reshape(1,-1)				


			except Exception as e:
				exit(e)
				
			# fix data tyles
			markups[:,0] = markups[:,0].astype(str)
			markups[:,1] = markups[:,1].astype(int)

			# copy swcs
			for markup in markups:

				swc_id, label = markup
				prev_label = None

				if gcut:

					# check if swc already exists in other folders
					if exists(join(self.region_path, sub_region_name, 'normalized_oblique_gcut_0', swc_id + '.swc')):
						prev_label = 0
					if exists(join(self.region_path, sub_region_name, 'normalized_oblique_gcut_1', swc_id + '.swc')):
						prev_label = 1
					if exists(join(self.region_path, sub_region_name, 'normalized_oblique_gcut_2', swc_id + '.swc')):
						prev_label = 2
					if exists(join(self.region_path, sub_region_name, 'normalized_oblique_gcut_4', swc_id + '.swc')):
						prev_label = 4

					# same label
					if prev_label == label:
						continue

					print(swc_id)

					if prev_label is None:
						print('Copying to markup folder ' + str(label))
						source = join(self.region_path, sub_region_name, 'normalized_oblique_gcut', swc_id + '.swc')
						destination = join(self.region_path, sub_region_name ,'normalized_oblique_gcut_' + str(label), swc_id + '.swc')
						copyfile(source, destination)
					else:
						print('Moving from '+ str(prev_label) + ' to ' + str(label))
						source = join(self.region_path, sub_region_name,'normalized_oblique_gcut_' + str(prev_label), swc_id + '.swc')
						destination = join(self.region_path, sub_region_name,'normalized_oblique_gcut_' + str(label), swc_id+ '.swc')
						move(source, destination)

				else:

					# check if swc already exists in other folders
					if exists(join(self.region_path,sub_region_name,'normalized_oblique_0', swc_id + '.swc')):
						prev_label = 0
					if exists(join(self.region_path,sub_region_name,'normalized_oblique_1', swc_id + '.swc')):
						prev_label = 1
					if exists(join(self.region_path,sub_region_name,'normalized_oblique_2', swc_id + '.swc')):
						prev_label = 2
					if exists(join(self.region_path,sub_region_name,'normalized_oblique_3', swc_id + '.swc')):
						prev_label = 3
					if exists(join(self.region_path,sub_region_name,'normalized_oblique_4', swc_id + '.swc')):
						prev_label = 4
			
					# same label
					if prev_label == label:
						continue

					# keep track of swcs to gcut
					if label == 3:
						swcs_to_gcut[sub_region_name].append(swc_id)
					
					print(swc_id)

					if prev_label is None:
						print('Copying to markup folder ' + str(label))
						source = join(self.region_path, sub_region_name,  'normalized_oblique', swc_id + '.swc')
						destination = join(self.region_path, sub_region_name,'normalized_oblique_' + str(label), swc_id + '.swc')
						copyfile(source, destination)
					else:
						print('Moving from ' + str(prev_label)  + ' to ' + str(label))
						source = join(self.region_path, sub_region_name, 'normalized_oblique_' + str(prev_label), swc_id + '.swc')
						destination = join(self.region_path, sub_region_name, 'normalized_oblique_' + str(label), swc_id + '.swc')
						move(source, destination) 

			print()

		# print out which swcs need to be gcut
		if not gcut:
			print()
			print('############')
			print('SWCs to gcut')
			print('############')
			print()

			for sub_region_name in self.sub_region_list:
			
				print('Sub Region:', sub_region_name)
				for swc_id in swcs_to_gcut[sub_region_name]:
					print(swc_id)

				print()


	def copy_pruned_swcs(self):
		
		"""
		copy swcs after pruning

		"""

		print()
		print('#####################')
		print('Copy Pruned SWCs')
		print('#####################')
		print()

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			print(sub_region_name)

			pruned_path = join(params.ARUN_BASE_PATH, self.region, sub_region_name, 'manual_pruning')
			swc_names = [f for f in listdir(pruned_path) if f.endswith('.swc')]

			norm_oblique_files = listdir(join(self.region_path, sub_region_name, 'normalized_oblique_2'))
			norm_oblique_gcut_files = listdir(join(self.region_path, sub_region_name , 'normalized_oblique_gcut_2'))

			out_path = join(self.region_path, sub_region_name ,'normalized_oblique_2_pruned')
			if not exists(out_path):
				mkdir(out_path)

			out_path_gcut = join(self.region_path, sub_region_name ,'normalized_oblique_gcut_2_pruned')
			if not exists(out_path_gcut):
				mkdir(out_path_gcut)

			for swc_name in swc_names:
			
				if swc_name in norm_oblique_files:

					source = join(pruned_path, swc_name)
					destination = join(out_path, swc_name)
					copyfile(source, destination)

				elif swc_name in norm_oblique_gcut_files:

					source = join(pruned_path, swc_name)
					destination = join(out_path_gcut, swc_name)
					copyfile(source, destination)

				else:
					print('Error: ' + swc_name + ' not in normalized_oblique folders')


	def merge_processed_swcs(self):

		"""
		merge processed arrays after initial markup

		"""

		print()
		print('#####################')
		print('Merge Processed SWCs')
		print('#####################')
		print()

		dirs = ['normalized_oblique_1','normalized_oblique_2_pruned','normalized_oblique_gcut_1','normalized_oblique_gcut_2_pruned' ]	
	
		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			out_path = join(self.region_path, sub_region_name , 'normalized_oblique_merged')
			if not exists(out_path):
				mkdir(out_path)

			for directory in dirs:
			
				swc_names = listdir(join(self.region_path, sub_region_name, directory))

				for swc_name in swc_names:
		
					source = join(self.region_path, sub_region_name, directory, swc_name)
					destination = join(out_path, swc_name)
					copyfile(source, destination)


	def collapse_soma(self):

		"""
		collapse soma

		"""

		print()
		print('#############')
		print('Collapse Soma')
		print('#############')
		print()

		in_path = join(self.region_path)
		in_dir = 'normalized_oblique_merged'
		out_dir = 'normalized_oblique_collapsed_soma'
		swcs_dict = self.load_swcs(in_path, in_dir=in_dir, out_dir=out_dir)


		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			swcs = swcs_dict[sub_region_name]

			out_path = join(self.region_path, sub_region_name, out_dir)

			for swc in swcs:

				swc.collapse_soma()
				swc.save_swc(out_path)

	def create_skeletons_for_dendrite_markup(self):

		"""
		Creates skeletons and tiffs for dendrite markup
	
		"""

		print()
		print('#####################################')
		print('Create Skeletons for Dendrite Markup')
		print('#####################################')
		print()

		initial_swc_base_path = join(self.ORIGINAL_SWC_PATH, self.region)
		cropped_image_path = join(self.SWC_CROPPED_IMAGE_PATH , self.region)

		in_path = self.region_path
		in_dir = 'normalized_oblique_collapsed_soma'
		out_dir = 'dendrite_markup_oblique'

		# load SWCs
		swcs_dict = self.load_swcs(in_path, in_dir=in_dir, out_dir=out_dir)

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			swcs = swcs_dict[sub_region_name]

			out_path = join(self.region_path, sub_region_name, out_dir)

			for swc in swcs:

				print(swc.get_swc_name())
			
				# get soma from jason swc to unnormalize
				# check for swc in gcut folder
				initial_swc_path = join(initial_swc_base_path, sub_region_name, 'gcut_neurons', swc.get_swc_name(with_ext=True))
	
				if not exists(initial_swc_path):
					initial_swc_path = join(initial_swc_base_path, sub_region_name, swc.get_swc_name(with_ext=True))
	
				original_swc = SWC(initial_swc_path)
				original_soma = original_swc.node_soma.get_coords()

				# get cropped image path
				brain_id = swc.get_swc_info()['brain']
				soma_id = swc.get_swc_info()['id']
				cropped_image_soma_path = join(cropped_image_path, sub_region_name, brain_id + '_reconstructed', soma_id)

				raw_cropped_files = [f for f in listdir(cropped_image_soma_path) if 'raw_cropped' in f]
				if len(raw_cropped_files) != 1:
					exit('Error: Amount of raw_cropped files is ' + str(len(raw_cropped_files)) + ' but needs to be 1')
				cropped_path = join(cropped_image_soma_path, raw_cropped_files[0])

				# get cropped image size 
				cropped_shape = imread(cropped_path).shape

				# convert to cropped coordinates
				swc.normalized_to_cropped('oblique', original_soma, cropped_shape)
		
				# save stem info	
				swc.save_branch_info_for_structure_label_markup(out_path)

				# save skeleton			
				swc.save_skeleton(join(out_path, swc.get_swc_name() + '_dilated_skeleton.tif'), image_shape=cropped_shape, compress_ids=True, dilate=True, dilate_width=3)

	def oblique_to_coronal_for_dendrite_markup(self):

		"""
		Converts skeletons and branches from oblique to coronal
	
		"""

		print()
		print('###########################################')
		print('Create Coronal Skeletons for Dendrite Markup')
		print('###########################################')
		print()

		
		cropped_image_path = join(self.SWC_CROPPED_IMAGE_PATH , self.region)


		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			in_path = join(self.region_path, sub_region_name, 'dendrite_markup_oblique')
			out_path = join(self.region_path, sub_region_name, 'dendrite_markup_coronal')

			if not exists(out_path):
				mkdir(out_path)

			already_processed = [f[:-4] for f in listdir(out_path) if f.endswith('.txt')]

			swc_names = [f[:-4] for f in listdir(in_path) if f.endswith('.txt') and f[:-4] not in already_processed]
		
			print('# SWCs:', len(swc_names))

			for swc_name in sorted(swc_names):

				print(swc_name)

				# convert skeleton and raw to coronal

				# skeleton path
				skeleton_path = join(in_path, swc_name + '_dilated_skeleton.tif')

				# get cropped image path
				brain_id = swc_name.split('_')[0]
				soma_id = swc_name.split('_')[-1]
				cropped_image_soma_path = join(cropped_image_path, sub_region_name, brain_id + '_reconstructed',soma_id)

				raw_cropped_files = [f for f in listdir(cropped_image_soma_path) if 'raw_cropped' in f]
				if len(raw_cropped_files) != 1:
					exit('Error: Amount of raw_cropped files is ' + str(len(raw_cropped_files)) + ' but needs to be 1')
				cropped_path = join(cropped_image_soma_path, raw_cropped_files[0])

				args = skeleton_path + '?' + cropped_path + '?' + out_path + '?' + swc_name
				exec_string = params.IMAGEJ_PATH +  ' --console -macro ' + params.OBLIQUE_TO_CORONAL_MACRO_PATH + ' ' + args

				# run imagej from command line
				system(exec_string)

				# convert branch file to coronal
				branch_file_path = join(in_path, swc_name + '.txt')
				coords, ids = utils.parseBranchMarkupFile(branch_file_path)

				print('Load tif...')
				oblique_tif_path = join(in_path, swc_name + '_dilated_skeleton.tif')
				oblique_tif = imread(oblique_tif_path)
				x_len, y_len ,z_len = oblique_tif.shape[::-1]
	
				print('Transform Coordinate from Oblique to Coronal Orientation...')
				coronal_coords = []
				fused_dims = {'x':{'length':x_len,'min':0},'y':{'length':y_len,'min':0},'z':{'length':z_len,'min':0}}
				for coord in coords:
					coronal_coords.append(utils.obliqueToCoronal(coord, fused_dims, params.SHEAR_FACTOR_ANISOTROPIC))

				print('Saving Coronal Branch Coordinates...')
				utils.writeBranchMarkupFile(join(out_path, swc_name + '.txt'), coronal_coords, ids)

		
	def set_structure_labels(self):

		"""
		Sets structure labels from markup

		"""

		print()
		print('#####################')
		print('Set Structure Labels')
		print('#####################')
		print()

		in_path = self.region_path
		in_dir = 'normalized_oblique_collapsed_soma'
		out_dir = 'normalized_oblique_collapsed_soma_with_structure_labels'
		swcs_dict = self.load_swcs(in_path, in_dir=in_dir, out_dir=out_dir)

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			swcs = swcs_dict[sub_region_name]

			out_path = join(self.region_path, sub_region_name, out_dir)

			if not exists(out_path):
				mkdir(out_path)

			markup_path = join(params.ARUN_BASE_PATH, self.region, sub_region_name, 'dendrite_markup_coronal')

			files_not_exist = []

			for swc in swcs:
			
				swc_name = swc.get_swc_name()

				print(swc_name)

				# load markup info
				markup_path_local = join(markup_path, swc_name + '.txt')
				label_info = utils.read_branch_markups(markup_path_local)

				if label_info:
		
					# set structure labels
					swc.set_structure_labels(label_info)

					# save SWC
					swc.save_swc(out_path)

				else:
					print('Structure label file doesnt exist...')	
					files_not_exist.append(swc_name)
		
			print('Files that dont exist:', files_not_exist)		


	def smooth(self, sub_region=None):

		"""
		smooth
	
		"""

		print()
		print('#######')
		print('Smooth')
		print('#######')
		print()

		in_path = self.region_path
		in_dir = 'normalized_oblique_collapsed_soma_with_structure_labels'
		out_dir = 'normalized_oblique_collapsed_soma_with_structure_labels_smoothed'
		swcs_dict = self.load_swcs(in_path, in_dir=in_dir, out_dir=out_dir, sub_region=sub_region)

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			swcs = swcs_dict[sub_region_name]

			out_path = join(self.region_path, sub_region_name, out_dir)


			if not exists(out_path):
				mkdir(out_path)
	
			for swc in swcs:

				swc.smooth(orientation='oblique')
				swc.save_swc(out_path)

	def oblique_to_coronal(self, sub_region=None):

		"""
		oblique to coronal
		
		"""

		print()
		print('##################')
		print('Oblique to Coronal')
		print('###################')
		print()


		in_path = self.region_path
		in_dir = 'normalized_oblique_collapsed_soma_with_structure_labels_smoothed'
		out_dir = 'coronal_full_size'
		swcs_dict = self.load_swcs(in_path, in_dir=in_dir, out_dir=out_dir, sub_region=sub_region)

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			swcs = swcs_dict[sub_region_name]

			out_path = join(self.region_path, sub_region_name, out_dir)
			outpath_normalized = join(self.region_path, sub_region_name, 'coronal_normalized')

			# create output paths
			if not exists(out_path):
				mkdir(out_path)
			if not exists(outpath_normalized):
				mkdir(outpath_normalized)
	
			# iterate through swcs
			for swc in swcs:
	
				print(swc.get_swc_name())
				swc_id = int(swc.get_swc_info()['id'])
	
				# get soma coordinate from soma detection
				brain_path = join(self.SOMA_BASE_PATH, swc.get_swc_info()['brain'])
				soma_path = join(brain_path, 'regions', sub_region_name + '.csv')
				somas = np.genfromtxt(soma_path, delimiter=',', dtype=str)
				volume = somas[swc_id-1,4]
				soma_volume_coord = [int(x) for x in somas[swc_id-1,5:8]]
	
				# get stitching matrices and full size image lengths
				xml_path = join(brain_path,'stitching_parameters.xml')
				stitching_parameters = StitchingXML(xml_path, verbose=False)

				# get matrices and dimensions
				stitching_matrices = stitching_parameters.getStitchingMatrices(volume)
				fused_dims = stitching_parameters.getFusedDimensions()

				# convert from microns to pixels
				swc.microns_to_pixels(orientation='oblique', soma_center=list(soma_volume_coord))

				# apply stitching matrices
				swc.apply_stitching_matrices(stitching_matrices)

				# convert from oblique to coronal
				swc.oblique_to_coronal(fused_dims, params.SHEAR_FACTOR_ANISOTROPIC)

				# save full size data		
				swc.save_swc(out_path)

				# normalize
				swc.pixels_to_microns(orientation= 'coronal', center_around_soma=True)

				# save normalized data		
				swc.save_swc(outpath_normalized)


	def register(self, sub_region=None):

		"""
		register

		"""
		
		print()
		print('########')
		print('Register')
		print('########')
		print()


		in_path = self.region_path
		in_dir = 'coronal_full_size'
		out_dir = 'registered'
		swcs_dict = self.load_swcs(in_path, in_dir=in_dir, out_dir=out_dir, sub_region=sub_region)

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:

			swcs = swcs_dict[sub_region_name]

			out_path = join(self.region_path, sub_region_name, out_dir)
			outpath_normalized = join(self.region_path, sub_region_name, 'registered_normalized')

			# create output directories
			if not exists(out_path):
				mkdir(out_path)
			if not exists(outpath_normalized):
				mkdir(outpath_normalized)

			for swc in swcs:

				print(swc.get_swc_name())

				# register SWC
				swc.register()

				# save swc
				swc.save_swc(out_path)

				swc.center_around_soma()
				swc.save_swc(outpath_normalized)

	def remove_zero_swcs(self):
		
		"""
		SWCs marked up as 0 in later steps need to be removed from further processing

		"""

		print()
		print('##########################')
		print('Remove SWCs marked up as 0')
		print('##########################')
		print()

		# iterate through sub regions
		for sub_region_name in self.sub_region_list:
	
			sub_region_path = join(self.region_path, sub_region_name)
			
			# get all swcs marked up as 0
			normalized_oblique_0_path = join(sub_region_path,'normalized_oblique_0')
			normalized_oblique_gcut_0_path = join(sub_region_path,'normalized_oblique_gcut_0')
			zero_swcs = [f[:-4] for f in listdir(normalized_oblique_0_path)] + [f[:-4] for f in listdir(normalized_oblique_gcut_0_path)]
		
			for swc_name in zero_swcs:

				# remove from normalized_oblique_2_pruned
				path = join(sub_region_path, 'normalized_oblique_2_pruned', swc_name + '.swc')
				if exists(path):
					print('Removing:', path)
					remove(path)

				# remove from normalized_oblique_gcut_2_pruned
				path = join(sub_region_path, 'normalized_oblique_gcut_2_pruned', swc_name + '.swc')
				if exists(path):
					print('Removing:', path)
					remove(path)

				# remove from merged
				path = join(sub_region_path,'normalized_oblique_merged', swc_name + '.swc')
				if exists(path):
					print('Removing:', path)
					remove(path)

				# remove from collapsed soma
				path = join(sub_region_path, 'normalized_oblique_collapsed_soma', swc_name + '.swc')
				if exists(path):
					print('Removing:', path)
					remove(path)

				# remove from dendrite_markup_oblique
				path = join(sub_region_path,'dendrite_markup_oblique', swc_name + '.tif')
				if exists(path):
					print('Removing:', path)
					remove(path)
				path = join(sub_region_path,'dendrite_markup_oblique', swc_name + '.txt')
				if exists(path):
					print('Removing:', path)
					remove(path)
				path = join(sub_region_path,'dendrite_markup_oblique', swc_name + '_dilated_skeleton.tif')
				if exists(path):
					print('Removing:', path)
					remove(path)

				# remove from dendrite_markup_coronal
				path = join(sub_region_path,'dendrite_markup_coronal', swc_name + '.txt')
				if exists(path):
					print('Removing:', path)
					remove(path)
				path = join(sub_region_path,'dendrite_markup_coronal', swc_name + '_raw.tif')
				if exists(path):
					print('Removing:', path)
					remove(path)
				path = join(sub_region_path,'dendrite_markup_coronal', swc_name + '_skeleton.tif')
				if exists(path):
					print('Removing:', path)
					remove(path)
				path = join(sub_region_path,'dendrite_markup_coronal', swc_name + '_skeleton_dilated.tif')
				if exists(path):
					print('Removing:', path)
					remove(path)
	
	def createSkeleton(self, inpath, outpath, inres, outres, imageshape, max_proj):

		swcs = self.loadSWCs(inpath)

		if max_proj:
			outimage = np.zeros(shape=imageshape[1:], dtype=np.uint8)
		else:
			outimage =  np.zeros(shape=imageshape, dtype=np.uint8)
		downsampling = inres/outres
	
		for id_, swc in enumerate(swcs, start=1):
			
			print(swc)
	
			swcArray = (swc.generateSWCArray(onlyCoords=True)*downsampling).round().astype(int)

			if max_proj:
				outimage[swcArray[:,1],swcArray[:,0]] = id_
			else:
				outimage[swcArray[:,2],swcArray[:,1],swcArray[:,0]] = id_



		# save
		tif.imsave(outpath, outimage, imagej=True)

	def add_layer_info_to_names(self, sub_region=None):


		for sub_region_name in self.sub_region_list:

			if sub_region_name != sub_region:
				continue

			in_path = join(self.region_path, sub_region_name, 'registered')
			swc_names = listdir(in_path)
			layer_counts = [0]*6

			for swc_name in swc_names:
			
				swc_info = utils.get_swc_info(swc_name)
				layer = swc_info['layer']

				if layer is None:

					print(swc_name)
				
					# load somas
					somas_path = join(self.SOMA_BASE_PATH, swc_info['brain'] , 'regions', sub_region_name + '.csv')
					somas = np.genfromtxt(somas_path, delimiter=',')
				
					# get layer info
					layer = int(somas[int(swc_info['id'])-1,0])

					# correction for layer 1			
					if (self.region == 'barrel_cortex' or self.region == 'motor_cortex') and layer == 1:
						layer = 2
					if self.region == 'pfc' and (layer == 1 or layer == 2 or layer == 3):
						layer = 2
				
					# add layer to file name
					source_path = join(in_path, swc_name)
					dest_path = join(in_path, str(layer) + '_' + swc_name)
					move(source_path, dest_path)

				layer_counts[layer-1] += 1

			if self.region == 'motor_cortex':
				mop_path = '/data/elowsky/OLST/swc_analysis/automatically_traced/motor_cortex/first_round_of_reconstructions/final_swcs_1_19_2021/'
				layerCounts[1] += len(listdir(join(mop_path, 'layer_'+str(2))))
				layerCounts[4] += len(listdir(join(mop_path, 'layer_'+str(5))))
				layerCounts[5] += len(listdir(join(mop_path, 'layer_'+str(6))))


			print('Layer Counts: ', layer_counts)

	def getLayerInfo(self, directory):

		inpath = join(self.regionPath, directory)
		swcNames = listdir(inpath)
		layerCounts = [0]*6

		for swcName in swcNames:
	
			swcInfo = utils.getSWCInfo(swcName)
			layer = swcInfo['layer']

			if layer is None:

				print(swcName)
				
				# load somas
				somasPath = join(params.SOMA_BASE_PATH, swcInfo['brain'],'regions', self.region + '.csv')
				somas = np.genfromtxt(somasPath, delimiter=',')
				
				# get layer info
				layer = int(somas[int(swcInfo['id'])-1,0])

				# correction for layer 1			
				if (self.region == 'barrel_cortex' or self.region == 'motor_cortex') and layer == 1:
					layer = 2

				if self.region == 'pfc' and (layer == 1 or layer == 2 or layer == 3):
					layer = 2
				

			layerCounts[layer-1] += 1

		print('Layer Counts: ', layerCounts)
				
		
if __name__ == '__main__':


	REGION = 'retrosplenial_cortex'
	BRAIN = None
	swcs = SWC_Processing(REGION, brain=BRAIN)
	
	# reorder and invert y
	#swcs.reorder_and_invert_y_from_vaa3D(gcut=False)

	# normalize 
	#swcs.normalize(gcut=False)

	# Tell Arun to do initial markup

	# copy to markup folders after initial markup
	#swcs.copy_swcs_after_initial_markup(gcut=False)

	# 0 - no good
	# 1 - good without fixing
	# 2 - maybe good, needs fixing possibly
	# 3 - gcut
	# 4 - potentially good, but not going deal with right now

	# send jason gcut list 

	# reorder and invert y for gcut swcs
	#swcs.reorder_and_invert_y_from_vaa3D(gcut=True)

	# normalize for gcut swcs
	#swcs.normalize(gcut=True)

	# Tell Arun to do initial markup for gcut

	# copy to markup folders after initial markup for gcut swcs
	#swcs.copy_swcs_after_initial_markup(gcut=True)

	# Tell Arun to do pruning

	# copy pruned swcs
	#swcs.copy_pruned_swcs()

	# merge acceptable swcs
	#swcs.merge_processed_swcs()

	# collapse soma
	#swcs.collapse_soma()

	# create oblique skeletons and txt files for dendrite markup
	#swcs.create_skeletons_for_dendrite_markup()

	# convert coronal skeletons and raw image for dendrite markup
	#swcs.oblique_to_coronal_for_dendrite_markup()

	# go through all that are marked up as 0 and remove files from further processing steps
	#swcs.remove_zero_swcs()

	# set strucure labels
	#swcs.set_structure_labels()

	# smooth
	#swcs.smooth()

	# oblique to coronal
	#swcs.oblique_to_coronal()

	# register
	#swcs.register()

	# change swc names to include layer information
	#swcs.add_layer_info_to_names()


	######## END  ##########

	# not part of pipeline
	#swcs.getLayerInfo('normalized_oblique_merged')


	# skeleton
	#region = 'barrel_cortex'
	#swcs = SWC_Processing(region)
	#out_res = 10
	#inpath = '/data/elowsky/OLST/swc_analysis/automatically_traced/'+region+'/registered/'
	#outpath = '/data/elowsky/OLST/swc_analysis/automatically_traced/'+region+'/barrel_cortex_' +str(out_res) + 'um.tif'
	#inres = np.array([1, 1, 1])
	#outres = np.array([out_res, out_res, out_res])
	#outshape_25um = [528, 320, 456]

	#outshape = [int(x*(25/out_res)) for x in outshape_25um]

	#swcs.createSkeleton(inpath, outpath, inres, outres, outshape, max_proj=False)
	
	



