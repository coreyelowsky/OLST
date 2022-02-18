import numpy as np
import math
import tifffile as tif
import copy
import matplotlib.pyplot as plt
import params
import utils
import SWC_Array

from scipy.ndimage import gaussian_filter
from Branch import Branch
from Node import Node
from os.path import exists, join
from os import listdir, stat
from sys import exit
from mpl_toolkits.mplot3d import Axes3D
from skimage.io import imread

"""
SWC object
	- swc is internally represented as a tree of nodes
	- implements functions for processing of swcs
"""

class SWC:

	def __init__(self, swc_input, build_tree=True):

		"""
		Params:	swcInput - can either be a path or a 2d Numpy Array
		
		This constructor calls build_swc_tree() to create the tree structure
		that will be used to represent the SWC

		"""		

		if swc_input is not None:

			# make sure input is either a string (path) or a numpy array
			if isinstance(swc_input, str):
				self.swc_path = swc_input
			elif not isinstance(swc_input, np.ndarray):
				exit('Error: swcInput is of type ' + str(type(swc_input)) + ' but must be either a string or numpy array')


			# if input is string then load SWC from file
			if isinstance(swc_input, str): 
				# if pathe exists
				if exists(self.swc_path):
					# if file is empty
					if stat(self.swc_path).st_size == 0:
						swc_array = np.array([[1,1,0,0,0,0,-1]])
					else:
						swc_array = np.genfromtxt(self.swc_path)
				else:
					exit('ERROR: SWC Path Does Not Exist!!!')
			else:
				swc_array = swc_input
		
			# if 1d array make 3d
			if swc_array.ndim == 1:
				swc_array = swc_array.reshape(1,-1)

			# dimension errors
			if swc_array.ndim != 2:
				exit('Error: # Dimensions of SWC Array is ' + str(swc_array.ndim) + ' and must be 2!')

			if swc_array.shape[1] != len(params.SWC_INDICES):
				exit('Error: # Columns of SWC Array is ' + str(swc_array.shape[1]) + ' and must be '+str(len(params.SWC_INDICES)))

			# modify data types
			swc_array = SWC_Array.set_data_types(swc_array)

			self.swc_array = swc_array

			if build_tree:
				self.build_swc_tree()

	@classmethod
	def swc_from_node(cls, node):

		"""
		Alternative constructor that takes a node to be the root

		"""
		
		swc = cls(None)
		node.parent_node = None
		swc.node_soma = node		

		return swc

	@classmethod
	def swc_from_branch(cls, branch):

		"""
		Alternative constructor that takes a node to be the branch

		"""
		
		swc = cls(None)
		node.parent_node = None
		swc.node_soma = node		

		return swc

	def get_swc_info(self):

		"""
		calls get_swc_info() from utils

		"""
	
		return utils.get_swc_info(self.swc_path)


	def get_soma_coords(self):

		"""
		Returns soma coordinates

		"""
		
		return self.node_soma.get_coords()


	def get_swc_name(self, with_ext=False):
	
		"""
		returns name of swc with or without extension

		"""

		if with_ext:
			return self.swc_path.split('/')[-1]
		else:
			return self.swc_path.split('/')[-1].split('.')[0]
		

	def build_swc_tree(self):

		"""
		builds a tree structure where each coordinate in the SWC is
		represented by a Node in the tree

		Params: swc_input

		"""
		swc_info = self.get_swc_info()
		print('Building SWC Tree:', swc_info['name'])

		
		# get soma row
		soma_row = SWC_Array.get_soma_row_from_swc_array(self.swc_array)

		# create soma node
		self.node_soma = Node.node_from_array(soma_row)

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()

			# get all children
			child_rows = SWC_Array.get_child_rows_from_swc_array(self.swc_array, node.id_)

			# build nodes and push children on stack
			for child_row in child_rows:
				node_child = Node.node_from_array(child_row, parent_node=node)
				node.child_nodes.insert(0, node_child)
				stack.append(node_child)

	def contains_node(self, other_node, node_root=None):

		"""
		Check if SWC contains a node

		"""

		if node_root is None:
			node_root = self.node_soma

		stack = [node_root]

		while len(stack):

			node = stack.pop()

			if node == other_node:
				return True

			for child_node in node.child_nodes:
				stack.append(child_node)
				

		return False


	def get_nodes_in_tree(self, node_root=None):

		"""
		Return: list of all nodes in tree from root

		"""
			
		if node_root is None:
			node_root = self.node_soma

		nodes = []
		stack = [node_root]

		while len(stack):

			node = stack.pop()

			nodes.append(node)

			for child_node in node.child_nodes:
				stack.append(child_node)
				
		return nodes

	def translate_coords(self, shifts, node_root=None):

		"""
		Params: shifts, list of 3 floats [x shift, y shift, z shift]

		Translates (shifts) all coordinates

		"""

		if node_root is None:
			node_root = self.node_soma

		stack = [node_root]

		while len(stack):

			node = stack.pop()

			node.translate(shifts)

			for child_node in node.child_nodes:
				stack.append(child_node)

	def scale_coords(self, scale_factors, node_root=None, from_array=False):

		"""
		Params: scaleFactors, list of 3 floats [x scale, y scale, z scale] or just scalar 

		Scales all coordinates

		"""

		if node_root == None:
			node_root = self.node_soma

		stack = [node_root]

		while len(stack):

			node = stack.pop()
			
			node.scale(scale_factors)

			for child_node in node.child_nodes:
				stack.append(child_node)


	def invert_y_pixel_coords(self, image_y_size, node_root=None):

		"""
		Inverts all y coordinates in pixel coordinates from shape of image
		
		"""
		self.invert_y()
		self.translate_coords([0,image_y_size-1,0], node_root=node_root)		

	def invert_y(self, node_root=None):

		"""
		Inverts all y coordinates
		Assumes soma is centered already
		
		"""

		self.scale_coords(scale_factors=[1,-1,1], node_root=node_root)

	def pixels_to_microns(self, orientation, *, center_around_soma=False,  node_root=None):

		"""
		Params: orientation - string, orientation of coordinates (coronal or oblique)

		Converts all coordinates from pixels to microns
	
		"""
		
		if orientation == 'oblique':
			res_factors = params.RES_OBLIQUE
		elif orientation == 'coronal':
			res_factors = params.RES_CORONAL
		else:
			exit('Orientation Variable Not Valid: ' + orientation)

		if center_around_soma:
			self.center_around_soma()

		self.scale_coords(scale_factors=res_factors, node_root=node_root)

	def microns_to_pixels(self, orientation, soma_center=None, node_root=None):

		"""
		Params: orientation - string, orientation of coordinates (coronal or oblique)

		Converts all coordinates from microns to pixels
	
		"""
		
		if orientation == 'oblique':
			res_factors = params.RES_OBLIQUE
		elif orientation == 'coronal':
			res_factors = params.RES_CORONAL
		else:
			exit('Orientation Variable Not Valid: ' + orientation)

		self.scale_coords(scale_factors=list(1/np.array(res_factors)), node_root=node_root)

		if soma_center != None:
			self.translate_coords(soma_center, node_root=node_root)


	def normalized_to_cropped(self, orientation, soma_center, cropped_shape):

		"""
		Converts normalized coords to cropped image coords

		"""

		self.invert_y()
		self.microns_to_pixels(orientation='oblique', soma_center=soma_center)

		# accounts for flipped y
		self.invert_y()
		self.translate_coords([0,cropped_shape[1]-1,0])


	def apply_stitching_matrices(self, stitching_matrices):

		"""
		Applys Stitching Matrices 		

		"""

		self.transform_coords(stitching_matrices['calibration'])
		self.transform_coords(stitching_matrices['Translation to Regular Grid'])
		self.transform_coords(stitching_matrices['Stitching Transform'])
		self.transform_coords(np.linalg.inv(stitching_matrices['calibration']))


	def transform_coords(self, transformation_matrix):

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()
			node_coords = transformation_matrix @ np.array([node.x,node.y,node.z,1]).reshape(-1,1)
			node.x = node_coords.flatten()[0]
			node.y = node_coords.flatten()[1]
			node.z = node_coords.flatten()[2]

			for child_node in node.child_nodes:
				stack.append(child_node)


	def extract_type_tree(self, structures):

		"""
		Params: structure - string, name of structure (basal dendrite, apical dendrite, etc...)
		
		Removes all branches that are not of the same type as structure. For example, if only 
		want to do analysis on basal tree, make struture = 'basal dendrite'

		"""
	
		if not type(structures) is list:
			exit('Error: structures must be a list!')		
	
		for structure in structures:
			if not structure in params.SID_DICT:
				exit('Structure Not In Structure Dictionary: ' + str(structure) +'\nMust Be In: ' + str(list(params.SID_DICT.keys())))
	

		struture_ids = [params.SID_DICT[structure] for structure in structures]

		new_child_list = [node for node in self.node_soma.child_nodes if node.structure_id in struture_ids]
		self.node_soma.child_nodes = new_child_list


	def oblique_to_coronal(self, fused_dimensions, shear_factor):

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()

			coronal_coords = utils.oblique_to_coronal([node.x,node.y,node.z], fused_dimensions, shear_factor)
			node.x, node.y, node.z = coronal_coords

			for child_node in node.child_nodes:
				stack.append(child_node)


	def center_around_soma(self):

		"""
		Translates all coordinates so soma is at (0,0,0)

		"""

		self.translate_coords(shifts=[-self.node_soma.x, -self.node_soma.y ,-self.node_soma.z])

	def num_nodes(self, node_root=None):
		
		"""
		Return: int, number of nodes in SWC 		
	
		"""

		if node_root is None:
			node_root = self.node_soma

		return len(self.get_nodes_in_tree(node_root=node_root))


	def num_stems(self):

		"""
		Return: int, number of stems (number of children of soma)		
	
		"""
		
		return self.node_soma.num_children()
		

	def set_radii(self, soma_radius, radius):

		"""
		Params: somaRadius - radius of soma node 
			radius - radius of all other node besides soma

		Sets the radius of all nodes in the tree

		"""

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()

			if node is self.node_soma:
				node.radius = soma_radius
			else:
				node.radius = radius

			for child_node in node.child_nodes:
				stack.append(child_node)

	def set_branch_ids_same(self):

		"""
		Sets all ids to ids of branch (id of child of soma

		"""
		
		for child in self.node_soma.child_nodes:

			branch_id = child.id_
			
			stack = [child]

			while len(stack):

				node = stack.pop()

				node.id_ = branch_id

				for child_node in node.child_nodes:
					stack.append(child_node)

	def save_skeleton(self, out_path, *, image_shape=None, compress_ids=False, downsampling=1, pixel_connected=False, gaussian_blur=False, dilate=False, dilate_width=3):

		""" 
		Params:
		

		Saves tiff file that is skeleton
		"""
		
		dilate_width_half = int(dilate_width/2)

		# if image shape is not provided will use bounds of swc
		mins = np.array([0,0,0])


		if image_shape is None:
			mins = np.min(self.generate_swc_array(only_coords=True)/downsampling,axis=0).astype(int)
			maxs = np.max(self.generate_swc_array(only_coords=True)/downsampling,axis=0).astype(int)
			image_shape = (maxs-mins)[::-1]+1
	
		if compress_ids:
			skeleton = np.zeros(shape=image_shape, dtype=np.uint8)
			print('Saving Tiff as 8 bit (Branch Ids Same)...')
		else:
			# calculate how many nodes
			num_nodes = self.num_nodes()

			if num_nodes <= 2**8 - 1:
				print('Saving Tiff as 8 bit...')
				skeleton = np.zeros(shape=image_shape, dtype=np.uint8)
			elif num_nodes <= 2**16 - 1:
				print('Saving Tiff as 16 bit...')
				skeleton = np.zeros(shape=image_shape, dtype=np.uint16)
			elif num_nodes <= 2**32 - 1:
				print('Saving Tiff as 32 bit...')
				skeleton = np.zeros(shape=image_shape, dtype=np.uint32)
			elif num_nodes <= 2**64 - 1:
				print('Saving Tiff as 64 bit...')
				skeleton = np.zeros(shape=image_shape, dtype=np.uint64)


		# set soma intensity to 1
		soma = (np.array(self.node_soma.get_coords())/downsampling).astype(int) - mins
		skeleton[soma[2],soma[1],soma[0]] = 1
		 

		id_ = 1
		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()

			coords = (np.array(node.get_coords())/downsampling).astype(int) - mins			
			
			if compress_ids:
				if node.parent_node is self.node_soma:
					id_ += 1

				if dilate:
					
					skeleton[coords[2]-dilate_width_half:coords[2]+dilate_width_half+1,coords[1]-dilate_width_half:coords[1]+dilate_width_half+1,coords[0]-dilate_width_half:coords[0]+dilate_width_half+1] = id_					
				else:
					skeleton[coords[2],coords[1],coords[0]] = id_

			else:
				skeleton[coords[2],coords[1],coords[0]] = node.id_


			for child_node in node.child_nodes:
				stack.append(child_node)	

				# if skeleton needs to be pixel connected
				if pixel_connected:
					child_coords = (np.array(child_node.get_coords())/downsampling).astype(int) - mins
					connecting_coords = utils.Bresenham3D(coords, child_coords)
					for coords in connecting_coords:
						if compress_ids:

							if node.parent_node is self.node_soma:
								id_ += 1

							skeleton[coords[2],coords[1],coords[0]] = id_

						else:
							skeleton[coords[2],coords[1],coords[0]] = node.id_
				
		if out_path.endswith('tif'):
			tif.imsave(out_path,skeleton)
		

			if gaussian_blur:
				print('Gaussian Blurring (This may take a moment...)')
				blurred_skeleton = gaussian_filter(skeleton,sigma=1)
				tif.imsave(out_path.split('.')[0] + '_gaussianBlurred.tif',blurred_Skeleton)
				print('Done Blurring...')
		else:
			tif.imsave(join(out_path,self.get_swc_name() + '.tif'),skeleton)

	def save_branch_info_for_structure_label_markup(self, out_path, downsampling=1, normalize=False):

		"""
		Params:

		Saves branch info to text file for label markup

		"""
		out_file = join(out_path, self.get_swc_name() + '.txt')

		if normalize:
			mins = np.min(self.generate_swc_array(only_coords=True)/downsampling,axis=0).astype(int)
		else:
			mins = np.array([0,0,0])

		with open(out_file, 'w') as fp:

			# soma info
			soma = (np.array(self.node_soma.get_coords())/downsampling).astype(int)-mins
			fp.write('Soma: ' + str(soma[0]) + ',' + str(soma[1]) + ',' + str(soma[2]) + '\n\n')
			
			# reverse so in increasing order of id
			self.node_soma.child_nodes.reverse()

			for i,child in enumerate(self.node_soma.child_nodes, start=2):
				
				# branch info
				coord = (np.array(child.get_coords())/downsampling).astype(int)-mins
				fp.write(str(i) + ' : ' + str(int(child.id_)) + ' : ' + str(coord[0]) + ' ' + str(coord[1]) + ' ' + str(coord[2]))

				if i < len(self.node_soma.child_nodes) + 1:
					fp.write('\n') 

			# undo reversal
			self.node_soma.child_nodes.reverse()


	def generate_swc_array(self, only_coords=False, *, preserve_sid=True, preserve_radius=True):

		"""
		Return: 2D array, SWC in 2D array
	
		"""

		out_swc = []

		stack = [(self.node_soma, params.SOMA_PID)]

		while len(stack):
		
			node_tup = stack.pop()
			node, pid = node_tup[0], node_tup[1]
			
			out_swc.append([node.id_ , node.structure_id, node.x, node.y, node.z, node.radius, pid])			

			for c in node.child_nodes:
				stack.append((c,node.id_))
	
		out_swc = np.array(out_swc)

		if not preserve_sid:
			out_swc[:,params.SWC_INDICES['sid']] = 0
			
		if not preserve_radius:
			out_swc[:,params.SWC_INDICES['radius']] = 0			

		if only_coords:
			return out_swc[:,[params.SWC_INDICES['x'], params.SWC_INDICES['y'], params.SWC_INDICES['z']]]
		else:
			return out_swc

	def save_swc(self, out_path, fmt=params.SWC_FORMAT_FLOAT, *, preserve_sid=True, preserve_radius=True, from_array=False):

		"""
		Params: outPath - string, path to save SWC
			fmt - array of strings, datatype specifier for each column in SWC
			      ex: ['%d','%d','%d','%d','%d','%d','%d']

		Generates and Saves SWC

		"""

		# create path
		if not '.swc' in out_path:
			out_path = join(out_path, self.get_swc_name(with_ext=True)) 

		# dont overwrite original
		if out_path == self.swc_path:
			exit('Error: Do not overwrite SWC')

		print('Saving SWC:',out_path)

		if from_array:
			swc_array = self.swc_array
		else:
			# remove duplicates
			dups = self.check_for_duplicates()
			if len(dups):
				print('# Duplicates:',len(dups))
			self.remove_duplicates()

			# reset ids
			self.reset_node_ids()

			# generate SWC Array
			swc_array = self.generate_swc_array(preserve_sid=preserve_sid, preserve_radius=preserve_radius)

		np.savetxt(out_path, swc_array, delimiter=params.SWC_DELIMITER, fmt=fmt)


	def get_bif_nodes(self, include_soma=True, node_root=None):
	
		if node_root is None:
			nodeRoot = self.node_soma

		if include_soma:
			return [node for node in self.get_nodes_in_tree(node_root=node_root) if node.is_bifurcation()]
		else:
			return [node for node in self.get_nodes_in_tree(node_root=node_root) if node.is_bifurcation() and not node.is_soma()]

	
	def get_terminal_nodes(self,node_root=None):

		if node_root is None:
			node_root = self.node_soma

		terminal_nodes = []
		stack = [node_root]

		while len(stack):

			node = stack.pop()
			
			if node.is_terminal():
				terminal_nodes.append(node)

			for child_node in node.child_nodes:
				stack.append(child_node)

		return terminal_nodes
	
	def num_terminal_nodes(self, node_root=None):

		"""
		Return: int, number of terminal nodes

		"""

		return len(self.get_terminal_nodes(node_root=node_root))


	def collapse_soma(self,collapse_radius=params.SOMA_COLLAPSE_RADIUS):

		"""
		Params: collapseRadius, float - remove all nodes within collapseRadius from soma

		"""
		
		print()
		print('Collapse Soma: ' + str(collapse_radius) + ' (um)')

		inner_nodes, outer_nodes = self.get_inner_outer_nodes(collapse_radius)
		inner_nodes = [inner_node for inner_node in inner_nodes if inner_node != self.node_soma]
		outer_to_inner_nodes = [inner_node.parent_node for inner_node in inner_nodes if inner_node.parent_node in outer_nodes]

		print('# Nodes (before):',self.num_nodes())
		print('# Inner Nodes:',len(inner_nodes))
		print('# Outer Nodes:',len(outer_nodes))
		print('# Outer to Inner Nodes:',len(outer_to_inner_nodes))

		if len(inner_nodes) == 0:
			print('No Inner Nodes - Do not modify SWC...')
			return

		# get all terminal nodes that are outside radius
		terminal_nodes = [terminal_node for terminal_node in self.get_terminal_nodes() if self.node_soma.distance(terminal_node) >= collapse_radius]

		# track backwards from all terminal nodes until reaches node in radius
		new_soma_child_nodes = []

		for terminal_node in terminal_nodes:
			
			pointer = terminal_node

			while pointer:

				# parent reference might already be removed
				if pointer.parent_node is None:
					break

				if self.node_soma.distance(pointer.parent_node) < collapse_radius:
					pointer.parent_node = None
					new_soma_child_nodes.append(pointer)
					break
		
				pointer = pointer.parent_node

		# reset somas child nodes
		self.node_soma.child_nodes = new_soma_child_nodes
		for child in self.node_soma.child_nodes:
			child.parent_node = self.node_soma


		# remove all children of outer to inner nodes
		# this handles case where branch goes back into soma and then goes back out
		for outer_to_inner_node in outer_to_inner_nodes:
			outer_to_inner_node.child_nodes = []
				
		# Error Checks
		inner_nodes = [node for node in self.get_nodes_in_tree(self.node_soma) if node.distance(self.node_soma) < collapse_radius and node != self.node_soma]
		if len(inner_nodes) > 0:	
			exit('Error: Node found < ' + str(collapse_radius) + ' microns from soma')

		print('# Nodes (after):', self.num_nodes())


	def remove_duplicates(self):

		"""
		Removes duplicate nodes in treee
		Nodes are considered duplicates if same coordinates and are in parent-child relationship
		

		"""

		stack = [self.node_soma]
		num_dup = 0
		
		while len(stack):

			node = stack.pop()

			if node is not self.node_soma:
				if node == node.parent_node:
					
					num_dup += 1

					# get children of parent node excluding node
					parent_node_children = [n for n in node.parent_node.child_nodes if node.id_ != n.id_]
					node.parent_node.child_nodes = node.child_nodes + parent_node_children

					# set childrens parent pointers to grandparent
					for child_node in node.parent_node.child_nodes:
						child_node.parent_node = node.parent_node

			for child_node in node.child_nodes:
				stack.append(child_node)

		if num_dup > 0:
			print('REMOVED DUPLICATES:',num_dup)



	def check_for_duplicates(self):

		"""
		Return: list, contains ids of duplcate nodes

		"""

		coords_visited, duplicate_ids = [], []

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()

			if node is not self.node_soma:
				if node == node.parent_node:
					duplicate_ids.append(node.id_)

			for child_node in node.child_nodes:
				stack.append(child_node)

		return duplicate_ids


	def smooth(self, orientation):

		"""
		Params: orientation, string, 'corona'l or 'oblique', oblique smooths z, coronal smooths y
			

		- smooths coordinates to  (non biological) sharp edges
		- doesnt change soma, children of soma, terminal nodes, bifurcation nodes, or children of bifurcation nodes
		
		"""

		if orientation != 'coronal' and orientation != 'oblique':
			exit('Error: Invalid Orientation!')

		diff = self.smoothing_diffs(orientation)

		diffChange, iters = np.Inf, 0

		while diffChange > params.SMOOTHING_DIFF_THRESHOLD or iters < params.SMOOTHING_ITERS_THRESHOLD:

			# smooth
			stack = [self.node_soma]

			while len(stack):

				node = stack.pop()
	
				# make sure not soma and has 1 child (discludes bifurcation and terminal nodes
				if node is not self.node_soma and len(node.child_nodes) == 1 and len(node.parent_node.child_nodes) == 1:

					if orientation == 'coronal':
						node.tmp = (node.parent_node.y + node.y + node.child_nodes[0].y)/3
					elif orientation == 'oblique':
						node.tmp = (node.parent_node.z + node.z + node.child_nodes[0].z)/3

				for child_node in node.child_nodes:
					stack.append(child_node)

			# update coordinates
			for node in self.get_nodes_in_tree():
				if node.tmp:
					if orientation == 'coronal':
						node.y = node.tmp
					elif orientation == 'oblique':
						node.z = node.tmp


			# recaulate differences
			diffSmoothed = self.smoothing_diffs(orientation)
			diffChange = abs(diffSmoothed - diff)
			diff = diffSmoothed
			iters += 1

		print('# iters:',iters)


	def smoothing_diffs(self, orientation):

		stack = [self.node_soma]

		diff = 0

		while len(stack):

			node = stack.pop()

			for child_node in node.child_nodes:
				if node != self.node_soma:
					if orientation == 'coronal':
						diff += abs(node.y - child_node.y)
					elif orientation == 'oblique':
						diff += abs(node.z - child_node.z)
				stack.append(child_node)

		return diff

	
	def get_inner_outer_nodes(self, collapse_radius):

		inner_nodes, outer_nodes = [], []

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()
			
			if node.distance(self.node_soma) < collapse_radius:
				inner_nodes.append(node)
			else:
				outer_nodes.append(node)

			for child_node in node.child_nodes:
				stack.append(child_node)

		if len(inner_nodes) + len(outer_nodes) != self.num_nodes():
			sys.exit('Error: Inner nodes + Outer nodes != All nodes')

		return inner_nodes, outer_nodes
		
	def set_structure_labels(self, label_info):
		
		self.node_soma.structure_id = 1

		for child_node in self.node_soma.child_nodes:
			
			for label_tup in label_info:
				s_id = None
				if label_tup[0] == child_node.id_:
					s_id = label_tup[1]
					break

			if s_id is None:
				exit('Error: Soma Child with Label Not Found!')

			stack = [child_node]

			while len(stack):

				node = stack.pop()
			
				node.structure_id = s_id

				for child_node in node.child_nodes:
					stack.append(child_node)
		
		# Error Check
		nodes = self.get_nodes_in_tree()
		if sum([node for node in nodes if node.structure_id == 0]) > 0:
			exit('Some nodes do not have label!')

	def register(self):

		"""

		Assumes SWC is in coronal orientation

		"""
		
		swc_info = self.get_swc_info()
		brain = swc_info['brain']

		cropDict = utils.loadCropInfo()

		# apply cropping shift
		shift = cropDict[brain]
		self.translate_coords([0,-shift*params.CROP_FACTOR,0])
				
		# get coordinates
		swcCoords = self.generate_swc_array(only_coords=True)

		brainPath = join(params.REGISTRATION_BASE_PATH, brain)

		# write coords to file
		outPath = join(brainPath,'swc_in.txt')
		utils.writeCoordsForRegistration(outPath, swcCoords)

		# run transformix
		tpPath = join(brainPath,'TransformParameters_labels.1_full_scale.txt')
		inPath = join(brainPath,'swc_in.txt')
		utils.runTransformix(brainPath, tpPath, inPath)

		# read in registered coordinates
		transformedCoordsPath = join(brainPath,'outputpoints.txt')
		registeredCoords = utils.parseTransformixOutput(transformedCoordsPath)

		swc_array = self.generate_swc_array()
		swc_array[:,[params.SWC_INDICES['x'], params.SWC_INDICES['y'], params.SWC_INDICES['z']]] = registeredCoords

		# re build SWC Tree
		self.swc_array = swc_array
		self.build_swc_tree()

		self.scale_coords(params.REGISTRATION_RES_MICRONS)


	def overlayOnReferenceSpace(self, outPath, * , downsampling=25):

		"""
		Params: outPath - string, path to save output image
			downsampling - int, resolution downsampling

		Assumes SWC coordinates are in microns 
		"""

		# read in reference space image
		referenceImage = imread(params.MOP_REFERENCE_PATHS[downsampling])



		# convert to swc array
		swcArray = (self.generateSWCArray(onlyCoords=True)/downsampling).round().astype(int)

		# overlay on reference image
		referenceImage[swcArray[:,2],swcArray[:,1],swcArray[:,0]] = params.OVERLAY_INTENSITY

		# save
		tif.imsave(outPath, referenceImage)


	def get_mopul_layer(self, downsampling=25):

		# read in reference space image
		referenceImage = imread(params.MOP_REFERENCE_PATHS[downsampling])

		# get soma
		soma = self.node_soma.getCoords()

		# get layer
		layer = referenceImage[int(soma[2]/downsampling),int(soma[1]/downsampling),int(soma[0]/downsampling)]

		return layer


	def plot(self, *, dim=2, save=False, out_path=None, x_lim=None, y_lim=None, from_array=False):

		"""
		Params: dim - int, number of dimensions to plot int (2 or 3)

		Plots SWC in 2d or 3d. If 2d, plots x and y
		"""
	
		# covert SWC to array
		if from_array:
			swc_array = self.swc_array
		else:
			swc_array = self.generate_swc_array()

		# get all unique sIDs in SWC		
		unique_sid = np.unique(swc_array[:,params.SWC_INDICES['sid']]).astype(int)		

		if dim == 2:

			# plot 2D
			fig, ax = plt.subplots()

			if x_lim is not None:
				ax.set_xlim(x_lim)
			if y_lim is not None:
				ax.set_ylim(y_lim)
		
			for sid in unique_sid:
				sid_coords = swc_array[swc_array[:,params.SWC_INDICES['sid']] == sid]
				ax.scatter(sid_coords[:,params.SWC_INDICES['x']],-sid_coords[:,params.SWC_INDICES['y']], color=params.SID_COLORS[sid], s=params.SID_PLOTTING_RADII[sid])
		elif dim == 3:

			# plot 3D

			fig = plt.figure()
			ax = fig.add_subplot(111,projection='3d')
			for sid in unique_sid:
				sidCoords = swc_array[swc_array[:,params.SWC_INDICES['sid']] == sid]
				ax.scatter(sidCoords[:,params.SWC_INDICES['x']],-sidCoords[:,params.SWC_INDICES['y']],sidCoords[:,params.SWC_INDICES['z']], color=params.SID_COLORS[sid], s=params.SID_PLOTTING_RADII[sid])			
		else:
			exit('Error: # Dimensions must be 2 or 3')
	
		if save:
			plt.savefig(out_path)
			plt.close('all')
		else:
			plt.show()


	def get_node_from_coords(self, coords):

		""" 
		return node given coords

		"""

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()
			
			if coords == node.get_coords():
				return node

			for child_node in node.child_nodes:
				stack.append(child_node)

		return False
	


	def get_node_from_id(self, node_id):

		""" 
		return node given an ID

		"""

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()
			
			if node.id_ == node_id:
				return node

			for child_node in node.child_nodes:
				stack.append(child_node)
	
		return False


	def prune_from_id(self, node_id, *, prune_backwards=False):
		
		"""
		Removes all nodes that are descendents
		If prune Backwards is True, will also remove all ancestors up to bifurcation point

		"""

		# get node
		node = self.get_node_from_id(node_id)
		
		# remove all descendents
		node.remove_children()

		# remove all ancestors up to bifurcation
		if prune_backwards:
			ancestor_bif = node.get_ancestor_after_bif()
			ancestor_bif.parent_node.child_nodes.remove(ancestor_bif)

		# reset node ids
		self.reset_node_ids()

	def reset_node_ids(self):

		"""
		Resets nodes ID's

		"""

		id_ = 1
		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()

			node.id_ = id_
			id_ += 1
		
			for child_node in node.child_nodes:
				stack.append(child_node)

	def get_branches(self):

		"""
		Return list of Branch objects

		"""

		branch_list = []
		
		stack = [(self.node_soma, self.node_somas)]

		while len(stack):

			node_b, node_a = stack.pop()		

			if (node_b.is_bifurcation() or node_b.is_terminal()) and not node_b.is_soma():
				branch_list.append(Branch(node_a, node_b))
				node_a = node_b


			for child_node in node_b.child_nodes:
				stack.append((child_node, node_a))
			

		return branch_list
			
		
	def contains_any_node_in_branch(self, branch, check_node_a=True):

		pointer = branch.nodeB
		while pointer != branch.nodeA:

			if self.containsNode(pointer):
				return True

			pointer = pointer.parentNode
		

		if checkNodeA and self.containsNode(branch.nodeA):
			return True

		return False

	def contains_nan_coords(self):

		coords = self.swc_array[:,2:5]
		num_nan = np.isnan(coords).any(axis=1).sum()
		
		return np.isnan(coords).any(), num_nan

	def check_nan_coords(self):

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()
			
			if node.coords_are_nan():
				if node.is_bifurcation():
					print(node.id_)

			for child_node in node.child_nodes:
				stack.append(child_node)
	
		return False

	def remove_nan_nodes(self):

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()
			
			# if node is nan and not bif
			# connect parent to child
			if node.coords_are_nan():
	
				#if not node.is_bifurcation():
				#	if not node.is_terminal():
				#		node.child_nodes[0].parent_node = node.parent_node
				#		node.parent_node.child_nodes.append(node.child_nodes[0])
				#	node.parent_node.child_nodes.remove(node)
				#else:
				#	print(node.id_)

				node.parent_node.child_nodes = []
				node.parent_node = None
				
			for child_node in node.child_nodes:
				stack.append(child_node)


	def __str__(self):

		return self.swc_path


if __name__ == '__main__':

	
	#path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/all_structures/2_171012_13.swc'
	#swc = SWC(path)

	#swc.plot(dim=3)


	#path = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/coronal_full_size/190327_68.swc'
	#swc = SWC(path)
	#swc.register()
	#swc.center_around_soma()
	#swc.save_swc('/data/elowsky/')

	#path = '/data/elowsky/OLST/swc_analysis/automatically_traced/pfc/limbic_cortex/oblique_reordered_y_inverted_gcut/181004_5.swc'
	#out_path = '/data/elowsky/OLST/swc_analysis/automatically_traced/pfc/limbic_cortex/oblique_reordered_y_inverted_gcut/'
	#swc = SWC(path)
	#swc.save_skeleton(out_path,compress_ids=True, image_shape=(225,1490,919) )
	#swc.register()
	#swc.center_around_soma()
	#swc.save_swc('/data/elowsky/')

	path = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/coronal_normalized/190123_52.swc'

	swc = SWC(path)
	swc.plot(dim=2)
	#print(swc.check_nan_coords())
	#swc.remove_nan_nodes()
	#print(swc.check_nan_coords())


	







	

	

	







