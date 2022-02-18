from SWC import SWC
import numpy as np
import json

"""
SWC object for recreating Lmeasure functions

	- this functionality could have been included directly in SWC
	  but chose to keep separate as there are some peculiarities 
	  that lmeasure has when only one node is a soma
	
	- the functions in this class are implemented to match
	  lmeasure (regardless of whether it is correct or not)

	- many of the standard deviations do not match those of 
	  lmeasure (not sure what lmeasure is doing differently)

"""

class Lmeasure_SWC(SWC):
	
	def N_stems(self):

		n_stems = self.node_soma.num_children()

		return {'total':n_stems}

	def N_bifs(self, include_soma=True):

		n_bifs = len(self.get_bif_nodes())

		if not include_soma:
			n_bifs -= 1

		return {'total':n_bifs}

	def N_branch(self):

		n_branches = self.N_bifs(include_soma=False)['total'] + self.num_terminal_nodes()
	
		# correction for lmeasure single soma issue
		n_branches += 2

		return {'total':n_branches}

	def Width(self):

		coords = [node.x for node in self.get_nodes_in_tree()]

		coords = np.sort(np.array(coords))
		
		#coords = coords[int(len(coords)*.025):len(coords)-int(len(coords)*.025)]

		width = coords[-1] - coords[0]	

		return {'total':width}

	def Height(self):

		coords = np.array([node.y for node in self.get_nodes_in_tree()])

		coords = np.sort(coords)
		
		#coords = coords[int(len(coords)*.025):len(coords)-int(len(coords)*.025)]

		height = coords[-1] - coords[0]	

		return {'total':height}

	def Depth(self):

		coords = np.array([node.z for node in self.get_nodes_in_tree()])
		coords = np.sort(coords)
		
		#coords = coords[int(len(coords)*.025):len(coords)-int(len(coords)*.025)]

		depth = coords[-1] - coords[0]	

		return {'total':depth}

	def Length(self):

		lengths = []

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()

			for child_node in node.child_nodes:
				lengths.append(node.distance(child_node))
				stack.append(child_node)

		return {'total':np.sum(lengths),'min':np.min(lengths),'avg':np.mean(lengths),'max': np.max(lengths),'std':np.std(lengths)}

	def EucDistance(self):
		
		distances = []

		stack = [self.node_soma]

		while len(stack):

			node = stack.pop()
			distances.append(node.distance(self.node_soma))

			for child_node in node.child_nodes:
				stack.append(child_node)
		
		return {'total':np.sum(distances),'min':np.min(distances),'avg':np.mean(distances),'max': np.max(distances),'std':np.std(distances)}



	def PathDistance(self):

		distances = []

		stack = [(self.node_soma,0)]

		while len(stack):

			node_tup = stack.pop()
			node = node_tup[0]
			dist = node_tup[1]

			distances.append(dist)

			for child_node in node.child_nodes:
				stack.append((child_node,dist + node.distance(child_node)))

		return {'total':np.sum(distances),'min':np.min(distances),'avg':np.mean(distances),'max': np.max(distances),'std':np.std(distances)}

	def Branch_Order(self):

		orders = []

		stack = [(self.node_soma,0)]

		while len(stack):

			node_tup = stack.pop()
			node = node_tup[0]
			order = node_tup[1]

			orders.append(order)

			for child_node in node.child_nodes:

				if node.is_bifurcation() and not node.is_soma():
					stack.append((child_node, order+1))
				else:
					stack.append((child_node, order))
		
		# correction for lmeasure single soma issue
		orders.extend([0,0])

		return {'total':np.sum(orders),'min':np.min(orders),'avg':np.mean(orders),'max':np.max(orders),'std':np.std(orders)}

	def Terminal_degree(self):

		terminal_degrees = []

		all_nodes = self.get_nodes_in_tree()

		for node in all_nodes:
			terminal_degrees.append(self.num_terminal_nodes(node_root=node))

		# correction for lmeasure single soma issue
		terminal_degrees[0] += 2
		terminal_degrees.extend([1,1])
			
		return {'total':np.sum(terminal_degrees),'min':np.min(terminal_degrees),'avg':np.mean(terminal_degrees),'max':np.max(terminal_degrees),'std':np.std(terminal_degrees)}

	def TerminalSegment(self):

		terminal_segments = 0

		terminal_nodes = self.get_terminal_nodes()

		for terminal_node in terminal_nodes:
		
			pointer = terminal_node
			
			while pointer is not None and not pointer.is_bifurcation():
				terminal_segments += 1
				pointer = pointer.parent_node

		# correction for lmeasure single soma issue
		terminal_segments += 2		

		return {'total':np.sum(terminal_segments)}

	def Branch_pathLength(self):

		lengths = []

		stack = [(self.node_soma,0)]

		while len(stack):

			node_tup = stack.pop()
			node, length = node_tup

			if not node.is_soma() and (node.is_bifurcation() or node.is_terminal()):
				lengths.append(length)

			for child_node in node.child_nodes:
				if node.num_children() == 1:
					stack.append((child_node,length+node.distance(child_node)))
				else:
					stack.append((child_node,node.distance(child_node)))
		
		return {'total':np.sum(lengths),'min':np.min(lengths),'avg':np.mean(lengths),'max': np.max(lengths),'std':np.std(lengths)}


	def Contraction(self):

		contractions = []

		stack = [(self.node_soma,0,None)]

		while len(stack):

			node_tup = stack.pop()
			node, length, parent_bif_node = node_tup

			if parent_bif_node is not None and not node.is_soma() and (node.is_bifurcation() or node.is_terminal()):
				euc_dist = node.distance(parent_bif_node)
				contractions.append(euc_dist/length)

			for child_node in node.child_nodes:
				if node.num_children() == 1:
					stack.append((child_node,length+node.distance(child_node),parent_bif_node))
				else:
					stack.append((child_node,node.distance(child_node),node))
		
		if len(contractions):
			return {'total':np.sum(contractions),'min':np.min(contractions),'avg':np.mean(contractions),'max': np.max(contractions),'std':np.std(contractions)}
		else:
			return {'total':0,'min':0,'avg':0,'max':0,'std':0}

	def Partition_asymmetry(self):

		values = []
		
		bifurcation_nodes = self.get_bif_nodes(include_soma=False)

		for bif_node in bifurcation_nodes:
			
			n1 = self.num_terminal_nodes(node_root=bif_node.child_nodes[0])
			n2 = self.num_terminal_nodes(node_root=bif_node.child_nodes[1])
	
			if n1 + n2 <= 2:
				values.append(0)
			else:
				values.append(abs(n1-n2) / (n1 + n2 - 2))

		# correction for lmeasure single soma issue
		values.append(1)

		return {'total':np.sum(values),'min':np.min(values),'avg':np.mean(values),'max': np.max(values),'std':np.std(values)}

	def Bif_ampl_remote(self):

		angles = []
		
		bifurcation_nodes = self.get_bif_nodes(include_soma=False)


		for bif_node in bifurcation_nodes:
			
			# descendant nodes
			pointer_a = bif_node.child_nodes[0].get_bif_or_term_descendent()
			pointer_b = bif_node.child_nodes[1].get_bif_or_term_descendent()

			# get vectors
			vec_a = np.array(pointer_a.get_coords()) - np.array(bif_node.get_coords())
			vec_b = np.array(pointer_b.get_coords()) - np.array(bif_node.get_coords())

			# angle
			angle = np.arccos(np.dot(vec_a/np.linalg.norm(vec_a),vec_b/np.linalg.norm(vec_b)))*(180/np.pi)
			angles.append(angle)
		
		if len(angles) == 0:
			angles = [0]		

		return {'total':np.sum(angles),'min':np.min(angles),'avg':np.mean(angles),'max': np.max(angles),'std':np.std(angles)}

	def Bif_tilt_remote(self):

		angles = []
		
		bifurcation_nodes = self.get_bif_nodes(include_soma=False)

		for bif_node in bifurcation_nodes:

			# get ancestor bif
			ancestor_bif = bif_node.get_ancestor_bif()
			if not ancestor_bif:
				angles.append(0)
				continue
			
			# descendant nodes
			pointer_a = bif_node.child_nodes[0].get_bif_or_term_descendent()
			pointer_b = bif_node.child_nodes[1].get_bif_or_term_descendent()
			
			# get vectors
			vec_ancestor = np.array(ancestor_bif.get_coords()) - np.array(bif_node.get_coords())
			vec_a = np.array(pointer_a.get_coords()) - np.array(bif_node.get_coords())
			vec_b = np.array(pointer_b.get_coords()) - np.array(bif_node.get_coords())

			# angles
			angle_a = np.arccos(np.dot(vec_a/np.linalg.norm(vec_a),vec_ancestor/np.linalg.norm(vec_ancestor)))*(180/np.pi)
			angle_b = np.arccos(np.dot(vec_b/np.linalg.norm(vec_b),vec_ancestor/np.linalg.norm(vec_ancestor)))*(180/np.pi)

			angles.append(min([angle_a,angle_b]))

		if len(angles) == 0:
			angles = [0]

		return {'total':np.sum(angles),'min':np.min(angles),'avg':np.mean(angles),'max': np.max(angles),'std':np.std(angles)}

	def sholl(self):
		pass

	def runAllFunctions(self):

		print('N_stems',self.N_stems())
		print('N_bifs', self.N_bifs())
		print('N_branches', self.N_branch())
		print('Length', self.Length())
		print('EucDistance', self.EucDistance())
		print('PathDistance', self.PathDistance())
		print('Branch_pathLength', self.Branch_pathLength())
		print('Branch_Order', self.Branch_Order())
		print('TerminalSegment', self.TerminalSegment())
		print('Terminal_degree', self.Terminal_degree())
		print('Contraction', self.Contraction())
		print('Partition Assymetry', self.Partition_asymmetry())
		print('Width', self.Width())
		print('Height', self.Height())
		print('Depth', self.Depth())
		print('Bif Ampl Remote', self.Bif_ampl_remote())
		print('Bif Tilt Remote', self.Bif_tilt_remote())

	def run_functions_from_json(self, json_path):

		with open(json_path,'r') as fp:
			morphometrics_dict = json.load(fp)

		MORPHOMETRIC_NAME_TO_FUNCTION = {
				'Bif_ampl_remote':self.Bif_ampl_remote,
				'Bif_tilt_remote':self.Bif_tilt_remote,
				'Branch_Order':self.Branch_Order,
				'Branch_pathlength':self.Branch_pathLength,
				'Contraction':self.Contraction,
				'Depth':self.Depth,
				'EucDistance':self.EucDistance,
				'Height':self.Height,
				'Length':self.Length,
				'N_bifs':self.N_bifs,
				'N_branch':self.N_branch,
				'N_stems':self.N_stems,
				'Partition_asymmetry':self.Partition_asymmetry,
				'PathDistance':self.PathDistance,
				'TerminalSegment':self.TerminalSegment,
				'Terminal_degree':self.Terminal_degree,
				'Width':self.Width,
		}

		MORPHOMETRIC_TYPE_CONVERSION = {
				'Average':'avg',
				'Maximum':'max',
				'Total_Sum':'total',
		}


		for morphometric_name in sorted(morphometrics_dict):
			
			metric_list = morphometrics_dict[morphometric_name]
			morpho_data = MORPHOMETRIC_NAME_TO_FUNCTION[morphometric_name]()


			for metric in metric_list:
				
				data = morpho_data[MORPHOMETRIC_TYPE_CONVERSION[metric]]
				print(morphometric_name, '-', metric, '-',np.round(data,4))
				
				

	



if __name__ == '__main__':

	swc_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/swcs/removed_nan_swcs/basal_apical/2_171012_9.swc'
	json_path = '/data/elowsky/OLST/swc_analysis/paper/analysis/pia_white_matter_normalized/lmeasure_morphometrics.json'
	swc = Lmeasure_SWC(swc_path)
	swc.run_functions_from_json(json_path)


