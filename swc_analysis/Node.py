import math
import params

from sys import exit


class Node:

	def __init__(self, id_, structure_id, x, y, z, radius, parent_node=None):

		self.id_ = id_
		self.structure_id = structure_id
		self.x, self.y, self.z = x,y,z
		self.radius = radius
		self.parent_node = parent_node
			
		self.child_nodes = []

		# temporary coord for smoothing
		self.tmp = None

	@classmethod
	def node_from_array(cls, swc_row, *, parent_node=None):

		"""
		Alternative constructor that takes a 1D array

		"""
		
		if len(swc_row) != len(params.SWC_INDICES):
			exit('ERROR: SWC Array has length ' + str(len(swc_row)) + ' but must be ' + str(len(params.SWC_INDICES)))

		node = cls(
			swc_row[params.SWC_INDICES['id']],
			swc_row[params.SWC_INDICES['sid']],
			swc_row[params.SWC_INDICES['x']],
			swc_row[params.SWC_INDICES['y']],
			swc_row[params.SWC_INDICES['z']],
			swc_row[params.SWC_INDICES['radius']],
			parent_node=parent_node,
		)

		return node

	def coords_are_nan(self):

		return math.isnan(self.x) or math.isnan(self.y) or math.isnan(self.z)


	def nodes_to_vector(self, other_node):
		
		"""
		Subtracts two nodes to create a vector

		"""

		return [other_node.x - self.x, other_node.y - self.y, other_node.z - self.z]

	def distance(self, other_node):
		
		"""
		Euclidean Distance between self and another node

		"""

		return math.sqrt((self.x-other_node.x)**2 + (self.y-other_node.y)**2 + (self.z-other_node.z)**2)

	def is_soma(self):

		"""
		Boolean whether node is soma

		"""

		return True if self.structure_id == 1 or self.parent_node is None else False

	def num_children(self):

		"""
		Number of child nodes

		"""

		return len(self.child_nodes)

	def is_bifurcation(self):
		
		"""
		Boolean whether node has more than one child
	
		"""

		return True if self.num_children() > 1 else False

	def is_terminal(self):
	
		"""
		Boolean whether node is a terminal node
	
		"""

		return True if self.num_children() == 0 else False

	def get_coords(self):

		"""
		Returns coordinates in a list

		"""

		return [self.x, self.y, self.z]

	def get_ancestor_bif(self):

		"""
		Returns first ancestor of node that is a bifurcation node

		"""
		if self.parent_node is None:
			return False

		pointer = self.parent_node

		while not pointer.is_bifurcation():
			pointer = pointer.parent_node

			if pointer is None:
				return False

		return pointer

	def get_ancestor_after_bif(self):

		"""
		Returns first ancestor of node that is a child of first bifurcation node

		"""

		pointer = self

		while True:
	
			if pointer.parent_node is None:
				return False

			if pointer.parent_node.is_bifurcation():
				return pointer

			pointer = pointer.parent_node

		return pointer

	def get_bif_or_term_descendent(self):

		"""
		Returns first descendant of node that is a bifurcation or terminal node
	
		"""

		pointer = self
		
		while not pointer.is_terminal() and not pointer.is_bifurcation():
			pointer = pointer.child_nodes[0]

		return pointer

	def translate(self, shifts):
		
		"""
		Translates coordinates by shifts
	
		"""

		if isinstance(shifts, float) or isinstance(shifts, int):
			shifts = [shifts]*3
		
		self.x += shifts[0]
		self.y += shifts[1]
		self.z += shifts[2]

	def scale(self, scale_factors):

		"""
		Scales coordinates by scaleFactors
	
		"""

		if isinstance(scale_factors, float) or isinstance(scale_factors, int):
			scale_factors = [scale_factors]*3
		
		self.x *= scale_factors[0]
		self.y *= scale_factors[1]
		self.z *= scale_factors[2]

	def remove_children(self):
	
		"""
		Removes all children of a given node
	
		"""

		self.child_nodes = []

	def get_all_descendants(self):

		"""
		returns list of all descendents

		"""

		nodes = []

		stack = [self]

		while len(stack):

			node = stack.pop()
			nodes.append(node)

			for child_node in node.child_nodes:
				stack.append(child_node)

		nodes.pop(0)

		return nodes


	def is_descendant(self, other_node):

		"""	
		checks if 'this' node is descendent of otherNode

		"""

		stack = [other_node]

		while len(stack):

			node = stack.pop()

			if node == self:
				return True		

			for child_node in node.child_nodes:
				stack.append(child_node)

		return False			 


	def __eq__(self, other_node):

		"""
		Nodes are considered equal if their coordinates are the same		
	
		"""

		return self.x == other_node.x and self.y == other_node.y and self.z == other_node.z

	def __repr__(self):

		return str(self.id_) + ': (' + str(self.x)+','+str(self.y)+','+str(self.z)+')'



if __name__ == '__main__':
	
	a = Node(1,0,1,1,1,0)
	b = Node(2,0,2,2,2,0,a)
	c = Node(3,0,3,3,3,0,b)
	d = Node(4,0,4,4,4,0,a)
	a.childNodes = [b,d]
	b.childNodes = [c]


	print(a.getAllDescendants())




