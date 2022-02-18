class Branch:

	def __init__(self, nodeA, nodeB):
		
		self.nodeA = nodeA
		self.nodeB = nodeB


	def getVectorForward(self, dist=None):

		if dist is None:
			self.nodeA.nodesToVector(self.nodeB)

		# figure out which child to traverse
		for child in self.nodeA.childNodes:

			if self.nodeB.isDescendant(child):
				localDist = -1
				pointer = child
				while pointer != self.nodeB:
					localDist = self.nodeA.distance(pointer)
					if localDist >= dist:
						break
					pointer = pointer.childNodes[0]
				
				if localDist < dist:
					pointer = self.nodeB

				return self.nodeA.nodesToVector(pointer)
	
		exit('Error: Node B is not in any subtree')		
			
	def getVectorBackward(self, dist=None):
		
		if dist is None:
			self.nodeB.nodesToVector(self.nodeA)

		pointer = self.nodeB
		
		while pointer != self.nodeA:
			localDist = self.nodeB.distance(pointer)
			if localDist >= dist:
				break
			pointer = pointer.parentNode

		if localDist < dist:
			pointer = self.nodeA

		return self.nodeB.nodesToVector(pointer) 

		

	def eucDistance(self):

		"""
		Euclidean distsance between starting and eneding node of branch		

		"""
		
		return self.nodeA.distance(self.nodeB)

	def pathDistance(self):

		"""
		Path distsance between starting and eneding node of branch		

		"""

		pointer = self.nodeB
		dist = 0

		while pointer != self.nodeA:

			dist += pointer.distance(pointer.parentNode)
			pointer = pointer.parentNode

		return dist

	def numNodes(self):

		return len(self.getNodesInBranch())

	def getNodesInBranch(self):

		"""
		Returns list of all nodes in branch

		"""
		
		nodes = []
		pointer = self.nodeB

		while pointer != self.nodeA:

			nodes.append(pointer)
			pointer = pointer.parentNode

		# append node A
		nodes.append(pointer)

		# reverse order so node A is firt
		nodes.reverse()

		return nodes

	def numTerminalNodes(self):

		stack = [self.nodeA]
		numTerminalNodes = 0

		while len(stack):

			node = stack.pop()
			
			if node.isTerminal():
				numTerminalNodes += 1

			for childNode in node.childNodes:
				stack.append(childNode)

		return numTerminalNodes
		

	def __eq__(self, otherBranch):

		return self.nodeA == otherBranch.nodeA and self.nodeB == otherBranch.nodeB

	def __repr__(self):

		return str(self.nodeA)+' - '+str(self.nodeB)


if __name__ == '__main__':
	
	pass
