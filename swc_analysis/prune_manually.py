from os import listdir, mkdir
from os.path import isfile, join, exists
from sys import exit
from SWC import SWC
from mpl_toolkits.mplot3d import Axes3D
from skimage import io

import params
import shutil
import numpy as np
import matplotlib.pyplot as plt

########################################################################################
REGION = 'barrel_cortex'
SWC_INPATH = '/data/elowsky/OLST/test/'
SWC_OUTPATH = '/data/elowsky/OLST/test_out/'
CROPPED_DIRECTORY = join('/data/palmer/data/reconstructed_brains/', REGION)
ORIGINAL_SWC_PATH = join('/data/palmer/data/swc_files/', REGION)
########################################################################################

def update_annot(ind):
	pos = sc.get_offsets()[ind["ind"][0]]
	annot.xy = pos
	annot.set_text(int(swcArray[ind["ind"][0],0]))

def press(event):
	
	if event.key == 'enter':

		# save swc 

		swc.saveSWC(SWC_OUTPATH)
		plt.close()
	elif event.key == 'escape':

		# exit program

		exit()
	elif event.key == 'n':

		# go to next swc
	
		plt.close()

	elif event.key == 'x':
		
		# connect soma to node

		nodeID = input('Enter node ID to connect to soma: ')

		while not nodeID.isdigit():
			print('Please enter an integer...')
			nodeID = input('Enter node ID to connect to soma: ')

		nodeID = int(nodeID)

		node = swc.getNodeFromID(nodeID)

		# remove from parents children
		node.parentNode.childNodes.remove(node)

		# make parent soma
		node.parentNode = swc.nodeSoma

		# append to soma children
		swc.nodeSoma.childNodes.append(node)
		
		# reset node ids
		swc.resetNodeIds()

		# update plot	
		update_plot_and_save_skeleton()

		print()

	elif event.key == 'c':
		
		# press 'c' to connect two nodes
		# first node cant be a bifurcation
		# second node must be a descendant of the first

		firstID = input('Enter first ID (parent) to connect: ')

		while not firstID.isdigit():
			print('Please enter an integer...')
			firstID = input('Enter first ID (parent) to connect: ')

		firstID = int(firstID)

		# dont let first node be bifurcation
		while swc.getNodeFromID(firstID).isBifurcation():
			print('You chose a bifurcartion node, please choose a different node...')
			firstID = int(input('Enter first ID (parent) to connect: '))

		# get first node
		firstNode = swc.getNodeFromID(firstID)
			
		secondID = input('Enter second ID (child) to connect: ')
	
		while not secondID.isdigit():
			print('Please enter an integer...')
			secondID = input('Enter second ID (child) to connect: ')

		secondID = int(secondID)
		
		# get second node
		secondNode = swc.getNodeFromID(secondID)

		# make sure second node is descendant of first node
		while not secondNode.isDescendant(firstNode):
			print('Second Node must be a descendant of the first, please choose a different node...')
			secondID =int(input('Enter second ID (child) to connect: '))
			secondNode = swc.getNodeFromID(secondID)

		# remove all children of start node
		firstNode.removeChildren()

		# connect parent to child
		firstNode.childNodes = [secondNode]
		secondNode.parentNode = firstNode

		# reset node ids
		swc.resetNodeIds()

		# update plot	
		update_plot_and_save_skeleton()

		print()
		
	elif event.key == 'r':

		# press 'r' in to remove a branch
		# option to remove only descendents of node chosen or entire branch up to first bifurcaiton ancestor

		idToRemove = input('Enter ID of Branch to Remove: ')

		while not idToRemove.isdigit():
			print('Please Enter an Integer...')
			idToRemove = input('Enter ID of Branch to Remove: ')
		
		idToRemove = int(idToRemove)
		pruneBackwards = bool(int(input('Should Anscestors be Pruned (Enter 1 or 0): ')))

		print('ID Selected:',idToRemove)
		print('Removing Ancestors:',pruneBackwards)

		# prune
		swc.pruneFromID(idToRemove, pruneBackwards=pruneBackwards)

		# update plot
		update_plot_and_save_skeleton()

		print()


def update_plot_and_save_skeleton():

	global sc, annot, swcArray

	# convert to cropped coordinates
	swc.normalizedToCropped('oblique', originalSomaCoords, croppedShape)

	# write skeleton 
	skeletonOutPath = join(skeletonPath, str(swc.numNodes()) + '.tif')
	swc.saveSkeleton(skeletonOutPath, imageShape=croppedShape, gaussianBlur=True)

	# clear plot
	plt.cla()
	
	# normalize
	swc.pixelsToMicrons(orientation='oblique', centerAroundSoma=True)

	# get coordinates
	swcArray = swc.generateSWCArray()
	
	# plot soma
	ax.scatter(swcArray[0,params.SWC_INDICES['x']],swcArray[0,params.SWC_INDICES['y']],swcArray[0,params.SWC_INDICES['z']],s=100)

	# plot tracing points
	sc = ax.scatter(swcArray[:,params.SWC_INDICES['x']],swcArray[:,params.SWC_INDICES['y']],swcArray[:,params.SWC_INDICES['z']],s=.2)

	# plot lines connecting soma to children
	for child in swc.nodeSoma.childNodes:
		ax.plot([0,child.getCoords()[0]],[0,child.getCoords()[1]],[0,child.getCoords()[2]],color='black')

	annot = ax.annotate("",xy=(0,0),xytext=(20,20),textcoords="offset points",bbox=dict(boxstyle="round",fc="w"),arrowprops=dict(arrowstyle="->"))
	annot.set_visible(False)
	ax.grid(False)
	ax.axis('off')


def hover(event):
	vis = annot.get_visible()
	if event.inaxes == ax:
		cont, ind = sc.contains(event)

		if cont:
			update_annot(ind)
			annot.set_visible(True)
			fig.canvas.draw_idle()
		else:
			if vis:
				annot.set_visible(False)
				fig.canvas.draw_idle()



if __name__ == '__main__':

	print()
	print('##################')
	print('SWC Branch Removal')
	print('##################')
	print()
	
	# get list of swcs that need to be pruned
	completed = [f for f in listdir(SWC_INPATH) if isfile(join(SWC_OUTPATH, f)) and f.endswith('.swc')]

	# get list of swc for branch removal
	swcFiles = [f for f in listdir(SWC_INPATH) if isfile(join(SWC_INPATH, f)) and f.endswith('.swc') and f not in completed]

	print('# SWCs Completed:', len(completed))
	print('# SWCs Remaining:',len(swcFiles))
	print()


	# iterate through SWCs
	for swcFile in sorted(swcFiles):

		print('###########')
		print(swcFile[:-4])
		print('###########')

		# load SWC
		swcPath = join(SWC_INPATH, swcFile)
		swc = SWC(swcPath)
		swcInfo = swc.getSWCInfo()
	
		# create output diretory for tif skeletons
		skeletonPath = join(SWC_OUTPATH, swc.getSWCName() + '_skeletons')
		if exists(skeletonPath):
			shutil.rmtree(skeletonPath)
		mkdir(skeletonPath)
		
		# get dimensions from cropped image
		croppedPath = join(CROPPED_DIRECTORY,swcInfo['brain'] + '_reconstructed', str(swcInfo['id']))
		rawCroppedFiles = [f for f in listdir(croppedPath) if 'raw_cropped' in f]
		if len(rawCroppedFiles) != 1:
			exit('Error: Amount of raw_cropped files is ' + str(len(rawCroppedFiles)) + ' but needs to be 1')
		croppedPath = join(croppedPath,rawCroppedFiles[0])
		print('Raw Cropped Path:', croppedPath)
		print('Loading raw image (this may take a moment)...')
		croppedRaw = io.imread(croppedPath)
		croppedShape = croppedRaw.shape
			
		# get soma from original SWC
		# check if exits in gcut folder first
		if exists(join(ORIGINAL_SWC_PATH,'gcut_neurons',swc.getSWCName(wExt=True))):
			originalSWCpath = join(ORIGINAL_SWC_PATH,'gcut_neurons',swc.getSWCName(wExt=True))
		else:
			originalSWCpath = join(ORIGINAL_SWC_PATH, swc.getSWCName(wExt=True))
		print('Original SWC Path:', originalSWCpath)
		swcOriginal = SWC(originalSWCpath)
		originalSomaCoords = swcOriginal.getSomaCoords()

		# set up swc plot
		fig = plt.figure()
		fig.canvas.set_window_title(swcInfo['brainSWCID'])
		ax = fig.add_subplot(111, projection='3d')

		# plot swc
		update_plot_and_save_skeleton()

		# set up event listeners
		fig.canvas.mpl_connect('motion_notify_event', hover)
		fig.canvas.mpl_connect('key_press_event', press)

		# plot
		plt.show()

		print()


