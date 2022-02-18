from os import listdir, remove
from os.path import join, isfile
from sys import exit
from SWC import SWC
from Lmeasure_SWC import Lmeasure_SWC 
from scipy.stats import zscore
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.decomposition import PCA
from utils import angleBetweenVectors
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.setrecursionlimit(100000)

#######################################################
BASE_PATH = '/data/elowsky/OLST/swc_analysis/pruning/'
LOCAL_ANGLE_DIST = 10
#######################################################

def compileBranchFeatures(outpath):

	beforePruningPath = join(BASE_PATH,'before_pruning')
	afterPruningPath = join(BASE_PATH,'after_pruning')

	# get filenames
	beforePruningFiles = listdir(beforePruningPath)
	afterPruningFiles = listdir(afterPruningPath)

	# remove before files that are not included in after
	for f in beforePruningFiles:
	
		if not isfile(join(BASE_PATH, 'after_pruning', f)):
			print('Removing File (Not used in pipeline):',f)
			remove(join(BASE_PATH, 'before_pruning', f))

	# just making sure all after pruning files have <= num nodes
	# plus some preliminary data exploration

	afterPruningFiles.sort()

	numNotPruned = 0
	numPruned = 0

	for swcName in afterPruningFiles:

		beforeSize = len(np.genfromtxt(join(beforePruningPath,swcName)))
		afterSize = len(np.genfromtxt(join(afterPruningPath,swcName)))

		if afterSize > beforeSize:
			exit('Error: After Pruning size is greater than before pruning size for ' + swcName)
		elif afterSize == beforeSize:
			numNotPruned += 1
		else:
			numPruned += 1

	print()
	print('########################')
	print('# SWCs:', len(afterPruningFiles))
	print('# Pruned:', numPruned)
	print('# Not Pruned:', numNotPruned)
	print('########################')
	print()

	# for all swcs (before), extract metrics, and also get label
	data = []
	labelsDist = [0]*6

	for swcName in afterPruningFiles:
		

		print()
		print('----------')
		print(swcName.split('.')[0])
		print('----------')

		swcBefore = SWC(join(beforePruningPath,swcName))
		swcAfter = SWC(join(afterPruningPath, swcName))

		print()
		print('# Nodes (before):', swcBefore.numNodes())
		print('# Nodes (after):', swcAfter.numNodes())

		branchesBefore = swcBefore.getBranches()
		branchesAfter = swcAfter.getBranches()

		print()
		print('# Branches (before):',len(branchesBefore))
		print('# Branches (after):',len(branchesAfter))
		print()

		# classify every before branch
		#	0 - branch start node is soma
		#	1 - branch pruned at bifurcation
		#	2 - branch pruned not at bifurcation
		#	3 - branch is descendant of pruned branch 
		#	4 - branch fully in after
		#	5 - branch in after but other branch was removed at bifurcation

	
		branchesBeforeLabels = []
		for i, branchBefore in enumerate(branchesBefore):
	
			label = None

			# case 0 - start node of branch is soma (regardless of whether it is removed or not, assumption of algorithm will be to not remove these)
			if branchBefore.nodeA.isSoma():
				label = 0
				branchesBeforeLabels.append((branchBefore,label))
				continue

			# case 1- branch pruned at bifurcation
			if swcAfter.containsNode(branchBefore.nodeA) and not swcAfter.containsNode(branchBefore.nodeB) and branchBefore.nodeA.isBifurcation():

				# make sure no nodes of branch (besides node A) are in swc after
				if not swcAfter.containsAnyNodeInBranch(branchBefore, checkNodeA=False):
					label = 1
					branchesBeforeLabels.append((branchBefore,label))
					continue

			# case 2 - branch pruned not at bifurcation
			if swcAfter.containsNode(branchBefore.nodeA) and not swcAfter.containsNode(branchBefore.nodeB):				
				label = 2
				branchesBeforeLabels.append((branchBefore,label))
				continue

			# case 3 - branch is descendant of pruned
			if not swcAfter.containsNode(branchBefore.nodeA) and not swcAfter.containsNode(branchBefore.nodeB):
				label = 3
				branchesBeforeLabels.append((branchBefore,label))
				continue

			# case 4 - branch is fully in after
			for branchAfter in branchesAfter:
				if branchBefore == branchAfter:
					label = 4
					branchesBeforeLabels.append((branchBefore,label))
					break
		
			if not label is None:
				continue

			# case 5 - in after but other branch was removed at bifurcation
			if swcAfter.containsNode(branchBefore.nodeA) and swcAfter.containsNode(branchBefore.nodeB):
				label = 5
				branchesBeforeLabels.append((branchBefore,label))
				continue
		
			exit('Error: Should never get here')


		
	
		# get distribution of labels		
		labels = [b[1] for b in branchesBeforeLabels]
		print(labels)


		for i in range(6):
			labelsDist[i] += labels.count(i)

		# get features
		for branchBeforeLabel in branchesBeforeLabels:
		
			branch, label = branchBeforeLabel

			# list to store features for branch
			branchFeatures = [swcBefore.getSWCName(), label]

			# local distance features
			eucDistLocal = branch.eucDistance()
			branchFeatures.append(eucDistLocal)	
			pathDistLocal = branch.pathDistance()
			branchFeatures.append(pathDistLocal)

			# local angles

			# find first node that is descendent of A whos euc distance is >= dist	
			vectorForward = branch.getVectorForward(dist=LOCAL_ANGLE_DIST)


			# find all other vectors for all other descendant branches of nodeA
			forwardAngles = []
			for branchOther in branchesBefore:
				if branchOther.nodeA == branch.nodeA and branchOther.nodeB != branch.nodeB:
					vectorOther = branchOther.getVectorForward(dist=LOCAL_ANGLE_DIST)
					angle = angleBetweenVectors(vectorForward, vectorOther)
					forwardAngles.append(angle)
			if len(forwardAngles):
				avgDescAngle = np.mean(forwardAngles)
			else:
				avgDescAngle = 0
			branchFeatures.append(avgDescAngle)


			# get parent branch
			parentBranch = None
			for branchOther in branchesBefore:
				if branchOther.nodeB == branch.nodeA:
					parentBranch = branchOther
					break
			
			if parentBranch:
				vectorBackward = parentBranch.getVectorBackward(dist=LOCAL_ANGLE_DIST)
				anscAngle = angleBetweenVectors(vectorForward, vectorBackward)
			else:
				anscAngle = 0
			branchFeatures.append(anscAngle)


			# global features

			# number of terminal nodes
			numTerminalNodes = branch.numTerminalNodes()
			branchFeatures.append(numTerminalNodes)

			# lmeasure features for tree that follows branch
			
			# re read in swc				
			lmeasureSWC = Lmeasure_SWC(join(beforePruningPath,swcName))
			nodeA = lmeasureSWC.getNodeFromCoords(branch.nodeA.getCoords())
			nodeA.parentNode = None

			for child in nodeA.childNodes:
				if branch.nodeB.isDescendant(child):
					nodeA.childNodes = [child]
					break
			lmeasureSWC.nodeSoma = nodeA
			
			numBifs = lmeasureSWC.N_bifs()['total']
			branchFeatures.append(numBifs)

			numBranches = lmeasureSWC.N_branch()['total']
			branchFeatures.append(numBranches)

			width = lmeasureSWC.Width()['total']
			branchFeatures.append(width)

			height = lmeasureSWC.Height()['total']
			branchFeatures.append(height)

			depth = lmeasureSWC.Depth()['total']
			branchFeatures.append(depth)

			length = lmeasureSWC.Length()['total']
			branchFeatures.append(length )

			eucDistance = lmeasureSWC.EucDistance()['max']
			branchFeatures.append(eucDistance)

			pathDistanceAvg = lmeasureSWC.PathDistance()['avg']
			branchFeatures.append(pathDistanceAvg)

			pathDistanceMax = lmeasureSWC.PathDistance()['max']
			branchFeatures.append(pathDistanceMax)

			branchOrder = lmeasureSWC.Branch_Order()['max']
			branchFeatures.append(branchOrder)

			terminalDegree = lmeasureSWC.Terminal_degree()['max']
			branchFeatures.append(terminalDegree)

			terminalSegment = lmeasureSWC.TerminalSegment()['total']
			branchFeatures.append(terminalSegment)

			branchPathLength = lmeasureSWC.Branch_pathLength()['avg']
			branchFeatures.append(branchPathLength)

			contraction = lmeasureSWC.Contraction()['avg']
			branchFeatures.append(contraction)

			partitionAssymetry = lmeasureSWC.Partition_asymmetry()['avg']
			branchFeatures.append(partitionAssymetry)

			bifAmpleRemote = lmeasureSWC.Bif_ample_remote()['avg']
			branchFeatures.append(bifAmpleRemote)

			bifTiltRemote = lmeasureSWC.Bif_tilt_remote()['avg']
			branchFeatures.append(bifTiltRemote)

			# append to data
			data.append(branchFeatures)



	print(labelsDist)
	data = np.array(data)

	print('Branch Label Distribution:',labelsDist)

	np.savetxt(outpath, data,fmt='%s')


def crossValidate(inpath, model):
	
	# load data
	print()
	print('Inpath:',inpath)
	print()

	data = np.genfromtxt(inpath,dtype=object)
	data[:,0] = data[:,0].astype(str)
	data[:,1] = data[:,1].astype(int)
	data[:,2:] = data[:,2:].astype(float)
	
	print('# SWCs',len(np.unique(data[:,0])))
	print('# Branches:', len(data))
	print('# Features:', data.shape[1]-2)
	print()
	
	# remove all branches with label of 0
	# these start from soma and we are making assumption to never prune these
	# this may create some false negatives but we are trying to constrain the problem
	rowsLabeled0 = np.where(data[:,1] == 0)[0]
	data = np.delete(data,rowsLabeled0,axis=0)
	print('# Branches after removing label 0:',len(data))
	print()

	# remove all branches with label of 3
	# these are descendant branches of pruned branches
	# we do not care whether these are kept or removed
	# obviously if the target branch isnt pruned, it will matter
	# but this is another assumption to constrain the problem
	rowsLabeled3 = np.where(data[:,1] == 3)[0]
	data = np.delete(data,rowsLabeled3,axis=0)
	print('# Branches after removing label 3:',len(data))
	print()

	# Z score all features within each SWC
	# change all nans to 0
	swcs = np.unique(data[:,0])
	for swc in swcs:
		swcData = data[data[:,0] == swc]
		numericalData = swcData[:,2:].astype(float)
		numericalData[np.isnan(numericalData)] = 0
		normalizedData = zscore(numericalData)
		data[data[:,0] == swc,2:] = normalizedData

	# PCA
	pca = PCA(n_components=2)
	pca.fit(data[:,2:])
	#print(pca.explained_variance_ratio_)
	labels = data[:,1]
	#print(labels)

	transformed = pca.transform(data[:,2:])
	
	fig,ax = plt.subplots()

	colors = {1:'black', 2:'yellow', 4:'yellow', 5:'yellow'}
	for l in np.unique(labels):
		labeledData = transformed[labels==l]
		ax.scatter(labeledData[:,0],labeledData[:,1],s=5, c=colors[l], label=l)
	ax.legend()

	plt.show()

	# train and test a model for each SWC
	swcs = np.unique(data[:,0])

	trueNegativeRates = []
	truePositiveRates = []

	for swc in swcs:

		# split into train and test
		trainingRows = data[data[:,0] != swc ]
		testingRows = data[data[:,0] == swc ]

		X_train = trainingRows[:,2:]
		y_train = trainingRows[:,1]

		X_test = testingRows[:,2:]
		y_test = testingRows[:,1]

		# change labels of all except 1 to 0 (dont remove)
		y_train = (y_train == 1)
		y_test = (y_test == 1)

		if model == 'logistic regression':
			clf = LogisticRegression(solver='lbfgs').fit(X_train, y_train)
		elif model == 'random forest':
			clf = RandomForestClassifier(n_estimators=100).fit(X_train, y_train)
		elif model == 'svm':
			clf = svm.SVC().fit(X_train,y_train)
		elif model == 'knn':
			clf = KNeighborsClassifier(n_neighbors=3).fit(X_train, y_train)

		trainPredictions = clf.predict(X_train)
		testPredictions = clf.predict(X_test)
	
		print(clf.score(X_train,y_train))
		print(clf.score(X_test,y_test))

		P = y_test.sum()
		N = (y_test==0).sum()
		
		fp = sum([1 for y,y_hat in zip(y_test, testPredictions) if y == 0 and y_hat == 1])
		fn = sum([1 for y,y_hat in zip(y_test, testPredictions) if y == 1 and y_hat == 0])
		tp = sum([1 for y,y_hat in zip(y_test, testPredictions) if y == 1 and y_hat == 1])
		tn = sum([1 for y,y_hat in zip(y_test, testPredictions) if y == 0 and y_hat == 0])

		print('# Positive (Training):',y_train.sum())
		print('# Negatives (Testing):',(y_train==0).sum())
		print('# Positive:',P)
		print('# Negatives:',N)
		print('fp:',fp)
		print('fn:',fn)
		print('tp:',tp)
		print('tn:',tn)

		trueNegativeRate = tn/N
		trueNegativeRates.append(trueNegativeRate)

		if P > 0:
			truePositiveRate = tp/P
			truePositiveRates.append(truePositiveRate)

	print('True Negative Rates:',trueNegativeRates)
	print('True Positive Rates:',truePositiveRates)
			



if __name__ == '__main__':


	#outpath = join(BASE_PATH,'data_local_global.txt')
	#compileBranchFeatures(outpath)

	inpath = '/data/elowsky/OLST/swc_analysis/pruning/data_local_global.txt'
	crossValidate(inpath, model='knn')




	
