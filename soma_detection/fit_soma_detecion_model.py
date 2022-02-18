import numpy as np
from sklearn import linear_model
from sklearn import metrics
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from utils import *
from sklearn import linear_model
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier


def fit_model(X,Y,model):
	if model == 'log reg':
		clf = linear_model.LogisticRegression(C=.01,solver='newton-cg').fit(X,Y)
	elif model == 'svm':
		clf = SVC(gamma='scale',C=10).fit(X,Y)
	elif model == 'naive bayes':
		clf = GaussianNB().fit(X,Y)
	elif model == 'knn':
		clf = KNeighborsClassifier(n_neighbors=6).fit(X,Y)
	elif model == 'decision tree':
		clf = tree.DecisionTreeClassifier().fit(X,Y)
	elif model == 'random forests':
		clf = RandomForestClassifier(n_estimators=100).fit(X,Y)
	else:
		print('Error: Invalid Model Name')
	
	return clf


print()

###########################
model ='svm'
#model ='log reg'
#model = 'naive bayes'
#model = 'knn'
#model = 'decision tree'
#model = 'random forests'

print('Model: ', model)
print()
###########################

################
NORMALIZE = True
################

###################
BRAIN_ID = 180926
THRESHOLD = 1000
###################

#########################################################################################################
brain_path = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '_Emxcre_Reconstruction/'
model_path = brain_path + 'soma_detection/model/' + str(THRESHOLD) + '/'
train_test_path = model_path + 'train_test_set/'
#########################################################################################################

########### Load data ##########################
print('Loading training data...')
x_train_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/training/X_train_features.npy'
X_train = np.load(x_train_path,allow_pickle=True)

y_train_path =  '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/training/y_train.npy'
y_train = np.load(y_train_path,allow_pickle=True)
##################################################


# Only need columns with relevant features (0 to 6 are coordinate info)
X_train = X_train[:,7:].astype(np.int32)

if NORMALIZE:
	X_train = (X_train - np.mean(X_train,axis=0)) / np.std(X_train,axis=0)


print('X train: ', X_train.shape)
print('y train: ', y_train.shape)


########################################################
# K FOLD Cross Validation #

K = 10

kf = KFold(n_splits=K)

ppvs = []
npvs = []
tprs = []
tnrs = []
acc = []
fscores = []

if True:
	for train_index, test_index in kf.split(X_train):

		X_fold_train = X_train[train_index]
		y_fold_train = y_train[train_index]

		X_fold_test = X_train[test_index]
		y_fold_test = y_train[test_index]

		# Train Model
		clf = fit_model(X_fold_train,y_fold_train,model)

		# make predictions
		predictions = clf.predict(X_fold_test)

		tp = [1 for a,b in zip(predictions,y_fold_test) if a==1 and b == 1]
		fp = [1 for a,b in zip(predictions,y_fold_test) if a==1 and b == 0]
		fn = [1 for a,b in zip(predictions,y_fold_test) if a==0 and b == 1]
		tn = [1 for a,b in zip(predictions,y_fold_test) if a==0 and b == 0]

		precision = sum(tp) / (sum(tp) + sum(fp))
		recall = sum(tp) / (sum(tp) + sum(fn))
		
		ppvs.append(precision)
		npvs.append(sum(tn) / (sum(tn) + sum(fn)))
		tprs.append(recall)
		tnrs.append(sum(tn) / (sum(fp) + sum(tn)))


		fscores.append(2*(precision * recall)/(precision + recall))

		acc.append(sum(np.equal(predictions,y_fold_test))/len(predictions))
	
	print()
	n_test = len(X_fold_test)
	print('X Fold Test: ', X_fold_test.shape)
	print()

	print('PPV (Mean): ' , int(100*np.mean(ppvs)),'%')
	print('NPV (Mean): ' , int(100*np.mean(npvs)),'%')
	print('TPR (Mean): ' ,  int(100*np.mean(tprs)),'%')
	print('TNR (Mean): ' , int(100*np.mean(tnrs)),'%')

	print('Accuracy (Mean): ', int(100*np.mean(acc)),'%')
	print('F Score (Mean): ', int(100*np.mean(fscores)),'%')



#################################################################################

if False:
	# Train on entire train set
	clf = fit_model(X_train,y_train,model)

	# load test set
	print('Loading training data...')
	x_test_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/training/X_train_features.npy'
	X_test = np.load(x_test_path,allow_pickle=True)
	X_test_orig = X_test

	y_test_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/1000/train_test_set/training/y_train.npy'
	y_test = np.load(y_test_path,allow_pickle=True)

	# normalize
	X_test = X_test[:,7:].astype(np.int32)

	if NORMALIZE:
		X_test = (X_test - np.mean(X_test,axis=0)) / np.std(X_test,axis=0)

	# Predict on test set
	print()
	print('Test Set: ')

	predictions = clf.predict(X_test)

	tp = [1 for a,b in zip(predictions,y_test) if a==1 and b == 1]
	fp = [1 for a,b in zip(predictions,y_test) if a==1 and b == 0]
	fn = [1 for a,b in zip(predictions,y_test) if a==0 and b == 1]
	tn = [1 for a,b in zip(predictions,y_test) if a==0 and b == 0]

	print()
	print('True Positives: ' + str(sum(tp)))
	print('False Positives: ' + str(sum(fp)))
	print('False Negatives: ' + str(sum(fn)))
	print('True Negatives: ' + str(sum(tn)))

	ppv = sum(tp) / (sum(tp) + sum(fp))
	npv = sum(tn) / (sum(tn) + sum(fn))
	tpr = sum(tp) / (sum(tp) + sum(fn))
	tnr = sum(tn) / (sum(fp) + sum(tn))

	f_score = 2*(ppv * tpr)/(ppv + tpr)

	print()
	print('Testing (PPV): ', round(ppv*100),'%')
	print('Testing (NPV): ', round(npv*100),'%')
	print('Testing (TPR): ', round(tpr*100),'%')
	print('Testing (TNR): ', round(tnr*100),'%')
	print('Testing (F Score): ', round(f_score*100),'%')

	print()

	tp = [c for a,b,c in zip(predictions,y_test,X_test_orig) if a==1 and b == 1]
	fp = [c for a,b,c in zip(predictions,y_test,X_test_orig) if a==1 and b == 0]
	fn = [c for a,b,c in zip(predictions,y_test,X_test_orig) if a==0 and b == 1]
	tn = [c for a,b,c in zip(predictions,y_test,X_test_orig) if a==0 and b == 0]


	print()
	print('PPV: given what model calls somas, % that actually are somas')
	print('NPV: given what model calls not somas, % that actually are not somas')
	print('TPR: given actual somas, % that model says are somas')
	print('TNR: given actual not somas, % that model says  are not somas')
	print()

	#tp = np.array(tp)
	#fp = np.array(fp)
	#fn = np.array(fn)
	#tn = np.array(tn)	

	#np.save(train_test_path + 'tp',tp)
	#np.save(train_test_path + 'fp',fp)
	#np.save(train_test_path + 'tn',tn)
	#np.save(train_test_path + 'fn',fn)	
	























