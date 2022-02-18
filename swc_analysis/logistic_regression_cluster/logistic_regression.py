import sys
import numpy as np

from sklearn.linear_model import LogisticRegression
from os.path import exists

def load_clusters(path):
	
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

def bootstrap_logistic_regression(
		X, 
		y, 			
		bootstrap_samples=1000, 
		c_precision=2, 
		max_iter=1000, 
		feature_index=None):

	# output array		
	coeffs = np.zeros((bootstrap_samples, len(np.unique(y))))
	
	# iterate through bootstrap runs
	for i in range(bootstrap_samples):

		sys.stdout.flush()

		print(f'i: {i}')

		# get features data
		feature_data = X[:,feature_index]

		# randomly sample with replacement
		sampled_data = np.random.choice(feature_data, size=len(feature_data), replace=True)

		# place sampled data in new array
		X_sampled = np.copy(X)
		X_sampled[:,feature_index] = sampled_data

		# find optimal c value
		c = logistic_regression_find_optimal_c(
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
		feature_coeffs = model.coef_[:,feature_index]
		
		# place in output array
		coeffs[i,:] = feature_coeffs

	return coeffs

def logistic_regression_find_optimal_c( 
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



				
if __name__ == '__main__':


	print()
	print("#############################")
	print("Logistic Regression Bootstrap")
	print("#############################")
	print()

	# parse arguments
	feature_index = int(sys.argv[1])
	data_path = sys.argv[2]
	morphometrics_path = sys.argv[3]
	clusters_path = sys.argv[4]
	bootstrap_samples = int(sys.argv[5])
	out_path = sys.argv[6]

	print("Feature Index:", feature_index)
	print("Data Path:", data_path)
	print("Morphometrics Path:", morphometrics_path)
	print("Clusters Path:", clusters_path)
	print("Bootstrap Samples:", bootstrap_samples)
	print("Out Path:", out_path)

	# load clusters
	clusters = load_clusters(clusters_path)

	# load data
	data_matrix = np.genfromtxt(morphometrics_path, delimiter=' ', dtype=object)
	feature_names = data_matrix[0, 2:].astype(str)
	X = data_matrix[1:,2:].astype(float)

	# get clusters	
	y = clusters[:,2].astype(int)

	
	# run boostrap logistic regression
	# on indivudual feature							
	bootstrap_coeffs = bootstrap_logistic_regression(
					X, 
					y, 
					c_precision=2, 
					max_iter=1000, 
					bootstrap_samples=bootstrap_samples,
					feature_index=feature_index)


	# save coeffs
	np.save(out_path,bootstrap_coeffs)

