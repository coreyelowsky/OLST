import sys
import numpy as np

from os.path import join
from os import remove

print()
print("###############################")
print("Merge Logistic Regression Coeffs")
print("###############################")
print()

# parse arguments
data_path = sys.argv[1]
num_features = int(sys.argv[2])
bootstrap_samples = int(sys.argv[3])

print("Data Path:", data_path)
print("# Features:", num_features)
print("Bootstrap Samples:", bootstrap_samples)

# iterate through features
for i in range(num_features):

	feature_coeffs = np.load(join(data_path, f'coeffs_{i}.npy'))
	
	if i == 0:
		coeffs = np.zeros((num_features, bootstrap_samples, feature_coeffs.shape[1]))
		coeffs[i,:,:] = feature_coeffs
	else:
		coeffs[i,:,:] = feature_coeffs

	# delete feature coeffs file
	remove(join(data_path, f'coeffs_{i}.npy'))

np.save(join(data_path, 'bootstrap_coeffs.npy'), coeffs)
		

