import numpy as np
import os
from utils import *
from timeit import default_timer as timer

# Clustering Algorithm to Detect Somas in Volume

print()
print('##########################')
print('DETECT SOMAS (CLUSTERING)')
print('##########################')
print()

###################
INTENSITY_THRESH = 1000
RADIUS_THRESH = 30 # microns
BRAIN = 190327
###################

print('Intensity Threshold:',INTENSITY_THRESH)
print('Radius Threshold:',RADIUS_THRESH,'(um)')
print()

####################################################################
SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
####################################################################

# initilialize paths
brain_path = os.path.join(SOMA_BASE_PATH,str(BRAIN))
xml_full_path = os.path.join(brain_path,'stitching_parameters.xml')
output_path = os.path.join(SOMA_BASE_PATH,str(BRAIN))

# initialize output arrays
global_soma_centroids_detailed = np.empty(shape=(0,7))
global_soma_centroids = np.empty(shape=(0,4))

# Extract information from xml
files, registrations, stitchings = extract_stitching_parameters(xml_full_path)

# iterate through all volumes
for file_dict in files:

	volume_full_path = file_dict['file name']
	volume_id = volume_full_path[-11:-4]

	if os.path.exists(volume_full_path):
		# load volume
		print("Loading Volume:", volume_full_path)
		load_start = timer()
		raw_volume = tif_to_numpyArray(volume_full_path)
		load_end = timer()
		print(str(round(load_end - load_start)) + " seconds")
		print()
	else:
		print('VOLUME DOES NOT EXIST:',volume_full_path)
		continue

	z_size, y_size, x_size = raw_volume.shape

	# get coordinates of all voxels greter than threshold
	print('Get all voxels > threshold...')
	start = timer()
	soma_coordinates = np.where(raw_volume > INTENSITY_THRESH)
	end = timer()
	print(str(round(end - start)) + " seconds")
	print()
	
	# sort voxels in descending order of intensity
	print('Sort Voxels...')
	start = timer()
	intensities = raw_volume[soma_coordinates]
	order = np.argsort(-intensities)
	(soma_x, soma_y, soma_z) = (soma_coordinates[2][order],soma_coordinates[1][order],soma_coordinates[0][order])
	intensities = intensities[order]
	print(str(round(end - start)) + " seconds")
	print()

	# find centroids of all somas
	# each element in list is [(np array of centroid), num voxels, average intensity, max intensity]
	soma_centroids = np.empty(shape=(0,6))

	print('Cluster Somas...')
	start = timer()
	# iterate through all soma voxels in volume
	for x,y,z,intensity in zip(soma_x,soma_y,soma_z,intensities):
		# if there are no centroids (first soma voxel become centroid)
		if len(soma_centroids) == 0:
			centroid = [x,y,z,1,intensity,intensity]
			soma_centroids = np.vstack((soma_centroids,centroid))
		else:

			# find centroid with closest distance
			min_dist = np.Inf
			min_index = -1
			for i,centroid in enumerate(soma_centroids):
				coord = centroid[0:3]
				dist = distance((x,y,z),coord,'voxels','microns')
				if dist < min_dist:
					min_dist = dist
					min_index = i
			
			# if not close enough to any centroids, make new centroid
			# else update centroid
			if min_dist > RADIUS_THRESH:
				centroid = [x,y,z,1,intensity,intensity]
				soma_centroids = np.vstack((soma_centroids,centroid))
			else:	
				centroid_coord = soma_centroids[min_index][0:3]			
				centroid_count = soma_centroids[min_index][3]
				centroid_average_intensity = soma_centroids[min_index][4]

				soma_centroids[min_index][0:3] = (np.asarray(centroid_coord)*centroid_count+np.asarray((x,y,z)))/(centroid_count+1)
				soma_centroids[min_index][3] += 1
				soma_centroids[min_index][4] = (centroid_average_intensity*centroid_count + intensity) / (centroid_count + 1)
	
	soma_centroids = np.round(soma_centroids,1).astype(np.int32)
	end = timer()
	print(str(round(end - start)) + " seconds")
	print()

	# concat soma centroids for volume
	soma_centroids_detailed = np.hstack((np.full((len(soma_centroids),1),volume_id,dtype=object),soma_centroids) )
	soma_centroids = soma_centroids_detailed[:,0:4]

	# concat soma centroids for brain
	global_soma_centroids_detailed = np.vstack((global_soma_centroids_detailed,soma_centroids_detailed))
	global_soma_centroids = np.vstack((global_soma_centroids,soma_centroids))

	print("# Somas Detected (Volume):",len(soma_centroids))
	print("# Somas Detected (Whole Brain):",len(global_soma_centroids))
	print()


	if len(soma_centroids) > 0:
		np.save(os.path.join(output_path,'detected_somas_clustering'),global_soma_centroids)
		np.save(os.path.join(output_path,'detected_somas_clustering_detailed'),global_soma_centroids_detailed)
		np.savetxt(os.path.join(output_path,'detected_somas_clustering_detailed.csv'),global_soma_centroids_detailed,delimiter=",",fmt="%s")
		np.savetxt(os.path.join(output_path,'detected_somas_clustering.csv'),global_soma_centroids,delimiter=",",fmt="%s")






