import csv
from utils import *

####################################################################
# Reads in soma locations (truth) from downloaded google drive csv #
# Gets stitching coordinates for somas and saves soma info         # 
####################################################################

###################
BRAIN_ID = 180206 #
###################

################################################################################
raw_path = '/mnt/brainstore8/anarasim/data/180206_Emxcre_Reconstruction/100pc/'
################################################################################

# load raw truth data
soma_locations_csv = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '_Emxcre_Reconstruction/soma_detection/truth/somas_unparsed.csv'

# Stitching Information
brain_path = '/data/elowsky/OLST/reconstruction/'+ str(BRAIN_ID) +'_Emxcre_Reconstruction/'
brain_full_file = brain_path + 'dataset_raw.xml'
files, registrations, stitchings = extract_stitching_parameters(brain_full_file)


# array to save truth
soma_array_truth = np.zeros((0,7),dtype=object)

with open(soma_locations_csv) as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	volume_z = ''
	volume_y = ''
	for row in readCSV:
		# ignore blank rows
		if all(x == '' for x in row):
			continue
		if row[4] == '-1':
			continue
		
		# parse volume name
		if row[0] != '':
			volume_z = row[0]
		if row[1] != '':
			volume_y = row[1]


		volume = volume_z + '_' + volume_y

		# parse local volume coordinates
		(local_x,local_y,local_z) = (int(row[2]),int(row[3]),int(row[4]))
		if local_x < 0 or local_x >= 2048 or local_y < 0 or local_y >= 1024 or local_z < 0 or local_z >= 4000:
			print('Error: local coordinate out of bounds')
			print(volume)
			print(local_x,local_y,local_z)

		# get stitched coordinates
		volume_full_file = raw_path + volume + '.tif'			
		translation_to_grid_matrix, stitching_matrix = get_stitching_matrices(files,registrations,volume_full_file)
		soma_coord = np.array([local_x,local_y,local_z,1])
		transformed_soma_coord = np.matmul(translation_to_grid_matrix,soma_coord)
		transformed_soma_coord = np.matmul(stitching_matrix,transformed_soma_coord)
		(stitched_x,stitched_y,stitched_z) = (int(transformed_soma_coord[0]),int(transformed_soma_coord[1]),int(transformed_soma_coord[2]))
	
		soma = np.array([int(stitched_x),int(stitched_y),int(stitched_z),volume,int(local_x),int(local_y),int(local_z)],dtype=object)
		soma_array_truth = np.vstack((soma_array_truth,soma))

# save
np.save(brain_path + 'soma_detection/truth/detected_somas_truth.npy',soma_array_truth)
np.savetxt(brain_path + 'soma_detection/truth/detected_somas_truth.csv',soma_array_truth,delimiter=",",fmt="%s")










