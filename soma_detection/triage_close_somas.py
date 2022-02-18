import numpy as np
import csv
import os
from scipy.spatial.distance import cdist

###########################################################################################################################
BRAINS = [170329, 171012, 180206, 180523, 180606, 180614, 180926, 181004 ,181115, 190123, 190306, 190327, 190416, 190522] #
MICRON_THRESH = 200                                                                                                       #
###########################################################################################################################

BRAINS = [170329]

num_cells = 0

for BRAIN_ID in BRAINS:
	
	out_path = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '/'
	soma_path = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '/triaged_somas_duplicates.npy'

	###########################################################
	if not os.path.exists(soma_path):
		continue

	if os.path.exists(out_path + '/triaged_somas_close_' + str(MICRON_THRESH) + '.npy'):
		continue	

	print('Brain: ', BRAIN_ID)
	###########################################################

	somas = np.load(soma_path, allow_pickle=True)
	average_somas = np.empty(shape=(len(somas),3), dtype=object)	

	for i,soma in enumerate(somas):

		num_dup = int(len([x for x in soma if x != None])/ 7)

		new_x = 0
		new_y = 0
		new_z = 0

		for j in range(num_dup):
			new_x += soma[j*7]
			new_y += soma[1 + j*7]
			new_z += soma[2 + j*7]
					
		new_x = int(round(new_x/num_dup*.406))
		new_y = int(round(new_y/num_dup*.406)) 
		new_z = int(round(new_z/num_dup*2.5))
		
		average_somas[i] = np.array([new_x,new_y,new_z])
	
	pairwise_dist = cdist(average_somas,average_somas)
	np.fill_diagonal(pairwise_dist,np.Inf)
	
	remove_indices = np.where(np.any(pairwise_dist < MICRON_THRESH,axis=1))[0]
	
	somas_out = np.delete(somas,remove_indices,axis=0)
	print(len(somas_out))
	num_cells += len(somas_out)

	np.save(out_path + '/triaged_somas_close_' + str(MICRON_THRESH) + '.npy',somas_out)
	np.savetxt(out_path + '/triaged_somas_close_' + str(MICRON_THRESH) + '.csv',somas_out,delimiter=",",fmt="%s")

print('Total:',num_cells)
  
