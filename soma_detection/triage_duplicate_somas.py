import numpy as np
import os
import copy

###########################################################################################################################
BRAINS = [170329, 171012, 180206, 180523, 180606, 180614, 180926, 181004 ,181115, 190123, 190306, 190327, 190416, 190522] #
###########################################################################################################################

BRAINS = [170329]

############ Helper Funtion #################
def string_soma_to_array(str_soma):

	soma_array = str_soma[1:-1].split()
	soma_array[0] = int(soma_array[0])
	soma_array[1] = int(soma_array[1])
	soma_array[2] = int(soma_array[2])
	soma_array[3] = soma_array[3][1:-1]
	soma_array[4] = int(soma_array[4])
	soma_array[5] = int(soma_array[5])
	soma_array[6] = int(soma_array[6])

	return np.array(soma_array,dtype=object)
##############################################

RES_FACTORS = np.array([.406,.406,2.5])

for BRAIN_ID in BRAINS:

	soma_path = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '/updated_somas_after_cnn_2.npy'
	somas = np.load(soma_path, allow_pickle=True)


	###########################################################
	if not os.path.exists(soma_path):
		continue

	if os.path.exists('/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) +'/triaged_somas_duplicates.npy'):
		continue	

	print(BRAIN_ID)
	###########################################################


	# For each soma find its closest somas thats in a different volume

	out = np.empty((len(somas),7),dtype=object)

	for i in range(len(somas)):
	
		if i % 100 == 0:
			print(i)	

		soma_a = somas[i]

		soma_a_stitched_coords = soma_a[0:3]
		soma_a_volume = soma_a[3]
		soma_a_volume_coords = soma_a[3:6]

		min_dist = np.Inf
		min_soma_number = -1

		soma_a_z = soma_a_volume[1:3]
		soma_a_y = soma_a_volume[5:7]

		for j in range(len(somas)):

			if i != j:
				soma_b = somas[j]

				soma_b_stitched_coords = soma_b[0:3]
				soma_b_volume = soma_b[3]
				soma_b_volume_coords = soma_b[3:6]

				if soma_a_volume != soma_b_volume:
					dist = np.sqrt(np.sum(np.square((soma_a_stitched_coords - soma_b_stitched_coords)*RES_FACTORS)))

					if dist < min_dist:
						min_dist = dist
						min_soma_number = j+1
						out[i,0] = int(round(min_dist))
						out[i,1] = int(round(abs(soma_a_stitched_coords[0]-soma_b_stitched_coords[0]))*RES_FACTORS[0])
						out[i,2] = int(round(abs(soma_a_stitched_coords[1]-soma_b_stitched_coords[1]))*RES_FACTORS[1])
						out[i,3] = int(round(abs(soma_a_stitched_coords[2]-soma_b_stitched_coords[2]))*RES_FACTORS[2])
						out[i,4] = str(soma_a)
						out[i,5] = str(soma_b)
				
						soma_b_z = soma_b_volume[1:3]
						soma_b_y = soma_b_volume[5:7]

						if soma_a_z != soma_b_z and soma_a_y != soma_b_y:
							out[i,6] = 'zy'
						elif soma_a_z != soma_b_z:
							out[i,6] = 'z'
						elif soma_a_y != soma_b_y:
							out[i,6] = 'y'


	#np.savetxt('/home/elowsky/Desktop/test.csv',out,delimiter=",",fmt="%s")	

	duplicate_array = out
	duplicate_array = duplicate_array[duplicate_array[:,0].argsort()]
	duplicate_array_copy = copy.deepcopy(duplicate_array)
	predictions = [0]*len(duplicate_array)
	# sort duplicate array by distances

	duplicates = []


	while len(duplicate_array) > 0:

		soma_pair = duplicate_array[0]	

		dup_tuple_soma_array = ()
		dup_tuple_duplicate_array = ()	
		direction_tuple = ()

		dist = soma_pair[0]
		dist_x = soma_pair[1]
		dist_y = soma_pair[2]
		dist_z = soma_pair[3]
	
		soma_1 = soma_pair[4]
		soma_1_array = string_soma_to_array(soma_1)
		soma_1_index_in_somas = np.where(np.all(somas == soma_1_array,axis=1))[0][0]
		soma_1_index_in_duplicates = np.where(duplicate_array[:,4] == soma_1)[0][0]

		soma_2 = soma_pair[5]
		soma_2_array = string_soma_to_array(soma_2)
		
		# if soma_2 is no longer in duplicate array
		if not np.any(duplicate_array[:,4] == soma_2):
			duplicate_array = np.delete(duplicate_array,soma_1_index_in_duplicates,0)
			continue

		soma_2_index_in_somas = np.where(np.all(somas == soma_2_array,axis=1))[0][0]
		soma_2_index_in_duplicates = np.where(duplicate_array[:,4] == soma_2)[0][0]

		direction = soma_pair[6]

		if dist <= 100:	
			if duplicate_array[soma_2_index_in_duplicates,5] == soma_1:

					# prediction list update
					predictions[np.where(soma_1 == duplicate_array_copy[:,4])[0][0]] = 1
					predictions[np.where(soma_2 == duplicate_array_copy[:,4])[0][0]] = 1

					# add to duplicate tuple
					dup_tuple_soma_array = dup_tuple_soma_array + (soma_1_index_in_somas,soma_2_index_in_somas)
					dup_tuple_duplicate_array = dup_tuple_duplicate_array + (soma_1_index_in_duplicates,soma_2_index_in_duplicates)
					direction_tuple = direction_tuple + (direction,)

					# find all other somas that are partners with 1 and 2
					indices_of_partners_with_1 = list(np.where(duplicate_array[:,5] == soma_1)[0])
					indices_of_partners_with_2 = list(np.where(duplicate_array[:,5] == soma_2)[0])

					# remove partner indices
					indices_of_partners_with_1.remove(soma_2_index_in_duplicates)
					indices_of_partners_with_2.remove(soma_1_index_in_duplicates)
		
					for ind in indices_of_partners_with_1:
						soma_pair = duplicate_array[ind,:]
						dist = soma_pair[0]
						soma = soma_pair[4]
	
						direction = soma_pair[6]
						soma_array = string_soma_to_array(soma)
						soma_index_in_somas = np.where(np.all(somas == soma_array,axis=1))[0][0]
						if dist < 100:
							if direction not in direction_tuple or ('zy' in direction_tuple and direction_tuple.count(direction)<2) :
								dup_tuple_soma_array = dup_tuple_soma_array + (soma_index_in_somas,)
								dup_tuple_duplicate_array = dup_tuple_duplicate_array + (ind,)
								direction_tuple = direction_tuple + (direction,)
								# prediction list update
								predictions[np.where(soma == duplicate_array_copy[:,4])[0][0]] = 1


					for ind in indices_of_partners_with_2:
						soma_pair = duplicate_array[ind,:]
						dist = soma_pair[0]
						soma = soma_pair[4]

						direction = soma_pair[6]
						soma_array = string_soma_to_array(soma)
						soma_index_in_somas = np.where(np.all(somas == soma_array,axis=1))[0][0]
						if dist < 100:
							if direction not in direction_tuple or ('zy' in direction_tuple and direction_tuple.count(direction)<2) :
								dup_tuple_soma_array = dup_tuple_soma_array + (soma_index_in_somas,)
								dup_tuple_duplicate_array = dup_tuple_duplicate_array + (ind,)
								direction_tuple = direction_tuple + (direction,)
								# prediction list update
								predictions[np.where(soma == duplicate_array_copy[:,4])[0][0]] = 1


					# remove somas that are grouped from duplicates array
					duplicate_array = np.delete(duplicate_array,dup_tuple_duplicate_array,0)

					# add to final duplcate list
					duplicates.append(dup_tuple_soma_array)
			else:
				duplicate_array = np.delete(duplicate_array,soma_1_index_in_duplicates,0)
			

		else:
			# Distance is > 100, delete
			duplicate_array = np.delete(duplicate_array,[soma_1_index_in_duplicates,soma_2_index_in_duplicates],0)

	soma_count = 0
	size_2 = 0
	size_3 = 0
	size_4 = 0

	for dup in duplicates:
		if len(dup) == 2:
			size_2 += 1
		if len(dup) == 3:
			size_3 += 1
		if len(dup) == 4:
			size_4 += 1
		soma_count += len(dup)

	print()
	print('# Somas:',len(somas))
	print()
	print('# Total Duplicate Groups:', len(duplicates))
	print('# Total Somas in Groups:', soma_count)
	print('# 2-Groups:', size_2)
	print('# 3-Groups:', size_3)
	print('# 4-Groups:', size_4)
	print()


	# create new soma list

	somas_output = np.empty(shape=(len(somas),len(somas[0])*4),dtype=object)
	somas_output[:,0:len(somas[0])] = copy.deepcopy(somas)

	for dup in duplicates:
		dup = sorted(dup)

		first_ind = dup.pop(0)
		first_soma = somas[first_ind]
		first_soma_ind_in_output = np.where(np.all(somas_output[:,0:len(first_soma)] == first_soma,axis=1))[0][0]
	
		dup_in_output = []
		for ind in dup:
			soma = somas[ind]
			dup_soma_ind_in_output = np.where(np.all(somas_output[:,0:len(soma)] == soma,axis=1))[0][0]
			dup_in_output.append(dup_soma_ind_in_output)	

			first_soma = np.hstack((first_soma,soma))

		

		somas_output[first_soma_ind_in_output,0:len(first_soma)] = first_soma

		# delete from somas output
		somas_output = np.delete(somas_output,dup_in_output,0)

	print(BRAIN_ID)
	print(len(somas_output))
	np.save('/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) +'/triaged_somas_duplicates.npy',somas_output)
	np.savetxt('/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) +'/triaged_somas_duplicates.csv',somas_output,delimiter=",",fmt="%s")	



