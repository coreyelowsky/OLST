import numpy as np
import csv
import matplotlib.pyplot as plt
import copy


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



soma_path = '/data/elowsky/OLST/reconstruction/180926/updated_somas_after_cnn_2.npy'
somas = np.load(soma_path ,allow_pickle=True)

# read csv
csv_path = '/data/elowsky/OLST/reconstruction/triage/triage_duplicates_truth_180926.csv'

duplicate_array = np.empty(shape=(len(somas),8),dtype=object)

with open(csv_path) as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	for i,row in enumerate(readCSV):

		dist = int(row[0])
		(dist_x, dist_y, dist_z) = (int(row[1]),int(row[2]),int(row[3]))
		soma_1 = row[4]
		soma_2 = row[5]
		direction = row[6]
		duplicate = int(row[7])

		duplicate_array[i,:] = np.array([dist,dist_x,dist_y,dist_z,soma_1,soma_2,direction,duplicate],dtype=object)

duplicate_array_copy = copy.deepcopy(duplicate_array)

duplicates = []

predictions = [0]*len(duplicate_array)
truth = list(duplicate_array[:,-1])

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

			print('not partners')

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



tp = sum([1 for x in zip(predictions,truth) if x[0] == 1 and x[1]==1])
tn = sum([1 for x in zip(predictions,truth) if x[0] == 0 and x[1]==0])
fp = sum([1 for x in zip(predictions,truth) if x[0] == 1 and x[1]==0])
fn = sum([1 for x in zip(predictions,truth) if x[0] == 0 and x[1]==1])
acc = round(sum([1 for x in zip(predictions,truth) if x[0] == x[1]])/len(truth)*100,1)

print('Accuracy:',acc,'%')
print('# tp:',tp)
print('# tn:',tn)
print('# fp:',fp)
print('# fn:',fn)

#for i,x in enumerate(zip(predictions,truth)):
#	if x[0] == 0 and x[1] == 1:
#		print(i+1)
	






		


		



		

