from utils import *
import matplotlib.pyplot as plt

# Given truth soma location and those from a model
# finds true positives, false positive, true negative

###############################################################################################
BRAIN_ID = 180206
brain_path = '/data/elowsky/OLST/reconstruction/'+ str(BRAIN_ID) +'_Emxcre_Reconstruction/'
out_path = brain_path + 'soma_detection/model/1000/'
################################################################################################

#############################################
MICRON_THRESH = -1
PIXEL_THRESH = 50
RES_FACTORS = np.array([0.406,0.406,2.5])
DIST_THRESH = PIXEL_THRESH
##############################################

########################################################################################################
# load truth somas
soma_array_truth = np.load(brain_path + 'soma_detection/truth/detected_somas_truth.npy',allow_pickle=True)
num_somas_truth = len(soma_array_truth)
print('# Somas (truth): ' + str(num_somas_truth))
####################################################################################################

###########################################################################################################
# load model somas
soma_array_model = np.load(brain_path + 'soma_detection/model/1000/detected_somas_model.npy',allow_pickle=True)
num_somas_model = len(soma_array_model)
print('# Somas (model): ' + str(num_somas_model))
###################################################################################################


# index 8: soma index of matched soma in truth
# index 9: microns from matched soma in truth
tp_array = np.zeros((0,9),dtype=object)

fp_array = np.zeros((0,7),dtype=object)

# index 8: soma index of closest to in false positives
# index 9: microns from closest to in false positives
fn_array = np.zeros((0,9),dtype=object)


volumes_truth = np.unique(soma_array_truth[:,3])

# loop through all volumes in truth
for i in range(len(volumes_truth)):
	volume = volumes_truth[i]

	z_in = int(volume[1:3])
	y_in = int(volume[5:7])

	if z_in == 9 or z_in == 13 or z_in == 17 or z_in == 20  or z_in == 25 or z_in == 30 or z_in == 34 or z_in == 38:
	


		# only somas from specific volume
		soma_array_truth_volume = soma_array_truth[soma_array_truth[:,3] == volume]
		soma_array_model_volume = soma_array_model[soma_array_model[:,3] == volume]

		# matrix to hold pairwise distances
		dist_matrix = np.full((len(soma_array_truth_volume),len(soma_array_model_volume)), np.Inf)
	
		if len(soma_array_model_volume) == 0:
			# false negatives
			fn_array = np.vstack((fn_array,np.hstack((soma_array_truth_volume[:,0:7],np.full((len(soma_array_truth_volume),2),None)))))
		else:
			# popluate distances matrix for volume
			for j in range(len(soma_array_truth_volume)):
				soma_truth = soma_array_truth_volume[j]
				(truth_stitched,truth_local) = (soma_truth[0:3],soma_truth[4:7])
	
				for k in range(len(soma_array_model_volume)):
					soma_model = soma_array_model_volume[k]
					(model_stitched,model_local) = (soma_model[0:3],soma_model[4:7])
			
					# MICRON DISTANCE
					#dist_matrix[j,k] = np.sqrt(np.sum(np.square(np.multiply(model_stitched-truth_stitched,RES_FACTORS))))
					# PIXEL DISTANCE
					dist_matrix[j,k] = np.sqrt(np.sum(np.square(model_stitched-truth_stitched)))

			# arrays to keep track of true indicies in volume
			truth_index_list = list(range(0,len(soma_array_truth_volume))) 
			model_index_list = list(range(0,len(soma_array_model_volume))) 		

			# algorithm to parse through distance matrix
			while dist_matrix.size > 0:

				# find index of minimum distance in matrix
				(truth_index, model_index) = np.unravel_index(np.argmin(dist_matrix,axis=None),dist_matrix.shape)
				(soma_truth, soma_model) = (soma_array_truth_volume[truth_index_list[truth_index],:],soma_array_model_volume[model_index_list[model_index],:])
				min_dist = dist_matrix[truth_index,model_index]

				if min_dist < DIST_THRESH:
					# true positive
					index_in_truth = np.where((soma_truth == soma_array_truth).all(axis=1))[0][0] + 1
					tp_array = np.vstack((tp_array,np.hstack((soma_model[0:7],index_in_truth,int(round(min_dist))))))
				
					# delete col
					dist_matrix = np.delete(dist_matrix,model_index,1)
				
					# remove from model_array
					t_index = np.where((soma_model == soma_array_model).all(axis=1))[0][0]
					soma_array_model = np.delete(soma_array_model,t_index,0)
					model_index_list.remove(model_index_list[model_index])
				else:
					# false negative
					fn_array = np.vstack((fn_array,np.hstack((soma_truth[0:7],np.full(2,None)))))

				# delete row
				dist_matrix = np.delete(dist_matrix,truth_index,0)
				truth_index_list.remove(truth_index_list[truth_index])
		
				# if truth somas left, but none in model, all left must be false negatives
				if dist_matrix.shape[0] > 0 and dist_matrix.shape[1] == 0:		
					for ind in range(len(dist_matrix)):
						local_soma = soma_array_truth_volume[truth_index_list[ind],:]
						fn_array = np.vstack((fn_array,np.hstack((local_soma[0:7],np.full(2,None)))))


fp_array = soma_array_model[:,0:7]

# since not all volumes were marked up
fp_array_new = np.zeros((0,7),dtype=object)
for row in fp_array:
	z = row[3][0:3]

	if z == 'Z09' or z == 'Z13' or z == 'Z17' or z == 'Z20' or z == 'Z25' or z == 'Z30' or z == 'Z34' or z == 'Z38':
		fp_array_new = np.vstack((fp_array_new,row))

fp_array = fp_array_new
#########################

print()
print('True Positives: ' + str(len(tp_array)))	
print('False Positives: ' + str(len(fp_array)))	
print('False Negatives: ' + str(len(fn_array)))
print()


np.save(out_path +'tp.npy',tp_array[:,0:7])
np.save(out_path +'fp.npy',fp_array[:,0:7])
np.savetxt(out_path + 'tp.csv',tp_array[:,0:7],delimiter=",",fmt="%s")
np.savetxt(out_path +'fp.csv',fp_array[:,0:7],delimiter=",",fmt="%s")

# for all false negatives, find closest in false positive and plot histogram
if True:
	dists = []
	fp_indices = []
	fp_somas = []
	fn_somas = []
	volumes = np.unique(fn_array[:,3])
	for i in range(len(volumes)):
		volume = volumes[i]

		# only somas from specific volume
		fn_volume = fn_array[fn_array[:,3] == volume]
		fp_volume = fp_array[fp_array[:,3] == volume]
	
		for j in range(len(fn_volume)):
			fn_soma = fn_volume[j] 
			fn_coords = fn_soma[4:7]
			(min_dist,min_fp_soma) = (np.Inf,None)
			for k in range(len(fp_volume)):
				fp_soma = fp_volume[k]
				fp_coords = fp_soma[4:7]
				#MICRONS
				#dist = np.sqrt(np.sum(np.square(np.multiply(fn_coords-fp_coords,RES_FACTORS))))
				#PIXELS
				dist = np.sqrt(np.sum(np.square(fn_coords-fp_coords)))
				if dist < min_dist:
					min_dist = dist
					min_fp_soma = fp_soma
			if min_dist != np.Inf:
				ind_fn = np.where((fn_soma == fn_array).all(axis=1))[0][0]
				ind_fp = np.where((min_fp_soma == fp_array).all(axis=1))[0][0]
				fn_array[ind_fn,7] = ind_fp
				fn_array[ind_fn,8] = int(round(min_dist))


np.save(out_path +'fn.npy',fn_array)
np.savetxt(out_path +'fn.csv',fn_array,delimiter=",",fmt="%s")	


plt.hist([a for a in list(fn_array[:,8]) if a != None],bins=200)
plt.show()





