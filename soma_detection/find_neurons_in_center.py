import numpy as np

BRAIN_ID = 180606

#####################################################################################################################
path = '/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '/motor_cortex_somas.npy'
somas = np.load(path, allow_pickle=True)
#####################################################################################################################

###########################################
center = np.array([1024,512,2000])
CENTER_THRESHOLD = 200
###########################################

volumes = np.unique(somas[:,4])

output_somas = []

for volume in volumes:

	y = int(volume[-2:])

	volume_array = somas[np.where(somas[:,4] == volume)]

	distances_x = abs(volume_array[:,5] - center[0])
	distances_y = abs(volume_array[:,6] - center[1])

	
	zipped = [(a,b) for a,b in zip(distances_x,distances_y)]

	if np.any([x[0] < CENTER_THRESHOLD and x[1] < CENTER_THRESHOLD for x in zipped]):
		
		poss_somas = volume_array[np.where([x[0] < CENTER_THRESHOLD and x[1] < CENTER_THRESHOLD for x in zipped])]

		for soma in poss_somas:
			output_somas.append(soma)



output_array = np.asarray(output_somas)

np.save('/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '/motor_cortex_somas_center', output_array)
np.savetxt('/data/elowsky/OLST/reconstruction/' + str(BRAIN_ID) + '/motor_cortex_somas_center.csv',output_array,delimiter=",",fmt="%s")

