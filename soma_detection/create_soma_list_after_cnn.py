import numpy as np

###################
brain_id = 170329 #
###################

old_soma_path = '/data/anarasim/data/soma_detections/' + str(brain_id) + '_Emxcre/soma_detection/model/1000/detected_somas_model.npy'
old_somas = np.load(old_soma_path,allow_pickle=True)

output_path = '/data/elowsky/OLST/reconstruction/' + str(brain_id) + '/'

for i in range(2):
	
	if i == 0:

		cnn_predictions_path = '/data/elowsky/OLST/reconstruction/' + str(brain_id) + '/cnn_soma_predictions.npy'
		cnn_predictions = np.load(cnn_predictions_path)
		cnn_predictions = np.array(cnn_predictions, dtype=np.bool)

		new_somas = old_somas[cnn_predictions]
		new_somas = new_somas[:,0:7]

		print(len(cnn_predictions))
		print(sum(cnn_predictions))

		np.save(output_path + 'updated_somas_after_cnn.npy',new_somas)
		np.savetxt(output_path + 'updated_somas_after_cnn.csv',new_somas,delimiter=",",fmt="%s")

	elif i == 1:

		cnn_predictions_path = '/data/elowsky/OLST/reconstruction/' + str(brain_id) + '/cnn_soma_predictions_cnn_2.npy'
		cnn_predictions = np.load(cnn_predictions_path)
		cnn_predictions = np.array(cnn_predictions, dtype=np.bool)

		new_somas = old_somas[cnn_predictions]
		new_somas = new_somas[:,0:7]

		print(len(cnn_predictions))
		print(sum(cnn_predictions))

		np.save(output_path + 'updated_somas_after_cnn_2.npy',new_somas)
		np.savetxt(output_path + 'updated_somas_after_cnn_2.csv',new_somas,delimiter=",",fmt="%s")
