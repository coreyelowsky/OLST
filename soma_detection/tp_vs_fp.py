import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------#
truth_somas_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/truth/detected_somas_truth.npy'
truth_somas = np.load(truth_somas_path,allow_pickle=True)

model_somas_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/2000/detected_somas_model.npy'
model_somas = np.load(model_somas_path,allow_pickle=True)
#---------------------------------------#

#---------------------------------------#
fn_somas_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/2000/fn/fn.npy'
fn_somas = np.load(fn_somas_path,allow_pickle=True)

tp_somas_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/2000/tp/tp.npy'
tp_somas = np.load(tp_somas_path,allow_pickle=True)

fp_somas_path = '/data/elowsky/OLST/reconstruction/180926_Emxcre_Reconstruction/soma_detection/model/2000/fp/fp.npy'
fp_somas = np.load(fp_somas_path,allow_pickle=True)
#---------------------------------------#

#-------------------#
RAW_MEAN_INDEX = 7
RAW_STD_INDEX = 8
RAW_MAX_INDEX = 9
SEG_MEAN_1_INDEX = 10
SEG_STD_1_INDEX = 11
SEG_MEAN_2_INDEX = 12
SEG_STD_2_INDEX = 13
SEG_MAX_INDEX = 14
SEG_COUNT_INDEX = 15
Z_INDEX = 16
Y_INDEX = 17
#--------------------#

INDEX = 12

tp = [a for a in list(tp_somas[:,INDEX]) if a != None]
fp = [a for a in list(fp_somas[:,INDEX]) if a != None]


plt.hist([tp,fp],bins=500,label=['tp','fp'])
plt.legend(loc='upper right')


plt.show()
