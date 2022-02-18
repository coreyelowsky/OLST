from os.path import join
import numpy as np

BRAIN = '190522'
SOMAS_BASE_PATH = '/data/elowsky/OLST/soma_detection/'

mop_somas_path = join(SOMAS_BASE_PATH, BRAIN,'motor_cortex_somas.npy')
mop_somas = np.load(mop_somas_path, allow_pickle=True)

cluster_somas_path = join(SOMAS_BASE_PATH, BRAIN,'detected_somas_clustering.npy')
cluster_somas = np.load(cluster_somas_path, allow_pickle=True)

ids = []

for soma in mop_somas:
	soma_id = np.where(np.all(cluster_somas[:,3:7] == soma[4:8],axis=1))[0][0]+1
	ids.append(soma_id)

ids = np.array(ids).reshape(-1,1)

np.savetxt(join(SOMAS_BASE_PATH,BRAIN,'cluster_ids_of_mop.txt'),ids,fmt=['%d'])



