import params
from os.path import join, exists
from os import rename, remove, mkdir, walk, listdir
import shutil
import numpy as np
from sys import exit
import csv


# rename files
def rename_files(source, target):

	for brain in params.BRAINS:

		brain_path = join(params.BASE_PATH, brain)
		source_path = join(brain_path, source)
		target_path = join(brain_path, target)
	
		if exists(source_path):
			rename(source_path, target_path)


# delete files
def delete(source):

	for brain in params.BRAINS:

		brain_path = join(params.BASE_PATH, brain)
		source_path = join(brain_path, source)
	
		if exists(source_path):
			remove(source_path)

# make directory
def make_directory(dir_name):

	for brain in params.BRAINS:

		brain_path = join(params.BASE_PATH, brain)
		source_path = join(brain_path,dir_name)
		mkdir(source_path)

# move files
def move(source, target):

	for brain in params.BRAINS:

		brain_path = join(params.BASE_PATH, brain)

		source_path = join(brain_path, source)
		target_path = join(brain_path, target, source)
	
		
		if exists(source_path):
			shutil.move(source_path, target_path)

def barrel_cortex_changes():

	for brain in params.BRAINS:

		source = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex_old_annotation/'
		target = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/'

		brain_path = join(params.BASE_PATH, brain)

		old_somas_path = join(brain_path,'regions', 'barrel_cortex_old_annotation.csv')
		somas_path =  join(brain_path,'regions', 'barrel_cortex.csv')
		
		if not exists(old_somas_path) or not exists(somas_path):
			continue

		print()
		print(brain)	

		old_somas = np.genfromtxt(old_somas_path, delimiter=',',dtype=str)
		somas = np.genfromtxt(somas_path, delimiter=',',dtype=str)
		
		print('# somas old:', len(old_somas))
		print('# somas new:', len(somas))

		old_ids_to_remove = []
		old_ids_to_move = []
		new_ids_to_reconstruct = []
		
		for soma_id, soma in enumerate(somas, start=1):

			# check if its in old exactly
			if (soma == old_somas).all(axis=1).any():
				soma_id_old = np.where((soma == old_somas).all(axis=1))[0][0] + 1
				old_ids_to_move.append((soma_id_old, soma_id))
			elif (soma[1:] == old_somas[:,1:]).all(axis=1).any():
				soma_id_old = np.where((soma[1:] == old_somas[:,1:]).all(axis=1))[0][0] + 1
				old_ids_to_move.append((soma_id_old, soma_id))
			else:
				new_ids_to_reconstruct.append(soma_id)
				

		for i in range(1, len(old_somas)+1):
			if i not in [x[0] for x in old_ids_to_move]:
				old_ids_to_remove.append(i)

		if True:
			for soma_id_tuple in old_ids_to_move:
			
				soma_id_old, soma_id_new = soma_id_tuple

				dirs = [x[0] for x in walk(source)][1:]
				for dir_name in dirs:
				
					files = [x for x in listdir(dir_name) if brain + '_' + str(soma_id_old) + '.' in x or brain + '_' + str(soma_id_old) + '_' in x ]
					for f in files:
						dir_n = dir_name.split('/')[-1]

						if '.swc' in f:
							out_name = brain + '_' + str(soma_id_new) + '.swc'
						elif '.txt' in f:
							out_name = brain + '_' + str(soma_id_new) + '.txt'
						elif 'skeleton' in f:
							out_name = brain + '_'+ str(soma_id_new) + '_dilated_skeleton.tif'
						elif 'raw' in f:
							out_name = brain + '_' +str(soma_id_new) + '_raw.tif'
						else:
							exit('Error')

					
						shutil.copy(join(source,dir_n,f), join(target,dir_n, out_name))
		if False:
			if brain == '180606':
				somas = somas.tolist()

				for soma_tup in old_ids_to_move:
			
					_, new_soma_id = soma_tup
			
					somas[new_soma_id-1].append('RECONSTRUCTED')



			
	
				out_name = join(brain_path,'regions', 'barrel_cortex_for_update.csv')


				with open(out_name, "w") as f:
				    writer = csv.writer(f)
				    writer.writerows(somas)



				print('to remove:',old_ids_to_remove)
				print('to move:', old_ids_to_move)
				print('to reconstruct:', new_ids_to_reconstruct)


		if False:
			info_dict = {'change':old_ids_to_move, 'remove':old_ids_to_remove, 'reconstruct':new_ids_to_reconstruct}
			np.save(join('/data/elowsky/OLST/soma_detection/',brain,'regions','barrel_cortex_update'), info_dict) 


def copy_dirs():

	source = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex_old_annotation/'
	target = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/'

	dirs = [x[0] for x in walk(source)][1:]

	for dir_name in dirs:
		
		mkdir(join(target, dir_name.split('/')[-1]))




if __name__ == '__main__':
	

	#read_dictionary = np.load('/data/elowsky/OLST/soma_detection/170329/regions/barrel_cortex_update.npy',allow_pickle='TRUE').item()
	#print(read_dictionary)
	barrel_cortex_changes()
	#copy_dirs()


	

	


	

	

