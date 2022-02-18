import os
import sys
import BRAINS

REGISTRATION_BASE_PATH = '/data/elowsky/OLST/registration/'

whole_brain_labels_path = os.path.join(REGISTRATION_BASE_PATH,'annotation_025_coronal.nrrd')
rodrigo_registration_path = os.path.join(REGISTRATION_BASE_PATH,'rodrigo_registration_files','MOp_registrateredFiles')

for brain in BRAINS.BRAINS:
	
	print()
	print('Brain:',brain)

	brain_path = os.path.join(REGISTRATION_BASE_PATH,brain)
	if not os.path.exists(brain_path):
		os.mkdir(brain_path)

	if brain == '170329_500':
		brain = '170329'
	tp_brain_path = os.path.join(rodrigo_registration_path,'elastixOutput_'+brain+'_Emxcre_25um_isotropic_removedStripes_crop')
	tp_0_path = os.path.join(tp_brain_path,'TransformParameters_labels.0.txt')
	tp_1_path = os.path.join(tp_brain_path,'TransformParameters_labels.1.txt')

	if not os.path.exists(tp_1_path):
		print('REGISTRATION FILES DONT EXIST')
		continue	

	# change path in transform parameters file
	with open(tp_1_path,'r') as fp:
		lines = fp.readlines()
		for i,line in enumerate(lines):
			if line.startswith('(InitialTransformParametersFileName'):
				lines[i] = '(InitialTransformParametersFileName "' + tp_0_path + '")\n'
	with open(tp_1_path,'w') as fp:
		fp.writelines(lines)

	# run transformix
	transformix_command = 'transformix -out ' + brain_path + ' -tp ' + tp_1_path  + ' -in ' + whole_brain_labels_path
	os.system(transformix_command)

	# change name
	os.rename(os.path.join(brain_path,'result.nrrd'),os.path.join(brain_path,'whole_brain_labels_registered.nrrd'))	






