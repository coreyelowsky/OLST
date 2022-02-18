import os
import sys
import BRAINS

###########################################################
REGISTRATION_BASE_PATH = '/data/elowsky/OLST/registration/'
###########################################################

mop_labels_path = os.path.join(REGISTRATION_BASE_PATH,'MOpul_Layers_528.tif')
rodrigo_registration_path = os.path.join(REGISTRATION_BASE_PATH,'rodrigo_registration_files','MOp_registrateredFiles')

for brain in BRAINS.BRAINS:
	
	print()
	print('Brain:',brain)

	# output brain path
	brain_path = os.path.join(REGISTRATION_BASE_PATH,brain)

	# change brain id fo 170329
	if brain == '170329_500':
		brain = '170329'

	# transform parameter paths for specific brain
	tp_brain_path = os.path.join(rodrigo_registration_path,'elastixOutput_'+brain+'_Emxcre_25um_isotropic_removedStripes_crop')
	tp_0_path = os.path.join(tp_brain_path,'TransformParameters_labels.0.txt')
	tp_1_path = os.path.join(tp_brain_path,'TransformParameters_labels.1.txt')

	# if registration path doesnt exist
	if not os.path.exists(tp_1_path):
		print('REGISTRATION FILES DONT EXIST')
		continue	

	##### change path in transform parameters file #####
	#with open(tp_1_path,'r') as fp:
	#	lines = fp.readlines()
	#	for i,line in enumerate(lines):
	#		if line.startswith('(InitialTransformParametersFileName'):
	#			lines[i] = '(InitialTransformParametersFileName "' + tp_0_path + '")\n'
	#with open(tp_1_path,'w') as fp:
	#	fp.writelines(lines)
	####################################################

	# run transformix to register mop labels
	transformix_command = 'transformix -out ' + brain_path + ' -tp ' + tp_1_path  + ' -in ' + mop_labels_path
	os.system(transformix_command)

	# change name
	os.rename(os.path.join(brain_path,'result.nrrd'),os.path.join(brain_path,'mop_labels_registered.nrrd'))






