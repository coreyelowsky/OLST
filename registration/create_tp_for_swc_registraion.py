from BRAINS import *
import os
from shutil import copyfile

REGISTRATION_PATH = '/data/elowsky/OLST/registration/'
RODRIGO_REGISTRATION_PATH = '/data/elowsky/OLST/registration/rodrigo_registration_files/MOp_registrateredFiles'

for brain in BRAINS:
	
	if brain == '190306':
		continue
	print(brain)

	if brain == '170329_500':
		brain_short = '170329'
	else:
		brain_short = brain
	rod_reg_brain_path = os.path.join(RODRIGO_REGISTRATION_PATH,'elastixOutput_' + brain_short + '_Emxcre_25um_isotropic_removedStripes_crop')
	tp_0_path = os.path.join(rod_reg_brain_path,'TransformParameters_labels.0.txt')
	tp_1_path = os.path.join(rod_reg_brain_path,'TransformParameters_labels.1.txt')

	reg_brain_path = os.path.join(REGISTRATION_PATH,brain)
	tp_0_path_out = os.path.join(reg_brain_path,'TransformParameters_labels.0.txt')
	tp_1_path_out = os.path.join(reg_brain_path,'TransformParameters_labels.1.txt')
	tp_0_path_out_fs = os.path.join(reg_brain_path,'TransformParameters_labels.0_full_scale.txt')
	tp_1_path_out_fs = os.path.join(reg_brain_path,'TransformParameters_labels.1_full_scale.txt')

	# copy tp files
	copyfile(tp_0_path,tp_0_path_out)
	copyfile(tp_1_path,tp_1_path_out)
	copyfile(tp_0_path,tp_0_path_out_fs)
	copyfile(tp_1_path,tp_1_path_out_fs)

	# change path and spacing in tp files
	spacing = '.01624 .1 .01624'	
	
	# spacing for 0
	command = 'sed -i \'s/.*(Spacing.*/(Spacing ' + spacing  +')/g\' ' + tp_0_path_out_fs
	os.system(command)

	# spacing for 1
	command = 'sed -i \'s/.*(Spacing.*/(Spacing ' + spacing  +')/g\' ' + tp_1_path_out_fs
	os.system(command)
	
	# path for 1
	tp_0_path_out_fs = tp_0_path_out_fs.replace('/','\/')
	command = 'sed -i \'s/.*(InitialTransformParametersFileName.*/(InitialTransformParametersFileName "' + tp_0_path_out_fs + '")/g\' '+ tp_1_path_out_fs
	os.system(command)

	
	


	
