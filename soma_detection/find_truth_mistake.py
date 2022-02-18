import os
import numpy as np
import tifffile as tif
from timeit import default_timer as timer
from utils import *

path = '/data/anarasim/data/180926_Emxcre_Reconstruction/100pc/'

x_c = [1270]

y_c = [296]

z_c = [3298]


SOMA_RADIUS = 20 #microns
RES_FACTORS = np.array([0.406,0.406,2.5])
x_pix_radius = int(SOMA_RADIUS/RES_FACTORS[0])
y_pix_radius = int(SOMA_RADIUS/RES_FACTORS[1])
z_pix_radius = int(SOMA_RADIUS/RES_FACTORS[2])

data = np.zeros((0,4))
for f in os.listdir(path):

	# load volume
	print("Loading Volume: " + path + f)
	load_start = timer()
	raw_volume = tif_to_numpyArray(path + f)
	load_end = timer()
	print(str(round(load_end - load_start)) + " seconds")

	for i in range(len(x_c)):
		x = int(x_c[i])
		y = int(y_c[i])
		z = int(z_c[i])

		if z-z_pix_radius < 0:
			min_z = 0
		else:
			min_z = z-z_pix_radius

		if y-y_pix_radius < 0:
			min_y = 0
		else:
			min_y = y-y_pix_radius

		if x-x_pix_radius < 0:
			min_x = 0
		else:
			min_x = x-x_pix_radius
	
		cropped_volume = raw_volume[min_z:z+z_pix_radius,min_y:y+y_pix_radius,min_x:x+x_pix_radius]
		mean = np.mean(cropped_volume)
		variance = np.var(cropped_volume)
		max_intensity = np.amax(cropped_volume)
		a = np.array([f,int(mean),int(variance),int(max_intensity)])
		data = np.vstack((data,a))

np.savetxt('/home/elowsky/Desktop/find_prob_z49_y10.csv',data,delimiter=",",fmt="%s")
		

