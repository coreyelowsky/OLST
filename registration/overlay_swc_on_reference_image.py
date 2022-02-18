import os
import numpy as np
from skimage.io import imread
import tifffile as tif

####################################################################################
IMAGE_PATH = '/data/elowsky/OLST/registration/MOpul_Layers_528_10x10x10.tif'
SWC_PATH = '/data/elowsky/OLST/swc_analysis/automatically_traced/flagship/layer_2/registered/171012_13.swc'
OUT_PATH = '/data/elowsky/OLST/swc_analysis/paper/swc_overlays/'
SWC_INTENSITY = 7
SWC_SCALING_FACTOR = 2.5
####################################################################################

# load image
print('Loading Image...')
image = imread(IMAGE_PATH)

# load SWC
print('Loading SWC...')
swc = np.genfromtxt(SWC_PATH)

# scale swc coordinates from 25um to 10um
coords = np.round(swc[:,2:5]*SWC_SCALING_FACTOR).astype(np.uint64)

# overlay swc coords on image
print('Overlay SWC Coordinate...')
image[coords[:,2],coords[:,1],coords[:,0]] = SWC_INTENSITY

# save image
print('Writing Image...')
out_path = os.path.join(OUT_PATH,SWC_PATH.split('/')[-1][:-4]+'.tif')
tif.imwrite(out_path,image)


