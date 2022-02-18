library(nat)
library(stringr)

LAYER <- 5
BASE_PATH <- "/data/elowsky/OLST/swc_analysis/automatically_traced/flagship"

layer_path = file.path(BASE_PATH,paste('layer_',LAYER,sep=''))

out_path = file.path(layer_path,'swc_images')

path_ones = file.path(layer_path,'normalized_oblique_1')
path_twos = file.path(layer_path,'normalized_oblique_2_branches_removed')

swcs_ones = list.files(path_ones,full.names=TRUE)
swcs_twos = list.files(path_twos,full.names=TRUE)

swcs = c(swcs_ones,swcs_twos)

open3d()
for (swc in swcs){

	neuron = read.neuron(swc)
	plot3d(neuron,soma=40,color='black')

	swc_id = str_replace(tail(unlist(strsplit(swc,'/')),n=1),".swc",".pdf")

	rgl.postscript(file.path(out_path,swc_id),"pdf")

	clear3d()

}


