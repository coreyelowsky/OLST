library(RcppCNPy)
library(dynamicTreeCut)
library(RColorBrewer)
library(dendextend)
library(ggdendro)
library(ggplot2)


print("")
print("########################")
print("Hierarchical Clustering")
print("########################")
print("")

# base path of data
base_path = "/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/analysis/removed_nan_nodes/scaled_1000/soma_centered/"

# whether to show clusters below plot
CLUSTERS = TRUE

#data_types = c("morphometrics", "arbor_density", "persistent_homology")
data_types = c("morphometrics")
#layers = c("2", "5", "6", "all")
layers = c("all")
#structure_types = c("basal", "apical", "basal_apical", "basal_apical_concat")
structure_types = c("basal_apical")
#norm_types = c("normalized", "pca")
norm_types = c("normalized")

colors_layer = brewer.pal(3, "Spectral" )

# iterate through feature representations
for(data_type in data_types){

	print("")
	print(paste("Feature Type: ", data_type))
	print("")

	data_type_path = paste(base_path, data_type,"/",  sep="")

	# iterate through lauyers
	for (layer in layers){

		print("")
		print(paste("Layer: ", layer, sep=""))
		print("")

		layer_path = paste(data_type_path, "layer_", layer,"/",  sep="")
		out_path_layer = paste(layer_path, "hierarchical_clustering/", sep="")

		# create dir
		dir.create(out_path_layer, showWarnings=FALSE)

		# iterate through normalization types
		for (norm_type in norm_types){

			print("")
			print(paste("Normalization Type: ",norm_type, sep=""))
			print("")

			norm_path_out = file.path(out_path_layer, norm_type)
			dir.create(norm_path_out, showWarnings=FALSE)
		
			# iterate through structure types
			for (structure in structure_types){

				print("")
				print(paste("Structure: ", structure, sep=""))
				print("")

				data_path = file.path(layer_path,norm_type, paste(structure, ".txt", sep=""))
	
				# make sure data path exists
				if(file.exists(data_path)){

					# load data
					data = read.table(data_path, header=TRUE, sep = " ")

					# print number of swcs
					print(paste("# SWCs: ",nrow(data),sep=""))

					# load layer info
					layer_info = data[,2]

					# num data
					num_data = data[,3:ncol(data)]

					# get pairwise distances
					distances = dist(num_data, method="euclidean")

					# create dendrogram
					dendro = hclust(distances, method="ward.D2")

					# cut dendrogram
					clusters = cutreeHybrid(dendro, distM=as.matrix(distances))
					cluster_labels = clusters$labels
					
					print(cluster_labels)
		
			


					## function to set label color
					labelCol <- function(x) {
						if (is.leaf(x)) {

							label <- attr(x, "label") 
			
							if (layer_info[label] == 2){
								color = "blue"
							} else if(layer_info[label] == 5){
								color = "red"
							} else if(layer_info[label] == 6){
								color = "green"

							}
							attr(x, "nodePar") <- list(lab.col=color)

						}
						return(x)
					}

					## apply color to number labels
					dendro <- dendrapply(as.dendrogram(dendro), labelCol)
	
					colors = cbind(cluster_labels, layer_info)

					# iterate through layers and set color
					for(j in 1:nrow(colors)){
						if(colors[j,2] == 2){
							colors[j,2] = colors_layer[1]
						}else if(colors[j,2] == 5){
							colors[j,2] = colors_layer[2]
						}else if(colors[j,2] == 6){
							colors[j,2] = colors_layer[3]
	
						}
		
					}



					if(CLUSTERS){
						out_path = "test"
					}else{
						out_path = file.path(norm_path_out, paste(structure, ".png", sep=""))
					}
			

					#png(file=out_path, width=1800, height=1500)
					plot(dendro, hang=-1, axes=F, xlab="", ylab="", sub="", main="")



					if(CLUSTERS){
						colored_bars(colors=colors, dend=dendro, rowLabels=c("cluster", "layer"), cex.rowLabels=1, text_shift=0, palette="Spectral")
					}

					#dev.off()
				}
			}
}
	}


}








