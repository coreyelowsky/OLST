from SWC_Analysis import SWC_Analysis
from os.path import join, exists




def hierarchical_clustering_batch(min_cluster_size=20):


	base_path = '/data/elowsky/OLST/swc_analysis/analysis/mop/pia_white_matter_normalized/analysis/removed_nan_nodes/scaled_1000/soma_centered/'
	data_types = ['morphometrics']
	layers = ['6']
	normalize_types = ['pca']
	structure_types = ['basal_apical']

	for data_type in data_types:
		for layer in layers:
			for normalize_type in normalize_types:
				for structure_type in structure_types:
				
					data_path = join(base_path, data_type,'layer_' + layer, normalize_type, structure_type + '.txt')
					out_path = join(base_path, data_type,'layer_' + layer,'hierarchical_clustering', normalize_type, structure_type + '.png')

					if exists(data_path):
						SWC_Analysis.hierarchical_clustering(data_path, min_cluster_size=min_cluster_size, plot=True, save=False, out_path=out_path)
			
	
	


if __name__ == '__main__':

	min_cluster_size = 3
	
	hierarchical_clustering_batch(min_cluster_size)
	

