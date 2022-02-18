# Paths
LABEL_MARKUP_PATH = '/data/anarasim/data/mop/'
INITIAL_SWC_BASE_PATH = '/data/palmer/data/reconstructed_brains/'
SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
SWC_BASE_PATH = '/data/elowsky/OLST/swc_analysis/automatically_traced/'
REGISTRATION_BASE_PATH = '/data/elowsky/OLST/registration/'
SWC_CROPPED_IMAGE_PATH = '/data/palmer/data/reconstructed_brains/'
ARUN_BASE_PATH = '/data/anarasim/data/'
JASON_SWC_PATH = '/data/palmer/data/swc_files/'
IMAGEJ_PATH = '/data/elowsky/applications/Fiji.app/ImageJ-linux64'
OBLIQUE_TO_CORONAL_MACRO_PATH = '/data/elowsky/OLST/swc_analysis/oblique_to_coronal_for_label_markup.ijm'

# Oblique Full Scale Resolution
X_RES_OBLIQUE = .406
Y_RES_OBLIQUE  = .406
Z_RES_OBLIQUE  = 2.5
RES_OBLIQUE = [X_RES_OBLIQUE,Y_RES_OBLIQUE,Z_RES_OBLIQUE]

# Coronal Full Scale Resolution
X_RES_CORONAL = .406
Y_RES_CORONAL = 2.5
Z_RES_CORONAL = .406
RES_CORONAL = [X_RES_CORONAL,Y_RES_CORONAL,Z_RES_CORONAL]

# Anisotropy factor
ANISOTROPY_FACTOR = Z_RES_OBLIQUE/X_RES_OBLIQUE

# Shearing Factor for Oblique -> Coronal
SHEAR_FACTOR_ISOTROPIC = -.7 
SHEAR_FACTOR_ANISOTROPIC = SHEAR_FACTOR_ISOTROPIC / ANISOTROPY_FACTOR

# All Brains
BRAINS = ['170329_500','171012','180206','180523','180606','180614','180926','181004','181115','190123','190306','190327','190416','190522']

# Volume Sizes
X_SIZE = 2048
Y_SIZE = 1024
Z_SIZE = {'170329_500':4000,'171012':4000, '180206':4000, '180523':4000, '180606':4000, '180614':4000, '180926':4000, '181004':4000, '190306':4000, '190327':4000, '190416':4000, '190522':4000,'181115':3600, '190123':3600}

# SWC Data Type Formats
SWC_FORMAT_INT = ['%d','%d','%d','%d','%d','%d','%d']
SWC_FORMAT_FLOAT = ['%d','%d','%g','%g','%g','%d','%d']

# SWC Indices
SWC_INDICES = {'id':0 , 'sid':1 , 'x':2 , 'y':3 ,'z':4 ,'radius':5 ,'pid':6 }

# Soma Parent ID
SOMA_PID = -1

# SWC delimiter
SWC_DELIMITER = ' '

# Structure ID Dictionary
SID_DICT = {'soma': 1, 'axon': 2, 'basal dendrite': 3, 'apical dendrite': 4}

# Structure Color Dictionary
SID_COLORS = {0: 'green',1: 'black', 2: 'yellow', 3: 'red', 4: 'blue',5:'purple'}
SID_PLOTTING_RADII = {0:1,1:100,2:1,3:1,4:1,5:1}

# soma collapse radius
SOMA_COLLAPSE_RADIUS = 30


METRIC_DISPLAY_DICT = {'Width':'Width','Terminal_degree':'Terminal Degree', 'TerminalSegment':'Terminal Segment','Partition_asymmetry':'Partition Asymmetry','Depth':'Depth','Contraction':'Contraction','Branch_Order':'Branch Order','Bif_tilt_remote':'Bifurcation Tilt Remote','Bif_ampl_remote':'Bifurcation Amplitude Remote','PathDistance':'Path Distance', 'Branch_pathlength':'Branch Path Length', 'N_stems':'# Stems', 'Length':'Length', 'EucDistance':'EucDistance', 'Height':'Height'}

# Smoothing Parameters
SMOOTHING_DIFF_THRESHOLD = .075
SMOOTHING_ITERS_THRESHOLD = 10

# Reigstration Parameters
CROP_INFO_PATH = '/data/elowsky/OLST/registration/crop_info.txt'
CROP_FACTOR = 10
REGISTRATION_RES_MICRONS = 25
MOP_REFERENCE_PATHS = {10:'/data/elowsky/OLST/registration/MOpul_Layers_528_10x10x10.tif',25:'/data/elowsky/OLST/registration/MOpul_Layers_528.tif'}
OVERLAY_INTENSITY = 100
CCF_TIF_SHAPE = [13200, 8000, 11400]

############
# Analysis #
############

SWC_STRUCTURE_TYPES = ['basal', 'apical', 'basal_apical']
SWC_STRUCTURE_DISPLAY = {'basal':'Basal', 'apical':'Apical', 'basal_apical':'Basal + Apical', 'basal_apical_concat':'Basal + Apical (Concatenated Features)'}

LAYERS = [2,5,6]
LAYERS_DISPLAY = {2:'Layer 2/3', 5:'Layer 5', 6:'Layer 6'}

BOXPLOT_PSHIFT = {'':0,'*':.025, '**':.04, '***':.06, '****':.082, '*****':.095}

BOX_PLOT_METRIC_DISPLAY_DICT = {'Width':{'Average':'Average Width (um)'},'Terminal_degree':{'Maximum':'Maximum Terminal Degree'}, 'TerminalSegment':{'Total_Sum':'Terminal Segment'},'Partition_asymmetry':{'Average':'Average Partition Asymmetry'},'Depth':{'Average':'Average Depth (um)'},'Contraction':{'Average':'Average Contraction'},'Branch_Order':{'Maximum':'Maximum Branch Order'},'Bif_tilt_remote':{'Average':' Average Bifurcation Tilt Remote'},'Bif_ampl_remote':{'Average':'Average Bifurcation Amplitude Remote'},'PathDistance':{'Maximum':'Maximum Path Distance (um)', 'Average':'Average Path Distance (um'}, 'Branch_pathlength':{'Average':'Average Branch Path Length (um)'}, 'N_stems':{'Total_Sum':'# Stems'}, 'Length':{'Total_Sum':'Total Length (um)'}, 'EucDistance':{'Maximum':'Maximum Euclidean Distance (um)'}, 'Height':{'Average':'Average Height (um)', 'Maximum':'Maximum Height (um)'},'N_bifs':{'Total_Sum':'# Bifurcations'},'N_branch':{'Total_Sum':'# Branches'}}

# LMeasure
LMEASURE_EXE = '/data/elowsky/applications/Lmv5.3_64bit/lmeasure'
LMEASURE_FUNCTION_IDS  = {'N_stems':1,'N_bifs':2,'N_branch':3,'Width':5,'Height':6,'Depth':7, 'Length':11, 'EucDistance':15, 'PathDistance':16, 'Branch_Order':18, 'Terminal_degree':19,  'TerminalSegment':20, 'Branch_pathlength':23,'Contraction':24, 'Partition_asymmetry':28, 'Bif_ampl_remote':34, 'Bif_tilt_remote':36}
LMEASURE_INDICES = {'Total_Sum':2,'Minimum':5,'Average':6,'Maximum':7,'S.D.':8}

# sholl
SHOLL_INDEPENDENT_VARIABLES = ['PathDistance', 'EucDistance']
SHOLL_DEPENDENT_VARIABLES = ['Length']

LAYERS_COLORS_DICT = {2:'tab:orange', 5:'tab:blue', 6:'tab:brown'}

# clustering for logistic regression
L1_CLUSTERING_SOMA_CENTERED = {2:((14,19),(3,4,11,17),(5,6,9,18),(1,2,10,16),(7,8,12,13,15,20,21))}
L1_CLUSTERING = {2:((1,2,4,10,11,14,16,17,18,19,20),(3,5,6,8,13,21),(7,9,12,15)),5:((3,5,6,9,11,15,16,18),(4,7,8,21,23),(1,2,10,12,13,14,17,19,20,22)),6:((3,5,7,8,9,11,12,16,17),(13,18,20),(1,2,4,6,10,14,15,19,21,22))}
L1_CLUSTERING_ALL_LAYERS_BASAL = ((3, 5, 7, 12, 13, 15, 21),(1, 2, 4, 6, 8, 9, 10, 14, 17, 18, 19, 20),(11, 16, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45),(47, 49, 51, 52, 53, 55, 56, 57, 60, 61, 62, 63, 64),(46, 48, 50, 54, 58, 59, 65, 66))

MAX_C_L1 = 10000

CCF_WIDTH = 11400







