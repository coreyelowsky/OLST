#########################
swc_path = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/normalized_oblique/'
markup_path = '/data/anarasim/data/barrelcortex/markups/swc_markups.csv'
ROTATION_PERIOD = 20
SOMA_SIZE = 20
#######################

# load nblast libarary
library(nat)

# function to save swc_id and label
save_data <- function(swc_id,label) {
	data = cbind(swc_id,label)
 	write.table(data,file=markup_path, sep=",",  col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
}

# if markups file doesnt exist create it
if (!file.exists(markup_path)){
	file.create(markup_path)
}

# all swcs in folder
swc_files = list.files(swc_path)

# read in markups
markups <- try(read.table(file = markup_path,header=FALSE,sep=",",stringsAsFactors=FALSE))
if (!inherits(markups, 'try-error')){
	swcs_completed = markups[,1]
	labels_completed = markups[,2]
}else{
	swcs_completed = c()
	labels_completed = c()
}

# open 3d window
open3d()


print("#######################")
print(paste('# of Neurons Left:',length(swc_files)-length(swcs_completed)))
print("#######################")

# iterate through all swcs in folder
for (swc_file in swc_files) {

	# remove .swc from file name
	swc_id = substr(swc_file,0,nchar(swc_file)-4)

	print(swc_id)

	# if neuron already marked up, skip it
	if ((swc_id %in% swcs_completed)){next}


	# read swc
	swc_nblast = read.neuron(paste(swc_path,swc_file,sep=""))
		
	execute = TRUE
	manual = FALSE
	while(execute){

		# clear 3D window
		clear3d()

		# plot neuron in 3D
		plot3d(swc_nblast,soma=SOMA_SIZE)

		# decide whether to allow manual manipulation or automatic rotating
		if (manual){
			Sys.sleep(ROTATION_PERIOD)
		}else{
			# play rotatingvideo
			spriteid <- NULL
			spin1 <- spin3d(rpm = 60/ROTATION_PERIOD ) # the scene spinner
			spin2 <- spin3d(rpm = 1 ) # the sprite spinner
			f <- function(time) {
			    par3d(skipRedraw = TRUE) # stops intermediate redraws
			    on.exit(par3d(skipRedraw = FALSE)) # redraw at the end

			    rgl.pop(id = spriteid) # delete the old sprite
			    cubeid <- shade3d(cube3d(), col = "red")
			    spriteid <<- sprites3d(0:1, 0:1, 0:1, shape = cubeid,
			                   userMatrix = spin2(time, 
			                     base = spin1(time)$userMatrix)$userMatrix)
			    spin1(time)
			}
	  		play3d(f, duration = ROTATION_PERIOD)

		}

		# get input from user
		cat("Label: ")
		label <- readLines("stdin",n=1)

		# execute different inputs
		if (label == "manual"){
			manual=TRUE
			execute=TRUE
		} else if (label == "quit"){
			quit()
		} else if (label == "0" || label == "1" || label == "2" || label == "3"){
			save_data(swc_id,label)
			execute=FALSE
		} else {
			execute=TRUE
		}
	}
}








