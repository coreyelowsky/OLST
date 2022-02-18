#########################
swc_path = '/data/elowsky/OLST/swc_analysis/automatically_traced/barrel_cortex/coronal_normalized/180606_8.swc'
ROTATION_PERIOD = 10
SOMA_SIZE = 20
#######################

# load nblast libarary
library(nat)


# open 3d window
open3d()


# read swc
swc_nblast = read.neuron(swc_path)


# plot neuron in 3D
plot3d(swc_nblast,soma=SOMA_SIZE)


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

movie3d(f, duration = ROTATION_PERIOD, dir='/data/elowsky/science_cafe/')






