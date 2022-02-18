setBatchMode(true);

// parse arguments
args = getArgument();
args = split(args, "?");

skeletonPath = args[0];
rawPath = args[1];
outPath = args[2];
swcName = args[3];

print("");
print("SWC: " + swcName);
print("Raw Path: " + rawPath);
print("Skeleton Path: " + skeletonPath);
print("Out Path: " + outPath);
print("");

// transform raw from oblique to coronal
print("Open Raw...");
run("Bio-Formats (Windowless)", "open="+rawPath);

//print("Reslice...");
run("Reslice [/]...", "output=1.000 start=Left flip avoid");

//print("Save Reslice...");
reslicePath = outPath + "/" +  swcName + "_reslice_raw.tif";
saveAs("Tiff", reslicePath);

//print("Shear...");
shearPath =  outPath + "/" + swcName + "_shear_raw.tif";
run("shear ","inputpath=" + reslicePath + " shearfile=/data/elowsky/OLST/registration/shear_anisotropic outputpath=" + shearPath);

//print("Load Sheared...");
open(shearPath);

//print("Reslice 2...");
run("Reslice [/]...", "output=1.000 start=Left avoid");

//print("Rotate...");
run("Rotate 90 Degrees Right");

//print("Save...")
outPathRaw = outPath + "/" + swcName + "_raw.tif";
saveAs("Tiff", outPathRaw);

// delete temp files
File.delete(reslicePath);
File.delete(shearPath);

// transform skeleton from oblique to coronal

print("Open Skeleton...");
open(skeletonPath);

print("Reslice...");
run("Reslice [/]...", "output=1.000 start=Left flip avoid");

print("Save Reslice...");
reslicePath = outPath + "/" + swcName + "_reslice_skeleton.tif";
saveAs("Tiff", reslicePath);

print("Shear...");
shearPath =  outPath + "/" + swcName + "_shear_skeleton.tif";
run("shear ","inputpath=" + reslicePath + " shearfile=/data/elowsky/OLST/registration/shear_anisotropic outputpath=" + shearPath);

print("Load Sheared...");
open(shearPath);

print("Reslice 2...");
run("Reslice [/]...", "output=1.000 start=Left avoid");

print("Rotate...");
run("Rotate 90 Degrees Right");

print("Save...")
outPathSkeleton = outPath + "/" + swcName + "_skeleton_dilated.tif";
saveAs("Tiff", outPathSkeleton);

// delete temp files
File.delete(reslicePath);
File.delete(shearPath);

run("Quit");

