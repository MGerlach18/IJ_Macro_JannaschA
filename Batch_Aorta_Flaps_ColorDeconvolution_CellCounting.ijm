//getting input paramaters
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ Integer (label="Pyramid Level to do the StarDist analyis", style="slider", min=1, max=8, value= 2, stepSize=1, persist=true) level

//Preparing Stage
run("Bio-Formats Macro Extensions");
setOption("BlackBackground", true);
run("Set Measurements...", "area redirect=None decimal=0");
print("\\Clear");
run("Clear Results");
close("*");
setBatchMode(false);

if (roiManager("count")>0) {
roiManager("Deselect");
roiManager("Delete");
}

File.makeDirectory(output + "\\Manual_ROIS");
File.makeDirectory(output + "\\Nuclear_ROIS");

//Processing the folder to create previes on lower scale
processFolder_regionSelect(input);
waitForUser("Preselection Phase Finished", "Preselection loop has been completed, processing will continue in background now");
setBatchMode(true);
Table.create("Summary");

//Processing the folder again to do the complete Processing
processFolder_regionProcess(input);

//save Results as .csv-file and clean up
selectWindow("Summary");
saveAs("Results", output + "\\Results_Level_" + level + ".csv");
close("*");
print("Batch processing completed");

print("            { ");
print("            { ");
print("         {   } ");
print("          }_{ { ");
print("       .-{   }   }-. ");
print("      (   }     {   ) ");
print("   |`-.._____..-'| ");
print("   |                    ;--. ");
print("   |                   (__     | ");
print("   |                    |    )  ) ");
print("   |                    |  /  / ");
print("   |                     /  /    ");
print("   |                   (  / ");
print("   |                    y ");
print("    `-.._____..-'; ");

//End of Macro




//Definition of functions below


// function to scan folders/subfolders/files to find files with correct suffix and process the region creation on a lower scale
function processFolder_regionSelect(input) {
	list = getFileList(input);
	A=list.length;
	print("0 of "  + A + " files processed");
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder_regionSelect(input + File.separator + list[i]);
		if(endsWith(list[i], ".czi"))
			prepareRois(input, output, list[i]);
	}
}
// function to scan folders/subfolders/files to find files with correct suffix and do the final processing
function processFolder_regionProcess(input) {
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			rocessFolder_regionProcess(input + File.separator + list[i]);
		if(endsWith(list[i], ".czi"))
			processFile(input, output, list[i]);
	}
}

//function to annotate ROIs in low-magnification images which load faster 
function prepareRois(input, output, file) {
//Import
run("Bio-Formats Importer", "open=[" + input + "\\" + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_3");
title=getTitle();
array1=split(title," ");
title=array1[0];

//check wheter the presegmentation for this file already exists
if (File.exists(output + "\\Manual_ROIS\\" + title  +".zip")) {
	close("*");
} else {

//preparing file
rename(title);
Stack.getDimensions(width, height, channels, slices, frames);

//adjusting to widescreen for improved segmentation properties
if (width<height) {
	run("Rotate 90 Degrees Right");
}

//Detecting total slice area
run("Duplicate...", "title=Mask_Aorta duplicate");
run("RGB Color");
run("8-bit");
getStatistics(area, mean, min, max, std, histogram);
run("Duplicate...", "title=Mask_Aorta_backfill duplicate");

//Compensate missing tiles
if (min==0) {
setThreshold(0, 1);
//
run("Convert to Mask");
run("Create Selection");
selectWindow("Mask_Aorta (RGB)");
run("Restore Selection");
run("Fill", "slice");
run("Select None");
}

//Compensate White Tiles
selectWindow("Mask_Aorta (RGB)");
run("Duplicate...", "title=Mask_Aorta_frontfill duplicate");
setThreshold(255, 255);
run("Convert to Mask");
run("Create Selection");
run("Make Inverse");

//Tissue Detection
selectWindow("Mask_Aorta (RGB)");
run("Restore Selection");
setAutoThreshold("Huang");
run("Convert to Mask");
run("Remove Outliers...", "radius=3 threshold=50 which=Bright");
run("Remove Outliers...", "radius=3 threshold=50 which=Dark");
run("Create Selection");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "Aorta_Flap_Total_Area");

close("Mask*");

//create manual ROIs for analysis
selectWindow(title);
setTool(3);
waitForUser("Create Total Area");
roiManager("add");
roiManager("Select", 1);
roiManager("rename", "Total_Area");
run("Select None");

setTool(3);
waitForUser("Create Smooth Area");
roiManager("add");
roiManager("Select", 2);
roiManager("rename", "Smooth_Area");
run("Select None");

setTool(3);
waitForUser("Create Center Area");
roiManager("add");
roiManager("Select", 3);
roiManager("rename", "Center_Area");
run("Select None");

setTool(3);
waitForUser("Create Rough Area");
roiManager("add");
roiManager("Select", 4);
roiManager("rename", "Rough_Area");
run("Select None");

// save the RoiSet for further processing steps
roiManager("save", output + "\\Manual_ROIS\\" + title  +".zip");
roiManager("deselect");
roiManager("delete");
run("Close All");
}
}

// function to do the actual processing (nuclear counbting with StarDist)
function processFile(input, output, file) {

// Bio-Formats Importer opens files in specified pyramid stage and gets metadata
//getting the Pyramid stage scale
run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_2");
height2=getHeight();
close();
run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_3");
height3=getHeight();
close();
base=round(height2/height3);

run("Bio-Formats Importer", "open=[" + input + "\\" + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+ level);
title=getTitle();
array1=split(title," ");
title=array1[0];
rename(title);
Stack.getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);


//correction for widescreen 
if (width<height) {
	run("Rotate 90 Degrees Right");
}

//Correct for wrong scaling in Pyramid formats
pixelWidth2=(pixelWidth*pow(base, level-1));
pixelHeight2=(pixelHeight*pow(base, level-1));
run("Properties...", "channels="+channels+" slices="+slices+" frames="+frames+" pixel_width="+pixelWidth2+" pixel_height="+pixelHeight2+" voxel_depth=3.0000000");

//Preparing&Getting statistics

//Creating Measurements tabel
selectWindow("Summary");
Table.set("File Name", i, title);

//Adapting ROIs to actual pyramid size and measure total area of ROIs
roiManager("open", output + "\\Manual_ROIs\\" + title + ".zip");
for (o=0; o<roiManager("count"); ++o) {
	roiManager("Select", o);
	ROI=Roi.getName;
	// Scale ROI
	run("Scale... ", "x="+ pow(base, 3-level)+" y="+ pow(base, 3-level));

	// Replace old ROI with scaled one
	roiManager("update");
	
	//measure ROI areas
	run("Measure");
	area=getResult("Area", o);
	selectWindow("Summary");
	Table.set(ROI + "_Total area [µm²]", i, area);
}

Table.update;
run("Clear Results");

//Color Deconvolution 2 - Preset H&E 2
selectWindow(title);
run("RGB Color");
rename(title);
run("Colour Deconvolution2", "vectors=[H&E 2] output=32bit_Absorbance simulated cross hide");
print("\\Clear");
selectWindow(title + "-(Colour_1)A");
run("8-bit");
close(title + "-(Colour_2)A");
close(title + "-(Colour_3)A");


//Preparing measurement sub-images for nucleus counting
selectWindow(title + "-(Colour_1)A");
roiManager("select", 2);
getSelectionBounds(x_smooth, y_smooth, w, h);
run("Duplicate...", "title=Smooth_Area");
run("Clear Outside");
run("Select None");

selectWindow(title + "-(Colour_1)A");
roiManager("select", 3);
getSelectionBounds(x_center, y_center, w, h);
run("Duplicate...", "title=Center_Area");
run("Clear Outside");
run("Select None");

selectWindow(title + "-(Colour_1)A");
roiManager("select", 4);
getSelectionBounds(x_rough, y_rough, w, h);
run("Duplicate...", "title=Rough_Area");
run("Clear Outside");
run("Select None");

//applying StarDist to the presegmented images and count
roiManager("Deselect");
roiManager("Delete");

selectWindow("Smooth_Area");
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Smooth_Area', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'25.0', 'percentileTop':'100.0', 'probThresh':'0.4', 'nmsThresh':'0.4', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");

// count detected nuclei (ROIs) and put it into results table
count=roiManager("count");
selectWindow("Summary");
Table.set("Smooth_Area_CellCount", i, count);
Table.update;
roiManager("Deselect");

// Shift ROIs to their original position in the image for verification purposes
for (o=0; o<roiManager("count"); ++o) {
	roiManager("Select", o);
	run("Translate... ", "x="+x_smooth+" y="+y_smooth);
	roiManager("update");
}
//save RoiSet controls for review
roiManager("Save", output + "\\Nuclear_ROIs\\" + title + "_Smooth_Area.zip");
roiManager("delete");

// count detected nuclei (ROIs) and put it into results table
selectWindow("Center_Area");
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Center_Area', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'25.0', 'percentileTop':'100.0', 'probThresh':'0.4', 'nmsThresh':'0.4', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
count=roiManager("count");
selectWindow("Summary");
Table.set("Center_Area_CellCount", i, count);
Table.update;
roiManager("Deselect");
// Shift ROIs to their original position in the image  for verification purposes
for (o=0; o<roiManager("count"); ++o) {
	roiManager("Select", o);
	run("Translate... ", "x="+x_center+" y="+y_center);
	roiManager("update");
}

//save RoiSet controls for review
roiManager("Save", output + "\\Nuclear_ROIs\\" + title + "_Center_Area.zip");
roiManager("delete");

// count detected nuclei (ROIs) and put it into results table
selectWindow("Rough_Area");
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Rough_Area', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'25.0', 'percentileTop':'100.0', 'probThresh':'0.4', 'nmsThresh':'0.4', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
count=roiManager("count");
selectWindow("Summary");
Table.set("Rough_Area_CellCount", i, count);
Table.update;
roiManager("Deselect");
// Shift ROIs to their original position in the image for verification purposes
for (o=0; o<roiManager("count"); ++o) {
	roiManager("Select", o);
	run("Translate... ", "x="+x_rough+" y="+y_rough);
	roiManager("update");
}
//save RoiSet controls for review

roiManager("Save", output + "\\Nuclear_ROIs\\" + title + "_Rough_Area.zip");
roiManager("delete");

//Clean up

roiManager("Deselect");
roiManager("Delete");
close("*");
}
