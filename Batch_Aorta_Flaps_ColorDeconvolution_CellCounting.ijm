#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ Integer (label="Pyramid Level to analyze", style="slider", min=1, max=8, value= 2, stepSize=1, persist=true) level

//Obsolete, since Color Deconvolution 1 has been discontinued?
//#@ String(label="Deconvolution Algorithm", choices={"ColorDeconvolution 1", "ColorDeconvolution 2"}, value="ColorDeconvolution 1", style="radioButtonHorizontal") T

run("Bio-Formats Macro Extensions");

//Preparing Stage
//setBatchMode(true);
setOption("BlackBackground", true);
print("\\Clear");
close("*");
if (roiManager("count")>0) {

roiManager("Deselect");
roiManager("Delete");

}
Table.create("Summary");
File.makeDirectory(output + "\\Manual_ROIS");
File.makeDirectory(output + "\\Nuclear_ROIS");

//Processing the folder 
processFolder(input);

//save Results as .csv-file and clean up
selectWindow("Summary");
saveAs("Results", output + "\\Results_Level_" + level + ".csv");
close("*");
print("Batch processing completed");

print("            { ");
print("            { ");
print("         {   } ");
print("          }_{ __{ ");
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

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	A=list.length;
	print("0 of "  + A + " files processed");
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], ".czi"))
			processFile(input, output, list[i]);
	}
}

// function to do the actual processing
function processFile(input, output, file) {

// Bio-Formats Importer opens files in specified pyramid stage and gets metadata

run("Bio-Formats Importer", "open=[" + input + "\\" + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+ level);
title=getTitle();
Stack.getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);

//Correct for wrong scaling in Pyramid formats
pixelWidth2=(pixelWidth*pow(2, level-1));
pixelHeight2=(pixelHeight*pow(2, level-1));
run("Properties...", "channels="+channels+" slices="+slices+" frames="+frames+" pixel_width="+pixelWidth2+" pixel_height="+pixelHeight2+" voxel_depth=3.0000000");

//Preparing&Getting statistics
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
roiManager("save", output + "\\Manual_ROIS" + title  +".zip");

//Creating Measurements 
selectWindow("Summary");
Table.set("Slice", i, title);

//Measuring Area for all Rois in ROI Manager and
//put them in a "Summary" table
for (a = 0; a < roiManager("count"); a++) {
roiManager("select", a);
ROI=Roi.getName;
Table.set(ROI + "_Total area [µm²]", i, (Roi.size*pixelWidth2*pixelHeight2));
}
Table.update;

//Color Deconvolution - Preset H&E DAB
selectWindow(title);
run("RGB Color");
rename(title);
run("Colour Deconvolution", "vectors=H&E");
print("\\Clear");
close(title + "-(Colour_2)");
close(title + "-(Colour_3)");


//Preparing measuremnt images for nucleus counting
selectWindow(title + "-(Colour_1)");
run("Duplicate...", "title=Smooth_Area");
roiManager("select", 2);
run("Clear Outside");
run("Select None");
run("Invert");
run("Grays");

selectWindow(title + "-(Colour_1)");
run("Duplicate...", "title=Center_Area");
roiManager("select", 2);
run("Clear Outside");
run("Select None");
run("Invert");
run("Grays");

selectWindow(title + "-(Colour_1)");
run("Duplicate...", "title=Rough_Area");
roiManager("select", 2);
run("Clear Outside");
run("Select None");
run("Invert");
run("Grays");

//applying StarDist to the presegmented images and count
roiManager("Deselect");
roiManager("Delete");

selectWindow("Smooth_Area");
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Smooth_Area', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'0.0', 'percentileTop':'100.0', 'probThresh':'0.4', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'20', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
count=roiManager("count");
selectWindow("Summary");
Table.set("Smooth_Area_CellCount", i, count);
Table.update;
roiManager("Deselect");
roiManager("Save", output + "\\Nuclear_ROIs\\" + title + "Smooth_Area.zip");
roiManager("delete");

selectWindow("Center_Area");
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Center_Area', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'0.0', 'percentileTop':'100.0', 'probThresh':'0.4', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'20', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
count=roiManager("count");
selectWindow("Summary");
Table.set("Center_Area_CellCount", i, count);
Table.update;
roiManager("Deselect");
roiManager("Save", output + "\\Nuclear_ROIs\\" + title + "Center_Area.zip");
roiManager("delete");

selectWindow("Rough_Area");
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Rough_Area', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'0.0', 'percentileTop':'100.0', 'probThresh':'0.4', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'20', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
count=roiManager("count");
selectWindow("Summary");
Table.set("Rough_Area_CellCount", i, count);
Table.update;
roiManager("Deselect");
roiManager("Save", output + "\\Nuclear_ROIs\\" + title + "Center_Area.zip");
roiManager("delete");

//Clean up

roiManager("Deselect");
roiManager("Delete");
close("*");
print((i+1)  + " of "  + A + " files processed");
}
