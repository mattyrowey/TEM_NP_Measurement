//This script will measure the average diameter of nanoparticles from Transmission Electron Microscopy images.
//This script was written on ImageJ version 1.54f
//The following plugins are required:
	//BIOP
//A cellpose python environment must be installed onto the computer so that the BIOP plugin can direct the Cellpose macro.
//To do this, a custom Cellpose 2.0 model has been trained.
//This script was written by Matthew C. Rowe

//Analysis functions
purge();
initialise();
analyse();
exportData();

//Close any open images and clear existing ROIs from the ROI manager
function purge(){
	if(nImages>=1){
		run("Close All");
}
	roiManager("reset");
}

//Opens a dialog box to define key model files and analysis parameters.
//The folder containing all the images to analyse is selected
//The Cellpose Model file that will be used to segment the nanoparticles must be selected
//The "Average NP Diameter" is used to estimate the Cellpose model prediction
//"Pixel Size (A)" is used to define the number of Angstroms per pixel, calculated during image acquisition
//The save directories are susequently defined
function initialise(){
	Dialog.create("Nanoparticle Measurement");
	Dialog.addDirectory("Select Folder", "");
	Dialog.addFile("Select Cellpose Model","");
	Dialog.addNumber("Average NP Diameter (px)", 120);
	Dialog.addNumber("Pixel Size (A)", 1.55);
	Dialog.show();
	
	dir=Dialog.getString();
	roiDir=dir+fs+"ROIs"+fs;
	File.makeDirectory(roiDir);
	outputDir=dir+fs+"Output"+fs;
	File.makeDirectory(outputDir);
	cellposeModel=Dialog.getString();
	NPDiameter=Dialog.getNumber();
	//The pixel size in Angstroms is converted into nanometers
	pixelSize=Dialog.getNumber()/10;
}

//This function will perform the analysis and Loop through all images in the folder
function analyse(){
	list=getFileList(dir);
	for(i=0;i<list.length;i++){
	image=dir+list[i];
	if(endsWith(image, ".tif") || endsWith(image, ".tiff")){
		numImages=numImages+1;
	}
	}
	for(i=0;i<list.length;i++){
	image=dir+list[i];
	print("\\Clear");
	if(endsWith(image, ".tif") || endsWith(image, ".tiff")){
		open(image);
		imageName=getTitle();
		renameA=indexOf(imageName, ".tiff");
		justName=substring(imageName, 0, renameA);
		rename(justName);
		imageName=getTitle();
		print("\\Update:Analysing Image "+(i+1)+" of "+numImages+" - "+imageName);
		getDimensions(xWidth, xHeight, channels, slices, frames);
		//Each image is calibrated with the defined pixel units
		Stack.setXUnit("nm");
		run("Properties...", "channels=1 slices=1 frames=1 pixel_width="+pixelSize+" pixel_height="+pixelSize+" voxel_depth=1");
		run("Duplicate...", " ");
		selectWindow(imageName+"-1");
		run("Scale Bar...", "width=50 height=2 thickness=50 font=100 color=Black bold overlay");
		run("Flatten");
		saveAs("Tiff", outputDir+imageName+" - 10nm Scale Bar");
		close(imageName+" - 10nm Scale Bar.tif");
		close(imageName+"-1");
		selectWindow(imageName);
		Stack.setXUnit("nm");
		run("Properties...", "channels=1 slices=1 frames=1 pixel_width="+pixelSize+" pixel_height="+pixelSize+" voxel_depth=1");
		//Cellpose segments a binary mask of all nanoparticles using a custom trained model
		run("Cellpose ...", "env_path=["+cellPoseEnvironment+"] env_type=conda model= model_path=["+cellposeModel+"] diameter="+NPDiameter+" ch1=0 ch2=-1 additional_flags=[--flow_threshold, "+cellFlow+", --cellprob_threshold, "+cellProbability+", --use_gpu]");
		selectWindow(imageName+"-cellpose");
		imageNameCellpose=getTitle();
		setBatchMode("exit and display");
		selectWindow(imageNameCellpose);
		run("Label image to ROIs", "rm=[RoiManager[size=10000, visible=true]]");
		numNPs=roiManager("count");
		//The ROIs of every nanoparticle will be shrunk by 1 pixel to prevent overlap of nanoparticle measurements.
		for(n=0;n<numNPs;n++){
			roiManager("select",n);
			run("Enlarge...", "enlarge=-1");
			roiManager("Update");
}
	close(imageNameCellpose);
	newImage("NP Mask", "16-bit black", xWidth, xHeight, 1);
	Stack.setXUnit("nm");
	run("Properties...", "channels=1 slices=1 frames=1 pixel_width="+pixelSize+" pixel_height="+pixelSize+" voxel_depth=1");
	NPMask=getTitle();
	setForegroundColor(255, 255, 255);
	roiManager("deselect");
	roiManager("combine");
	roiManager("fill");
	roiManager("reset");
	//The diameter, area and circularity is measured for all nanoparticles.
	run("Set Measurements...", "area limit feret's redirect=None decimal=3");
	setThreshold(1.0000, 1000000000000000000000000000000.0000);
	//Particles smaller than 20pixels^2 or particles with a circularity measurement less than 0.60 were excluded from the final results
	run("Analyze Particles...", "size=20-Infinity circularity=0.60-1.00 show=Masks display exclude clear summarize add");
	numNPs=roiManager("count");
	numNPsArray=Array.concat(numNPsArray,roiManager("count"));
	if(numNPs!=0){
		roiManager("save", roiDir+imageName+" - NP ROIs.zip");
}
	//The measured results of individual nanoparticles are output into multiple arrays
	for(p=0;p<numNPs;p++){
		NPAreaArray=Array.concat(NPAreaArray,getResult("Area", numNPs-p-1));
		NPDiameterFeretArray=Array.concat(NPDiameterFeretArray,getResult("Feret", numNPs-p-1));
		imageNameArray=Array.concat(imageNameArray,imageName);
}
	//The size distribution of nanoparticles is plotted and the values are painted onto a new image
	run("Distribution...", "parameter=Feret or=30 and=0-60");
	selectWindow("Feret Distribution");
	Plot.getValues(bins, counts);
	if(isOpen("Diameter Distribution")){
}
	else{
		newImage("Diameter Distribution", "16-bit black", 30, 2, 1);
		distributionImage=getTitle();
		for(x=0;x<30;x++){
		makeRectangle(x, 0, 1, 1);
		run("Add...", "value="+bins[x]);
}
}
	selectWindow(distributionImage);
	for(x=0;x<30;x++){
		makeRectangle(x, 1, 1, 1);
		run("Add...", "value="+counts[x]);
}
	close("Results");
	Table.rename("Summary", "Results");
	//The average of all measured nanoparticles is output into multiple arrays
	NPAverageAreaArray=Array.concat(NPAverageAreaArray,getResult("Average Size"));
	NPAverageDiameterFeretArray=Array.concat(NPAverageDiameterFeretArray,getResult("Feret"));
	imageNameShortArray=Array.concat(imageNameShortArray,imageName);
	selectWindow("Mask of "+NPMask);
	run("Invert LUT");
	saveAs("Tiff", outputDir+imageName+" - NP Segemented Mask");
	close();
	close(NPMask);
	
	close("Results");
	close(imageName);
	roiManager("reset");
	close("Feret Distribution");
}
}
}

//The data is printed into a .CSV file to be used for subsequent GraphPad Prism Analysis
function exportData(){
	selectWindow(distributionImage);
	run("Set Measurements...", "mean redirect=None decimal=3");
	for(x=0;x<30;x++){
		run("Clear Results");
		makeRectangle(x, 0, 1, 1);
		run("Measure");
		binSizeArray=Array.concat(binSizeArray,getResult("Mean"));
		run("Clear Results");
		makeRectangle(x, 1, 1, 1);
		run("Measure");
		totalCountsArray=Array.concat(totalCountsArray,getResult("Mean"));
}
	selectWindow(distributionImage);
	saveAs("Tiff", outputDir+fs+"Distribution Image");
	close();
	
	//A .CSV file of binned nanoparticle size distrubutions is saved
	tableName="Distribution Data";
	tableName2="["+tableName+"]";
	run("Table...","name="+tableName2+" width=800 height=250");
	print(tableName2,"\\Headings:Bin Start (nm)\tTotal Counts");
	for(i=0;i<30;i++){
	print(tableName2, binSizeArray[i]+"\t"+totalCountsArray[i]+"\n");
}
	saveAs("Text", outputDir+fs+"Distribution Data.csv");
	run("Close");
	
	//A .CSV file of the averaged results of all images analysed combined is saved
	tableName="Summary Data";
	tableName2="["+tableName+"]";
	run("Table...","name="+tableName2+" width=800 height=250");
	print(tableName2,"\\Headings:Image Name\tNumber of NPs\tNP Average Area (nm^2)\tNP Average Diameter (nm)");
	for(i=0;i<NPAverageAreaArray.length;i++){
	print(tableName2, imageNameShortArray[i]+"\t"+numNPsArray[i]+"\t"+NPAverageAreaArray[i]+"\t"+NPAverageDiameterFeretArray[i]+"\n");
}
	saveAs("Text", outputDir+fs+"Summary Data.csv");
	run("Close");
	
	//A .CSV file of every individual nanoparticle measurement is saved
	tableName="All NP Data";
	tableName2="["+tableName+"]";
	run("Table...","name="+tableName2+" width=800 height=250");
	print(tableName2,"\\Headings:Image Name\tParticle Number\tNP Area (nm^2)\tNP Diameter (nm)");
	for(i=0;i<NPAreaArray.length;i++){
	print(tableName2, imageNameArray[i]+"\t"+(i+1)+"\t"+NPAreaArray[i]+"\t"+NPDiameterFeretArray[i]+"\n");
}
	saveAs("Text", outputDir+fs+"All NP Data.csv");
	run("Close");
}
////////VARIABLES//////////////
var fs=File.separator;
var dir="";
var list="";
var roiDir="";
var outputDir="";
var cellposeModel="";
var NPDiameter=0;
var pixelSize=0;
var xWidth=0;
var xHeight=0;
var imageName="";
var imageNameCellpose="";
var NPMask="";
var numNPs=0;
//cellProbaility and cellFlow are fixed values that were assigned during the training of the model
var cellProbability=0;
var cellFlow=0.4;
//cellposeEnvironment is a system dependent directory according to where the Cellpose package is installed
var cellPoseEnvironment="C:\\Anaconda3\\envs\\CellPose";
var numImages=0;
var distributionImage="";

var numNPs=0;
var numNPsArray=newArray();
var NPAverageAreaArray=newArray();
var NPAverageDiameterFeretArray=newArray();
var NPAreaArray=newArray();
var NPDiameterFeretArray=newArray();
var imageNameArray=newArray();
var imageNameShortArray=newArray();
var bins=newArray();
var counts=newArray();
var binSizeArray=newArray();
var totalCountsArray=newArray();