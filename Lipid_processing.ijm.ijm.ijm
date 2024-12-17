
// @String (label="Directory with all tiff images", style="directory") AllDir

//speed up processing；
setBatchMode(true);

// Result file Creation

resultsPath = "/Users/zhenghonghui/Downloads/lipid_results.csv";
File.saveString("sample_id, GfapArea, Tuj1Area, DapiArea, lipid_count, total_lipid_area, avg_lipid_size, Gfap_lipid, Gfap_lipid_count, TuJ1_lipid, TuJ1_lipid_count\n", resultsPath);

// Paredirectory to start the processing
processParentDirectory(AllDir);


setBatchMode(false);
	
	
function processParentDirectory(parentDirectory) {
    subDirectories = getFileList(parentDirectory);

    for (i = 0; i < subDirectories.length; i++) {
        subDir = parentDirectory + "/" + subDirectories[i];

        if (File.isDirectory(subDir)) {
            // Get the list of TIFF images in the current subdirectory
            imageList = getFileList(subDir);
            print(subDir);
            print(subDir.indexOf("40x"));
           	if (subDir.indexOf("40x") != -1) {
            
            // Process each TIFF image in the subdirectory
	            for (j = 0; j < imageList.length; j++) {
	                currentImage = subDir + "/" + imageList[j];
	                print(currentImage);
	                // Call a function to process the current TIFF image
	                processImage_binary(currentImage);
            
           		}
             	lipid_parameter_cal(subDir);
           	}
           
            print("finished one folder");    
            //pass the subfolder name into final csv 

           
            
            close("*");
        }
    }
}

function processImage_binary(imagePath) {
	
	// Check if the file name ends with "tif" and did not contains "rgb" "binary"
    if (endsWith(imagePath, ".tif") && imagePath.indexOf("RGB") == -1 
    && imagePath.indexOf("binary") == -1) {
        // If yes, perform the Z projection
        open(imagePath);
        run("Set Scale...", "distance=3.2181 known=1 unit=µm global");
   
  		getDimensions(w, h, channels, slices, frames);
  	
        if (nSlices() > 1 && channels == 1) {
            // Perform the Z projection
            run("Z Project...", "projection=[Max Intensity]");
            // Get the original file name without extension
        	originalFileName = File.getNameWithoutExtension(imagePath);
			// TUJ1 DAPI OSTU, GFAP lipidspot488 MAxentropy       
			print(imagePath);
        	if (endsWith(imagePath, "C3.tif")) {
        		setAutoThreshold("Otsu dark");
        		rename("TuJ1");
        		
        	}
        	
        	if (endsWith(imagePath, "C1.tif") ) {
        		setAutoThreshold("Otsu dark");
        		rename("Dapi");
        		
        	}
        	
        	if (endsWith(imagePath, "C2.tif")) {
        		run("Subtract Background...", "rolling=10");
        		setAutoThreshold("MaxEntropy dark");
        		rename("Lipidspot");
        		
        	}
        	
        	if (endsWith(imagePath, "C4.tif")) {
        		setAutoThreshold("MaxEntropy dark");
        		rename("GFAP");
    
        	}
        	
        	run("Convert to Mask");
        	
        	// Save the result with the new name in the same directory
        	currentTitle = getTitle(); // keep the mask name
        	outputImagePath = File.getDirectory(imagePath) + "/binary_" + originalFileName + ".tif";
        	saveAs("Tiff", outputImagePath);
        	rename(currentTitle); // keep the mask name
        	
        	print("Processed and saved: " + outputImagePath);
        	
     
        } else {
        print("Skipped processing for: " + imagePath); 
        }
        
    } else {
        print("Skipped processing for: " + imagePath);
    }
}

function lipid_parameter_cal(subDir) {
	

	selectWindow("TuJ1");
	setAutoThreshold("Default");
	run("Measure");
	Tuj1Area = getResult("Area", nResults-1);

	selectWindow("Dapi");
	setAutoThreshold("Default");
	run("Measure");
	DapiArea = getResult("Area", nResults-1);

	selectWindow("GFAP");
	setAutoThreshold("Default");
	run("Measure");
	GfapArea = getResult("Area", nResults-1);

	// clear the roi manager in case the formmer rois are calculated;
	if (roiManager("count") > 0){
		roiManager("Delete");
	}
	
	if (roiManager("count") > 0){
		roiManager("Delete");
	}
	
	selectWindow("Lipidspot");
	run("Analyze Particles...", "size=10-Infinity pixel show=Masks display summarize add");
	selectWindow("Summary");
	IJ.renameResults("Summary","Results");
	lipid_count = getResult("Count", 0);
	total_lipid_area = getResult("Total Area", 0);
	avg_lipid_size = getResult("Average Size", 0);			

	run("Clear Results");
	selectImage("Mask of Lipidspot");
	imageCalculator("AND create", "TuJ1", "Mask of Lipidspot");
	selectWindow("Result of TuJ1");
	setAutoThreshold("Default");
	run("Measure");
	TuJ1_lipid = getResult("Area", nResults-1);

	run("Clear Results");
	TuJ1_lipid_count = 0;
	for (i = 0; i < roiManager("count"); i++)
	{ 
		selectWindow("Result of TuJ1");
		setAutoThreshold("Default");
		roiManager("Select", i); // count positive roi
		run("Measure");
		roi_area= getResult("Area", nResults-1);
		if (roi_area >0) {
			TuJ1_lipid_count = TuJ1_lipid_count + 1;
		}
	}

	run("Clear Results");
	
	imageCalculator("AND create", "GFAP", "Mask of Lipidspot");
	selectWindow("Result of GFAP");
	setAutoThreshold("Default");
	run("Measure");
	Gfap_lipid = getResult("Area", nResults-1);

	run("Clear Results");
	Gfap_lipid_count = 0;
	for (i = 0; i < roiManager("count"); i++)
	{ 
		selectWindow("Result of GFAP");
		setAutoThreshold("Default");
		roiManager("Select", i); // count positive roi
		run("Measure");
		roi_area= getResult("Area", nResults-1);
		if (roi_area >0) {
			Gfap_lipid_count = Gfap_lipid_count + 1;
		}
	}
	
	if (roiManager("count") > 0){
		roiManager("Delete");
	}
	run("Clear Results");
	
	
	s = split(subDir, '/');
	sample_id =  toString(s[lengthOf(s)-1]);

    // append all the results into csv file;	
	File.append(sample_id + ',' + toString(GfapArea) + ',' + toString(Tuj1Area) + ',' + toString(DapiArea) + ',' + toString(lipid_count) 
	+ ',' + toString(total_lipid_area) + ',' + toString(avg_lipid_size) + ',' + toString(Gfap_lipid) + ',' + toString(Gfap_lipid_count) 
	+ ',' + toString(TuJ1_lipid) + ',' + toString(TuJ1_lipid_count)+ "\n", resultsPath);
}







