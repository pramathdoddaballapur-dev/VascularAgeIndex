dirpath = getDirectory ("Choose Source Directory ");
dirpath_out = dirpath + "/Collagen/"; 
File.makeDirectory(dirpath_out); 

//fnames = newArray(1,2,3,4R,5,6,7,8,9,10,11,12,13,14,15,16,17,18); ' list of stacks to process
fnames = newArray(100, 80); // list of stacks to process


for(i=0;i<=0;i++)
{
	stack = fnames[i];
	//open(dirpath + stack + "_edit/Transform/" + stack + "_Composite_rotated.tif");
	open(dirpath + stack +"_edit/" + stack + "_Composite.tif");
	fpath = dirpath_out + fnames[i];
	rename("s.tif");
	run("Duplicate...", "duplicate channels=1-1");
	selectWindow("s.tif");
	close();
	selectWindow("s-1.tif");

	// Collagen straightness
	makeRectangle(160, 160, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	roiManager("Show All");
	makeRectangle(420, 420, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	makeRectangle(680, 160, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	makeRectangle(160, 680, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	makeRectangle(680, 680, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	
	run("Set Measurements...", "area mean min stack redirect=None decimal=3");
	setTool("freeline");
	waitForUser("Draw the profile of 5 fibers, and then their edge-to-edge distance. Repeat for 5 subvolumes. Click OK when you have selected 25 fibers."); //wait for user action
	roiManager("Measure");
	saveAs("Table", fpath + "_straightness.csv");
	waitForUser("Copied?"); //wait for user action
	roiManager("Select All");
	roiManager("Delete");
	Table.deleteRows(0, 1000);

	// Collagen bundle width
	
	makeRectangle(160, 160, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	roiManager("Show All");
	makeRectangle(420, 420, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	makeRectangle(680, 160, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	makeRectangle(160, 680, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	makeRectangle(680, 680, 208, 209);
	Roi.setPosition(1);
	roiManager("Add");
	setTool("line");
	waitForUser("Select 5 transversal section of collagen bundles"); //wait for user action
	roiManager("Measure");
	saveAs("Table", fpath + "_bundle_width2.csv");
	waitForUser("Copied?"); //wait for user action
	Table.deleteRows(0, 1000);
	roiManager("Select All");
	roiManager("Delete");	
	run("Close");
	
}