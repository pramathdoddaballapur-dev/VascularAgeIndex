// CELL COUNTER
// - ADVENTITIAL, SM AND ENDOTH CELLS CHOOSEN BY SHAPE
// - SELECT 6 VOLUMES CONSISTENT WITH THICKNESS (AND VOLUME FRACTION) 

  
dirpath = getDirectory ("Choose Source Directory ");

//fnames = newArray(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18); // list of stacks to process
fnames = newArray(100, 80); // list of stacks to process.

for(i=0;i<=0;i++)
{

	stack = fnames[i];
	//open(dirpath + stack + "_edit/Transform/" + stack + "_Composite.tif");  //Alexia 
	//open(dirpath + "/Composites/"  + stack + "_Composite.tif");  			  //New
	open(dirpath + stack +"_edit/" + stack + "_Composite.tif");  // change based on your way to save Composite
	rename("s");
	run("Split Channels");
	close("C1-s");
	close("C2-s");
	selectWindow("C3-s");
	
for(k=0;k<=2;k++)
{
	
	for(g=0;g<=1;g++)
	{ 
			makeRectangle(300*k+120, 300*g+250, 207, 207);
			run("Grays");
			run("Brightness/Contrast...");
			run("Duplicate...", "title=s duplicate");
			setTool("multipoint");
			waitForUser("Pick adventitial cells. Write number and press Enter"); //wait for user action
			run("Select None");
			waitForUser("Pick SMCs. Write number and press Enter"); //wait for user action
			run("Select None");
			waitForUser("Pick Endothelial cells. Write number and press Enter"); //wait for user action
			close("s");
	}
}
run("Close All");
}



