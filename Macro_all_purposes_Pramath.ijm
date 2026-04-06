
'------------------------------------------------------------------------------------------------------------------------------
'------------------------------------------    STEP 1 = CLEANING      ---------------------------------------------------------
'------------------------------------------------------------------------------------------------------------------------------
//RUN THE PARTS SEPARATED BY "//COMMENT LINE WITH INSTRUCTIONS" INDIVIDUALLY 
// SELECT COMMANDS AND PRESS CTRL+SHIFT+R

// OPTION1 = Inport input image. Insert right file name and directory first
// OPTION2 = FASTER, Drag and drop all folder containing images (example: "15")
run("Image Sequence...", "open=[H:/Cristina_Projects/2P/Ed_Mice/20220502_RPA_2Y/RPA_23_1/15/14-09-18_PMT - PMT [Blue] _C00_z-Stepper Z0000.ome.tif]");


// From stack to hyperstack
	run("Stack to Hyperstack...", "order=xyztc channels=3 slices="+nSlices/3+" frames=1 display=Color");
	Stack.setXUnit("um");
	run("Properties...", "channels=3 frames=1 pixel_width=0.48 pixel_height=0.48 voxel_depth=2");
	title = getTitle();

// NEED TO ROTATE FIGURE? APPLIES ONLY IF THE SPECIMEN IS CLEARLY TILTED AT THE LEFT OR RIGHT END. Edit y-angle only. Repeat until you obtain the optimal orientation. 
	run("TransformJ Rotate", "z-angle=0.0 y-angle=5 x-angle=0.0 interpolation=Linear background=0.0 adjust");
	makeRectangle(0, 0, 1042, 1042);
	waitForUser("Remove black borders"); //wait for user action		
	run("Crop");


// CH 1 - Clean artefacts
	title = getTitle();
	run("Split Channels");
	run("Grays");
	imageCalculator("Subtract create stack", "C1-"+title,"C2-"+title);
	run("Z Project...", "projection=[Max Intensity]");
	getRawStatistics(count, mean, min, max, std);
	waitForUser("Check if MAX is an ARTEFACT"); //wait for user action		
	close()
	setSlice(24);
	run("Color Balance...");
	setMinAndMax(0, 14875);
	run("Apply LUT", "stack");
	
	run("Reverse"); // run this line only if you do not need to reverse (reverse when the shape of the vessel surface is convex) 



// CH2 - Clean artefacts
	run("Remove Outliers...", "radius=50 threshold=14000 which=Bright stack");
	run("Remove Outliers...", "radius=50 threshold=8000 which=Bright stack");
	run("Remove Outliers...", "radius=50 threshold=5200 which=Bright stack");
	run("Remove Outliers...", "radius=50 threshold=2500 which=Bright stack");
	run("Remove Outliers...", "radius=50 threshold=1000 which=Bright stack"); 
	// then clean manually using the rectangle tool 
	// done? you can run the following commands 
	setSlice(24);
	run("Color Balance...");
	setMinAndMax(0, 1000); //constant ->65000
	run("Apply LUT", "stack");

	run("Reverse"); // run this line only if you do not need to reverse (reverse when the shape of the vessel surface is convex) 


// CH3 - Clean artefacts
	// clean manually using the rectangle tool 
	// done? you can run the following commands 
	
	setSlice(24);
	run("Color Balance...");
	setMinAndMax(0, 13000); //constant ->65000
	run("Apply LUT", "stack");

	run("Reverse"); // run this line only if you do not need to reverse (reverse when the shape of the vessel surface is convex) 


//CHANGE STACK NAME AND CH AND SAVE
	//create subdir only once for each stack
	splitDir= "H:/Cristina_Projects/2P/Ed_Mice/20220502_RPA_2Y/"; //edit  
	File.makeDirectory(splitDir); 
	splitDir= "F:/_F27mo/240826_12-56-28_F27moRPA_4_15mmHg"+ "/15_edit/"; //edit
	//print(splitDir);
	File.makeDirectory(splitDir); 
	saveAs("F:/_F27mo/240826_12-56-28_F27moRPA_4_15mmHg"+ "/15_edit/C3-15.tif"); //edit C1 - C2 - C3 + stack number


// RUN MATLAB CODE Unwrap.m


'------------------------------------------------------------------------------------------------------------------------------
'---------------------------------------   STEP 2 = POST-TRANSFORMATION   -----------------------------------------------------
'------------------------------------------------------------------------------------------------------------------------------


// FROM YZ TO XY AFTER TRASFORMATION WITH ROTATION A POSTERIORI

// 1.Reslice 
	dir = "F:/_M27mo/202300605_M27moCTRL/230605_M27moCTRL_3_15mmHg/"; //edit
	stack = 15;
	File.delete(dir+stack+"_edit/"+stack+".mat");
	//print(dir+stack+"_edit/");
	//File.openSequence(dir+stack+"_edit/" , "filter=_C");
	File.openSequence(dir+stack+"_edit/", "filter=_C");
	run("8-bit");
	run("Stack to Hyperstack...", "order=xyztc channels=3 slices="+nSlices/3+" frames=1 display=Color");
	makeRectangle(8, 63, 1042, 119);
	waitForUser("Resize rectangle"); //wait for user action
	run("Crop");
	Stack.setXUnit("um");
	run("Properties...", "channels=3 slices=1042 frames=1 pixel_width=0.48 pixel_height=1 voxel_depth=0.48"); //change value 1042 if needed after rotation
	run("Reslice [/]...", "output=1 start=Top avoid");
	run("Rotate 90 Degrees Left");

// 2. final rotation in Y or X ONLY IF NEEDED, Close result and repeat the process if needed
	run("TransformJ Rotate", "z-angle=0.0 y-angle=25.0 x-angle=0.0 interpolation=Linear background=0.0 adjust");
	makeRectangle(0, 0, 1042, 1042);
	waitForUser("Remove black borders"); //wait for user action		
	run("Crop");
	
	
// 3. Save Composite and delete profile images
	//change this 
	dir = "F:/_M27mo/202300605_M27moCTRL/230605_M27moCTRL_3_15mmHg/"; //edit
	stack = 15;
	saveAs("Tiff", dir+stack+"_edit/"+stack+"_Composite.tif");
	run("Close All");

	setBatchMode(true);
	for(i=1;i<=1042;i++)
	{
		File.delete(dir+stack+"_edit/"+stack+"_C1_"+i+".tiff");
		File.delete(dir+stack+"_edit/"+stack+"_C2_"+i+".tiff");
		File.delete(dir+stack+"_edit/"+stack+"_C3_"+i+".tiff");
	}