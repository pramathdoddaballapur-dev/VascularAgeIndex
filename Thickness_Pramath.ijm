

extension = ".tif";
dirpath = getDirectory ("Choose Source Directory ");
dirpath_out = dirpath + "/Thickness/"; 
File.makeDirectory(dirpath_out); 

fnames = newArray(5,5); //edit list of stacks to process
//open(dirpath + "25" +"_edit/" + "25" + "_Composite.tif");
for(i=0;i<=0;i++) //edit index i
{

	stack = fnames[i];
	print(stack);
	//open(dirpath + stack + "_edit/Transform/" + stack + "_Composite_rotated.tif");
	print(dirpath + stack +"_edit/" + stack + "_Composite.tif");
	open(dirpath + stack +"_edit/" + stack + "_Composite.tif");
	//open(dirpath +"100_edit/100_Composite.tif");


	fpath = dirpath_out + fnames[i];
	Table.create("Table_1");
	Table.create("Table_2");
	Table.create("Table_3");	


for(k=0;k<=2;k++)
{
	
	for(g=0;g<=1;g++)
	{
			makeRectangle(300*k+120, 300*g+250, 209, 208);
			run("Duplicate...", "duplicate");
			run("Split Channels");

  			run("Plot Z-axis Profile");

  			xpoints = newArray ();
  			ypoints = newArray ();
  			Plot.getValues (xpoints, ypoints);
 	 		run ("Clear Results");

	  		updateResults ();
	  	    selectWindow("Table_3");
	  	    Table.setColumn("y_"+k+"_"+g, ypoints);
	
 	 		close();
 	 		run("Close");
 	 		run("Plot Z-axis Profile");

 	 		xpoints = newArray ();
 	 		ypoints = newArray ();
	  		Plot.getValues (xpoints, ypoints);
 	 		run ("Clear Results");
 
  			updateResults ();
	  	    selectWindow("Table_2");
	  	    Table.setColumn("y_"+k+"_"+g, ypoints);
  			
			close();
 	 		run("Close");
 	 		run("Plot Z-axis Profile");

  			xpoints = newArray ();
 			ypoints = newArray ();
 			Plot.getValues (xpoints, ypoints);
  			run ("Clear Results");
  		
  			updateResults ();
  			selectWindow("Table_1");
	  	    Table.setColumn("y_"+k+"_"+g, ypoints);

	  	    close();
 	 		run("Close");

			}
		}

	selectWindow("Table_3");
	print("Results", fpath + "_3.csv");
	saveAs("Results", fpath + "_3.csv");
//	saveAs("Table_3", fpath + "_3.csv");
	selectWindow("Table_2");
	saveAs("Results", fpath + "_2.csv");
//	saveAs("Table_2", fpath + "_2.csv");
	selectWindow("Table_1");
	saveAs("Results", fpath + "_1.csv");
//	saveAs("Table_1", fpath + "_1.csv");
	
}

run("Close All");
close("*")

