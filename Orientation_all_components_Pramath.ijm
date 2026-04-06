// ORIENTATION DISTRIBUTION in 3 channel STACKs (COLLAGEN-ELASTIN-CELL NUCLEI)

g = newArray(2,3,4); //edit SAMPLE IDs
for(f=0;f<=0;f++)   
{
	/*if(f==0)
	{
		path_in = "D:/Dr. Cavinato/Pramath/20240515_Dar_Ao_arches/20240515_Dar_Ao_arches/Arch1/240515_11-34-43_Arch1_2/";
	}
	else if(f==1)
	{
		path_in = "D:/Dr. Cavinato/Pramath/20240515_Dar_Ao_arches/20240515_Dar_Ao_arches/Arch1/240515_11-44-33_Arch1_3/";
	}
	else 
	{
		path_in = "D:/Dr. Cavinato/Pramath/20240515_Dar_Ao_arches/20240515_Dar_Ao_arches/Arch1/240515_11-55-46_Arch1_4/";
	}*/
	

path_in = "F:/_F18mo/20201207_RPA_F_18m_0_3/unloaded/5/";

stack_list = newArray(5,5); //edit list of stacks to process

for(kk=0;kk<=0;kk++) //edit index
{
stack = stack_list[kk];

print(path_in + stack +"_edit/"+ stack + "_Composite.tif");
open(path_in + stack +"_edit/"+ stack + "_Composite.tif");
//open(path_in + stack + "_Composite.tif");

selectWindow(stack +"_Composite.tif");
makeRectangle(4, 5, 1020, 1010);
run("Crop");
rename(stack);
run("Split Channels");
splitDir= path_in + "/Orientation/"; 
File.makeDirectory(splitDir); 


//____________CHANNEL 3________________

st=d2s(stack,0);
selectWindow("C3-"+st);
depth = nSlices;

for(ii=1;ii<=nSlices;ii++)
{
	e=d2s(ii,0);
	
	run("Duplicate...", "title=s_"+e+".tif duplicate range="+e+"-"+e);

	run("OrientationJ Distribution", "log=0.0 tensor=3.0 gradient=4 min-coherency=10.0 min-energy=10.0 harris-index=on s-distribution=on hue=Gradient-X sat=Gradient-X bri=Gradient-X ");
	selectWindow("OJ-Histogram-1-slice-1");
	Plot.getValues(x, y);
	wait(10);
	selectWindow("s_"+e+".tif");
  	run("Close");
  	wait(10);
  	selectWindow("OJ-Histogram-1-slice-1");
  	run("Close");
  	wait(10);
  	selectWindow("OJ-Histogram-1-slice-1");
  	run("Close");
  
   	for (i=0; i<x.length; i++)
   	print(x[i], y[i]);
  
};
selectWindow("Log");
saveAs("Text", splitDir+"OJ_"+st+"_C3");
run("Close");
close();


//____________CHANNEL 1________________

wait(10);
st=d2s(stack,0);

selectWindow("C1-"+st);
depth = nSlices;

for(ii=1;ii<=nSlices;ii++)
{
	e=d2s(ii,0);
	run("Duplicate...", "title=s_"+e+".tif duplicate range="+e+"-"+e);

	run("OrientationJ Distribution", "log=0.0 tensor=2.0 gradient=4 min-coherency=20.0 min-energy=5.0 harris-index=on s-distribution=on hue=Gradient-X sat=Gradient-X bri=Gradient-X ");
	selectWindow("OJ-Histogram-1-slice-1");
	Plot.getValues(x, y);
	wait(10);
	selectWindow("s_"+e+".tif");
  	run("Close");
  	wait(10);
  	selectWindow("OJ-Histogram-1-slice-1");
  	run("Close");
  	wait(10);
  	selectWindow("OJ-Histogram-1-slice-1");
  	run("Close");
  
  
   	for (i=0; i<x.length; i++)
   	print(x[i], y[i]);
  
};
selectWindow("Log");
saveAs("Text", splitDir+"OJ_"+st+"_C1");
run("Close");
close();

//____________CHANNEL 2________________

wait(50);
st=d2s(stack,0);
selectWindow("C2-"+st);
depth = nSlices;

for(ii=1;ii<=nSlices;ii++)
{
	e=d2s(ii,0);
	run("Duplicate...", "title=s_"+e+".tif duplicate range="+e+"-"+e);

	run("OrientationJ Distribution", "log=0.0 tensor=2.0 gradient=4 min-coherency=1.0 min-energy=30.0 harris-index=on s-distribution=on hue=Gradient-X sat=Gradient-X bri=Gradient-X ");
	selectWindow("OJ-Histogram-1-slice-1");
	Plot.getValues(x, y);
	wait(10);
	selectWindow("s_"+e+".tif");
  	run("Close");
  	wait(10);
  	selectWindow("OJ-Histogram-1-slice-1");
  	run("Close");
  	wait(10);
  	selectWindow("OJ-Histogram-1-slice-1");
  	run("Close");
  
  
   	for (i=0; i<x.length; i++)
   	print(x[i], y[i]);
  
};
selectWindow("Log");
saveAs("Text", splitDir+"OJ_"+st+"_C2");
run("Close");
close();

}

close("*");

}



run("8-bit");
title = getTitle();
run("OrientationJ Distribution", "tensor=2.0 gradient=4 radian=on histogram=on table=on min-coherency=20.0 min-energy=5.0 ");
selectWindow("OJ-Distribution-1");
saveAs("Results", 'G:/Cristina_Projects/2P/Ed_H/Projections/'+ title+ '.txt');
close();
selectWindow(title);
selectWindow(title+'.txt');
run("Close" );
close();
close();

run("Rotate 90 Degrees Left");
