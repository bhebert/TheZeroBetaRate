# TheZeroBetaRate
Code for "The Zero Beta Rate" (Di Tella, Hebert, Kurlat, Wang)

The below instructions will replicate the results of the paper.

NOTE: THE RESULTS WILL NOT BE IDENTICAL DUE TO CRSP/COMPUSTAT CHANGES

	-The code will use your WRDS log-in to download CRSP and Compustat. These databases are updated periodically, and the version you download will not be identical to the version we downloaded when generating the results, although it should be very similar. As a result, the output will closely resemble, but not be exactly the same as, the numbers in our paper. If you need to perfectly replicate our numbers, please contact us.

Instructions:

1) Unzip archived data

 	a) cd to "TheZeroBetaRate" directory
 	b) cd to "Raw Data"
 	c) unzip Archive.zip
 
2) Run Python scripts to create dataset

	a) See README PYTHON.md for detailed instructions

3) Run MATLAB code to generate results

        a) Open MATLAB. tested with R2022a on an M1 Macbook pro and a Windows PC, known not to work on R2019b or earlier.
        
        b) cd to TheZeroBetaRate/MatlabCode
        
        c) run RunAll.m
        
        		i) this takes about 30m on an 2021 M1 Macbook Pro
        
        d) run MonetaryShock.m
        
        		i) this will run very quickly.
		ii) This will generate at Matlab warning about renaming column headers that can be ignored.
		
	e) cd to TheZeroBetaRate/Model
	
	f) run monetary_shock_model.m
	
		i) this will run very quickly.
		
4) Run STATA to generate tables

	a) open Stata. Tested with Stata MP 17.0 and 18.0 on an M1 Macbook Pro and a Windows PC.
	
	b) cd to TheZeroBetaRate/StataCode
	
	c) do MakeTables.do
	
		i) this will run very quickly.
		
5) Inspect output

	a) The output can be found in the "Output" folder
	
	b) To view the output as tables and figures, use LyX (tested with v2.3.6.2) to compile TablesAndFigures.lyx 



