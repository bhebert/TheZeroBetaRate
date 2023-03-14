# TheZeroBetaRate
Code for "The Zero Beta Rate" (Di Tella, Hebert, Kurlat, Wang)



Instructions:

1) Unzip archived data

 	a) cd to "TheZeroBetaRate" directory
 	b) cd to "Raw Data"
 	c) unzip Archive.zip
 
2) Run Python scripts to create dataset

	a) See README PYTHON.md for details

3) Run MATLAB code to generate results

        a) Open MATLAB. tested with R2022a on an M1 Macbook pro, known not to work on R2019b or earlier.
        
        b) cd to TheZeroBetaRate/MatlabCode
        
        c) run RunAll.m
        
        		i) this takes about 30m on an 2021 M1 Macbook Pro
        
        d) run MonetaryShock.m
        
        		i) this will run very quickly.
		
4) Run STATA to generate tables

	a) open Stata. Tested with Stata MP 17.0 on an M1 Macbook Pro.
	
	b) cd to TheZeroBetaRate/StataCode
	
	c) do MakeTables.do
	
		i) this will run very quickly.
		
5) Inspect output

	a) The output can be found in the "Output" folder



