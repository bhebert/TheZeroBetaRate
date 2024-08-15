# TheZeroBetaRate
Code for "The Zero Beta Rate" (Di Tella, Hebert, Kurlat, Wang)

The below instructions will replicate the results of the paper.

NOTE: THE RESULTS WILL NOT BE IDENTICAL IN THE EVENT OF CRSP/COMPUSTAT CHANGES

	-The code requires a Wharton Research Data Services (WRDS) account with access to CRSP, Compustat, "CCM" (which is CRSP/Compustat merged), and "Bond Returns by WRDS."

	-The code will use your WRDS log-in to download CRSP and Compustat. These databases are updated periodically, and the version you download will not be identical to the version we downloaded when generating the results, although it should be very similar. As a result, the output will closely resemble, but not be exactly the same as, the numbers in our paper. If you need to perfectly replicate our numbers, please contact us.

Instructions:

1) Unzip archived data

 	a) cd to "TheZeroBetaRate" directory
 	b) cd to "Raw Data"
 	c) unzip Archive.zip
		NOTE: type the above "unzip" a terminal/console. Don't double-click. You should end up with e.g. Raw Data/Factors/ff3.xlsx, not Raw Data/Archive/Factors/ff3.xlsx.
		If you end up with an "Archive" folder, you will need to manually move the contents out of that folder and into the Raw Data folder.
		
 
2) Run Python scripts to create dataset

	a) See README PYTHON.md for detailed instructions

3) Run MATLAB code to generate results. Required toolboxes: optimization and statistical learning

        a) Open MATLAB. Requires R2023b or later due to use of 'xscale' command in plots.
        
        b) cd to TheZeroBetaRate/MatlabCode
        
        c) run RunAll.m
        
        		i) this takes about 30m on an 2021 M1 Macbook Pro
		ii) when the parallel pool starts, you may be asked to allow Matlab to accept incoming connections. Allow this; the software should only use your local parallel pool unless you have specifically configured Matlab to do otherwise.
		
		iii) OPTIONAL: If you have the factors from the "Informative Factors" paper of  Amman, Hemauer, and Straumann, you can prepare that file and uncomment the associated code in RunAll.m. See the file for details.
        
        d) run MonetaryShock.m
        
        		i) this will run very quickly.
		ii) This will generate Matlab warnings about renaming column headers and extra legend entries that can be ignored.
		
	e) run EquityAnalysis.m
		
		i) this will run very quickly. 
		ii) the figure and most of the results will be generated by setting "spec=0"
		iii) the numbers quoted in the text are called "sdnpv" and "sdnpv_tbill"
		iii) one number quoted in the text is generated from running the script with "spec=2" and results but no specific numbers from "spec=1" are mentioned in the text
		
	f) cd to TheZeroBetaRate/Model
	
	g) run monetary_shock_model.m
	
		i) this will run very quickly.
		
4) Run STATA to generate tables. Requires mmerge, robreg, moremata, estout ("ssc install ..." if you need to install them)

	a) open Stata. Tested with Stata MP 17.0 and 18.0 on an M1 Macbook Pro and a Windows PC.
	
		i) depending on your Matlab and Stata versions and Mac/Windows, the dates in the intermediate .csv files might be "DMY" or "MDY" format
		   Try the "DMY" format on Mac and the "MDY" on Windows (lines 3-4 of MakeTables.do), or inspect ../Output/zero_beta_rateMain.csv to see what the format is 
	
	b) cd to TheZeroBetaRate/StataCode
	
	c) do MakeTables.do
	
		i) this will run very quickly.
		
		ii) OPTIONAL: if you ran the "Informative Factors" code in RunAll.m, include the "InfFactors" specification in the "local AltFactors ..." near the top of the .do file
		
5) Inspect output

	a) The output can be found in the "Output" folder
	
	b) To view the output as tables and figures, use LyX (tested with v2.3.6.2) to compile TablesAndFigures.lyx 



