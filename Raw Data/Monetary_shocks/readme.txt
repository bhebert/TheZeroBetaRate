The purpose of these files is to extend the Romer-Romer (2004) monetary policy shock series. Program and data are provided without any warranty. Please email jfwieland@ucsd.edu if you find any discrepancies.

When using this data in your work, cite Wieland-Yang (2019) along with Romer-Romer (2004) as a reference.

Thanks to:
1. Yeji Sung for noting a missing match in an earlier version of the code and for alerting me to Philadelphia Fed’s dataset, which includes data for August-October 1972 that was missing in the original Romer-Romer dataset. The values I use differ slightly from the Philadelphia Fed’s. This is because I do not compound annualized growth rates in order to be consistent with the original Greenbook forecasts. 
2. Pavel Kapinos for noting that the old target variable was shifted by one meeting from 2004-7.
3. Michael McMahon for noting that the May 15th 2001 meeting was incorrectly coded as May 18.

Contents:
1. RR_monetary_shock_update.do
This code file generates three Stata datasets, RR_monetary_shock_monthly.dta, RR_monetary_shock_quarterly.dta, and RR_monetary_shock_annual.dta. These correspond to the monetary shock series at monthly, quarterly, and annual frequency. Each Stata file contains four variables. The date variable "date", “resid” are the original Romer-Romer (2004) shocks, "resid_romer" are the monetary policy shocks based on the original Romer-Romer (2004) regression, and "resid_full" are the monetary policy shocks based on running the Romer-Romer (2004) regression on the full 1969-2007 sample.

2. RRimport.xls
This is the original Romer-Romer (2004) dataset of Greenbook forecasts updated to 2007. Also includes recently-published data for August-October 1972 that was missing in the original Romer-Romer dataset (these are marked in yellow).

3. RRshock_Quarterly_1.txt and RRshock_Quarterly_2.txt
These are Greenbook forecasts downloaded from the FRED database. They are used to check the entries in RRimport.xls.

4. ForecastRelease.xlsx
Contains forecast release dates for meetings where FRED does not have the data.

5. RRshock_xls folder
These are the digitized Greenbook forecast from the Philadelphia Fed website. These are used to cross-check the data from FRED and the entries in RRimport.xls.