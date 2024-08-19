
//uncomment one of these, depends on your machine type
local dateType DMY //Works on our Mac
//local dateType YMD //Works on our Windows

//make sure this lines up with the instruments used in "Main" specification
local insts L.rf L.ump L.ebp L.tsp L2.cpi_rolling 

//uncomment "informative factors" if you have that data
local allspecs AltData AltFactors AltInstruments AltPorts

local AltData Main NoDrop20 FF3Industry NDOnly Ridge
local AltFactors Main LinearCons WithConsSigma1 WithConsSigma5 WithConsSigma10 MktOnly FF3Only VaryingBetas InfFactors
local AltInstruments Main AltBAAS AltCAPE LagCons Shadow AltCPI AltSahm
local AltPorts Main CovSample CovPCA CovI

//old
//local specs Main NoDrop20 FF3Industry LinearCons WithConsSigma1 WithConsSigma5 MktOnly FF3Only AltBAAS AltCAPE LagCons Shadow VaryingBetas NDOnly NoCOVID //InformativeFactors

local nameMain RF UMP EBP TSP CPI_Rolling
local nameFF3Industry RF UMP EBP TSP CPI_Rolling
local nameLinearCons RF UMP EBP TSP CPI_Rolling
local nameWithConsSigma1 RF UMP EBP TSP CPI_Rolling
local nameWithConsSigma5 RF UMP EBP TSP CPI_Rolling
local nameWithConsSigma10 RF UMP EBP TSP CPI_Rolling
local nameMktOnly RF UMP EBP TSP CPI_Rolling
local nameFF3Only RF UMP EBP TSP CPI_Rolling
local nameAltBAAS RF UMP TSP CPI_Rolling BAAS
local nameAltCAPE RF UMP EBP CAPE TSP CPI_Rolling
local nameAltSahm RF SAHM_RT EBP TSP CPI_Rolling
local nameAltCPI RF UMP EBP TSP LagCPI
local nameLagCons RF UMP EBP TSP CPI_Rolling LagCons
local nameShadow RF UMP EBP TSP CPI_Rolling ShadowSpread
local nameVaryingBetas RF UMP EBP TSP CPI_Rolling
local nameNDOnly RF UMP EBP TSP CPI_Rolling
local nameRidge RF UMP EBP TSP CPI_Rolling
local nameNoCOVID RF UMP EBP TSP CPI_Rolling
local nameInfFactors RF UMP EBP TSP CPI_Rolling
local nameNoDrop20 RF UMP EBP TSP CPI_Rolling
local nameCovSample RF UMP EBP TSP CPI_Rolling
local nameCovPCA RF UMP EBP TSP CPI_Rolling
local nameCovI RF UMP EBP TSP CPI_Rolling

tempfile temp temp2

import delimited "../Input/Instruments.csv", clear



gen Date = date(date,"YMD")
gen month = mofd(Date)
drop date Date
format month %tm
order month
sort month


save `temp', replace




import delimited "../Output/zero_beta_rateMain.csv", clear

rename var1 date
rename var2 exprfRateReal
rename var3 zbRateReal
rename var4 portRetReal
rename var5 cons_g
rename var6 zbRateNom
rename var7 portRetNom

gen Date = date(date,"`dateType'")


drop if Date == .
gen month = mofd(Date)
drop date Date


format month %tm
order month
sort month

mmerge month using `temp', unmatched(both)

sort month
drop _merge


save `temp', replace

import delimited "../Output/zero_beta_rateMainNC.csv", clear

keep var1 var7 var6
rename var1 date
rename var6 zbRateNomNC
rename var7 portRetNomNC

gen Date = date(date,"`dateType'")


drop if Date == .
gen month = mofd(Date)
drop date Date
sort month

save `temp2', replace

use `temp', clear
mmerge month using `temp2', unmatched(master)
sort month
drop _merge

tsset month


local inames
local i = 0

foreach ins in `insts' {
	
	local i = `i'+1
	local nm = regexr("`ins'","\.","")
	
	local inames `inames' `nm'
	su `ins' if month >= ym(1973,3)
	gen ins_`nm' = `ins' - r(mean)
	
	su `ins' if month <= ym(2019,12) & month >= ym(1973,3)
	gen insNC_`nm' = `ins' - r(mean) if month <= ym(2019,12)
	local vname : word `i' of `nameMain'
	label var ins_`nm' `vname'
	label var insNC_`nm' `vname'
	egen z_`nm' = std(`ins') if month >= ym(1973,3)
}


label var exprfRateReal "Exp. Real Tsy. Ret."
label var zbRateReal "Real Zero-Beta Rate"
label var cons_g "Cons. Gr."


label var portRetReal "Real Z.B. Port. Ret"



save `temp', replace

//new code to merge in mkt return
import delimited "../Input/Factors_Nominal.csv", clear

keep date mkt
gen Date = date(date,"YMD")
gen month = mofd(Date)
drop date Date
format month %tm
order month
sort month

mmerge month using `temp', unmatched(using)
drop _merge
sort month

save `temp', replace

import delimited "../Output/ses_Main.csv", clear
mkmat var1, matrix(b_test)
mkmat var2*, matrix(V_test)

import delimited "../Output/ses_MainNC.csv", clear
mkmat var1, matrix(b_test2)
mkmat var2*, matrix(V_test2)


insheet using "../Output/BootstrapPvals.csv", clear

levelsof var2 if var1 == "Wald", local(pWald)
levelsof var2 if var1 == "WaldNC", local(pWaldNC)
levelsof var2 if var1 == "F-eff", local(pFeff)
levelsof var2 if var1 == "F-effNC", local(pFeffNC)

use `temp', clear

//changed to portfolio return instead of spread to be consistent with
//change in definition of gamma in paper
//gen spread = portRetReal - exprfRateReal
//label var spread "Spread"


reg portRetNom ins_*, robust

//check construction
predict zbrateTest, xb
gen errs = zbRateNom - zbrateTest
su errs

gen resids = portRetNom - zbRateNom

su resids

local residSD = `r(sd)'
su portRetNom

local portRetSD = `r(sd)'

local gmmR2 = 1 - (`residSD'*`residSD') / (`portRetSD'*`portRetSD')

gen residsNC = portRetNomNC - zbRateNomNC

su residsNC

local residSDNC = `r(sd)'
su portRetNomNC

local portRetSDNC = `r(sd)'

local gmmR2NC = 1 - (`residSDNC'*`residSDNC') / (`portRetSDNC'*`portRetSDNC')

local n_ins = e(df_m)


local names : colnames e(b)
local names : subinstr local names  "_cons" ""

local eN = e(N)

eststo simple_preg
estadd scalar F_p = Ftail(e(df_m), e(df_r), e(F))

//need to align variable names
rename ins_* insC_*
rename insNC_* ins_*
reg portRetNomNC ins_* if month <= ym(2019,12), robust
eststo simple_pregNC
estadd scalar F_p = Ftail(e(df_m), e(df_r), e(F))
local eN2 = e(N)

rename ins_* insNC_*
rename insC_* ins_*

//setup with COVID gmm results
matrix rownames b_test = _cons `names'
matrix rownames V_test = _cons `names'
matrix colnames V_test = _cons `names'

matrix b_test = b_test'
matrix list b_test
matrix list V_test

matrix nullh = J(1,`n_ins',0)
matrix colnames nullh = `names'

matrix list nullh


matrix nullh[1,colnumb(nullh,"ins_Lrf")] = 1

matrix list nullh

matrix bNC = b_test[1,2...] - nullh
matrix VNC = V_test[2...,2...]

matrix temp = bNC * invsym(VNC) * bNC'

local gmm_wald = temp[1,1]
local pval = 1 - chi2(`n_ins',`gmm_wald')

ereturn post b_test V_test, obs(`eN')

eststo gmm_preg
estadd scalar F = `gmm_wald'
estadd scalar F_p = `pval'
estadd scalar F_pb = `pWald'
estadd scalar r2 = `gmmR2'


//no covid GMM
matrix rownames b_test2 = _cons `names'
matrix rownames V_test2 = _cons `names'
matrix colnames V_test2 = _cons `names'

matrix b_test2 = b_test2'
matrix list b_test2
matrix list V_test2

matrix bNC2 = b_test2[1,2...] - nullh
matrix VNC2 = V_test2[2...,2...]

matrix temp = bNC2 * invsym(VNC2) * bNC2'

local gmm_wald2 = temp[1,1]
local pval2 = 1 - chi2(`n_ins',`gmm_wald2')


ereturn post b_test2 V_test2, obs(`eN2')
eststo gmm_preg2
estadd scalar F = `gmm_wald2'
estadd scalar F_p = `pval2'
estadd scalar F_pb = `pWaldNC'
estadd scalar r2 = `gmmR2NC'

esttab gmm_preg simple_preg gmm_preg2 simple_pregNC, scalars("F Wald/F" "F_p p-value (asymptotic)" "F_pb p-value (bootstrap)" "rmse RMSE") sfmt("%6.0g") se label obslast depvars nostar mtitles("GMM" "OLS (inf.)" "GMM" "OLS (inf.)" ) r2 mgroups("Including 2020" "Excluding 2020", pattern(1 0 1 0) span) nonumber 

esttab gmm_preg simple_preg gmm_preg2 simple_pregNC using ../Output/GMMReg.csv, scalars("F Wald/F" "F_p p-value (asymptotic)" "F_pb p-value (bootstrap)" "rmse RMSE") sfmt("%6.0g") se label obslast depvars nostar mtitles("GMM" "OLS (inf.)" "GMM" "OLS (inf.)" ) r2 replace mgroups("Including 2020" "Excluding 2020", pattern(1 0 1 0) span) nonumber 

esttab gmm_preg simple_preg gmm_preg2 simple_pregNC using ../Output/GMMReg.tex, scalars("F Wald/F" "F_p p-value (asymptotic)" "F_pb p-value (bootstrap)" "rmse RMSE") sfmt("%6.0g") se label obslast depvars nostar mtitles("GMM" "OLS (inf.)" "GMM" "OLS (inf.)" ) r2 replace mgroups("Including 2020" "Excluding 2020", pattern(1 0 1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber 




corr ins_*, cov
matrix zcov = r(C)


//new code for table with consumption
gen cons_gr_ann = 12*cons_g
gen covid = month > ym(2019,12)

reg cons_gr_ann ins_*, robust
eststo simple_cons

local T = e(N)
matrix best = e(b)
matrix temp = best[1,1..`n_ins'] * zcov * best[1,1..`n_ins']' 
matrix temp2 = e(V)
matrix temp3 = temp2[1..`n_ins',1..`n_ins'] * zcov
local trv = trace(temp3)
local effF = temp[1,1] / `trv'

estadd scalar F_eff = `effF'
estadd scalar F_effp = `pFeff'
//estadd scalar F_p = Ftail(e(df_m), e(df_r), e(F))


//feasible GLS with covid predictor
//not used in paper because it ends up being the same
//as just dropping covid
predict cresids, r
gen cresids2 = log(cresids * cresids)
reg cresids2 ins_* covid
predict pcresids, xb
gen iaw = exp(-pcresids)

reg cons_gr_ann ins_* [aweight=iaw], robust
eststo cons_fgls
estadd scalar F_p = Ftail(e(df_m), e(df_r), e(F))

//outlier robust regression with bi-weight
robreg m cons_gr_ann ins_*, biweight ftest
eststo cons_bi

matrix best = e(b)
matrix temp = best[1,1..`n_ins'] * zcov * best[1,1..`n_ins']' 
matrix temp2 = e(V)
matrix temp3 = temp2[1..`n_ins',1..`n_ins'] * zcov
local trv = trace(temp3)
local effF = temp[1,1] / `trv'

estadd scalar F_eff = `effF'

//display this for purpose of making matlab code match stata
local tuning = e(k)
disp "tuning constant: `tuning'"

//estadd scalar F = e(chi2)
estadd scalar F_p = e(p)
estadd scalar r2 = e(r2_p)

//mm-robust regression
robreg mm cons_gr_ann ins_*, ftest
eststo cons_mm
//estadd scalar F = e(chi2)
estadd scalar F_p = e(p)
estadd scalar r2 = e(r2_p)

matrix best = e(b)
matrix temp = best[1,1..`n_ins'] * zcov * best[1,1..`n_ins']' 
matrix temp2 = e(V)
matrix temp3 = temp2[1..`n_ins',1..`n_ins'] * zcov
local trv = trace(temp3)
local effF = temp[1,1] / `trv'

estadd scalar F_eff = `effF'

//adding mkt specification for comparison
reg mkt ins_*, robust
eststo simple_mkt
estadd scalar F_p = Ftail(e(df_m), e(df_r), e(F))

matrix best = e(b)
matrix temp = best[1,1..`n_ins'] * zcov * best[1,1..`n_ins']' 
matrix temp2 = e(V)
matrix temp3 = temp2[1..`n_ins',1..`n_ins'] * zcov
local trv = trace(temp3)
local effF = temp[1,1] / `trv'

estadd scalar F_eff = `effF'

//pre-covid
rename ins_* insC_*
rename insNC_* ins_*

corr ins_*, cov
matrix zcovNC = r(C)

reg cons_gr_ann ins_* if covid==0, robust
eststo simple_cons_pre
estadd scalar F_p = Ftail(e(df_m), e(df_r), e(F))

matrix best = e(b)
matrix temp = best[1,1..`n_ins'] * zcovNC * best[1,1..`n_ins']' 
matrix temp2 = e(V)
matrix temp3 = temp2[1..`n_ins',1..`n_ins'] * zcovNC
local trv = trace(temp3)
local effF = temp[1,1] / `trv'

estadd scalar F_eff = `effF'
estadd scalar F_effp = `pFeffNC'


//outlier robust regression with bi-weight
robreg m cons_gr_ann ins_* if covid==0, biweight ftest
eststo cons_biNC

matrix best = e(b)
matrix temp = best[1,1..`n_ins'] * zcovNC * best[1,1..`n_ins']' 
matrix temp2 = e(V)
matrix temp3 = temp2[1..`n_ins',1..`n_ins'] * zcovNC
local trv = trace(temp3)
local effF = temp[1,1] / `trv'

estadd scalar F_eff = `effF'
estadd scalar F_p = e(p)
estadd scalar r2 = e(r2_p)

reg mkt ins_* if month<=ym(2019,12), robust
eststo simple_mkt2
estadd scalar F_p = Ftail(e(df_m), e(df_r), e(F))

matrix best = e(b)
matrix temp = best[1,1..`n_ins'] * zcovNC * best[1,1..`n_ins']' 
matrix temp2 = e(V)
matrix temp3 = temp2[1..`n_ins',1..`n_ins'] * zcovNC
local trv = trace(temp3)
local effF = temp[1,1] / `trv'

estadd scalar F_eff = `effF'


rename ins_* insNC_*
rename insC_* ins_*

esttab simple_cons cons_bi simple_cons_pre cons_biNC simple_mkt simple_mkt2, scalars("r2 (Pseudo) R-squared" "F_eff F (effective)" "F_effp Bootstrap p-value") sfmt("%6.0g") se label obslast depvars nostar mtitles("OLS" "Robust (bi)" " OLS" "Robust (bi)" "OLS" "OLS") unstack  eqlabels(none) mgroups("Cons. incl. 2020" "Cons. excl. 2020" "Mkt. incl. 2020" "Mkt. excl. 2020", pattern(1 0 1 0 1 1) span ) nonumber


//drop(S: scale:)
esttab simple_cons cons_bi simple_cons_pre cons_biNC simple_mkt simple_mkt2 using ../Output/ConsTable.tex, scalars("r2 (Pseudo) R-squared" "F_eff F (effective)" "F_effp Bootstrap p-value") sfmt("%6.0g") se label obslast depvars nostar mtitles("OLS" "Robust (bi)" " OLS" "Robust (bi)" "OLS" "OLS") unstack eqlabels(none) replace mgroups("Cons. incl. 2020" "Cons. excl. 2020" "Mkt. incl. 2020" "Mkt. excl. 2020", pattern(1 0 1 0 1 1) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber 


//drop(S: scale:) 

foreach specgroup in `allspecs' {
disp "group: `specgroup'"

local specs ``specgroup''

disp "specs: `specs'"

local ests
local estsNC
local grp
local grpNC
local specNames
local specNamesNC
foreach spec in `specs' {
	disp "spec: `spec'"
	
	foreach i in 1 2 { 

	local basespec `spec'
	if `i' == 2 {
		local spec `spec'NC
		disp "specNC: `spec'"
	}
	
	import delimited "../Output/zero_beta_rate`spec'.csv", clear


	rename var2 exprfRateReal`spec'
	rename var3 zbRateReal`spec'
	rename var6 zbRateNom`spec'
	rename var7 portRetNom`spec'
	rename var8 p_cons`spec'
	
	//gen spread`spec' = zbRateReal`spec' - exprfRateReal`spec'
	
	su zbRateReal`spec'
	local mean`spec' = `r(mean)'
	local sd`spec' = `r(sd)'
	
	su portRetNom`spec'
	local sd_port`spec' = `r(sd)'
	
	gen resids`spec' = portRetNom`spec' - zbRateNom`spec'
	su resids`spec'
	local sd_resid`spec' = `r(sd)'
	
	local gmmR2`spec' = 1 - (`sd_resid`spec''*`sd_resid`spec'') / (`sd_port`spec''*`sd_port`spec'')
	
	corr zbRateReal`spec' p_cons`spec'
	local rho`spec' = `r(rho)'
	
	corr exprfRateReal`spec' p_cons`spec'
	local rhoC`spec' = `r(rho)'
	
	import delimited "../Output/ses_`spec'.csv", clear
	mkmat var1, matrix(b_`spec')
	mkmat var2*, matrix(V_`spec')
	
	matrix rownames b_`spec' = _cons `name`basespec''
	matrix rownames V_`spec' = _cons `name`basespec''
	matrix colnames V_`spec' = _cons `name`basespec''
	
	matrix b_`spec' = b_`spec''
	matrix list b_`spec'
	if ~regexm("`spec'","Ridge") {
		ereturn post b_`spec' V_`spec'
	}
	else {
		ereturn post b_`spec'
	}
	
	estadd scalar AnnMean = 12*`mean`spec''
	estadd scalar AnnSD = 12*`sd`spec''
	estadd scalar Corr = `rho`spec''
	estadd scalar CorrC = `rhoC`spec''
	estadd scalar AnnPortSD = sqrt(12)*`sd_port`spec''
	estadd scalar r2 = `gmmR2`spec''
	
	disp "post"
	eststo gmm_`spec'
	disp "stored"
	
	if `i' == 1 {
		if ~regexm("`spec'","InfFactors") {
			local ests `ests' gmm_`spec'
			local specNames `specNames' `spec'
		}
	}
	else {
		local estsNC `estsNC' gmm_`spec'
		local specNamesNC `specNamesNC' `basespec'
	}
	
	}
	
	if "`basespec'" == "Main" {
		local grp `grp' 1
		local grpNC `grpNC' 1
	}
	else {
		if ~regexm("`spec'","InfFactors") {
			local grp `grp' 0
		}
		local grpNC `grpNC' 0
	}
}
disp "`grp' `grpNC'"

esttab `ests' `estsNC', scalars("AnnMean Avg. Annualized Real ZB Rate" "AnnSD SD of Annualized Real ZB Rate" "Corr Corr(Real ZB Rate, Exp. Cons. Gr.)" "CorrC Corr(Exp. Real T-bill, Exp. Cons. Gr.)") sfmt("%6.0g") se label depvars nostar mtitles(`specNames' `specNamesNC') noobs r2 order(`nameMain')  mgroups("Including 2020" "Excluding 2020", pattern(`grp' `grpNC') span ) nonumber 
//prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) alignment(D{.}{.}{-1}) page(dcolumn)	booktabs
esttab `ests' `estsNC' using ../Output/Robustness_`specgroup'.tex, scalars("AnnMean Avg. Annualized Real ZB Rate" "AnnSD SD of Annualized Real ZB Rate" "Corr Corr(Real ZB Rate, Exp. Cons. Gr.)" "CorrC Corr(Exp. Real T-bill, Exp. Cons. Gr.)") sfmt("%6.0g") se label depvars nostar mtitles(`specNames' `specNamesNC') noobs replace r2 order(`nameMain') mgroups("Including 2020" "Excluding 2020", pattern(`grp' `grpNC') prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber //alignment(D{.}{.}{-1}) page(dcolumn)
	
}
	
use `temp', clear


import delimited ../Output/RRData_Main.csv, clear
gen dt = date(date,"DMY") //used to be MDY
gen month = mofd(dt)
sort month

tsset month

gen shock_rr = L.resid_full
label var shock_rr "RR_shock"

gen DZeroBeta = 12*(F4.zbr - L2.zbr)
gen DZeroBetaN = 12*(F4.zbn - L2.zbn)
gen DTBill = 12*(F4.rf - L2.rf)


gen Dp_cons = F4.pcons - L2.pcons
label var Dp_cons "Ex. C. Gr."
label var DTBill "Real Bill"
label var DZeroBeta "Real Z.B."

gen DConvenience = DZeroBeta - DTBill

label var DConvenience Conv

reg DTBill shock_rr, robust
eststo DTBillRR
reg DZeroBeta shock_rr, robust
eststo DZeroBetaRR

reg Dp_cons shock_rr, robust
eststo Dp_consRR

reg DZeroBetaN shock_rr, robust //DConvenience
eststo DConv

local nameMain RF UMP EBP TSP CPI_Rolling
local ests DTBillRR DZeroBetaRR Dp_consRR //DConv
local i = 1
foreach var in `nameMain' {
	gen g_`var' = 12*(F4.var4_`i' - L2.var4_`i')
	reg g_`var' shock_rr if DZeroBeta != ., robust
	eststo g_`var'
	local ests `ests' g_`var'
	local i = `i' + 1
}

esttab `ests', label nostar not gaps
esttab `ests' using ../Output/RREffects.tex, not label nostar replace gaps

esttab DTBillRR DZeroBetaRR Dp_consRR, not label nostar gaps

import delimited ../Output/NSData_Main.csv, clear
gen dt = date(date,"DMY") //used to be MDY
gen month = mofd(dt)
sort month

tsset month

gen shock = L.ns_shock
label var shock "NS_shock"
gen shock_rr = .
label var shock_rr "RR_shock"

gen DZeroBeta = 12*(F4.zbr - L2.zbr)
gen DZeroBetaN = 12*(F4.zbn - L2.zbn)
gen DTBill = 12*(F4.rf - L2.rf)


gen Dp_cons = F4.pcons - L2.pcons
label var Dp_cons "Ex. C. Gr."

label var DTBill "Real Bill"
label var DZeroBeta "Real Z.B."

gen DConvenience = DZeroBeta - DTBill

label var DConvenience Conv

reg DTBill shock, robust
eststo DTBill
reg DZeroBeta shock, robust
eststo DZeroBeta

reg Dp_cons shock, robust
eststo Dp_cons

reg DZeroBetaN shock, robust //DConvenience
eststo DConv


local nameMain RF UMP EBP TSP CPI_Rolling
local ests DTBill DZeroBeta Dp_cons //DConv
local i = 1
foreach var in `nameMain' {
	gen g_`var' = 12*(F4.var4_`i' - L2.var4_`i')
	reg g_`var' shock if DZeroBeta != ., robust
	eststo g_`var'
	local ests `ests' g_`var'
	local i = `i' + 1
}

esttab `ests', label nostar not gaps
esttab `ests' using ../Output/NSEffects.tex, not label nostar replace gaps

esttab DTBillRR DZeroBetaRR Dp_consRR DTBill DZeroBeta Dp_cons, not label nostar gaps
esttab DTBillRR DZeroBetaRR Dp_consRR DTBill DZeroBeta Dp_cons using ../Output/BothEffects.tex, not label nostar replace gaps
