
//make sure this lines up with the instruments used in "Main" specification
local insts L.rf L.ump L.ebp L.tsp L2.cpi_rolling 

//uncomment "informative factors" if you have that data
local specs Main NoDrop20 FF3Industry LinearCons WithConsSigma1 WithConsSigma5 MktOnly FF3Only AltBAAS AltCAPE LagCons Shadow VaryingBetas NDOnly NoCOVID //InformativeFactors

local nameMain RF UMP EBP TSP CPI_Rolling
local nameFF3Industry RF UMP EBP TSP CPI_Rolling
local nameLinearCons RF UMP EBP TSP CPI_Rolling
local nameWithConsSigma1 RF UMP EBP TSP CPI_Rolling
local nameWithConsSigma5 RF UMP EBP TSP CPI_Rolling
local nameMktOnly RF UMP EBP TSP CPI_Rolling
local nameFF3Only RF UMP EBP TSP CPI_Rolling
local nameAltBAAS RF UMP TSP CPI_Rolling BAAS
local nameAltCAPE RF UMP EBP CAPE TSP CPI_Rolling
local nameLagCons RF UMP EBP TSP CPI_Rolling LagCons
local nameShadow RF UMP EBP TSP CPI_Rolling ShadowSpread
local nameVaryingBetas RF UMP EBP TSP CPI_Rolling
local nameNDOnly RF UMP EBP TSP CPI_Rolling
local nameRidge RF UMP EBP TSP CPI_Rolling
local nameNoCOVID RF UMP EBP TSP CPI_Rolling
local nameInformativeFactors RF UMP EBP TSP CPI_Rolling
local nameNoDrop20 RF UMP EBP TSP CPI_Rolling

tempfile temp

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

gen Date = date(date,"DMY")
//gen Date = date(date,"YMD")
drop if Date == .
gen month = mofd(Date)
drop date Date


format month %tm
order month
sort month

mmerge month using `temp', unmatched(both)

sort month
drop _merge

tsset month


local inames
local i = 0

foreach ins in `insts' {
	
	local i = `i'+1
	local nm = regexr("`ins'","\.","")
	
	local inames `inames' `nm'
	su `ins'
	gen ins_`nm' = `ins' - r(mean)
	local vname : word `i' of `nameMain'
	label var ins_`nm' `vname'
	egen z_`nm' = std(`ins') if month >= ym(1973,3)
}


label var exprfRateReal "Exp. Real Tsy. Ret."
label var zbRateReal "Real Zero-Beta Rate"
label var cons_g "Cons. Gr."


label var portRetReal "Real Z.B. Port. Ret"



save `temp', replace

import delimited "../Output/ses_Main.csv", clear
mkmat var1, matrix(b_test)
mkmat var2*, matrix(V_test)

use `temp', clear

gen spread = portRetReal - exprfRateReal
label var spread "Convenience Spread"

//changed to portfolio return instead of spread to be consistent with
//change in definition of gamma in paper
reg portRetNom ins_*, robust

local n_ins = e(df_m)


local names : colnames e(b)
local names : subinstr local names  "_cons" ""

local eN = e(N)

eststo simple_preg
estadd scalar F_p = Ftail(e(df_m), e(df_r), e(F))


matrix rownames b_test = _cons `names'
matrix rownames V_test = _cons `names'
matrix colnames V_test = _cons `names'

matrix b_test = b_test'
matrix list b_test
matrix list V_test

matrix bNC = b_test[1,2...]
matrix VNC = V_test[2...,2...]

matrix temp = bNC * invsym(VNC) * bNC'

local gmm_wald = temp[1,1]
local pval = 1 - chi2(`n_ins',`gmm_wald')

ereturn post b_test V_test, obs(`eN')

eststo gmm_preg
estadd scalar F = `gmm_wald'
estadd scalar F_p = `pval'

esttab gmm_preg simple_preg, scalars("F Wald/F" "F_p p-value" "rmse RMSE") sfmt("%6.0g") se label obslast depvars nostar mtitles("GMM" "OLS (inf.)") 

esttab gmm_preg simple_preg using ../Output/GMMReg.csv, scalars("F Wald/F" "F_p p-value" "rmse RMSE") sfmt("%6.0g") se label obslast depvars nostar mtitles("GMM" "OLS (inf.)") replace

esttab gmm_preg simple_preg using ../Output/GMMReg.tex, scalars("F Wald/F" "F_p p-value" "rmse RMSE") sfmt("%6.0g") se label obslast depvars nostar mtitles("GMM" "OLS (inf.)") replace





local ests
foreach spec in `specs' {
	disp "spec: `spec'"
	import delimited "../Output/ses_`spec'.csv", clear
	mkmat var1, matrix(b_`spec')
	mkmat var2*, matrix(V_`spec')
	
	matrix rownames b_`spec' = _cons `name`spec''
	matrix rownames V_`spec' = _cons `name`spec''
	matrix colnames V_`spec' = _cons `name`spec''
	
	matrix b_`spec' = b_`spec''
	matrix list b_`spec'
	if ~regexm("`spec'","Ridge") {
		ereturn post b_`spec' V_`spec'
	}
	else {
		ereturn post b_`spec'
	}
	disp "post"
	eststo gmm_`spec'
	disp "stored"
	local ests `ests' gmm_`spec'
}

esttab `ests', se label depvars nostar mtitles(`specs') noobs
	
esttab `ests' using ../Output/Robustness.tex, se label depvars nostar mtitles(`specs') noobs replace
	
use `temp', clear


import delimited ../Output/RRData_Main.csv, clear
gen dt = date(date,"YMD")
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
gen dt = date(date,"YMD")
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
