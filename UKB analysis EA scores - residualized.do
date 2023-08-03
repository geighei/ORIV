*DO-FILE TO CHECK THE PREDICTIVE POWER OF THE VARIOUS POLYGENIC SCORES FOR EA IN THE SIBLING SAMPLE OF THE UKB
* NOW RESIDUALIZING THE OUTCOME
*HANS, 5 FEBRUARY 2021

capture log close
cd "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\"
*use "Birth rank & Genes\Analysis\Input\PGS_ldpred_UKB_EA_nosibsrel_splitsample_scores.dta", clear
*use "GEIGHEI\projects\PGS ranking\Analysis\Input\PGS_ldpred_plink_EA_height.dta", clear
use "GEIGHEI\projects\PGS ranking\Analysis\Input\PGS_ldpred_plink_EA_height_cvd.dta", clear
merge 1:1 ID using "GEIGHEI\projects\Siblings\Output\Relatedness_to_siblings_UKB.dta", gen(_merge2)
merge 1:1 ID using "Birth rank & Genes\Analysis\Input\ancestry.dta", gen(_merge3)
//Check for those who withrew consent by August 2020
merge 1:1 ID using "Birth rank & Genes\Analysis\Input\w41382_20210201_withdrew_consent.dta", gen(consent)
drop if consent>1
merge 1:1 ID using "Birth rank & Genes\Analysis\Input\w41382_20200820_withdrew_consent.dta", gen(consent2)
drop if consent2>1
*(98 observations deleted)
count 
// Keep Europeans 
keep if european==1

estimates drop _all

*KEEP SIBLING SAMPLE
keep if famid ~= .
bys famid: egen count = count(ea_ukb_23me_ld)
drop if count < 2
drop if ea_ukb_23me_ld == .

keep if EA_new ~= .
keep if height ~= .

drop EA
rename EA_new EA

*******************************************************************************
* REGULAR OLS
*******************************************************************************

*RESIDUALIZE THE OUTCOME
global controls sex i.YoB c.sex#i.YoB i.MoB e_PC*
reg EA $controls, cluster(famid)
predict EAres, res
replace EA=EAres

*STANDARDIZE VARIABLES
sum EA
replace EA = (EA - r(mean))/r(sd)
sum ea_ukb_23me_ld
replace ea_ukb_23me_ld = (ea_ukb_23me_ld-r(mean))/r(sd)
sum ea_ukb_ld
replace ea_ukb_ld = (ea_ukb_ld-r(mean))/r(sd)
sum ea_23me_ld
replace ea_23me_ld = (ea_23me_ld - r(mean))/r(sd)
sum ea_ukb_ld_0
replace ea_ukb_ld_0 = (ea_ukb_ld_0 - r(mean))/r(sd)
sum ea_ukb_ld_1
replace ea_ukb_ld_1 = (ea_ukb_ld_1 - r(mean))/r(sd)

*THIS IS USING THE INDIVIDUAL SCORES AS REGRESSORS
reg EA ea_ukb_ld, cluster(famid)
estimates store e_ukb
*local beta_Gstar = _b[ea_ukb_ld]
*local sdindv_ukb = sqrt((`beta_G'^2)/(`beta_Gstar'^2)-1)
nlcom _b[ea_ukb_ld]^2
sca r2_ukb =  r(b)[1,1]
sca r2_seukb = sqrt(r(V)[1,1])
sca coefbukb = _b[ea_ukb_ld]

reg EA ea_23me_ld, cluster(famid)
estimates store e_23me
*local beta_Gstar = _b[ea_23me_ld]
*local sdindv_23me = sqrt((`beta_G'^2)/(`beta_Gstar'^2)-1)
sca coefb23me = _b[ea_23me_ld]
nlcom _b[ea_23me_ld]^2
sca r2_23me =  r(b)[1,1]
sca r2_se23me = sqrt(r(V)[1,1])

reg EA ea_ukb_ld_0, cluster(famid)
estimates store e_ukb0
nlcom _b[ea_ukb_ld_0]^2
sca r2_ukb0 = r(b)[1,1]
sca r2_seukb0 = sqrt(r(V)[1,1])
sca coefukb0 = _b[ea_ukb_ld_0]

reg EA ea_ukb_ld_1, cluster(famid)
estimates store e_ukb1
nlcom _b[ea_ukb_ld_1]^2
sca r2_ukb1 = r(b)[1,1]
sca r2_seukb1 = sqrt(r(V)[1,1])
sca coefukb1 = _b[ea_ukb_ld_1]

*THIS IS THE META-ANALYSIS SCORE
reg EA ea_ukb_23me_ld, cluster(famid)
estimates store e_meta
*local beta_Gstar = _b[ea_ukb_23me_ld]
*local sdindv_meta = sqrt((`beta_G'^2)/(`beta_Gstar'^2)-1)
sca coefmeta = _b[ea_ukb_23me_ld]
nlcom _b[ea_ukb_23me_ld]^2
sca r2_meta =  r(b)[1,1]
sca r2_semeta = sqrt(r(V)[1,1])

*THIS IS THE PGI REPOSITORY CORRECTION PROCEDURE (BECKER ET AL.)
sca rho = sqrt(0.155/e(r2))
sca coefbecker=rho*_b[ea_ukb_23me_ld]
sca sebecker=rho*_se[ea_ukb_23me_ld]

*COMPUTE THE CORRELATION BETWEEN THE SCORES TO DETERMINE THE SCALING FACTOR
reg ea_23me_ld ea_ukb_ld
sca corrscores = _b[ea_ukb_ld]
reg ea_ukb_ld_0 ea_ukb_ld_1
sca corrscoresukb = _b[ea_ukb_ld_1]

*THIS IS APPLYING ORIV USING THE 23ANDME AND UKB AS IVS FOR EACH OTHER
preserve
qui replace ea_ukb_ld = ea_ukb_ld/sqrt(corrscores)
qui replace ea_23me_ld = ea_23me_ld/sqrt(corrscores)
expand 2, generate(replicant)
gen mainVar = ea_ukb_ld if replicant == 0
replace mainVar = ea_23me_ld if replicant == 1
gen instrument = ea_23me_ld if replicant == 0
replace instrument = ea_ukb_ld if replicant == 1
reghdfe mainVar instrument, absorb(replicant) cluster(famid ID)
testparm instrument
ivreghdfe EA (mainVar=instrument), absorb(replicant) cluster(famid ID)
estimates store e_iv
nlcom _b[mainVar]^2
sca r2_iv =  r(b)[1,1]
sca r2_seiv = sqrt(r(V)[1,1])
restore

*THIS IS APPLYING REGULAR IV TWO TIMES FOR UKB AND 23ME SCORES
preserve
qui replace ea_ukb_ld = ea_ukb_ld/sqrt(corrscores)
qui replace ea_23me_ld = ea_23me_ld/sqrt(corrscores)
ivregress 2sls EA (ea_ukb_ld = ea_23me_ld), cluster(famid)
estimates store e_iv1

ivregress 2sls EA (ea_23me_ld = ea_ukb_ld), cluster(famid)
estimates store e_iv2
restore

*THIS IS APPLYING ORIV USING THE SPLIT SAMPLE UKB SCORES AS IVS FOR EACH OTHER
preserve
qui replace ea_ukb_ld_0 = ea_ukb_ld_0/sqrt(corrscoresukb)
qui replace ea_ukb_ld_1 = ea_ukb_ld_1/sqrt(corrscoresukb)
expand 2, generate(replicant)
gen mainVar = ea_ukb_ld_0 if replicant == 0
replace mainVar = ea_ukb_ld_1 if replicant == 1
gen instrument = ea_ukb_ld_1 if replicant == 0
replace instrument = ea_ukb_ld_0 if replicant == 1
reghdfe mainVar instrument, absorb(replicant) cluster(famid ID)
testparm instrument
ivreghdfe EA (mainVar=instrument), absorb(replicant) cluster(famid ID)
estimates store e_ivukb
nlcom _b[mainVar]^2
sca r2_ivukb =  r(b)[1,1]
sca r2_seivukb = sqrt(r(V)[1,1])
restore

preserve
qui replace ea_ukb_ld_0 = ea_ukb_ld_0/sqrt(corrscoresukb)
qui replace ea_ukb_ld_1 = ea_ukb_ld_1/sqrt(corrscoresukb)
ivregress 2sls EA (ea_ukb_ld_0 = ea_ukb_ld_1), cluster(famid)
estimates store e_ivukb1

ivregress 2sls EA (ea_ukb_ld_1 = ea_ukb_ld_0), cluster(famid)
estimates store e_ivukb2
restore

set logtype text
log using "GEIGHEI\projects\Measurement error\tables-residualized.txt", replace
di "BETWEEN-FAMILY ESTIMATES"
estimates table e*,  drop(_cons) b(%9.3f) se(%9.3f) stats(r2 N)
scalar list r2_ukb r2_seukb r2_23me r2_se23me r2_ukb0 r2_seukb0 r2_ukb1 r2_seukb1 r2_meta r2_semeta r2_iv r2_seiv r2_ivukb r2_seivukb
scalar list coefbecker sebecker
log close

*******************************************************************************
* SIBLING FE
*******************************************************************************
xtset famid

*THIS IS USING THE INDIVIDUAL SCORES AS REGRESSORS
xtreg EA ea_ukb_ld, fe cluster(famid)
estimates store fe_ukb
nlcom _b[ea_ukb_ld]^2
sca r2_ukb =  r(b)[1,1]
sca r2_seukb = sqrt(r(V)[1,1])
xtreg EA ea_23me_ld, fe cluster(famid)
estimates store fe_23me
nlcom _b[ea_23me_ld]^2
sca r2_23me =  r(b)[1,1]
sca r2_se23me = sqrt(r(V)[1,1])

*THIS IS THE META-ANALYSIS SCORE
xtreg EA ea_ukb_23me_ld, fe cluster(famid)
estimates store fe_meta
nlcom _b[ea_ukb_23me_ld]^2
sca r2_meta =  r(b)[1,1]
sca r2_semeta = sqrt(r(V)[1,1])

*THIS IS APPLYING ORIV WITH USING THE 23ANDME AND UKB AS IVS FOR EACH OTHER
preserve
qui replace ea_ukb_ld = ea_ukb_ld/sqrt(corrscores)
qui replace ea_23me_ld = ea_23me_ld/sqrt(corrscores)
expand 2, generate(replicant)
egen newid = concat(famid replicant)
gen long newidreal = real(newid)
xtset newidreal
gen mainVar = ea_ukb_ld if replicant == 0
replace mainVar = ea_23me_ld if replicant == 1
gen instrument = ea_23me_ld if replicant == 0
replace instrument = ea_ukb_ld if replicant == 1
reghdfe mainVar instrument, absorb(newidreal) cluster(famid ID)
testparm instrument
*forvalues x=0/1 {
*	gen constant`x' = replicant == `x'
*	}
*xtivreg EA (mainVar = instrument) constant*, fe vce(cluster famid)
ivreghdfe EA (mainVar=instrument), absorb(newidreal) cluster(famid ID)
estimates store fe_iv
nlcom (_b[mainVar])^2
sca r2_iv =  r(b)[1,1]
sca r2_seiv = sqrt(r(V)[1,1])
restore

*THIS IS APPLYING ORIV WITH USING THE SPLIT SAMPLE UKB SCORES AS IVS FOR EACH OTHER
preserve
qui replace ea_ukb_ld_0 = ea_ukb_ld_0/sqrt(corrscoresukb)
qui replace ea_ukb_ld_1 = ea_ukb_ld_1/sqrt(corrscoresukb)
expand 2, generate(replicant)
egen newid = concat(famid replicant)
gen long newidreal = real(newid)
xtset newidreal
gen mainVar = ea_ukb_ld_0 if replicant == 0
replace mainVar = ea_ukb_ld_1 if replicant == 1
gen instrument = ea_ukb_ld_1 if replicant == 0
replace instrument = ea_ukb_ld_0 if replicant == 1
reghdfe mainVar instrument, absorb(newidreal) cluster(famid ID)
testparm instrument
ivreghdfe EA (mainVar=instrument), absorb(newidreal) cluster(famid ID)
estimates store fe_ivukb
nlcom _b[mainVar]^2
sca r2_ivukb =  r(b)[1,1]
sca r2_seivukb = sqrt(r(V)[1,1])
restore

log using "GEIGHEI\projects\Measurement error\tables-residualized.txt", append
di "WITHIN-FAMILY ESTIMATES"
estimates table fe*,  drop(_cons) b(%9.3f) se(%9.3f) stats(r2 N)
scalar list r2_ukb r2_seukb r2_23me r2_se23me r2_meta r2_semeta r2_iv r2_seiv r2_ivukb r2_seivukb
log close

*****************************************************************************************************************************************************************
* THIS IS MAKING THE FIGURE IN THE PAPER - WHERE I MANUALLY COPIED THE HERITABILITY AND STANDARD ERRORS FROM THE BETWEEN-FAMILY RESULTS IN A DATASET
*****************************************************************************************************************************************************************

use "GEIGHEI\projects\Measurement error\matrix-coefplot-revision.dta", clear
mkmat h2 lb ub, matrix(C)
matrix rownames C = "OLS 23andMe" "OLS UKB" "OLS Meta-analysis" "Split-sample ORIV" "2-sample ORIV" "PGI-RC (default)" "PGI-RC (GREML unc.)"
matrix transpose = C'
coefplot matrix(transpose), ci((2 3)) scheme(s1mono) xtitle("Estimated heritability") grid(none) xlabel(, format(%03.2f))
graph save "GEIGHEI\projects\Measurement error\eacoefplot.gph", replace
graph export "GEIGHEI\projects\Measurement error\eacoefplot.pdf", replace



