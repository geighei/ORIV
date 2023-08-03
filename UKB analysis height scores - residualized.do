*DO-FILE TO CHECK THE PREDICTIVE POWER OF THE VARIOUS POLYGENIC SCORES FOR HEIGHT IN THE SIBLING SAMPLE OF THE UKB
*	RESIDUALIZED THE OUTCOME
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
reg height $controls, cluster(famid)
predict heightres, res
replace height=heightres

*STANDARDIZE VARIABLES
sum height
replace height = (height - r(mean))/r(sd)
foreach i in height_ukb_ld_0 height_ukb_ld_1 height_ukb_giant_ld height_giant_ld height_ukb_ld {
    sum `i'
	replace `i' = (`i'-r(mean))/r(sd)
}

*THIS IS USING THE INDIVIDUAL SCORES AS REGRESSORS
reg height height_ukb_ld, cluster(famid)
estimates store e_ukb
nlcom _b[height_ukb_ld]^2
sca r2_ukb =  r(b)[1,1]
sca r2_seukb = sqrt(r(V)[1,1])

reg height height_giant_ld, cluster(famid)
estimates store e_giant
nlcom _b[height_giant_ld]^2
sca r2_giant =  r(b)[1,1]
sca r2_segiant = sqrt(r(V)[1,1])

*THIS IS THE META-ANALYSIS SCORE
reg height height_ukb_giant_ld, cluster(famid)
estimates store e_meta
nlcom _b[height_ukb_giant_ld]^2
sca r2_meta =  r(b)[1,1]
sca r2_semeta = sqrt(r(V)[1,1])

*THIS IS THE PGI REPOSITORY CORRECTION PROCEDURE (BECKER ET AL.)
sca rho = sqrt(0.530/e(r2))
sca coefbecker=rho*_b[height_ukb_giant_ld]
sca sebecker=rho*_se[height_ukb_giant_ld]

*COMPUTE THE CORRELATION BETWEEN THE SCORES TO DETERMINE THE SCALING FACTOR
reg height_giant_ld height_ukb_ld
sca corrscores = _b[height_ukb_ld]
reg height_ukb_ld_1 height_ukb_ld_0
sca corrscoresukb = _b[height_ukb_ld_0]

*THIS IS APPLYING ORIV WITH USING THE GIANT AND UKB AS IVS FOR EACH OTHER
preserve
qui replace height_ukb_ld = height_ukb_ld/sqrt(corrscores)
qui replace height_giant_ld = height_giant_ld/sqrt(corrscores)
expand 2, generate(replicant)
gen mainVar = height_ukb_ld if replicant == 0
replace mainVar = height_giant_ld if replicant == 1
gen instrument = height_giant_ld if replicant == 0
replace instrument = height_ukb_ld if replicant == 1
reghdfe mainVar instrument, absorb(replicant) cluster(famid ID)
testparm instrument
ivreghdfe height (mainVar=instrument), absorb(replicant) cluster(famid ID)
estimates store e_iv
nlcom (_b[mainVar])^2
sca r2_iv =  r(b)[1,1]
sca r2_seiv = sqrt(r(V)[1,1])
restore

*THIS IS APPLYING REGULAR IV TWO TIMES FOR UKB AND GIANT SCORES
preserve
qui replace height_ukb_ld = height_ukb_ld/sqrt(corrscores)
qui replace height_giant_ld = height_giant_ld/sqrt(corrscores)
ivregress 2sls height (height_ukb_ld = height_giant_ld), cluster(famid)
estimates store e_iv1

ivregress 2sls height (height_giant_ld = height_ukb_ld), cluster(famid)
estimates store e_iv2
restore

*THIS IS APPLYING ORIV WITH USING THE SPLIT SAMPLE UKB SCORES AS IVS FOR EACH OTHER
preserve
qui replace height_ukb_ld_0 = height_ukb_ld_0/sqrt(corrscoresukb)
qui replace height_ukb_ld_1 = height_ukb_ld_1/sqrt(corrscoresukb)
expand 2, generate(replicant)
gen mainVar = height_ukb_ld_0 if replicant == 0
replace mainVar = height_ukb_ld_1 if replicant == 1
gen instrument = height_ukb_ld_1 if replicant == 0
replace instrument = height_ukb_ld_0 if replicant == 1
reghdfe mainVar instrument, absorb(replicant) cluster(famid ID)
testparm instrument
ivreghdfe height (mainVar=instrument), absorb(replicant) cluster(famid ID)
estimates store e_ivukb
nlcom _b[mainVar]^2
sca r2_ivukb =  r(b)[1,1]
sca r2_seivukb = sqrt(r(V)[1,1])
restore

*THIS IS APPLYING REGULAR IV TWO TIMES FOR UKB SPLITSAMPLE
preserve
qui replace height_ukb_ld_0 = height_ukb_ld_0/sqrt(corrscoresukb)
qui replace height_ukb_ld_1 = height_ukb_ld_1/sqrt(corrscoresukb)
ivregress 2sls height (height_ukb_ld_0 = height_ukb_ld_1), cluster(famid)
estimates store e_ivukb1

ivregress 2sls height (height_ukb_ld_1 = height_ukb_ld_0), cluster(famid)
estimates store e_ivukb2
restore

set logtype text
log using "GEIGHEI\projects\Measurement error\tablesheight-residualized.txt", replace
di "BETWEEN-FAMILY ESTIMATES"
estimates table e*,  drop(_cons) b(%9.3f) se(%9.3f) stats(r2 N)
scalar list r2_ukb r2_seukb r2_giant r2_segiant r2_meta r2_semeta r2_iv r2_seiv r2_ivukb r2_seivukb
scalar list coefbecker sebecker
log close

*******************************************************************************
* SIBLING FE
*******************************************************************************
xtset famid

*THIS IS USING THE INDIVIDUAL SCORES AS REGRESSORS
xtreg height height_ukb_ld, fe cluster(famid)
estimates store fe_ukb
nlcom _b[height_ukb_ld]^2
sca r2_ukb =  r(b)[1,1]
sca r2_seukb = sqrt(r(V)[1,1])
xtreg height height_giant_ld, fe cluster(famid)
estimates store fe_giant
nlcom _b[height_giant_ld]^2
sca r2_giant =  r(b)[1,1]
sca r2_segiant = sqrt(r(V)[1,1])

*THIS IS THE META-ANALYSIS SCORE
xtreg height height_ukb_giant_ld, fe cluster(famid)
estimates store fe_meta
nlcom _b[height_ukb_giant_ld]^2
sca r2_meta =  r(b)[1,1]
sca r2_semeta = sqrt(r(V)[1,1])

*THIS IS APPLYING ORIV WITH USING THE GIANT AND UKB AS IVS FOR EACH OTHER
preserve
qui replace height_ukb_ld = height_ukb_ld/sqrt(corrscores)
qui replace height_giant_ld = height_giant_ld/sqrt(corrscores)
expand 2, generate(replicant)
egen newid = concat(famid replicant)
gen long newidreal = real(newid)
xtset newidreal
gen mainVar = height_ukb_ld if replicant == 0
replace mainVar = height_giant_ld if replicant == 1
gen instrument = height_giant_ld if replicant == 0
replace instrument = height_ukb_ld if replicant == 1
reghdfe mainVar instrument, absorb(newidreal) cluster(famid ID)
testparm instrument
ivreghdfe height (mainVar=instrument), absorb(newidreal) cluster(famid ID)
estimates store fe_iv
nlcom (_b[mainVar])^2
sca r2_iv =  r(b)[1,1]
sca r2_seiv = sqrt(r(V)[1,1])
restore

*THIS IS APPLYING ORIV WITH USING THE SPLIT SAMPLE UKB SCORES AS IVS FOR EACH OTHER
preserve
qui replace height_ukb_ld_0 = height_ukb_ld_0/sqrt(corrscoresukb)
qui replace height_ukb_ld_1 = height_ukb_ld_1/sqrt(corrscoresukb)
expand 2, generate(replicant)
egen newid = concat(famid replicant)
gen long newidreal = real(newid)
xtset newidreal
gen mainVar = height_ukb_ld_0 if replicant == 0
replace mainVar = height_ukb_ld_1 if replicant == 1
gen instrument = height_ukb_ld_1 if replicant == 0
replace instrument = height_ukb_ld_0 if replicant == 1
reghdfe mainVar instrument, absorb(newidreal) cluster(famid ID)
testparm instrument
ivreghdfe height (mainVar=instrument), absorb(newidreal) cluster(famid ID)
estimates store fe_ivukb
nlcom _b[mainVar]^2
sca r2_ivukb =  r(b)[1,1]
sca r2_seivukb = sqrt(r(V)[1,1])
restore

set logtype text
log using "GEIGHEI\projects\Measurement error\tablesheight-residualized.txt", append
di "WITHIN-FAMILY ESTIMATES"
estimates table fe*,  drop(_cons) b(%9.3f) se(%9.3f) stats(N r2)
scalar list r2_ukb r2_seukb r2_giant r2_segiant r2_meta r2_semeta r2_iv r2_seiv r2_ivukb r2_seivukb
log close

*****************************************************************************************************************************************************************
* THIS IS MAKING THE FIGURE IN THE PAPER - WHERE I MANUALLY COPIED THE HERITABILITY AND STANDARD ERRORS FROM THE BETWEEN-FAMILY RESULTS IN A DATASET
*****************************************************************************************************************************************************************

use "GEIGHEI\projects\Measurement error\matrix-coefplot-height-revision.dta", clear
mkmat h2 lb ub, matrix(C)
matrix rownames C = "OLS 23andMe" "OLS UKB" "OLS Meta-analysis" "Split-sample ORIV" "2-sample ORIV" "PGI-RC (default)" "PGI-RC (GREML unc.)"
matrix transpose = C'
coefplot matrix(transpose), ci((2 3)) scheme(s1mono) xtitle("Estimated heritability") grid(none) xlabel(, format(%03.2f))
graph save "GEIGHEI\projects\Measurement error\heightcoefplot.gph", replace
graph export "GEIGHEI\projects\Measurement error\heightcoefplot.pdf", replace
