*DO-FILE TO CHECK THE PREDICTIVE POWER OF THE VARIOUS POLYGENIC SCORES FOR DIASTOLIC BLOOD PRESSURE IN THE SIBLING SAMPLE OF THE UKB
* NOW RESIDUALIZING THE OUTCOME
*HANS, 31 AUGUST 2022

capture log close
cd "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\"
*use "Birth rank & Genes\Analysis\Input\PGS_ldpred_UKB_EA_nosibsrel_splitsample_scores.dta", clear
*use "GEIGHEI\projects\PGS ranking\Analysis\Input\PGS_ldpred_plink_EA_height.dta", clear
use "GEIGHEI\projects\PGS ranking\Analysis\Input\PGS_ldpred_plink_EA_height_cvd_bmi_dbp.dta", clear
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

bys famid: egen countdbp = count(dbp)
bys famid: egen total = count(ID)
drop if total ~= countdbp

*******************************************************************************
* REGULAR OLS
*******************************************************************************

*RESIDUALIZE THE OUTCOME
global controls sex i.YoB c.sex#i.YoB i.MoB e_PC*
reg dbp $controls, cluster(famid)
predict dbp_res, res
replace dbp=dbp_res

*STANDARDIZE VARIABLES
sum dbp
replace dbp = (dbp - r(mean))/r(sd)
sum dbp_ukb_ld
replace dbp_ukb_ld = (dbp_ukb_ld-r(mean))/r(sd)
sum dbp_ukb_ld_0
replace dbp_ukb_ld_0 = (dbp_ukb_ld_0-r(mean))/r(sd)
sum dbp_ukb_ld_1
replace dbp_ukb_ld_1 = (dbp_ukb_ld_1 - r(mean))/r(sd)

*THIS IS USING THE INDIVIDUAL SCORES AS REGRESSORS
reg dbp dbp_ukb_ld_0, cluster(famid)
estimates store e_ukb0
nlcom _b[dbp_ukb_ld_0]^2
sca r2_ukb0 =  r(b)[1,1]
sca r2_seukb0 = sqrt(r(V)[1,1])
sca coefbukb0 = _b[dbp_ukb_ld_0]

reg dbp dbp_ukb_ld_1, cluster(famid)
estimates store e_ukb1
nlcom _b[dbp_ukb_ld_1]^2
sca r2_ukb1 =  r(b)[1,1]
sca r2_seukb1 = sqrt(r(V)[1,1])
sca coefbukb1 = _b[dbp_ukb_ld_1]

*THIS IS THE META-ANALYSIS SCORE
reg dbp dbp_ukb_ld, cluster(famid)
estimates store e_meta
sca coefmeta = _b[dbp_ukb_ld]
nlcom _b[dbp_ukb_ld]^2
sca r2_meta =  r(b)[1,1]
sca r2_semeta = sqrt(r(V)[1,1])

*THIS IS THE PGI REPOSITORY CORRECTION PROCEDURE (BECKER ET AL.)
sca rho = sqrt(0.134/e(r2))
sca coefbecker=rho*_b[dbp_ukb_ld]
sca sebecker=rho*_se[dbp_ukb_ld]

*COMPUTE THE CORRELATION BETWEEN THE SCORES TO DETERMINE THE SCALING FACTOR
reg dbp_ukb_ld_0 dbp_ukb_ld_1
sca corrscores = _b[dbp_ukb_ld_1]

*THIS IS APPLYING ORIV USING THE SPLIT-SAMPLE UKB PGIs AS IVS FOR EACH OTHER
preserve
qui replace dbp_ukb_ld_0 = dbp_ukb_ld_0/sqrt(corrscores)
qui replace dbp_ukb_ld_1 = dbp_ukb_ld_1/sqrt(corrscores)
expand 2, generate(replicant)
gen mainVar = dbp_ukb_ld_0 if replicant == 0
replace mainVar = dbp_ukb_ld_1 if replicant == 1
gen instrument = dbp_ukb_ld_1 if replicant == 0
replace instrument = dbp_ukb_ld_0 if replicant == 1
reghdfe mainVar instrument, absorb(replicant) cluster(famid ID)
testparm instrument
ivreghdfe dbp (mainVar=instrument), absorb(replicant) cluster(famid ID)
estimates store e_iv
nlcom _b[mainVar]^2
sca r2_iv =  r(b)[1,1]
sca r2_seiv = sqrt(r(V)[1,1])
restore

set logtype text
log using "GEIGHEI\projects\Measurement error\tables-residualized-DBP.txt", replace
di "BETWEEN-FAMILY ESTIMATES"
estimates table e*,  drop(_cons) b(%9.3f) se(%9.3f) stats(r2 N)
scalar list r2_ukb0 r2_seukb0 r2_ukb1 r2_seukb1 r2_meta r2_semeta r2_iv r2_seiv 
scalar list coefbecker sebecker
log close

*******************************************************************************
* SIBLING FE
*******************************************************************************
xtset famid

*THIS IS USING THE INDIVIDUAL SCORES AS REGRESSORS
xtreg dbp dbp_ukb_ld_0, fe cluster(famid)
estimates store fe_ukb0
nlcom _b[dbp_ukb_ld_0]^2
sca r2_ukb0 =  r(b)[1,1]
sca r2_seukb0 = sqrt(r(V)[1,1])

xtreg dbp dbp_ukb_ld_1, fe cluster(famid)
estimates store fe_ukb1
nlcom _b[dbp_ukb_ld_1]^2
sca r2_ukb1 =  r(b)[1,1]
sca r2_seukb1 = sqrt(r(V)[1,1])

*THIS IS THE META-ANALYSIS SCORE
xtreg dbp dbp_ukb_ld, fe cluster(famid)
estimates store fe_meta
nlcom _b[dbp_ukb_ld]^2
sca r2_meta =  r(b)[1,1]
sca r2_semeta = sqrt(r(V)[1,1])

*THIS IS APPLYING ORIV WITH USING THE 23ANDME AND UKB AS IVS FOR EACH OTHER
preserve
qui replace dbp_ukb_ld_0 = dbp_ukb_ld_0/sqrt(corrscores)
qui replace dbp_ukb_ld_1 = dbp_ukb_ld_1/sqrt(corrscores)
expand 2, generate(replicant)
egen newid = concat(famid replicant)
gen long newidreal = real(newid)
xtset newidreal
gen mainVar = dbp_ukb_ld_0 if replicant == 0
replace mainVar = dbp_ukb_ld_1 if replicant == 1
gen instrument = dbp_ukb_ld_1 if replicant == 0
replace instrument = dbp_ukb_ld_0 if replicant == 1
reghdfe mainVar instrument, absorb(newidreal) cluster(famid ID)
testparm instrument
ivreghdfe dbp (mainVar=instrument), absorb(newidreal) cluster(famid ID)
estimates store fe_iv
nlcom (_b[mainVar])^2
sca r2_iv =  r(b)[1,1]
sca r2_seiv = sqrt(r(V)[1,1])
restore

log using "GEIGHEI\projects\Measurement error\tables-residualized-DBP.txt", append
di "WITHIN-FAMILY ESTIMATES"
estimates table fe*,  drop(_cons) b(%9.3f) se(%9.3f) stats(r2 N)
scalar list r2_ukb0 r2_seukb0 r2_ukb1 r2_seukb1 r2_meta r2_semeta r2_iv r2_seiv 
log close
