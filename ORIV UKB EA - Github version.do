* REPLICATION DO-FILE FOR EMPIRICAL ILLUSTRATION "STOP META-ANALYZING, START INSTRUMENTING"
* CHECK THE PREDICTIVE POWER OF THE VARIOUS POLYGENIC SCORES FOR EA IN THE SIBLING SAMPLE OF THE UKB
* HANS VAN KIPPERSLUIS, FEBRUARY 2021
* THIS PART COVERS THE IMPLEMENTATION OF ORIV AND OLS ON BASIS OF META-ANALYSIS SCORES
* THE SYNTAX FOR CONSTRUCTION OF POLYGENIC SCORES AND SAMPLE SELECTIONS ARE AVAILABLE UPON REQUEST
* COMMENTS, QUESTIONS TO hvankippersluis at ese.eur.nl

capture log close
estimates drop _all
cd "C:\Users\"
use "ukb.dta", clear

// Keep Europeans 
keep if european==1

*KEEP SIBLING SAMPLE
keep if famid ~= .
bys famid: egen count = count(ea_ukb_23me_ld)
drop if count < 2
drop if ea_ukb_23me_ld == .

keep if EA ~= .
keep if height ~= .

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
nlcom _b[ea_ukb_ld]^2
sca r2_ukb =  r(b)[1,1]
sca r2_seukb = sqrt(r(V)[1,1])
sca coefbukb = _b[ea_ukb_ld]

reg EA ea_23me_ld, cluster(famid)
estimates store e_23me
sca coefb23me = _b[ea_23me_ld]
nlcom _b[ea_23me_ld]^2
sca r2_23me =  r(b)[1,1]
sca r2_se23me = sqrt(r(V)[1,1])

*THIS IS THE META-ANALYSIS SCORE
reg EA ea_ukb_23me_ld, cluster(famid)
estimates store e_meta
nlcom _b[ea_ukb_23me_ld]^2
sca r2_meta =  r(b)[1,1]
sca r2_semeta = sqrt(r(V)[1,1])
sca coefmeta = _b[ea_ukb_23me_ld]

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

*THIS IS APPLYING ORIV WITH USING THE SPLIT SAMPLE UKB SCORES AS IVS FOR EACH OTHER
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

set logtype text
log using "\tables-residualized.txt", replace
di "BETWEEN-FAMILY ESTIMATES"
estimates table e*,  drop(_cons) b(%9.3f) se(%9.3f) stats(r2 N)
scalar list r2_ukb r2_seukb r2_23me r2_se23me r2_meta r2_semeta r2_iv r2_seiv r2_ivukb r2_seivukb
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

log using "tables-residualized.txt", append
di "WITHIN-FAMILY ESTIMATES"
estimates table fe*,  drop(_cons) b(%9.3f) se(%9.3f) stats(r2 N)
scalar list r2_ukb r2_seukb r2_23me r2_se23me r2_meta r2_semeta r2_iv r2_seiv r2_ivukb r2_seivukb
log close


