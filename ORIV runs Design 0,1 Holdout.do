*DO-FILE TO GENERATE RESULTS FROM DESIGNS 0 AND 1
*HANS, APRIL 2022

*cd "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\GEIGHEI\projects\Measurement error\"
*cd "C:\Users\Hans\Dropbox (Erasmus Universiteit Rotterdam)\GEIGHEI\projects\Measurement error\"
cd "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\Genoeconomics\Measurement error PGS\Runs"

capture log close
matrix drop _all
clear all

set seed 1

local repl = 100

mat samplesizesGWAS = (32000)
local rowsGWAS = rowsof(samplesizesGWAS)
mat samplesizesPRED = (1000 \ 2000 \ 4000 \ 8000 \ 16000)
local rowsPRED = rowsof(samplesizesPRED)

*INITIALIZE MATRICES
mat h2 = J(`rowsPRED',`repl',0)
mat h2se = J(`rowsPRED',`repl',0)

mat coefmeta = J(`rowsPRED',`repl',0)
mat coeforiv = J(`rowsPRED',`repl',0)
mat coefbecker = J(`rowsPRED',`repl',0)
mat semeta = J(`rowsPRED',`repl',0)
mat seoriv = J(`rowsPRED',`repl',0)
mat sebecker = J(`rowsPRED',`repl',0)
mat segreml = J(`rowsPRED',`repl',0)
mat biasmeta = J(`rowsPRED',`repl',0)
mat biasoriv = J(`rowsPRED',`repl',0)
mat biasbecker = J(`rowsPRED',`repl',0)
mat rmsemeta = J(`rowsPRED',`repl',0)
mat rmsemeta2 = J(`rowsPRED',`repl',0)
mat rmseoriv = J(`rowsPRED',`repl',0)
mat rmsebecker = J(`rowsPRED',`repl',0)
mat rmsegreml = J(`rowsPRED',`repl',0)
mat F = J(`rowsPRED',`repl',0)

mat coefmetafe = J(`rowsPRED',`repl',0)
mat coeforivfe = J(`rowsPRED',`repl',0)
mat semetafe = J(`rowsPRED',`repl',0)
mat seorivfe = J(`rowsPRED',`repl',0)
mat biasmetafe = J(`rowsPRED',`repl',0)
mat biasorivfe = J(`rowsPRED',`repl',0)
mat rmsemetafe = J(`rowsPRED',`repl',0)
mat rmseorivfe = J(`rowsPRED',`repl',0)
mat Ffe = J(`rowsPRED',`repl',0)

mat avg_coefmeta = J(`rowsPRED',1,0)
mat avg_coeforiv = J(`rowsPRED',1,0)
mat avg_coefbecker = J(`rowsPRED',1,0)
mat se_coefmeta = J(`rowsPRED',1,0)
mat se_coeforiv = J(`rowsPRED',1,0)
mat se_coefbecker = J(`rowsPRED',1,0)
mat se_coefgreml = J(`rowsPRED',1,0)
mat avg_biasmeta = J(`rowsPRED',1,0)
mat avg_biasoriv = J(`rowsPRED',1,0)
mat avg_biasbecker = J(`rowsPRED',1,0)
mat avg_rmsemeta = J(`rowsPRED',1,0)
mat avg_rmsemeta2 = J(`rowsPRED',1,0)
mat avg_rmseoriv = J(`rowsPRED',1,0)
mat avg_rmsebecker = J(`rowsPRED',1,0)
mat avg_rmsegreml = J(`rowsPRED',1,0)
mat avg_F = J(`rowsPRED',1,0)
mat avg_h2 = J(`rowsPRED',1,0)
mat avg_seh2 = J(`rowsPRED',1,0)

mat avg_coefmetafe = J(`rowsPRED',1,0)
mat avg_coeforivfe = J(`rowsPRED',1,0)
mat se_coeftruefe = J(`rowsPRED',1,0)
mat se_coefmetafe = J(`rowsPRED',1,0)
mat se_coeforivfe = J(`rowsPRED',1,0)
mat avg_biasmetafe = J(`rowsPRED',1,0)
mat avg_biasorivfe = J(`rowsPRED',1,0)
mat avg_rmsemetafe = J(`rowsPRED',1,0)
mat avg_rmseorivfe = J(`rowsPRED',1,0)
mat avg_Ffe = J(`rowsPRED',1,0)

forvalues h = 0/1 {
    
	if `h'==0 {
		local beta_G = sqrt(0.25)
		local beta_G_wf = sqrt(0.25)	
	}
	if `h'== 1 {
		local beta_G = sqrt(0.25)
		local beta_G_wf = sqrt(0.2)
	}
	di "design " `h'

forvalues j = 1/`rowsGWAS' {
	
	local s = samplesizesGWAS[`j',1]
	di "GWAS sample size " `s'
		
	forvalues k = 1/`rowsPRED' {
		
		di "Prediction sample size " samplesizesPRED[`k',1]
		if `k'==1 {
		    local t = 1000
		}
		if `k'==2 {
		    local t = 2000
		}
		if `k'==3 {
		    local t = 4000
		}
		if `k'==4 {
		    local t = 8000
		}
		if `k'==5 {
		    local t = 16000
		}
		
		forvalues i = 1/`repl' {
			
			if mod(`i',25)==0{ //show only every 25 iterations
			    di "Doing iteration `i' out of `repl'"
				}
			
			*LOAD HERITABILITY ESTIMATES (FROM DESIGN 3, SINCE HERE THE PREDICTION SAMPLE SIZE IS VARIED)
			qui import delimited "DESIGN.3.NHOLDOUT.`t'.RUN.`i'.HSq.out", clear
			qui sum heritability
			mat h2[`k',`i'] = r(mean)
			if r(mean) < 0.01 {
				mat h2[`k',`i']=.
			}
			qui sum standarderror
			mat h2se[`k',`i'] = r(mean)
						
			*LOAD DATA
			qui import delimited "DESIGN.`h'.NGWASPOOLED.`s'.RUN.`i'.pgi", clear

			*KEEP RANDOM SAMPLE TO CREATE DESIRED SAMPLE SIZE FOR THE HOLDOUT (PREDICTION) SAMPLE
			qui gen rnd = runiform()
			sort pid
			qui replace rnd = rnd[_n-1] if pid==pid[_n-1]
			sort rnd
			qui drop if _n > samplesizesPRED[`k',1]
			qui drop rnd
				
			*******************************************************************************
			* BETWEEN-FAMILY
			*******************************************************************************
			*THIS IS USING THE META-ANALYSIS PGI
			qui reg y pgigwaspooled, cluster(pid)
			estimates store e_3
			mat temp = r(table)
			mat coefmeta[`k',`i'] = temp[1,1]
			mat semeta[`k',`i'] = temp[2,1]
			mat biasmeta[`k',`i'] = (temp[1,1]-`beta_G')/`beta_G'
			mat rmsemeta[`k',`i'] = sqrt((temp[1,1]-`beta_G')^2)
			mat rmsemeta2[`k',`i'] = sqrt((`beta_G'*biasmeta[`k',`i'])^2 + semeta[`k',`i']^2)
			
			*BECKER ET AL. CORRECTION
			sca rho = sqrt(h2[`k',`i']/e(r2))
			mat coefbecker[`k',`i']=rho*temp[1,1]
			mat sebecker[`k',`i']=rho*temp[2,1]
			mat biasbecker[`k',`i'] = (coefbecker[`k',`i']-`beta_G')/`beta_G'
			mat rmsebecker[`k',`i'] = sqrt((coefbecker[`k',`i']-`beta_G')^2)
			
			*INCORPORATING GREML UNCERTAINTY
			mat segreml[`k',`i']=h2se[`k',`i']/(2*sqrt(h2[`k',`i']))
			mat rmsegreml[`k',`i']=sqrt(biasbecker[`k',`i']^2+segreml[`k',`i']^2)
						
			*COMPUTE THE CORRELATION BETWEEN THE SCORES TO DETERMINE THE SCALING FACTOR
			qui reg pgigwas1 pgigwas2, cluster(pid)
			sca corrscores = _b[pgigwas2]

			*THIS IS APPLYING ORIV WITH USING THE SCORES AS IVS FOR EACH OTHER
			preserve
			if corrscores > 0 { 
			qui replace pgigwas1 = pgigwas1/sqrt(corrscores) 
			qui replace pgigwas2 = pgigwas2/sqrt(corrscores) 
			quietly {
			expand 2, generate(replicant)
			}
			qui gen mainVar = pgigwas1 if replicant == 0
			qui replace mainVar = pgigwas2 if replicant == 1
			qui gen instrument = pgigwas2 if replicant == 0
			qui replace instrument = pgigwas1 if replicant == 1
			qui ivreghdfe y (mainVar=instrument), absorb(replicant) cluster(pid iid)
			mat coeforiv[`k',`i'] = _b[mainVar]
			mat seoriv[`k',`i'] = _se[mainVar]
			mat biasoriv[`k',`i'] = (_b[mainVar]-`beta_G')/`beta_G'
			mat rmseoriv[`k',`i'] = sqrt((_b[mainVar]-`beta_G')^2)
			mat F[`k',`i'] = e(widstat)
			}
			if corrscores < 0 {
				mat coeforiv[`k',`i'] = .
				mat seoriv[`k',`i'] = .
				mat biasoriv[`k',`i'] = .
				mat rmseoriv[`k',`i'] = .
				mat F[`k',`i'] = .
				di "corr scores < 0"
			}
			restore
			
			*******************************************************************************
			* WITHIN-FAMILY
			*******************************************************************************
			qui xtset pid
				
			*THIS IS USING THE META-ANALYSIS PGI, WITH FIXED EFFECTS
			qui xtreg y pgigwaspooled, fe cluster(pid)
			estimates store e_3
			mat temp = r(table)
			mat coefmetafe[`k',`i'] = temp[1,1]
			mat semetafe[`k',`i'] = temp[2,1]
			mat biasmetafe[`k',`i'] = (temp[1,1]-`beta_G_wf')/`beta_G_wf'
			mat rmsemetafe[`k',`i'] = sqrt((temp[1,1]-`beta_G_wf')^2)
			
			*ORIV WITH FIXED EFFECTS
			preserve
			if corrscores > 0 {
			qui replace pgigwas1 = pgigwas1/sqrt(corrscores) 
			qui replace pgigwas2 = pgigwas2/sqrt(corrscores) 
			quietly {
			expand 2, generate(replicant)
			}
			egen newid = concat(pid replicant)
			gen long newidreal = real(newid)
			qui xtset newidreal
			qui gen mainVar = pgigwas1 if replicant == 0
			qui replace mainVar = pgigwas2 if replicant == 1
			qui gen instrument = pgigwas2 if replicant == 0
			qui replace instrument = pgigwas1 if replicant == 1
			qui ivreghdfe y (mainVar=instrument), absorb(newidreal) cluster(pid iid)
			mat coeforivfe[`k',`i'] = _b[mainVar]
			mat seorivfe[`k',`i'] = _se[mainVar]
			mat biasorivfe[`k',`i'] = (_b[mainVar]-`beta_G_wf')/`beta_G_wf'
			mat rmseorivfe[`k',`i'] = sqrt((_b[mainVar]-`beta_G_wf')^2)
			mat Ffe[`k',`i'] = e(widstat)
			}
			if corrscores < 0 {
				mat coeforivfe[`k',`i'] = .
				mat seorivfe[`k',`i'] = .
				mat biasorivfe[`k',`i'] = .
				mat rmseorivfe[`k',`i'] = .
				mat Ffe[`k',`i'] = .
				di "corr scores FE < 0"
			}
			restore
			
			
			}
		}
	mata: st_matrix("avg_coef`j'meta", rowsum(st_matrix("coefmeta")/`repl'))
	mata: st_matrix("avg_coef`j'oriv", rowsum(st_matrix("coeforiv")/`repl'))
	mata: st_matrix("avg_coef`j'becker", rowsum(st_matrix("coefbecker")/`repl'))
	
	mata: st_matrix("se_coef`j'meta", rowsum(st_matrix("semeta")/`repl'))
	mata: st_matrix("se_coef`j'oriv", rowsum(st_matrix("seoriv")/`repl'))
	mata: st_matrix("se_coef`j'becker", rowsum(st_matrix("sebecker")/`repl'))
	mata: st_matrix("se_coef`j'greml", rowsum(st_matrix("segreml")/`repl'))
	
	mata: st_matrix("avg_bias`j'meta", rowsum(st_matrix("biasmeta")/`repl'))
	mata: st_matrix("avg_bias`j'oriv", rowsum(st_matrix("biasoriv")/`repl'))
	mata: st_matrix("avg_bias`j'becker", rowsum(st_matrix("biasbecker")/`repl'))	

	mata: st_matrix("avg_rmse`j'meta", rowsum(st_matrix("rmsemeta")/`repl'))
	mata: st_matrix("avg_rmse`j'meta2", rowsum(st_matrix("rmsemeta2")/`repl'))
	mata: st_matrix("avg_rmse`j'oriv", rowsum(st_matrix("rmseoriv")/`repl'))
	mata: st_matrix("avg_rmse`j'becker", rowsum(st_matrix("rmsebecker")/`repl'))
	mata: st_matrix("avg_rmse`j'greml", rowsum(st_matrix("rmsegreml")/`repl'))
	
	mata: st_matrix("avg`j'_F", rowsum(st_matrix("F")/`repl'))

	mata: st_matrix("avg_coef`j'metafe", rowsum(st_matrix("coefmetafe")/`repl'))
	mata: st_matrix("avg_coef`j'orivfe", rowsum(st_matrix("coeforivfe")/`repl'))

	mata: st_matrix("se_coef`j'metafe", rowsum(st_matrix("semetafe")/`repl'))
	mata: st_matrix("se_coef`j'orivfe", rowsum(st_matrix("seorivfe")/`repl'))

	mata: st_matrix("avg_bias`j'metafe", rowsum(st_matrix("biasmetafe")/`repl'))
	mata: st_matrix("avg_bias`j'orivfe", rowsum(st_matrix("biasorivfe")/`repl'))

	mata: st_matrix("avg_rmse`j'metafe", rowsum(st_matrix("rmsemetafe")/`repl'))
	mata: st_matrix("avg_rmse`j'orivfe", rowsum(st_matrix("rmseorivfe")/`repl'))
	
	mata: st_matrix("avg`j'_Ffe", rowsum(st_matrix("Ffe")/`repl'))
	
	mata: st_matrix("avg_h2`j'", rowsum(st_matrix("h2")/`repl'))
	mata: st_matrix("avg_seh2`j'", rowsum(st_matrix("h2se")/`repl'))
	

*STORE RESULTS IN A MATRIX AND PLOT THE RESULTS
mat repl`j' = J(`rowsPRED',1,`repl')

mat NGWAS`j' = J(`rowsPRED',1,samplesizesGWAS[`j',1])

mat bf`j' = (NGWAS`j', samplesizesPRED, avg_coef`j'meta, se_coef`j'meta, avg_coef`j'meta - 1.96*se_coef`j'meta, avg_coef`j'meta + 1.96*se_coef`j'meta, /*
	*/ avg_coef`j'oriv, se_coef`j'oriv, avg_coef`j'oriv - 1.96*se_coef`j'oriv, avg_coef`j'oriv + 1.96 * se_coef`j'oriv, /*
	*/ avg_coef`j'becker, se_coef`j'becker, avg_coef`j'becker - 1.96*se_coef`j'becker, avg_coef`j'becker + 1.96*se_coef`j'becker, /*
	*/ avg_coef`j'becker, se_coef`j'greml,  avg_coef`j'becker - 1.96*se_coef`j'greml,  avg_coef`j'becker + 1.96*se_coef`j'greml, repl`j')
matrix colnames bf`j' = GWAS sample_size coef_meta se_meta lb_meta ub_meta coef_oriv se_oriv lb_oriv ub_oriv coef_becker se_becker lb_becker ub_becker /*
	*/ coef_becker se_greml lb_greml ub_greml replications

mat overall_rmse`j' = (NGWAS`j', samplesizesPRED, avg_rmse`j'meta, avg_rmse`j'meta2, avg_rmse`j'oriv, avg_rmse`j'becker, avg_rmse`j'greml, repl`j')
matrix colnames overall_rmse`j' = GWAS sample_size rmse_meta rmse_meta2 rmse_oriv rmse_becker rmse_greml replications

mat wf`j' = (NGWAS`j', samplesizesPRED, avg_coef`j'metafe, se_coef`j'metafe, avg_coef`j'metafe - 1.96*se_coef`j'metafe, avg_coef`j'metafe + 1.96*se_coef`j'metafe, /*
	*/ avg_coef`j'orivfe, se_coef`j'orivfe, avg_coef`j'orivfe - 1.96*se_coef`j'orivfe, avg_coef`j'orivfe + 1.96 * se_coef`j'orivfe, repl`j')
matrix colnames wf`j' = GWAS sample_size coef_meta se_meta lb_meta ub_meta coef_oriv se_oriv lb_oriv ub_oriv replications

mat overall_rmse`j'fe = (NGWAS`j', samplesizesPRED, avg_rmse`j'metafe, avg_rmse`j'orivfe, repl`j')
matrix colnames overall_rmse`j'fe = GWAS sample_size rmse_meta rmse_oriv replications

mat overall_F`j' = (NGWAS`j', samplesizesPRED, avg`j'_F)
matrix colnames overall_F`j' = GWAS sample_size Fstat

mat overall_Ffe`j' = (NGWAS`j', samplesizesPRED, avg`j'_Ffe)
matrix colnames overall_Ffe`j' = GWAS sample_size Fstat

}

*WRITE THE RESULTS INTO EXCEL FILE

putexcel set "sim_results.xls", sheet("bf_design`h'") modify
putexcel A1=matrix(bf1), colnames
putexcel set "sim_results.xls", sheet("wf_design`h'") modify
putexcel A1=matrix(wf1), colnames
putexcel set "sim_results.xls", sheet("bf_rmse_design`h'") modify
putexcel A1=matrix(overall_rmse1), colnames
putexcel set "sim_results.xls", sheet("wf_rmse_design`h'") modify
putexcel A1=matrix(overall_rmse1fe), colnames
putexcel set "sim_results.xls", sheet("bf_F_design`h'") modify
putexcel A1=matrix(overall_F1), colnames
putexcel set "sim_results.xls", sheet("wf_F_design`h'") modify
putexcel A1=matrix(overall_Ffe1), colnames

*MAKE PLOTS
import excel "sim_results.xls", sheet("bf_design`h'") firstrow clear
drop if _n>5
mkmat coef_meta lb_meta ub_meta, matrix(meta) rownames(sample_size)
mkmat coef_becker lb_becker ub_becker, matrix(becker) rownames(sample_size)
mkmat coef_oriv lb_oriv ub_oriv, matrix(oriv) rownames(sample_size)
mkmat coef_becker lb_greml ub_greml, matrix(becker_greml) rownames(sample_size)
matrix transpose_meta = meta'
matrix transpose_becker = becker'
matrix transpose_oriv = oriv'
matrix transpose_greml = becker_greml'
coefplot matrix(transpose_meta) matrix(transpose_oriv) matrix(transpose_becker) matrix(transpose_greml), vertical grid(between) ylabel(0(0.2)0.85) ci((2 3)) scheme(s1mono) ytitle("Estimated coefficient") yline(`beta_G', lpattern(dash))  xtitle("Prediction sample size") legend(rows(1) lab(2 "Meta-analysis") lab(4 "ORIV") lab(6 "PGI-RC (default)") lab(8 "PGI-RC (GREML unc.)"))
graph export "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\GEIGHEI\projects\Measurement error\Figures revision\design`h'-bf-holdout-greml.pdf", replace

import excel "sim_results.xls", sheet("wf_design`h'") firstrow clear
drop if _n>5
mkmat coef_meta lb_meta ub_meta, matrix(meta) rownames(sample_size)
mkmat coef_oriv lb_oriv ub_oriv, matrix(oriv) rownames(sample_size)
matrix transpose_meta = meta'
matrix transpose_oriv = oriv'
coefplot matrix(transpose_meta) matrix(transpose_oriv), vertical grid(between) ci((2 3)) scheme(s1mono) yline(`beta_G_wf', lpattern(dash)) ytitle("Estimated coefficient") ylabel(0(0.2)0.85)  xtitle("Prediction sample size") legend(rows(1) lab(2 "Meta-analysis") lab(4 "ORIV"))
graph export "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\GEIGHEI\projects\Measurement error\Figures revision\design`h'-wf-holdout.pdf", replace

}


