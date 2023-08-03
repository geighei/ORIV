*IMPERFECT GENETIC CORRELATION
*DESIGN 7, 8
*HANS, 29 MAY 2022, 2 SEPTEMBER 2022

*cd "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\GEIGHEI\projects\Measurement error\"
*cd "C:\Users\Hans\Dropbox (Erasmus Universiteit Rotterdam)\GEIGHEI\projects\Measurement error\"
cd "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\Genoeconomics\Measurement error PGS\Runs"

capture log close

set seed 1

local repl = 100

foreach g in 7 8 {

mat rho = (0.00 \ 0.25 \ 0.50 \ 0.75 \ 1.00)
local rowsrho = rowsof(rho)
mat samplesizesPRED = (16000)
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

local beta_G = sqrt(0.25)
local beta_G_wf = sqrt(0.25)
	
forvalues j = 1/`rowsrho' {
	if `j' == 1 {
		local s = "0.00" 
	}
	if `j'==2 {
		local s = "0.25" 
	}
	if `j'==3 {
		local s = "0.50"
	}
	if `j'==4 { 
		local s = "0.75"
	}
	if `j'==5 {
		local s = "1.00"
	}
	di "Rho " `s'
	
	forvalues k = 1/`rowsPRED' {
				
		di "Prediction sample size " samplesizesPRED[`k',1]
		
		forvalues i = 1/`repl' {
			
			if mod(`i',25)==0{ //show only every 25 iterations
			    di "Doing iteration `i' out of `repl'"
				}

			qui import delimited "DESIGN.`g'.RHO_G.`s'.RUN.`i'.HSq.out", clear
			qui sum heritability
			mat h2[`k',`i'] = r(mean)
			qui sum standarderror
			mat h2se[`k',`i'] = r(mean)
			
			*LOAD DATA
			qui import delimited "DESIGN.`g'.RHO_G.`s'.RUN.`i'.pgi", clear
			
			*******************************************************************************
			* BETWEEN-FAMILY
			*******************************************************************************
			*THIS IS USING THE META-ANALYSIS PGI
			qui reg y pgigwaspooled, cluster(pid)
			estimates store e_3
			mat coefmeta[`k',`i'] = temp[1,1]
			mat semeta[`k',`i'] = temp[2,1]
			mat biasmeta[`k',`i'] = (temp[1,1]-`beta_G')/`beta_G'
			mat rmsemeta[`k',`i'] = sqrt((temp[1,1]-`beta_G')^2)
			
			*PGI-RC
			sca rhobecker = sqrt(h2[`k',`i']/e(r2))
			mat coefbecker[`k',`i']=rhobecker*temp[1,1]
			mat sebecker[`k',`i']=rhobecker*temp[2,1]
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
			qui replace pgigwas1 = pgigwas1/sqrt(corrscores) if corrscores > 0
			qui replace pgigwas2 = pgigwas2/sqrt(corrscores) if corrscores > 0
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
			qui replace pgigwas1 = pgigwas1/sqrt(corrscores) if corrscores > 0
			qui replace pgigwas2 = pgigwas2/sqrt(corrscores) if corrscores > 0
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

mat rho`j' = J(`rowsPRED',1,rho[`j',1])

mat bf`j' = (rho`j', samplesizesPRED, avg_coef`j'meta, se_coef`j'meta, avg_coef`j'meta - 1.96*se_coef`j'meta, avg_coef`j'meta + 1.96*se_coef`j'meta, /*
	*/ avg_coef`j'oriv, se_coef`j'oriv, avg_coef`j'oriv - 1.96*se_coef`j'oriv, avg_coef`j'oriv + 1.96 * se_coef`j'oriv, /*
	*/ avg_coef`j'becker, se_coef`j'becker, avg_coef`j'becker - 1.96*se_coef`j'becker, avg_coef`j'becker + 1.96*se_coef`j'becker, /*
	*/ avg_coef`j'becker, se_coef`j'greml,  avg_coef`j'becker - 1.96*se_coef`j'greml,  avg_coef`j'becker + 1.96*se_coef`j'greml, repl`j')
matrix colnames bf`j' = rho sample_size coef_meta se_meta lb_meta ub_meta coef_oriv se_oriv lb_oriv ub_oriv coef_becker se_becker lb_becker ub_becker /*
	*/ coef_becker se_greml lb_greml ub_greml replications

mat overall_rmse`j' = (rho`j', samplesizesPRED, avg_rmse`j'meta, avg_rmse`j'oriv, avg_rmse`j'becker, avg_rmse`j'greml, repl`j')
matrix colnames overall_rmse`j' = rho sample_size  rmse_meta rmse_oriv rmse_becker rmse_greml replications

mat wf`j' = (rho`j', samplesizesPRED, avg_coef`j'metafe, se_coef`j'metafe, avg_coef`j'metafe - 1.96*se_coef`j'metafe, avg_coef`j'metafe + 1.96*se_coef`j'metafe, /*
	*/ avg_coef`j'orivfe, se_coef`j'orivfe, avg_coef`j'orivfe - 1.96*se_coef`j'orivfe, avg_coef`j'orivfe + 1.96 * se_coef`j'orivfe, repl`j')
matrix colnames wf`j' = rho sample_size coef_meta se_meta lb_meta ub_meta coef_oriv se_oriv lb_oriv ub_oriv replications

mat overall_rmse`j'fe = (rho`j', samplesizesPRED, avg_rmse`j'metafe, avg_rmse`j'orivfe, repl`j')
matrix colnames overall_rmse`j'fe = rho sample_size rmse_meta rmse_oriv replications

mat overall_F`j' = (rho`j', samplesizesPRED, avg`j'_F)
matrix colnames overall_F`j' = rho sample_size Fstat

mat overall_Ffe`j' = (rho`j', samplesizesPRED, avg`j'_Ffe)
matrix colnames overall_Ffe`j' = rho sample_size Fstat

mat overall_h2`j' = (rho`j', samplesizesPRED, avg_h2`j')
matrix colnames overall_h2`j' = rho sample_size h2

mat overall_seh2`j' = (rho`j', samplesizesPRED, avg_seh2`j')
matrix colnames overall_seh2`j' = rho sample_size seh2

}

mat bf_rho = (bf1 \ bf2 \ bf3  \ bf4 \ bf5)
putexcel set "sim_results.xls", sheet("bf_design`g'") modify
putexcel A1=matrix(bf_rho), colnames
mat wf_rho = (wf1 \ wf2 \  wf3  \ wf4 \ wf5)
putexcel set "sim_results.xls", sheet("wf_design`g'") modify
putexcel A1=matrix(wf_rho), colnames
mat rmse_rho = (overall_rmse1 \ overall_rmse2 \ overall_rmse3 \ overall_rmse4 \ overall_rmse5)
putexcel set "sim_results.xls", sheet("bf_rmse_design`g'") modify
putexcel A1=matrix(rmse_rho), colnames
mat rmse_rho_fe = (overall_rmse1fe \ overall_rmse2fe \ overall_rmse3fe \ overall_rmse4fe \ overall_rmse5fe)
putexcel set "sim_results.xls", sheet("wf_rmse_design`g'") modify
putexcel A1=matrix(rmse_rho_fe), colnames
mat F_rho = (overall_F1 \ overall_F2 \ overall_F3 \ overall_F4 \ overall_F5)
putexcel set "sim_results.xls", sheet("bf_F_design`g'") modify
putexcel A1=matrix(F_rho), colnames
mat F_rho_fe = (overall_Ffe1 \ overall_Ffe2 \ overall_Ffe3 \ overall_Ffe4 \ overall_Ffe5)
putexcel set "sim_results.xls", sheet("wf_F_design`g'") modify
putexcel A1=matrix(F_rho_fe), colnames
mat h2 = (overall_h21 \ overall_h22 \ overall_h23 \ overall_h24 \ overall_h25)
putexcel set "sim_results.xls", sheet("h2_design`g'") modify
putexcel A1=matrix(h2), colnames
mat seh2 = (overall_seh21 \ overall_seh22 \ overall_seh23 \ overall_seh24 \ overall_seh25)
putexcel set "sim_results.xls", sheet("seh2_design`g'") modify
putexcel A1=matrix(seh2), colnames

}

local beta_G=sqrt(0.25)
import excel "sim_results.xls", sheet("bf_design7") firstrow clear
replace rho = rho*100
drop if rho==0
mkmat coef_meta lb_meta ub_meta, matrix(meta) 
mat rownames meta = 0.25 0.5 0.75 1
mkmat coef_becker lb_becker ub_becker, matrix(becker)
mat rownames becker = 0.25 0.5 0.75 1
mkmat coef_oriv lb_oriv ub_oriv, matrix(oriv)
mat rownames oriv = 0.25 0.5 0.75 1
mkmat coef_becker lb_greml ub_greml, matrix(becker_greml)
mat rownames becker_greml = 0.25 0.5 0.75 1
matrix transpose_meta = meta'
matrix transpose_becker = becker'
matrix transpose_oriv = oriv'
matrix transpose_greml = becker_greml'
coefplot matrix(transpose_meta) matrix(transpose_oriv) matrix(transpose_becker) matrix(transpose_greml), vertical grid(between) ylabel(0.1(0.1)0.8) yline(`beta_G', lpattern(dash)) /*
	*/ xtitle("Genetic correlation between discovery and prediction sample") /*
	*/ ci((2 3)) scheme(s1mono) ytitle("Estimated coefficient")  legend(rows(1) lab(2 "Meta-analysis") lab(4 "ORIV") lab(6 "PGI-RC (default)") lab(8 "PGI-RC (GREML unc.)"))
graph export "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\GEIGHEI\projects\Measurement error\Figures revision\design7-bf-greml.pdf", replace

import excel "sim_results.xls", sheet("bf_design8") firstrow clear
replace rho = rho*100
drop if rho==0
mkmat coef_meta lb_meta ub_meta, matrix(meta) 
mat rownames meta = 0.25 0.5 0.75 1
mkmat coef_becker lb_becker ub_becker, matrix(becker)
mat rownames becker = 0.25 0.5 0.75 1
mkmat coef_oriv lb_oriv ub_oriv, matrix(oriv)
mat rownames oriv = 0.25 0.5 0.75 1
mkmat coef_becker lb_greml ub_greml, matrix(becker_greml)
mat rownames becker_greml = 0.25 0.5 0.75 1
matrix transpose_meta = meta'
matrix transpose_becker = becker'
matrix transpose_oriv = oriv'
matrix transpose_greml = becker_greml'
coefplot matrix(transpose_meta) matrix(transpose_oriv) matrix(transpose_becker) matrix(transpose_greml), vertical grid(between) ci((2 3)) ylabel(0.1(0.1)0.8) yline(`beta_G', lpattern(dash)) /*
	*/ scheme(s1mono) ytitle("Estimated coefficient")  xtitle("Genetic correlation between the two GWAS samples") /*
	*/ legend(rows(1) lab(2 "Meta-analysis") lab(4 "ORIV") lab(6 "PGI-RC (default)") lab(8 "PGI-RC (GREML unc.)"))
graph export "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\GEIGHEI\projects\Measurement error\Figures revision\design8-bf-greml.pdf", replace
