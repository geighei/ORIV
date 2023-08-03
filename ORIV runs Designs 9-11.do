*DO-FILE TO GENERATE RESULTS FROM DESIGNS 9-11
*HANS, AUGUST 2022

cd "C:\Users\41153jvk\Dropbox (Erasmus Universiteit Rotterdam)\Genoeconomics\Measurement error PGS\Runs"

capture log close

set seed 1

local repl = 100

mat samplesizesGWAS9 = (6666 \ 50000)
mat samplesizesGWAS10 = (4000 \ 30000)
mat samplesizesGWAS11 = (2000 \ 15000)
mat samplesizesPRED = (1000 \ 4000 \ 16000)
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

forvalues h = 10/11 {
    
	local rowsGWAS = rowsof(samplesizesGWAS`h')	
	if `h' == 9 {
		local h2 = "0.15"
		local beta_G = sqrt(0.15)
		local beta_G_wf = sqrt(0.15)
	}
	if `h' == 10 {
		local h2 = "0.25"
		local beta_G = sqrt(0.25)
		local beta_G_wf = sqrt(0.25)
	}
	if `h' == 11 {
		local h2 = "0.5"
		local beta_G = sqrt(0.5)
		local beta_G_wf = sqrt(0.5)
	}
	
	di "design " `h'
	
forvalues j = 1/`rowsGWAS' {
	
	local s = samplesizesGWAS`h'[`j',1]
	di "GWAS sample size " `s'
		
	forvalues k = 1/`rowsPRED' {
		
		di "Prediction sample size " samplesizesPRED[`k',1]
		if `k'==1 {
		    local t = 1000
		}
		if `k'==2 {
		    local t = 4000
		}
		if `k'==3 {
		    local t = 16000
		}
		
		forvalues i = 1/`repl' {
			
			if mod(`i',25)==0{ //show only every 25 iterations
			    di "Doing iteration `i' out of `repl'"
				}
			
			*LOAD HERITABILITY ESTIMATES
			qui import delimited "DESIGN.`h'.HSQ.`h2'.NGWASPOOLED.`s'.NHOLDOUT.`t'.RUN.`i'.HSq.out", clear
			qui sum heritability
			mat h2[`k',`i'] = r(mean)
			if r(mean) < 0.01 {
				mat h2[`k',`i']=.
				di "too low heritability"
			}
			qui sum standarderror
			mat h2se[`k',`i'] = r(mean)

			
			*LOAD DATA
			qui import delimited "DESIGN.`h'.HSQ.`h2'.NGWASPOOLED.`s'.NHOLDOUT.`t'.RUN.`i'.pgi", clear
			
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

mat NGWAS`j' = J(`rowsPRED',1,samplesizesGWAS`h'[`j',1])

mat overall_coef`j' = (NGWAS`j', samplesizesPRED, avg_coef`j'meta, avg_coef`j'oriv, avg_coef`j'becker, repl`j')
matrix colnames overall_coef`j' = GWAS sample_size coef_meta coef_oriv coef_becker replications

mat overall_se`j' = (NGWAS`j', samplesizesPRED, se_coef`j'meta, se_coef`j'oriv, se_coef`j'becker, se_coef`j'greml, repl`j')
matrix colnames overall_se`j' = GWAS sample_size se_meta se_oriv se_becker se_greml replications

mat overall_bias`j' = (NGWAS`j', samplesizesPRED, avg_bias`j'meta, avg_bias`j'oriv, avg_bias`j'becker, repl`j')
matrix colnames overall_bias`j' = GWAS sample_size bias_meta bias_oriv bias_becker replications

mat overall_rmse`j' = (NGWAS`j', samplesizesPRED, avg_rmse`j'meta, avg_rmse`j'oriv, avg_rmse`j'becker, avg_rmse`j'greml, repl`j')
matrix colnames overall_rmse`j' = GWAS sample_size  rmse_meta rmse_oriv rmse_becker rmse_greml replications

mat overall_coef`j'fe = (NGWAS`j', samplesizesPRED, avg_coef`j'metafe, avg_coef`j'orivfe, repl`j')
matrix colnames overall_coef`j'fe = GWAS sample_size coef_meta coef_oriv replications

mat overall_se`j'fe = (NGWAS`j', samplesizesPRED, se_coef`j'metafe, se_coef`j'orivfe, repl`j')
matrix colnames overall_se`j'fe = GWAS sample_size se_meta se_oriv replications

mat overall_bias`j'fe = (NGWAS`j', samplesizesPRED, avg_bias`j'metafe, avg_bias`j'orivfe, repl`j')
matrix colnames overall_bias`j'fe = GWAS sample_size bias_meta bias_oriv replications

mat overall_rmse`j'fe = (NGWAS`j', samplesizesPRED, avg_rmse`j'metafe, avg_rmse`j'orivfe, repl`j')
matrix colnames overall_rmse`j'fe = GWAS sample_size rmse_meta rmse_oriv replications

mat overall_F`j' = (NGWAS`j', samplesizesPRED, avg`j'_F)
matrix colnames overall_F`j' = GWAS sample_size Fstat

mat overall_Ffe`j' = (NGWAS`j', samplesizesPRED, avg`j'_Ffe)
matrix colnames overall_Ffe`j' = GWAS sample_size Fstat

mat overall_h2`j' = (NGWAS`j', samplesizesPRED, avg_h2`j')
matrix colnames overall_h2`j' = GWAS sample_size h2

mat overall_seh2`j' = (NGWAS`j', samplesizesPRED, avg_seh2`j')
matrix colnames overall_seh2`j' = GWAS sample_size seh2

}

*WRITE THE RESULTS INTO EXCEL FILE

mat coef_GWAS = (overall_coef1, overall_coef2)
putexcel set "sim_results.xls", sheet("bf_coef_design`h'") modify
putexcel A1=matrix(coef_GWAS), colnames
mat se_GWAS = (overall_se1, overall_se2)
putexcel set "sim_results.xls", sheet("bf_se_design`h'") modify
putexcel A1=matrix(se_GWAS), colnames
mat coef_GWAS_fe = (overall_coef1fe, overall_coef2fe)
putexcel set "sim_results.xls", sheet("wf_coef_design`h'") modify
putexcel A1=matrix(coef_GWAS_fe), colnames
mat se_GWAS_fe = (overall_se1fe, overall_se2fe)
putexcel set "sim_results.xls", sheet("wf_se_design`h'") modify
putexcel A1=matrix(se_GWAS_fe), colnames
mat bias_GWAS = (overall_bias1, overall_bias2)
putexcel set "sim_results.xls", sheet("bf_bias_design`h'") modify
putexcel A1=matrix(bias_GWAS), colnames
mat bias_GWAS_fe = (overall_bias1fe, overall_bias2fe)
putexcel set "sim_results.xls", sheet("wf_bias_design`h'") modify
putexcel A1=matrix(bias_GWAS_fe), colnames
mat rmse_GWAS = (overall_rmse1, overall_rmse2)
putexcel set "sim_results.xls", sheet("bf_rmse_design`h'") modify
putexcel A1=matrix(rmse_GWAS), colnames
mat rmse_GWAS_fe = (overall_rmse1fe, overall_rmse2fe)
putexcel set "sim_results.xls", sheet("wf_rmse_design`h'") modify
putexcel A1=matrix(rmse_GWAS_fe), colnames
mat F_GWAS = (overall_F1, overall_F2)
putexcel set "sim_results.xls", sheet("bf_F_design`h'") modify
putexcel A1=matrix(F_GWAS), colnames
mat F_GWAS_fe = (overall_Ffe1, overall_Ffe2)
putexcel set "sim_results.xls", sheet("wf_F_design`h'") modify
putexcel A1=matrix(F_GWAS_fe), colnames

mat h2print = (overall_h21, overall_h22)
putexcel set "sim_results.xls", sheet("h2_design`h'") modify
putexcel A1=matrix(h2print), colnames
mat seh2print = (overall_seh21, overall_seh22)
putexcel set "sim_results.xls", sheet("seh2_design`h'") modify
putexcel A1=matrix(seh2print), colnames

}


