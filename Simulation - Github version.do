* REPLICATION DO-FILE FOR MAIN SIMULATIONS "STOP META-ANALYZING, START INSTRUMENTING" 
* DO POWER CALCULATIONS + RMSE FOR SAMPLE SIZE IN PRESENCE OF MEASUREMENT ERROR
* HANS VAN KIPPERSLUIS, 2020/2021
* THE SYNTAX FOR THE OTHER SIMULATIONS ARE AVAILABLE UPON REQUEST
* COMMENTS, QUESTIONS TO hvankippersluis at ese.eur.nl

capture clear all
set more off
set seed 31415

cd "/Users/`c(username)'/Downloads"
capture mkdir ORIV_simulations
cd "ORIV_simulations"

*SET NUMBER OF REPLICATIONS
local repl = 1000

*SET SAMPLE SIZES 
mat samplesizes = (1000 \ 2000 \ 3000 \ 4000 \ 5000 \ 10000 \ 20000)

*SET EFFECT SIZES G (SQUARED STANDARDIZED EFFECT SIZE IS (INCREMENTAL) R-SQUARED (SNP-BASED HERITABILITY))
mat beta_G = (0.5 \ 0.5 \ 0.5 \ 0.5) // this the true standardized effect size
mat beta_Gstar = (0.35 \ 0.35 \ 0.35 \ 0.35) // this is for the meta-analysis score - best predicting score on basis of EA3 including 23andMe, UKB, etc. (r-squared around 12.7%)
mat beta_Gstar2 = (0.25 \ 0.18 \ 0.14 \ 0.1) // this is for the independent PGSs, typical R-squared for 23andMe score is 7%, 2nd is the EA2 score (R-squared ~ 3.2%) and 3th is the EA1 score (R-squared ~2%), also consider score of 1%

*READ PARAMETER VALUES (STANDARDIZED GENOTYPE + OUTCOME)
local sd_G = 1
local sd_error = 1

*DO REPLICATIONS AND COMPUTE (I) BIAS, (II) POWER, AND (III) RMSE
local ss = rowsof(samplesizes)
local betas = rowsof(beta_G)
local rows = `ss'

*INITIALIZE MATRICES
mat pvalues3 = J(`rows',`repl',0)
mat pvaluesivsc = J(`rows',`repl',0)
mat sig3 = J(`rows',`repl',0)
mat sigivsc = J(`rows',`repl',0)
mat power3 = J(`rows',1,0)
mat powerivsc = J(`rows',1,0)
mat bias3 = J(`rows',`repl',0)
mat biasivsc = J(`rows',`repl',0)
mat matbias3 = J(`rows',1,0)
mat matbiasivsc = J(`rows',1,0)
mat rmse3 = J(`rows',`repl',0)
mat rmseivsc = J(`rows',`repl',0)
mat rrmse3 = J(`rows',1,0)
mat rrmseivsc = J(`rows',1,0)
mat coef3 = J(`rows',`repl',0)
mat coefivsc = J(`rows',`repl',0)
mat rcoef3 = J(`rows',1,0)
mat rcoefivsc = J(`rows',1,0)
mat se3 = J(`rows',`repl',0)
mat seivsc = J(`rows',`repl',0)
mat rse3 = J(`rows',1,0)
mat rseivsc = J(`rows',1,0)
mat r23 = J(`rows',`repl',0)
mat r2ivsc = J(`rows',`repl',0)
mat rr23 = J(`rows',1,0)
mat rr2ivsc = J(`rows',1,0)
mat Fivsc = J(`rows',`repl',0)

forvalues j = 1/`betas' {
	di "Going over beta `j' out of `betas'"
	forvalues k = 1/`ss' {
		di "Looping sample size `k' out of `ss'"
		sca nsig = 0
		forvalues i = 1/`repl' {
			capture clear
			if mod(`i',250)==0{ //show only every 250 iterations
			    di "Doing iteration `i' out of `repl'"
				}
			local beta_G = beta_G[`j',1]
			local beta_Gstar = beta_Gstar[`j',1]
			local beta_Gstar2 = beta_Gstar2[`j',1]
			local N = samplesizes[`k',1]
			qui set obs `N'
			qui gen id = _n
						
			*DRAW OUTCOME Y AND PGS GTRUE FROM A BIVARIATE NORMAL DISTRIBUTION WITH CORRELATION EQUAL TO THE STANDARDIZED EFFECT SIZE (=CORRELATION)
			matrix M = 0, 0
			matrix sd = `sd_G', `sd_error'
			matrix V = (1, `beta_G' \ `beta_G', 1)
			drawnorm Gtrue Y, corr(V) means(M) sds(sd)
	
			*DRAW INDEPENDENT PGSs (PGS1 AND PGS2) AND META-ANALYSIS PGS (PGS 3) FROM NORMAL DISTRIBUTION 
			local sdindv = sqrt((sd[1,1]^2)*((`beta_G'^2)/(`beta_Gstar2'^2)-1))
			qui gen PGS1 = Gtrue + rnormal(0,`sdindv')
			qui sum PGS1
			qui replace PGS1 = (PGS1 - r(mean))/r(sd)
			qui gen PGS2 = Gtrue + rnormal(0,`sdindv')
			qui sum PGS2
			qui replace PGS2 = (PGS2 - r(mean))/r(sd)
			local sdmeta = sqrt((sd[1,1]^2)*((`beta_G'^2)/(`beta_Gstar'^2) - 1))
			qui gen PGS3 = Gtrue + rnormal(0,`sdmeta')
			qui sum PGS3
			qui replace PGS3 = (PGS3 - r(mean))/r(sd)
			
			*REGRESS OUTCOME ON PGS3
			qui reg Y PGS3, robust
			mat temp = r(table)
			mat bias3[`k',`i'] = (temp[1,1]-`beta_G')/`beta_G'
			mat rmse3[`k',`i'] = sqrt((temp[1,1]-`beta_G')^2)
			mat pvalues3[`k',`i'] = temp[4,1]
			mat sig3[`k',`i'] = (pvalues3[`k',`i'] < 0.05)
			mat r23[`k',`i'] = e(r2)
			
			*ORIV
			
			*HERE I TAKE THE COVARIANCE BETWEEN THE TWO POLYGENIC SCORES, WHICH IS AN ESTIMATE OF THE VARIANCE OF THE TRUE LATENT G
			qui correlate PGS1 PGS2, cov
			local varG = r(cov_12)
			
			*RESCALE THE INDEPENDENT POLYGENIC SCORES
			qui correlate PGS1 PGS2
			sca corr12 = r(rho)
			qui gen PGS1_st = PGS1
			qui replace PGS1_st = PGS1/sqrt(corr12) if corr12>0
			qui gen PGS2_st = PGS2
			qui replace PGS2_st = PGS2/sqrt(corr12) if corr12>0
			
			*REPLICATE THE DATASET AND ALTERNATE MAIN INDEPENDENT VARIABLE AND INSTRUMENTAL VARIABLE
			qui expand 2, generate(replicant)
			forvalues x=0/1 {
				qui gen constant`x' = replicant == `x'
				}
			qui gen mainVar_st = PGS1_st if replicant == 0
			qui replace mainVar_st = PGS2_st if replicant == 1
			qui gen instrument_st = PGS2_st if replicant == 0
			qui replace instrument_st = PGS1_st if replicant == 1
			
			qui reg mainVar_st instrument_st constant*, cluster(id) nocons
			qui testparm instrument_st
			mat Fivsc[`k',`i'] = e(F)
			
			qui ivregress 2sls Y (mainVar_st = instrument_st) constant*, cluster(id) nocons
			mat temp = r(table)
			mat biasivsc[`k',`i'] = (temp[1,1]-`beta_G')/`beta_G'
			mat rmseivsc[`k',`i'] = sqrt((temp[1,1]-`beta_G')^2)
			mat pvaluesivsc[`k',`i'] = temp[4,1]
			mat sigivsc[`k',`i'] = (pvaluesivsc[`k',`i'] < 0.05)
			mat r2ivsc[`k',`i'] = temp[1,1]^2 
			}  // end of forvalues i
		mata: st_matrix("power`j'3", rowsum(st_matrix("sig3")/`repl'))
		mata: st_matrix("power`j'ivsc", rowsum(st_matrix("sigivsc")/`repl'))
		
		mata: st_matrix("matbias`j'3", rowsum(st_matrix("bias3")/`repl'))
		mata: st_matrix("matbias`j'ivsc", rowsum(st_matrix("biasivsc")/`repl'))
		
		mata: st_matrix("rrmse`j'3", rowsum(st_matrix("rmse3")/`repl'))
		mata: st_matrix("rrmse`j'ivsc", rowsum(st_matrix("rmseivsc")/`repl'))
		
		mata: st_matrix("rr2`j'3", rowsum(st_matrix("r23")/`repl'))
		mata: st_matrix("rr2`j'ivsc", rowsum(st_matrix("r2ivsc")/`repl'))
		
		mata: st_matrix("matF`j'ivsc", rowsum(st_matrix("Fivsc")/`repl'))
		}  // end of forvalues k
	mat repl`j' = J(`ss',1,`repl')
	
	mat overall_power`j' = (samplesizes, power`j'3, power`j'ivsc, repl`j')
	matrix colnames overall_power`j' = sample_size power`j'3 power`j'ivsc replications
	
	mat overall_bias`j' = (samplesizes, matbias`j'3, matbias`j'ivsc, repl`j')
	matrix colnames overall_bias`j' = sample_size bias`j'3 bias`j'ivsc replications
	
	mat overall_rmse`j' = (samplesizes, rrmse`j'3, rrmse`j'ivsc, repl`j')
	matrix colnames overall_rmse`j' = sample_size rmse`j'3 rmse`j'ivsc replications
	
	mat overall_r2`j' = (samplesizes, rr2`j'3, rr2`j'ivsc, repl`j')
	matrix colnames overall_r2`j' = sample_size r2`j'3 r2`j'ivsc replications
} // end of forvalues j


*STORE RESULTS IN A MATRIX AND PLOT THE RESULTS
forvalues j = 1/`betas' {
	svmat overall_power`j', names(matcol)
	svmat overall_bias`j', names(matcol)
	svmat overall_rmse`j', names(matcol)
	svmat overall_r2`j', names(matcol)
	svmat matF`j'ivsc
    local beta`j' = beta_G[`j',1]
}

twoway (line overall_power1power13 overall_power1sample_size, lpattern(solid)) /*
	*/ (line overall_power1power1ivsc overall_power1sample_size, lpattern(longdash))  (line overall_power2power2ivsc overall_power2sample_size, lpattern(dash)) /*
	*/ (line overall_power3power3ivsc overall_power3sample_size, lpattern(dash_dot)), /*
	*/ scheme(s1mono) ytitle("Power") xtitle("Prediction sample size") legend(rows(3) lab(1 "Meta-analysis - {it:R{sup:2}}=12% (EA3)") /*
	*/ lab(2 "ORIV - {it:{&rho}}=0.25 (UKB)") lab(3 "ORIV - {it:{&rho}}=0.13 (EA2)") lab(4 "ORIV - {it:{&rho}}=0.08 (EA1)"))	
graph save "power_EA1-3.gph", replace
graph export "power_EA1-3.pdf", replace

twoway (line overall_bias1bias13 overall_bias1sample_size, lpattern(solid)) /*
	*/ (line overall_bias1bias1ivsc overall_bias1sample_size, lpattern(longdash))  (line overall_bias2bias2ivsc overall_bias2sample_size, lpattern(dash)) /*
	*/ (line overall_bias3bias3ivsc overall_bias3sample_size, lpattern(dash_dot)), /*
	*/ scheme(s1mono) ytitle("Relative bias") xtitle("Prediction sample size") ylabel(-0.4(0.1)0.1)/*
	*/ legend(rows(3) lab(1 "Meta-analysis - {it:R{sup:2}}=12% (EA3)") /*
	*/ lab(2 "ORIV - {it:{&rho}}=0.25 (UKB)") lab(3 "ORIV - {it:{&rho}}=0.13 (EA2)") lab(4 "ORIV - {it:{&rho}}=0.08 (EA1)"))
graph save "bias_EA1-3.gph", replace
graph export "bias_EA1-3.pdf", replace

twoway (line overall_rmse1rmse13 overall_rmse1sample_size, lpattern(solid)) /*
	*/ (line overall_rmse1rmse1ivsc overall_rmse1sample_size, lpattern(longdash))  (line overall_rmse2rmse2ivsc overall_rmse2sample_size, lpattern(dash)) /*
	*/ (line overall_rmse3rmse3ivsc overall_rmse3sample_size, lpattern(dash_dot)), /*
	*/ scheme(s1mono) ytitle("RMSE") xtitle("Prediction sample size") legend(rows(3) lab(1 "Meta-analysis - {it:R{sup:2}}=12% (EA3)") /*
	*/ lab(2 "ORIV - {it:{&rho}}=0.25 (UKB)") lab(3 "ORIV - {it:{&rho}}=0.13 (EA2)") lab(4 "ORIV - {it:{&rho}}=0.08 (EA1)"))
graph save "rmse_EA1-3.gph", replace
graph export "rmse_EA1-3.pdf", replace

twoway (line matF1ivsc1 overall_bias1sample_size if matF1ivsc1 < 200, lpattern(longdash)) (line matF2ivsc1 overall_bias1sample_size, lpattern(dash)) /*
	*/ (line matF3ivsc1 overall_bias1sample_size, lpattern(dash_dot)) (function y = 10, range(0 20000) lwidth(medthick)), /*
	*/ scheme(s1mono) ytitle("{it:F}-statistic") xtitle("Prediction sample size") legend(rows(2) /*
	*/ lab(1 "ORIV - {it:{&rho}}=0.25 (UKB)") lab(2 "ORIV - {it:{&rho}}=0.13 (EA2)") lab(3 "ORIV - {it:{&rho}}=0.08 (EA1)") lab(4 "{it:F}-statistic=10"))
graph save "F_EA1-3.gph", replace
graph export "F_EA1-3.pdf", replace







	





		
	


