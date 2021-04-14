# ORIV
Stata syntax files to perform simulation and empirical analyses for the "Stop meta-analyzing, start instrumenting" paper, which can be found here:
https://www.biorxiv.org/content/10.1101/2021.04.09.439157v1

1. Simulation - Github version.do

This DO-file performs the simulations to create Figure 2 in the paper. 
It generates 3 polygenic scores: 1 meta-analysis polygenic score, and 2 independent polygenic scores used in ORIV. For varying sample sizes and degrees of measurement error it computes the relative bias, Root Mean Squared Error, and empirical power for (i) an OLS regression of the outcome on the meta-analysis based score; and (ii) an ORIV regression where the two independent polygenic scores are used as instruments for each other. 

Figures 3 and 5-8 are constructed in a very similar manner (changing a few parameters) and therefore the code is not reproduced here. It is of course, available upon request (hvankippersluis at ese.eur.nl). 

2. ORIV UKB EA - Github version.do

This DO-file performs the empirical analysis in the sibling subsample of the UKB, used to construct Figure 2 (left panel) and Table 1. Again it compares a meta-analysis based score (UKB & 23andMe) for Educational Attainment (EA) to ORIV where we use two independent polygenic scores (two-sample UKB and 23andMe, or split-sample UKB). It applies meta-analysis based scores and ORIV both between families and within families. 

Figure 2 (right panel) and Table 2 are constructed very similarly using height as the outcome and is therefore not reproduced here. The code is of course available upon request (hvankippersluis at ese.eur.nl). 

We are currently working on also releasing the syntax used to construct the UKB sample and the polygenic scores. This will be added later. 
