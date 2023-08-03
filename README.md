# ORIV
Stata syntax files to perform simulation and empirical analyses for the "Overcoming attenuation bias in regressions using polygenic indices" paper, which can be found here:
(https://www.nature.com/articles/s41467-023-40069-4). The main replication files are in Stata syntax; please also find an R package to implement ORIV, which is contributed by Hyeokmoon Kweon, h.kweon@vu.nl. 

1. ORIV runs Design X.do

These DO-files (with varying numbers X) use the simulated datasets generated in the GNAMES package (see https://github.com/devlaming/gnames) and perform the relevant analyses to create Figures 1-5 and Table 1 in the paper, and Supplementary Table 1-8 + Supplementary Figure 2 and 3 in the Supplementary Information.  

2. UKB analysis Y scores - residualized.do

These DO-files (with varying names Y) perform the empirical analyses in the sibling subsample of the UKB, used to construct Tables 2 and 3, and Figure 6 in the main text, and Supplementary Tables 10-12 in the Supplementary Information. How the UKB data is constructed and the genetic data is processed can be found on https://github.com/DilnozaM/Rank-Concordance-of-PGI/tree/main/GWAS%26PGS 
