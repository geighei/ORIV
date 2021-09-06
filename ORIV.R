# By Hyeokmoon Kweon, h.kweon@vu.nl

######################################################################################
# function GORIV: Runs ORIV 

    # Note: the function depends on "fixest", "data.table" packages

# Inputs:
    # fm: formula object to specify a model. Do not include a PGS here, but only covariates. 
    # PGS1: character object of first PGS name
    # PGS2: character object of second PGS name
    # data: data.frame or data.table 
    # IID: character object of indivdiual ID column name
    # FID: character object of family ID column name, Default to NULL. If provided, SE will be clustered at family, too. Has to be supplied for within-family estimation.  
    # resid: logical, indicating whether the target outcome should be residualized and standardized first. Default to TRUE
    # within: logical, indicating within-family estimation. Default to FALSE

# Output:
    # fixest object from "fixest" package
     # see https://www.rdocumentation.org/packages/fixest/versions/0.8.4/topics/feols

GORIV <- function(fm, PGS1, PGS2, data, IID, FID=NULL, resid=TRUE, within=FALSE){

library(data.table)        
data <- data.table(data)
data <- na.omit(data[, .SD, .SDcols=unique(c(all.vars(fm), IID, FID, PGS1, PGS2))])
if (within) data <- data[get(FID) %in% data[, .N, by=get(FID)][N>=2][[1]] ]

# Residualize / Set up model
if (resid==TRUE){
        data[, Y := lm(fm, data)$residuals]
        data[, Y := scale(Y)]

        if (within==TRUE){
                FM = paste0("Y ~ 1 | ", FID, "^rep | PGS_MAIN ~ PGS_IV")
        } else {
                FM = "Y ~ 1 | rep | PGS_MAIN ~ PGS_IV"
        }

} else if (resid==FALSE){
        fm = as.character(fm)        

        if (within==TRUE){
                FM = paste0(fm[2], "~", fm[3], "| ", FID, "^rep | PGS_MAIN ~ PGS_IV")
        } else {
                FM = paste0(fm[2], "~", fm[3], "| rep | PGS_MAIN ~ PGS_IV")
        }
}

# standardize observed PGS
data[, c(PGS1, PGS2) := lapply(.SD, scale), .SDcols=c(PGS1, PGS2)]

# compute scaling factor
if (within){
        data[, fam_n := .N, by=FID]
        data <- data[fam_n >= 2]
        data[, PGS1_dm := get(PGS1) - mean(get(PGS1)), by=FID]
        data[, PGS2_dm := get(PGS2) - mean(get(PGS2)), by=FID]
        R <- sqrt(data[, cor(PGS1_dm, PGS2_dm)])
} else {
        R <- sqrt(data[, cor(get(PGS1), get(PGS2))])
}

# scale PGS
data[, (PGS1) := get(PGS1) / R]
data[, (PGS2) := get(PGS2) / R]
cat("Scaling factor =", R, "\n")

# stack the data
N = nrow(data)
data = rbind(data, data)
data$rep = c(rep(0,N), rep(1,N))

data[, PGS_MAIN := ifelse(rep==0, get(PGS1), get(PGS2))]
data[, PGS_IV := ifelse(rep==0, get(PGS2), get(PGS1))]

# Run estimation
if (within){
        return(fixest::feols(as.formula(FM), data=data, se="twoway", cluster=c(FID, IID)))
} else if (!is.null(FID)){
        return(fixest::feols(as.formula(FM), data=data, se="twoway", cluster=c(FID, IID)))
} else {
        return(fixest::feols(as.formula(FM), data=data, se="cluster", cluster=IID))
}
}


### Example ###

# UKB <- fread("../TEMP/UKB_temp.csv")

# pc <- paste0("pc", 1:20)
# pc <- paste0(pc, collapse="+")
# fm = as.formula(paste0("edu ~ male*factor(yob) + geno + ", pc))

# GORIV(fm, "EA_UKB_PGS", "EA_23_PGS", data=UKB, IID="n_eid", FID="familyID")
# GORIV(fm, "EA_UKB_PGS", "EA_23_PGS", data=UKB, IID="n_eid", FID="familyID", resid=FALSE)

# GORIV(fm, "EA_UKB_PGS", "EA_23_PGS", data=UKB, IID="n_eid", FID="familyID", within=TRUE)
# GORIV(fm, "EA_UKB_PGS", "EA_23_PGS", data=UKB, IID="n_eid", FID="familyID", within=TRUE, resid=FALSE)







