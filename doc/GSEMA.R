### R code from vignette source 'GSEMA.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: GSEMA.Rnw:49-51
###################################################
library(GSEMA)
data("simulatedData")


###################################################
### code chunk number 3: GSEMA.Rnw:84-86
###################################################
#List of expression matrices
listMatrices <- list(study1Ex, study2Ex, study3Ex, study4Ex, study5Ex)


###################################################
### code chunk number 4: GSEMA.Rnw:91-93
###################################################
listPhenodata <- list(study1Pheno, study2Pheno, study3Pheno, study4Pheno,
    study5Pheno)


###################################################
### code chunk number 5: GSEMA.Rnw:108-109
###################################################
study1Pheno$Condition


###################################################
### code chunk number 6: GSEMA.Rnw:114-116
###################################################
head(names(GeneSets))
GeneSets[[1]]


###################################################
### code chunk number 7: GSEMA.Rnw:141-155
###################################################
listMatrices <- list(study1Ex, study2Ex, study3Ex, study4Ex, study5Ex)
listPhenodata <- list(study1Pheno, study2Pheno, study3Pheno, study4Pheno,
    study5Pheno)
phenoGroups <- c("Condition","Condition", "Condition", "Condition", "Condition")
phenoCases <- list("Case", "Case", "Case", "Case", "Case")
phenoControls <- list("Healthy", "Healthy", "Healthy", "Healthy", "Healthy")
objectMApathSim <- createObjectMApath(
    listEX = listMatrices,
    listPheno = listPhenodata, namePheno = phenoGroups,
    expGroups = phenoCases, refGroups = phenoControls,
    geneSets = GeneSets,
    pathMethod = "Zscore",
    n.cores = 1,
    internal.n.cores = 1, minSize = 7)


###################################################
### code chunk number 8: GSEMA.Rnw:173-174
###################################################
objectMApathSim <- filteringPaths(objectMApathSim, threshold = 0.65)


###################################################
### code chunk number 9: GSEMA.Rnw:208-210
###################################################
results <- metaAnalysisESpath(objectMApath = objectMApathSim,
    measure = "limma", typeMethod = "REM", missAllow = 0.3, proportionData = 1)


###################################################
### code chunk number 10: GSEMA.Rnw:215-217
###################################################
results[1,]
results[nrow(results),]


###################################################
### code chunk number 11: GSEMA.Rnw:257-260
###################################################
res <- heatmapPaths(objectMA=objectMApathSim, results,
            scaling = "zscor", regulation = "all",breaks=c(-2,2),
            fdrSig = 0.05, comES_Sig = 1, numSig=20, fontsize = 5)


###################################################
### code chunk number 12: GSEMA.Rnw:269-272
###################################################
Effects <- calculateESpath(objectMApath = objectMApathSim, measure = "limma")
Effects[[1]][1,]
Effects[[2]][1,]


###################################################
### code chunk number 13: GSEMA.Rnw:277-278
###################################################
sessionInfo()


