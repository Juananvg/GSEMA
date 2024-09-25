### R code from vignette source 'GSEMA.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: GSEMA.Rnw:50-52
###################################################
library(GSEMA)
data("simulatedData")


###################################################
### code chunk number 3: GSEMA.Rnw:85-87
###################################################
#List of expression matrices
listMatrices <- list(study1Ex, study2Ex, study3Ex, study4Ex, study5Ex)


###################################################
### code chunk number 4: GSEMA.Rnw:92-94
###################################################
listPhenodata <- list(study1Pheno, study2Pheno, study3Pheno, study4Pheno,
    study5Pheno)


###################################################
### code chunk number 5: GSEMA.Rnw:109-110
###################################################
study1Pheno$Condition


###################################################
### code chunk number 6: GSEMA.Rnw:115-117
###################################################
head(names(GeneSets))
GeneSets[[1]]


###################################################
### code chunk number 7: GSEMA.Rnw:142-156
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
### code chunk number 8: GSEMA.Rnw:174-175
###################################################
objectMApathSim <- filteringPaths(objectMApathSim, threshold = 0.65)


###################################################
### code chunk number 9: GSEMA.Rnw:209-211
###################################################
results <- metaAnalysisESpath(objectMApath = objectMApathSim,
    measure = "limma", typeMethod = "REM", missAllow = 0.3, proportionData = 1)


###################################################
### code chunk number 10: GSEMA.Rnw:216-218
###################################################
results[1,]
results[nrow(results),]


###################################################
### code chunk number 11: GSEMA.Rnw:258-261
###################################################
res <- heatmapPaths(objectMA=objectMApathSim, results,
            scaling = "zscor", regulation = "all",breaks=c(-2,2),
            fdrSig = 0.05, comES_Sig = 1, numSig=20, fontsize = 5)


###################################################
### code chunk number 12: GSEMA.Rnw:270-273
###################################################
Effects <- calculateESpath(objectMApath = objectMApathSim, measure = "limma")
Effects[[1]][1,]
Effects[[2]][1,]


###################################################
### code chunk number 13: GSEMA.Rnw:278-279
###################################################
sessionInfo()


