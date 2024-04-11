\name{simulatedData}
\alias{simulatedData}
\alias{GeneSets}
\alias{objectMApathSim}
\alias{study1Ex}
\alias{study2Ex}
\alias{study3Ex}
\alias{study4Ex}
\alias{study5Ex}
\alias{study1Pheno}
\alias{study2Pheno}
\alias{study3Pheno}
\alias{study4Pheno}
\alias{study5Pheno}
\docType{data}
\title{GSEMA synthetic data}
\description{
\itemize{
    \item study1Ex, study2Ex, study3Ex, study4Ex, study5Ex:
    five expression matrices.
    \item study1Pheno, study2Pheno, study3Pheno, study4Pheno, study5Pheno:
    five phenodata objects.
    \item GeneSets: a list of gene sets with each element are the genes that
    belong to a pathway.
    \item objectMApathSim: the meta-analysis object created from the different
    expression matrices and phenodatas.
}
}
\usage{data(simulatedData)}
\format{
    arrays (\code{study1Ex}, \code{study2Ex}, \code{study3Ex}, \code{study4Ex},
    \code{study5Ex}), data.frames (\code{study1Pheno}, \code{study2Pheno},
    \code{study3Pheno}, \code{study4Pheno}, \code{study5Pheno}),
    list vectors (\code{GeneSets}) and list of nested lists
    (\code{objectMApathSim}).
}
\source{
    study1Ex, study2Ex, study3Ex, study4Ex, study5Ex, study1Pheno, study2Pheno,
    study3Pheno, study4Pheno and study5Pheno are synthetic.
    GeneSets was created from the information of MsigDB  database by adding
    an invented pathway

    objectMApathSim was created after using createObjectMA GSEMA function to
    the different studies objects and with the information of GeneSets object

    maObjectDif was created after using createObjectMA DExMA function
    to the listExpressionSets object.
    maObject was obtained after using allSameID function to maObjectDif
    function.
}

\references{
Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P.
The Molecular Signatures Database (MSigDB) hallmark gene set collection.
Cell Syst. 2015;1: 417–425
}
\keyword{objectMApath}
