test_ES <- function() {
  library(GSEMA)
  checkEqualsNumeric(metaAnalysisDE(objectMApathSim, typeMethod = "REM", 
                    missAllow = 0.3, proportionData = 1)[1,2],
                    0.2329579, tolerance=1.0e-6)}
