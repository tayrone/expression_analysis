#---- This script generates quality reports for all data assessed 
# in this workflow ----

library(oligo)
library(arrayQualityMetrics)

input_data <- c("control", "case", "wnt", "shh", "group3", "group4") 

for(dataset in input_data){

  samples <- paste("./input_files/", dataset, sep = "") 
  
  celFiles <- list.celfiles(samples, full.names = TRUE)
  
  gexp <- read.celfiles(celFiles)
  
  gexp <- oligo::rma(gexp)
  
  arrayQualityMetrics(expressionset = gexp, force = TRUE,
                      outdir = paste("./quality_reports/", 
                                     dataset, "_quality_report", sep = ""))
}


