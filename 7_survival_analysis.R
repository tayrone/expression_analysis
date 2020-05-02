library("readODS")
library(tidyverse)

load(paste0("./rdata_files/network/", subgroup, "_rtn.RData"))


#---- Organizaing metadata makes everything easier later on ----

metadata <- read_ods("./input_files/data.ods", skip = 1)

colnames(metadata) <- tolower(colnames(metadata))

metadata$subgroup <- tolower(metadata$subgroup)

metadata$subgroup[metadata$subgroup == "group3"] <- "g3"
metadata$subgroup[metadata$subgroup == "group4"] <- "g4"


#---- Filtering and value setting makes graphs more readable----

#metadata <- metadata[!(metadata$study_id %in% bad_samples), ]

if(subgroup == "g34"){
  metadata <- metadata[(metadata$subgroup == "g3" | metadata$subgroup == "g4"), ]
}else{
  metadata <- metadata[metadata$subgroup == subgroup, ]
}

rownames(metadata) <- unlist(rtni@gexp)

filtered_metadata <- metadata %>% drop_na(age, dead, `os (years)`)

#selects event and time columns, which are required by pipeline
colnames(filtered_metadata)[6] <- "event"
colnames(filtered_metadata)[7] <- "time"

filtered_metadata <- select(filtered_metadata, event, time, everything())


#---- All preprocessed data is employed from now on ----

library(RTN)
library(RTNsurvival)

rtns <- tni2tnsPreprocess(rtni, survivalData = filtered_metadata, 
                          event = 1, time = 2)
rtns <- tnsGSEA2(rtns)

rtns <- tnsCox(rtns)
tnsPlotCox(rtns)

rtns <- tnsKM(rtns)

#tnsPlotGSEA2(rtns, "GSM2261539.CEL", regs = "AATF")


#---- Plotting a heatmap for all hazardous regulons might 
# provide interesting insights ----

library(pheatmap)

cox <- rtns@results$Cox$Table

hazardous_regulons <- cox[(cox$Lower95 > 1 | cox$Upper95 < 1), "Regulons"]

enrichmentScores <- tnsGet(rtns, "regulonActivity")
survival.data <- tnsGet(rtns, "survivalData")

annotation_col <- data.frame(annotation = survival.data[ , "Subtype"], 
                             row.names = rownames(survival.data))

pheatmap(t(enrichmentScores$dif), 
         annotation_col = annotation_col,
         show_colnames = FALSE, 
         annotation_legend = T)


#---- Kaplan-Meier plots for hazardous regulons are
# the principal goal of this script ----

tnsPlotKM(rtns, regs = hazardous_regulons, 
          fname = paste0("./survival_analysis_files/", subgroup, "_km_plots"), 
          fpath = ".", plotpdf = T, plotbatch = T)


save(rtns, hazardous_regulons, file = paste0("./survival_analysis_files/",
                                             subgroup, "_survival.RData"))

