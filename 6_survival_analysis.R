library(readODS)
library(tidyverse)
library(pheatmap)
library(RTNsurvival)
library(oligo)

load(paste0("./rdata_files/network/", subgroup, "_rtn.RData"))


#---- Organizaing metadata makes everything easier later on ----

metadata <- read_ods("./input_files/data.ods", skip = 1)

colnames(metadata) <- tolower(colnames(metadata))

metadata$subgroup <- tolower(metadata$subgroup)

metadata$subgroup[metadata$subgroup == "group3"] <- "g3"
metadata$subgroup[metadata$subgroup == "group4"] <- "g4"


#---- Filtering and value setting make graphs more readable----

if(subgroup == "g34"){
  metadata <- metadata[(metadata$subgroup == "g3" | metadata$subgroup == "g4"), ]
}else{
  metadata <- metadata[metadata$subgroup == subgroup, ]
}

samples <- paste0("./input_files/", subgroup)
subgroup_samples <- oligoClasses::list.celfiles(samples, full.names = FALSE)

rownames(metadata) <- subgroup_samples

bad_samples <- bad_samples_list[[subgroup]]

metadata <- metadata[!(subgroup_samples %in% bad_samples), ]

#metadata <- metadata %>% drop_na(age, dead, `os (years)`)

#"event" and "time" variables are required, so we make them easier to identify
metadata <- dplyr::rename(metadata, event = dead, time = `os (years)`)

metadata <- dplyr::select(metadata, time, event, everything())


#---- All preprocessed data is employed from now on ----

rtns <- tni2tnsPreprocess(rtni, survivalData = metadata, 
                          time = 1, event = 2)
rtns <- tnsGSEA2(rtns)

rtns <- tnsCox(rtns)
tnsPlotCox(rtns)

rtns <- tnsKM(rtns)

#tnsPlotGSEA2(rtns, "GSM2261539.CEL", regs = "AATF")


#---- Plotting a heatmap for all hazardous regulons might 
# provide interesting insights ----

cox <- rtns@results$Cox$Table

hazardous_regulons <- cox[(cox$Lower95 > 1 | cox$Upper95 < 1), "Regulons"]

enrichmentScores <- tnsGet(rtns, "regulonActivity")
survival_data <- tnsGet(rtns, "survivalData")

annotation_col <- data.frame(annotation = survival_data[, "subgroup"], 
                             row.names = rownames(survival_data))

colors <- list(subgroup = c(g3 = "#507512", g4 = "#15248b"))

graphics.off()

png(paste0("./survival_plots/", subgroup, "/heatmap.png"), 
           units = "in", width = 18, height = 10, res = 500)

pheatmap(t(enrichmentScores$dif), 
         annotation_col = annotation_col,
         show_colnames = FALSE,
         show_rownames = F,
         annotation_legend = T,
         treeheight_row = 0,
         treeheight_col = 10,
         annotation_colors = colors[1],
         annotation_names_col = F,
         annotation_names_row = F,
         drop_levels = F)

graphics.off()

#---- Kaplan-Meier plots for hazardous regulons are
# the principal goal of this script ----

if(length(hazardous_regulons) > 0){
  tnsPlotKM(rtns, regs = hazardous_regulons, 
            fname = paste0(subgroup, "_km"), 
            fpath = paste0("./survival_plots/", subgroup), 
            plotpdf = T, plotbatch = T, xlab = "Years")
}

save(rtns, hazardous_regulons, file = paste0("./rdata_files/survival/",
                                             subgroup, "_survival.RData"))

