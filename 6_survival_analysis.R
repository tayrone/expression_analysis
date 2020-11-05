library(readODS)
library(tidyverse)
library(pheatmap)
library(RTNsurvival)
library(oligo)

load(paste0("./rdata_files/network/", subgroup, "_rtn.RData"))

rm(rtna)

#---- Organizaing metadata makes everything easier later on ----

metadata <- read_ods("./input_files/data.ods", skip = 1)

colnames(metadata) <- tolower(colnames(metadata))

metadata$subgroup <- tolower(metadata$subgroup)

metadata$subgroup[metadata$subgroup == "group3"] <- "g3"
metadata$subgroup[metadata$subgroup == "group4"] <- "g4"


mb_gsm_dictionary <- read_delim("./input_files/mb_gsm.txt", skip = 1,
                                col_names = c("acession", "study_id"), 
                                delim = "\t")[seq(2, 1526, 2), 1:2]

metadata <- full_join(metadata, mb_gsm_dictionary, by = "study_id") %>% 
            select(acession, everything())


#---- Filtering and value setting make graphs more readable----

if(subgroup == "g34"){
  metadata <- metadata[(metadata$subgroup == "g3" | metadata$subgroup == "g4"), ]
}else if(subgroup == "no_wnt"){
  metadata <- metadata[(metadata$subgroup == "g3" | metadata$subgroup == "g4" | 
                        metadata$subgroup == "shh"), ]
}else{
  metadata <- metadata[metadata$subgroup == subgroup, ]
}


bad_samples <- str_remove(bad_samples_list[[subgroup]], ".CEL")

metadata <- filter(metadata, !(acession %in% bad_samples))


metadata <- metadata %>% 
  dplyr::rename(event = dead, time = `os (years)`) %>% 
  select(time, event, everything()) %>% 
  filter(!is.na(time) & !is.na(event)) %>% 
  tibble::column_to_rownames("acession")

rownames(metadata) <- paste0(rownames(metadata), ".CEL")


#---- All preprocessed data is employed from now on ----

#metadata$subgroup <- as.numeric(as.factor(metadata$subgroup))

rtns <- tni2tnsPreprocess(rtni, survivalData = metadata, 
                          time = 1, event = 2,
                          pAdjustMethod = "BH", keycovar = "subgroup")

rm(rtni)
#bonferroni results in a too-restricted output, compared to BH

rtns <- tnsGSEA2(rtns)

rtns <- tnsCox(rtns)
tnsPlotCox(rtns)

rtns <- tnsKM(rtns)

#tnsPlotGSEA2(rtns, "GSM2261539.CEL", regs = "AATF")


#---- Plotting a heatmap for all hazardous regulons might 
# provide interesting insights ----

cox <- rtns@results$Cox$Table

km <- rtns@results$KM$Table

hazardous_regulons <- cox[(cox$Lower95 > 1 | cox$Upper95 < 1) &
                          cox$Adjusted.Pvalue <= 0.06, ]

km <- km[km$Adjusted.Pvalue <= 0.06, ]

hazardous_regulons <- base::intersect(hazardous_regulons$Regulons, km$Regulons)

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
            fname = "hazardous_plots", 
            fpath = paste0("./survival_plots/", subgroup), 
            plotpdf = T, plotbatch = T, xlab = "Years")
}



save(rtns, hazardous_regulons, file = paste0("./rdata_files/survival/",
                                             subgroup, "_survival.RData"))

