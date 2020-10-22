library(gdata)
library(tidyverse)

#source("0a_quality_control.R")

rm(list = ls())

#source("0b_tfs.R")

rm(list = ls())


subgroups <- "no_wnt"

#---- 1_signature_preprocessing.R ----

signature_samples <- "./input_files/signature_complete_set"

signature_bad_samples <- c('GSM80618.CEL', 'GSM918580_mbt006-u133v2.CEL', 
                           'GSM918626_mbt127-u133v2.CEL', 'GSM918628_mbt136-u133v2.CEL',
                           'GSM918629_mbt140-u133v2.CEL', 'GSM918652_tbm111-u133v2.CEL', 
                           'GSM918653_tbm143-u133v2.CEL')

source("1_signature_preprocessing.R")

print("Script 1: Finished.")

#rm(list = ls())
gdata::keep(subgroups, sure = T)

#---- 2_signature_dictionary.R ----

preprocessed_signatures <- c("./rdata_files/signature/wnt_gexp",
                             "./rdata_files/signature/shh_gexp",
                             "./rdata_files/signature/g3_gexp",
                             "./rdata_files/signature/g4_gexp",
                             "./rdata_files/signature/g34_gexp",
                             "./rdata_files/signature/no_wnt_gexp")


for(rdata_file in preprocessed_signatures){
  
  load(paste0(rdata_file, ".RData"))
  subgroup_gexp <- eval(parse(text = ls(pattern = ".*gexp$")))
  
  source("2_signature_dictionary.R")
  
  gdata::keep(rdata_file, preprocessed_signatures, subgroups, sure = TRUE)
}

print("Script 2: Finished.")

#rm(list = ls())
gdata::keep(subgroups, sure = T)

#---- 3_signature.R ----

subgroups <- "no_wnt"

subgroup_labels <- list(wnt = c(rep("case", 7), rep("control", 8)),
                        shh = c(rep("case", 10), rep("control", 8)), 
                        g3 = c(rep("case", 15), rep("control", 8)), 
                        g4 = c(rep("case", 36), rep("control", 8)), 
                        g34 = c(rep("case", 51), rep("control", 8)),
                        no_wnt = c(rep("case", 61), rep("control", 8)))


for(subgroup in subgroups){
  
  source("3_signature.R")
  
  gdata::keep(subgroup_labels, subgroups, subgroup, sure = TRUE)
  
}

print("Script 3: Finished.")

gdata::keep(subgroups, sure = T)

#---- "4_network_dictionary.R" ----

bad_samples_list <- 
  list(wnt = c("GSM2261712.CEL", "GSM2262154.CEL", "GSM2262232.CEL"),
       shh = c("GSM2261629.CEL", "GSM2261633.CEL", "GSM2261689.CEL", 
               "GSM2261954.CEL", "GSM2261982.CEL", "GSM2262027.CEL",
               "GSM2262066.CEL", "GSM2262074.CEL", "GSM2262084.CEL",
               "GSM2262099.CEL", "GSM2262115.CEL", "GSM2262156.CEL",
               "GSM2262165.CEL", "GSM2262178.CEL", "GSM2262189.CEL",
               "GSM2262233.CEL", "GSM2262278.CEL", "GSM2262284.CEL",
               "GSM2262286.CEL", "GSM2262295.CEL"),
       g3 = c("GSM2261560.CEL", "GSM2261881.CEL", "GSM2262097.CEL", 
              "GSM2262155.CEL", "GSM2262229.CEL"),
       g4 = c("GSM2261610.CEL", "GSM2261640.CEL", "GSM2261831.CEL",
              "GSM2261884.CEL", "GSM2261985.CEL", "GSM2262013.CEL", 
              "GSM2262041.CEL", "GSM2262042.CEL", "GSM2262077.CEL",
              "GSM2262148.CEL", "GSM2262149.CEL", "GSM2262161.CEL",
              "GSM2262162.CEL", "GSM2262163.CEL", "GSM2262166.CEL",
              "GSM2262219.CEL", "GSM2262231.CEL", "GSM2262264.CEL",
              "GSM2262268.CEL", "GSM2262276.CEL", "GSM2262303.CEL",
              "GSM2262313.CEL"),
       g34 = c("GSM2261560.CEL", "GSM2261881.CEL", "GSM2262097.CEL", 
               "GSM2262155.CEL", "GSM2262229.CEL", "GSM2261610.CEL", 
               "GSM2261640.CEL", "GSM2261831.CEL", "GSM2262313.CEL",
               "GSM2261884.CEL", "GSM2261985.CEL", "GSM2262013.CEL", 
               "GSM2262041.CEL", "GSM2262042.CEL", "GSM2262077.CEL",
               "GSM2262148.CEL", "GSM2262149.CEL", "GSM2262161.CEL",
               "GSM2262162.CEL", "GSM2262163.CEL", "GSM2262166.CEL",
               "GSM2262219.CEL", "GSM2262231.CEL", "GSM2262264.CEL",
               "GSM2262268.CEL", "GSM2262276.CEL", "GSM2262303.CEL"),
       no_wnt = c("GSM2261629.CEL", "GSM2261633.CEL", "GSM2261689.CEL", 
                  "GSM2261954.CEL", "GSM2261982.CEL", "GSM2262027.CEL",
                  "GSM2262066.CEL", "GSM2262074.CEL", "GSM2262084.CEL",
                  "GSM2262099.CEL", "GSM2262115.CEL", "GSM2262156.CEL",
                  "GSM2262165.CEL", "GSM2262178.CEL", "GSM2262189.CEL",
                  "GSM2262233.CEL", "GSM2262278.CEL", "GSM2262284.CEL",
                  "GSM2262286.CEL", "GSM2262295.CEL",
                  "GSM2261560.CEL", "GSM2261881.CEL", "GSM2262097.CEL", 
                  "GSM2262155.CEL", "GSM2262229.CEL",
                  "GSM2261610.CEL", "GSM2261640.CEL", "GSM2261831.CEL",
                  "GSM2261884.CEL", "GSM2261985.CEL", "GSM2262013.CEL", 
                  "GSM2262041.CEL", "GSM2262042.CEL", "GSM2262077.CEL",
                  "GSM2262148.CEL", "GSM2262149.CEL", "GSM2262161.CEL",
                  "GSM2262162.CEL", "GSM2262163.CEL", "GSM2262166.CEL",
                  "GSM2262219.CEL", "GSM2262231.CEL", "GSM2262264.CEL",
                  "GSM2262268.CEL", "GSM2262276.CEL", "GSM2262303.CEL",
                  "GSM2262313.CEL"))


for(subgroup in subgroups){
  
  source("4_network_dictionary.R")
  
  gdata::keep(subgroup, subgroups, bad_samples_list, sure = TRUE)
  
}

print("Script 4: Finished.")

gdata::keep(subgroups, sure = T)

#---- "5_rtn.R" ----

for(subgroup in subgroups){
  
  source("5_rtn.R")

  gdata::keep(subgroup, subgroups, sure = TRUE)
  
  print(paste0("Script 5: Finished for", subgroup, "."))
  
}

gdata::keep(subgroups, sure = T)


#---- "6_survival_analysis.R" ----

#subgroups <- c("wnt", "shh", "g3", "g4", "g34")

bad_samples_list <- 
  list(wnt = c("GSM2261712.CEL", "GSM2262154.CEL", "GSM2262232.CEL"),
       shh = c("GSM2261629.CEL", "GSM2261633.CEL", "GSM2261689.CEL", 
               "GSM2261954.CEL", "GSM2261982.CEL", "GSM2262027.CEL",
               "GSM2262066.CEL", "GSM2262074.CEL", "GSM2262084.CEL",
               "GSM2262099.CEL", "GSM2262115.CEL", "GSM2262156.CEL",
               "GSM2262165.CEL", "GSM2262178.CEL", "GSM2262189.CEL",
               "GSM2262233.CEL", "GSM2262278.CEL", "GSM2262284.CEL",
               "GSM2262286.CEL", "GSM2262295.CEL"),
       g3 = c("GSM2261560.CEL", "GSM2261881.CEL", "GSM2262097.CEL", 
              "GSM2262155.CEL", "GSM2262229.CEL"),
       g4 = c("GSM2261610.CEL", "GSM2261640.CEL", "GSM2261831.CEL",
              "GSM2261884.CEL", "GSM2261985.CEL", "GSM2262013.CEL", 
              "GSM2262041.CEL", "GSM2262042.CEL", "GSM2262077.CEL",
              "GSM2262148.CEL", "GSM2262149.CEL", "GSM2262161.CEL",
              "GSM2262162.CEL", "GSM2262163.CEL", "GSM2262166.CEL",
              "GSM2262219.CEL", "GSM2262231.CEL", "GSM2262264.CEL",
              "GSM2262268.CEL", "GSM2262276.CEL", "GSM2262303.CEL",
              "GSM2262313.CEL"),
       g34 = c("GSM2261560.CEL", "GSM2261881.CEL", "GSM2262097.CEL", 
               "GSM2262155.CEL", "GSM2262229.CEL", "GSM2261610.CEL", 
               "GSM2261640.CEL", "GSM2261831.CEL", "GSM2262313.CEL",
               "GSM2261884.CEL", "GSM2261985.CEL", "GSM2262013.CEL", 
               "GSM2262041.CEL", "GSM2262042.CEL", "GSM2262077.CEL",
               "GSM2262148.CEL", "GSM2262149.CEL", "GSM2262161.CEL",
               "GSM2262162.CEL", "GSM2262163.CEL", "GSM2262166.CEL",
               "GSM2262219.CEL", "GSM2262231.CEL", "GSM2262264.CEL",
               "GSM2262268.CEL", "GSM2262276.CEL", "GSM2262303.CEL"))


for(subgroup in subgroups){
  
  source("6_survival_analysis.R")
  
  gdata::keep(subgroup, subgroups, bad_samples_list, sure = T)

  paste0("Finished for ", subgroup)
}

gdata::keep(subgroups, sure = T)

#---- "7_rtn_duals.R" ----

#subgroups <- c("wnt", "shh", "g3", "g4", "g34")

for(subgroup in subgroups){
  source("7_rtn_duals.R")
}

