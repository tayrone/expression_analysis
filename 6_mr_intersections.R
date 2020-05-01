#---- This script build a set intersection graph. It takes as
# input all transcription factors identified as master regulators,
# for each one of the analyses. -----

library(reshape2)
library(UpSetR)
library(tidyverse)


subgroups <- c("wnt", "shh", "g3", "g4", "g34")

load_mrs <- function(subgroup){
  
  load(paste0("./rdata_files/network/", subgroup, "_rtn.RData"))
  
  assign(paste0(subgroup, "_mrs"), 
         data.frame(tf = tna.get(rtna, what = "mra")[["Regulon"]], 
                    analysis = subgroup))
}

mrs_list <- lapply(subgroups, load_mrs)

rm(rtna, rtni, subgroups, load_mrs)

#---- Wide format is required by upset method ----

complete_map <- do.call("rbind", mrs_list)

complete_map_wide <- dcast(complete_map, tf~analysis, length, 
                           value.var = "analysis")

rownames(complete_map_wide) <- complete_map_wide$tf
complete_map_wide <- select(complete_map_wide, -tf)


upset(complete_map_wide, order.by = "freq")
 
