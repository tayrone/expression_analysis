# This scripts loads a list of human transcription factors made
# available by Lambert et al (2018).
# http://humantfs.ccbr.utoronto.ca/download.php

library(tidyverse)

tfs <- readr::read_csv("./input_files/tfs.csv")

tfs <- tfs[, c("HGNC symbol", "Is TF?")]

tfs <- tfs[tfs$`Is TF?` == "Yes", "HGNC symbol"]

tfs <- as.character(tfs$`HGNC symbol`)

save(tfs, file = "./rdata_files/tfs.RData")


#---- Lista de TFs de 2013 ----
# 
# library(Fletcher2013b)
# library(gdata)
# 
# data(miscellaneous)
# 
# gdata::keep(tfs, sure = T)
# 
# tfs <- names(tfs)
# 
# save(tfs, file = "./rdata_files/tfs.RData")


