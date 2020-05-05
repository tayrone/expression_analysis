library(RTNduals)

load(paste0("./rdata_files/network/", subgroup, "_rtn.RData"))
rm(rtna)

rmbr <- tni2mbrPreprocess(rtni)
rmbr <- mbrAssociation(rmbr)

mbrGet(rmbr, what = "summary")

overlap <- mbrGet(rmbr, what = "dualsOverlap")
correlation <- mbrGet(rmbr, what = "dualsCorrelation")

save(rmbr, overlap, correlation,
     file = paste0("./rdata_files/duals/", subgroup, "_duals.RData"))
