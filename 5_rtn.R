library(data.table)
library(RTN)

#subgroup <- "wnt"

load("./rdata_files/tfs.RData")
load(paste0("./rdata_files/signature/", subgroup, "_signature.RData"))
load(paste0("./rdata_files/network/", subgroup, "_network_dictionary.RData"))


#---- Check the amount of tfs that are recognized by dictionary object ----

tfs_ <- which(tfs %in% dictionary$gene_symbol)

tfs <- as.character(dictionary$gene_symbol[tfs_])

rm(tfs_)

#---- TNI processing ----

gexp <- as.matrix(gexp)

if(rownames(dictionary) == rownames(gexp))
  rownames(gexp) <- dictionary$gene_symbol

rtni <- tni.constructor(expData = gexp, regulatoryElements = tfs)

# Beware: this insctruction may take a really long time to process
rtni <- tni.permutation(rtni) 

rtni <- tni.bootstrap(rtni) # Unstable interactions are removed

# Removes the weakest links on TF-TF-target triangle interactions
rtni <- tni.dpi.filter(rtni) 

tni.get(rtni, what = "summary")
refnet <- tni.get(rtni, what = "refnet")
tnet <- tni.get(rtni, what = "tnet")

#g <- tni.graph(rtni) # Saves final network, ready to be plotted by RedeR package

#---- TNA processing ----

rtna <- tni2tna.preprocess(object = rtni, phenotype = pheno, hits = hits)

rtna <- tna.mra(rtna)

rtna <- tna.overlap(rtna)

rtna <- tna.gsea1(rtna, stepFilter = FALSE, nPermutations = 1000)

rtna <- tna.gsea2(rtna, tfs = names(tfs), nPermutations = 1000)

tna <- tna.get(rtna, what = "mra")

RegScores <- tni.gsea2(rtni)

save(rtni, rtna, file = paste("./rdata_files/network/", subgroup, "_rtn.RData"))


