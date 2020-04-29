library(oligo)
library(limma)
library(data.table)
library(plotly)
library(gdata)


#---- Generates expression matrix so the analysis can be carried on ----

signature_bad_samples <- paste0(signature_samples, "/", signature_bad_samples)

celFiles <- list.celfiles(signature_samples, full.names = TRUE)
celFiles <- celFiles[!(celFiles %in% signature_bad_samples)]

gexp <- read.celfiles(celFiles)

gexp <- oligo::rma(gexp)
gexp <- oligo::exprs(gexp)
gexp <- as.data.frame(gexp)

rm(signature_samples, celFiles, signature_bad_samples)


#---- Defines contrast conditions that will be used on model definition ----

conditions <- rep("control", 8)

gse37418_groups <- read.table(file = "./input_files/GSE37418_groups.txt", 
                              sep = "\t", fill = T, header = T)

gse37418_groups <- gse37418_groups[seq(2, 152, 2), c("Group", "Source.name")]
rownames(gse37418_groups) <- NULL

gse37418_groups <- gse37418_groups[-c(3, 49, 51, 52, 75, 76), ] 

conditions <- append(conditions, as.character(gse37418_groups$Source.name))

gexp <- gexp[ , - which(conditions == "U")]

conditions <- conditions[ - which(conditions == "U")]

rm(gse37418_groups)


#---- Principal component analysis allows checking how data is clustering ----

gexp2 <- t(as.matrix(gexp))

gexp.pca2 <- prcomp(gexp2, center = TRUE, scale. = TRUE)

pca <- gexp.pca2$x

pc1 <- pca[, 1]
pc2 <- pca[, 2]
pc3 <- pca[, 3]

plot_ly(x = pc1, y = pc2, type = "scatter", mode = "markers", 
        color = as.factor(conditions), 
        size = I(30),
        colors = c("#14b331", "#0050b3", "#b30024", "#ffbe0a", "#00d6d3"))

rm(gexp2, pca)


#---- Preprocessed signature objects are required by succeeding scripts ----

conditions <- tolower(conditions)
subgroups <- c("wnt", "shh", "g3", "g4", "g34")

for(subgroup in subgroups){
  
  if(subgroup == "g34")
    subgroup_gexp <- gexp[, which(conditions == "g3" | conditions == "g4")]
  else
    subgroup_gexp <- gexp[, which(conditions == subgroup)]
  
  subgroup_gexp <- cbind(subgroup_gexp, gexp[, which(conditions == "control")])
  
  assign(paste0(subgroup, "_gexp", sep = ""), subgroup_gexp)
  
  save(list = paste0(subgroup, "_gexp"), 
       file = paste0("./rdata_files/signature/", subgroup, "_gexp.RData"))
}



