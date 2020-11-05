library(dplyr)
library(oligo)
library(data.table)
library(matrixStats)

#---- Dictionary is built using probe ids and their respective gene symbols ----

platform_data <- 
  as.data.frame(fread(file  = "./input_files/GPL11532-32230.txt", 
                      header = TRUE, skip = 12)) 

dictionary <- as.data.frame(platform_data$ID, stringsAsFactors = FALSE)

gene_symbol <- platform_data$gene_assignment
gene_symbol <- strsplit(gene_symbol, " // ") 
gene_symbol <- as.data.frame(sapply(gene_symbol, `[`, 2))

dictionary <- cbind(dictionary, gene_symbol) 

colnames(dictionary) <- c("probe_id", "gene_symbol")

# Be aware some probes are not mapped to any gene
dictionary <- as.data.frame(dictionary[!is.na(dictionary$gene_symbol), ])

rownames(dictionary) <- dictionary$probe_id

rm(gene_symbol, platform_data) 


#---- Expression data is loaded and normalized, making it 
# possible to subsequently assemble the network ---- 

samples <- paste0("./input_files/", subgroup)

bad_samples <- unlist(bad_samples_list[subgroup])
bad_samples <- sapply(bad_samples, function(x) paste0(samples, "/", x))

celFiles <- list.celfiles(samples, full.names = TRUE)
celFiles <- celFiles[!(celFiles %in% bad_samples)]

gexp <- read.celfiles(celFiles)

gexp <- oligo::rma(gexp) 
gexp <- oligo::exprs(gexp)
gexp <- as.data.frame(gexp)

rm(samples, celFiles, bad_samples)


#---- Makes sure gexp and dictionary are both ordered in the same fashion.
# This makes all comparisons much more straightfoward ----

gexp <- gexp[rownames(gexp) %in% rownames(dictionary), ]

gexp <- gexp[sort(rownames(gexp)), ]

dictionary <- dictionary[sort(rownames(dictionary)), ]

all(rownames(gexp) == rownames(dictionary))


#--- Multiple probes which are mapped to the same symbol are 
# selected by greatest variance, so there are 
# no repeated gene symbol observation values in the end ----

duplicated_symbols <- 
  dictionary$gene_symbol[which(duplicated(dictionary$gene_symbol))]

duplicated_symbols <- unique(as.character(duplicated_symbols))

probes_from_duplicated_genes <- 
  dictionary[which(dictionary$gene_symbol %in% duplicated_symbols), ]

mtx <- as.matrix(gexp)
mtx <- mtx[as.character(probes_from_duplicated_genes$probe_id),  ]

probes_from_duplicated_genes <- cbind(probes_from_duplicated_genes, 
                                      variance = rowVars(mtx))

probes_with_greatest_variance <- 
  probes_from_duplicated_genes %>% 
  group_by(gene_symbol) %>%
  filter(variance == max(variance)) %>%
  # There are cases of several probes, mapped to the same gene, presenting
  # the same variance value. This makes the next line of code necessary.
  filter(row_number(gene_symbol) == 1) %>% 
  arrange(gene_symbol)

rm(mtx)


#---- After removing all duplicated gene symbols, dictionary is filled again
# with only one probe id for each symbol, so a unique expression value
# for each gene is obtained ----

repeated <- which(dictionary$gene_symbol %in% duplicated_symbols)
dictionary <- dictionary[-repeated, ]

probes_with_greatest_variance <- 
  select(probes_with_greatest_variance, "probe_id", "gene_symbol")

dictionary <- bind_rows(dictionary, probes_with_greatest_variance)

rownames(dictionary) <- dictionary$probe_id

rm(repeated)


#---- Sorts both gexp and dictionary to save them on RData file ----

gexp <- gexp[rownames(gexp) %in% rownames(dictionary), ]

gexp <- gexp[sort(rownames(gexp)), ]

dictionary <- dictionary[sort(rownames(dictionary)), ]

all(rownames(gexp) == rownames(dictionary))

save(dictionary, gexp, file = paste0("./rdata_files/network/",
                                     subgroup, "_network_dictionary.RData"))
