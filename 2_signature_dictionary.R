# Definition of dictionary for signature data.

library(dplyr)
library(matrixStats)

#---- Platform data is required for building probe-gene dictionary ----

dictionary <- read.table(file = "./input_files/GPL570-55999.txt", header = T, 
                         skip = 16, fill = T, sep = "\t", quote = "")

dictionary <- 
  dictionary %>%
    select(c("ID", "Gene.Symbol")) %>%
    filter(Gene.Symbol != "")


dictionary$Gene.Symbol <- 
  sapply(strsplit(as.character(dictionary$Gene.Symbol), " /// "), `[`, 1)


colnames(dictionary) <- c("probe_id", "gene")
rownames(dictionary) <- dictionary$probe_id


#---- Ordering dictionary and subgroup_gexp makes code easier to follow ----

subgroup_gexp <- 
  subgroup_gexp[rownames(subgroup_gexp) %in% rownames(dictionary), ]

subgroup_gexp <- subgroup_gexp[sort(rownames(subgroup_gexp)), ]

dictionary <- dictionary[sort(rownames(dictionary)), ]

all(rownames(subgroup_gexp) == rownames(dictionary))


#--- Multiple probes which are mapped to the same symbol are 
# selected by greatest variance, so there are 
# no repeated gene symbol observation values in the end ----

duplicated_symbols <- dictionary$gene[which(duplicated(dictionary$gene))]
duplicated_symbols <- unique(duplicated_symbols)

probes_from_duplicated_genes <- 
  dictionary[which(dictionary$gene %in% duplicated_symbols), ]

mtx <- as.matrix(subgroup_gexp)
mtx <- mtx[as.character(probes_from_duplicated_genes$probe_id),  ]

probes_from_duplicated_genes <- cbind(probes_from_duplicated_genes, 
                                      variance = rowVars(mtx))

probes_with_greatest_variance <- 
  probes_from_duplicated_genes %>% 
    group_by(gene) %>%
    filter(variance == max(variance)) %>%
    # There are cases of several probes, mapped to the same gene, 
    # with the same variance, so the next line of code is necessary.
    filter(row_number(gene) == 1) %>%
    arrange(gene)

rm(mtx)

#---- Final subgroup_gexp and dictionary objects are created for next script ----

signature_gexp <- 
  subset(subgroup_gexp, 
         !(rownames(subgroup_gexp) %in% probes_from_duplicated_genes$probe_id) |  
         (rownames(subgroup_gexp) %in% probes_with_greatest_variance$probe_id))

signature_dictionary <- subset(dictionary, 
                               rownames(dictionary) %in% rownames(signature_gexp))

save(signature_gexp, signature_dictionary, 
     file = paste0(rdata_file, "_dictionary.RData"))

