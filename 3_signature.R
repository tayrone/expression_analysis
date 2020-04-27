library(limma)

#---- Set objects required by statistical analysis on this script ----

names(labels) <- colnames(signature_gexp)

design <- model.matrix(~0+labels)
colnames(design) <- c("case", "control")


#---- Limma statical methods are used to do differential expression 
# analysis and, therefore, to build pheno, pheno_ids and hits objects, 
# which are mnecessary on master regulator analysis ----

fit <- lmFit(signature_gexp, design)
contrasts <- makeContrasts(case-control, levels = design)

# Gene ranking in by differential expression evidence
ct.fit <- eBayes(contrasts.fit(fit, contrasts), 0.05)

res.fit <- decideTests(ct.fit, method = "separate", 
                       adjust.method = "BH", p.value = 0.01, lfc = 1.5)


signature <- data.frame(logFC = ct.fit$coef, p.value = ct.fit$p.value, 
                        degenes = unclass(res.fit), stringsAsFactors = FALSE)


#---- Hits are genes that did not present a neutral regulation ----

signature_hits <- subset(signature, (signature$case...control.2 != 0))

hits <- as.character(rownames(signature_hits))

hits <- data.frame(probe_id = hits, stringsAsFactors = F)
hits <- merge(hits, signature_dictionary, by = "probe_id")
hits <- as.character(hits$gene)


#---- Pheno contains expression valur for all gene identifiers ----

pheno_ids <- data.frame(probe_id  = rownames(signature), 
                        exp_values = signature$case...control)

pheno_ids <- merge(pheno_ids, signature_dictionary, by = "probe_id")

pheno_ids <- pheno_ids[!duplicated(pheno_ids$probe_id), ]

pheno <- as.numeric(pheno_ids$exp_values)

signature_dictionary <- signature_dictionary[order(signature_dictionary$probe_id), ]
#all(names(pheno) == rownames(signature_dictionary))

names(pheno) <- signature_dictionary$gene


#---- Hits and pheno are mandatory on the workflow (pheno_ids is optional) ----

save(hits, pheno, pheno_ids, 
     file = paste0("./rdata_files/signature/", subgroup, "_signature.RData"))
