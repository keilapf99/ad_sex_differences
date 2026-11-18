# Title : Master Regulator Analysis - VIPER (molecular signature by deseq2)
# Author: Keila Adonai Pérez-Fuentes
# Date: May 14th 2026
# Objective: 


# Packages ---------------------------------------------------------------------

#BiocManager::install("viper", ask = FALSE, update = FALSE)
library(viper)


# 0. Data ----------------------------------------------------------------------
alzheimer <- read.table("/home/aperez/mra_alzheimer_sex/expression_alzheimer/ROSMAP_RNAseq_female_AD_counts_CRB.tsv",
                        header = TRUE, 
                        check.names = F)
# Convertir la columna feature en rownames
rownames(alzheimer) <- alzheimer$feature
# Eliminar la columna feature
alzheimer$feature <- NULL
# Convertir a matriz
alzheimer <- as.matrix(alzheimer)

expr_diff <- read.table("/home/aperez/mra_alzheimer_sex/deseq2/CRB_females_AD_vs_control_DEGs.csv",
                        header = TRUE,
                        sep = ",",
                        row.names = 1,
                        stringsAsFactors = FALSE)

# 1. Generate regulon ----------------------------------------------------------
regulon <- aracne2regulon("/home/aperez/alzheimer_sex/tf_networks_ad_sex/ROSMAP_RNAseq_female_AD_counts_CRB.network.txt", alzheimer)
# > regulon <- aracne2regulon("/home/aperez/alzheimer_sex/tf_networks_ad_sex/ROSMAP_RNAseq_female_AD_counts_CRB.network.txt", alzheimer)
# 
# Loading the dataset...
# Generating the regulon objects...
# number of iterations= 41 
# mu: -0.628521973755279, 0, 0.681995681306541. sigma: 0.13221927231519, 0.371083536314722, 0.118724733762793
# Hubo 50 o más avisos (use warnings() para ver los primeros 50)

# 2. Generate molecular signature ----------------------------------------------
signature <- expr_diff$stat
names(signature) <- rownames(expr_diff)

signature <- signature[!is.na(signature)]
signature <- signature[names(signature) != ""]
signature <- signature[!duplicated(names(signature))]
signature <- sort(signature, decreasing = TRUE)

# 3. Generate null model -------------------------------------------------------
nperm <- 1000
set.seed(1)

nullmodel <- sapply(seq_len(nperm), function(i) {
  perm <- signature
  names(perm) <- sample(names(signature))
  perm
})
colnames(nullmodel) <- paste0("perm_", seq_len(ncol(nullmodel)))

# 5. MRS with VIPER ------------------------------------------------------------
mrs <- msviper(signature, regulon, nullmodel, verbose = T,
               cores = 75)
summary(mrs, 20)


# 
# 
# # -----------------------------------------------------------------------------#
# # Leer tabla de jobs
# # -----------------------------------------------------------------------------#
# 
# jobs <- read.table(
#   "/home/aperez/alzheimer_sex/roots_job_mra.tsv",
#   header = TRUE,
#   sep = "\t",
#   stringsAsFactors = FALSE
# )
# 
# # Carpeta base donde se guardarán resultados
# base_outdir <- "/home/aperez/alzheimer_sex/MRA_results"
# dir.create(base_outdir, showWarnings = FALSE)
# 
# # -----------------------------------------------------------------------------#
# # Loop principal
# # -----------------------------------------------------------------------------#
# 
# for (i in seq_len(nrow(jobs))) {
#   
#   job <- jobs[i, ]
#   
#   outdir <- file.path(base_outdir, job$job_id)
#   
#   message("===================================================")
#   message("=== Corriendo job: ", job$job_id)
#   message("  AD      : ", job$ad)
#   message("  Health  : ", job$health)
#   message("  Regulon : ", job$regulon)
#   message("  Output  : ", outdir)
#   message("===================================================")
#   
#   # ---------------- Validaciones ---------------- #
#   
#   if (!file.exists(job$ad)) {
#     stop("❌ No existe AD: ", job$ad)
#   }
#   if (!file.exists(job$health)) {
#     stop("❌ No existe Health: ", job$health)
#   }
#   if (!file.exists(job$regulon)) {
#     stop("❌ No existe Regulon: ", job$regulon)
#   }
#   
#   dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#   
#   # ---------------- Leer expresión ---------------- #
#   
#   message("→ Leyendo expresión")
#   ad <- as.matrix(
#     read.table(job$ad, header = TRUE, row.names = 1, sep = "\t")
#   )
#   
#   health <- as.matrix(
#     read.table(job$health, header = TRUE, row.names = 1, sep = "\t")
#   )
#   
#   storage.mode(ad) <- "double"
#   storage.mode(health) <- "double"
#   
#   # ---------------- Regulon ---------------- #
#   
#   message("→ Generando regulon")
#   regulon <- aracne2regulon(job$regulon, ad)
#   
#   # ---------------- Signature ---------------- #
#   
#   message("→ Calculando signature")
#   signature <- rowTtest(ad, health)
#   signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) *
#                   sign(signature$statistic))[, 1]
#   
#   # ---------------- Null model ---------------- #
#   
#   message("→ Generando null model")
#   nullmodel <- ttestNull(ad, health, per = 1000, repos = TRUE)
#   
#   # ---------------- MRA ---------------- #
#   
#   message("→ Corriendo msVIPER")
#   mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)
#   
#   # ---------------- Guardar ---------------- #
#   
#   saveRDS(mrs, file = file.path(outdir, "mrs.rds"))
#   
#   write.table(
#     summary(mrs, 20),
#     file = file.path(outdir, "top20.txt"),
#     sep = "\t",
#     quote = FALSE,
#     row.names = TRUE
#   )
#   
#   # ---------------- PLOT (NUEVO) ---------------- #
#   
#   pdf(file.path(outdir, "mrs_plot.pdf"))
#   plot(mrs)
#   dev.off()
#   
#   message("✅ Job terminado: ", job$job_id)
# }
# 
# message("🎉 Todos los jobs finalizados")


save.image("/home/aperez/mra_alzheimer_sex/mra.RData")
