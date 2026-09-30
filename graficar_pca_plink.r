# =====================
# graficar_pca_plink.r
# 29 sep 2026
#
# ===================


library(ggplot2)


# ========
# 1. parámetros
# ========

umbral_r2 <- "0.1" #cambiar aquí por cada r2

ruta_carpeta_pca   <- paste0("/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/teocintle/25_sep/r2_", umbral_r2, "/plink/")
prefijo_pca        <- paste0(ruta_carpeta_pca, "pca_r2_", umbral_r2)
prefijo_rel        <- paste0(ruta_carpeta_pca, "rel_r2_", umbral_r2)
metadata     <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"


# ============================
# 2. leer eigenvec / eigenval
# ===========================
# el encabezado del .eigenvec empieza con "#FID", con comment.char="#" (el default de read.table) se pierde el encabezado
# y la primera fila de datos se toma como nombres de columna. Forzar comment.char="".

eigenvec <- read.table(paste0(prefijo_pca, ".eigenvec"), header = TRUE,
                       comment.char = "", stringsAsFactors = FALSE)
names(eigenvec)[1] <- "FID"  # el "#" queda pegado al nombre de la 1a columna

eigenval <- scan(paste0(prefijo_pca, ".eigenval"), quiet = TRUE)
n_pcs_reportados <- length(eigenval)
cat("PCs reportados por --pca:", n_pcs_reportados, "\n")


# ============================================================
# 3. varianza total real -- traza de la matriz de relacionalidad
# ============================================================
# rel.square es una matriz cuadrada de n_individuos x n_individuos, sin encabezado, 
# en el mismo orden que el .rel.id (que a su vez esta en el mismo orden que el .fam de entrada, no el mismo orden que el .eigenvec, que puede reordenar) 
# Solo necesitamos la diagonal, y la diagonal no depende del orden de filas/columnas mientras sea
# consistente (fila i y columna i son el mismo individuo)

rel <- as.matrix(read.table(paste0(prefijo_rel, ".rel"), header = FALSE))
varianza_total <- sum(diag(rel))
cat("Varianza total (traza de la GRM):", varianza_total, "\n")

pct_varianza <- eigenval / varianza_total * 100
cat("% varianza por PC:\n")
print(round(pct_varianza, 2))

pct_acumulado_reportados <- sum(pct_varianza)
cat("Los", n_pcs_reportados, "PCs reportados capturan en total",
    round(pct_acumulado_reportados, 1), "% de la varianza real\n")


# ===========================
# 4. unir con metadata taxon 
# ===========================

metadatos <- read.csv(metadata, stringsAsFactors = FALSE)

eigenvec_taxon <- merge(eigenvec, metadatos[, c("Sample_name", "Accession", "Taxon", "Race")],
                        by.x = "IID", by.y = "Sample_name", all.x = TRUE)

n_sin_taxon <- sum(is.na(eigenvec_taxon$Taxon))
if (n_sin_taxon > 0) {
  warning(n_sin_taxon, " individuos del .eigenvec no encontraron match de Taxon en ",
          "data_teosinte.csv -- revisar si el IID coincide exactamente con Sample_name.")
}

# ==========================
# 5. graficar PCA coloreado por Taxon
# ===================================

p <- ggplot(eigenvec_taxon, aes(x = PC1, y = PC2, color = Taxon)) +
  geom_point(alpha = 0.7, size = 1.5) +
  labs(
    title = paste0("PCA con podado a r2=", umbral_r2),
    x = paste0("PC1 (", round(pct_varianza[1], 2), "% de la varianza total)"),
    y = paste0("PC2 (", round(pct_varianza[2], 2), "% de la varianza total)")
  ) +
  theme_minimal(base_size = 12)

ggsave(paste0(ruta_carpeta_pca, "pca_taxon_r2_", umbral_r2, ".png"),
       p, width = 9, height = 7, dpi = 300)

write.csv(eigenvec_taxon, paste0(ruta_carpeta_pca, "pca_con_taxon_r2_", umbral_r2, ".csv"),
          row.names = FALSE)

cat("\nListo. Grafica y CSV guardados en:", ruta_carpeta_pca, "\n")


