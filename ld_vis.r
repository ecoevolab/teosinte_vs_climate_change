# ================
# ld_vis.r
# 25/28 sep 2026
#
# =============

library(snpStats)
library(LDheatmap)
library(ggplot2)


# ========
# 1. rutas
# ========

# para el de 1000, 25 sep
ruta_bfile <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_baby_subset_contiguo_1000"
ruta_carpeta_salida <- "/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/teocintle/25_sep/ld_viz_mejorada/"

# para el de 100, 28 sep
# ruta_bfile <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_baby_subset_contiguo_100"
# ruta_carpeta_salida <- "/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/teocintle/28_sep/ld_viz_mejorada/"

dir.create(ruta_carpeta_salida, recursive = TRUE, showWarnings = FALSE)


# =================
# 2. leer el bfile
# ================

datos_plink <- read.plink(paste0(ruta_bfile, ".bed"),
                          paste0(ruta_bfile, ".bim"),
                          paste0(ruta_bfile, ".fam"))

genotipos <- datos_plink$genotypes   # SnpMatrix, individuos x SNPs
mapa <- datos_plink$map              # columnas: chromosome, snp.name, cM, position, allele.1, allele.2

if (length(unique(mapa$chromosome)) > 1) {
  message("[AVISO] este bfile tiene mas de un cromosoma -- el heatmap y las distancias ",
          "fisicas no van a tener sentido biologico. Este script esta pensado para un ",
          "bfile de un solo cromosoma/tramo contiguo.")
}

cat("SNPs:", ncol(genotipos), "| Individuos:", nrow(genotipos), "\n")


# =================
# 3. heatmap de LD 
# =================

png(paste0(ruta_carpeta_salida, "ld_heatmap_paquete.png"), width = 5000, height = 5000,pointsize = 150)
resultado_heatmap <- LDheatmap(genotipos,
                               genetic.distances = mapa$position,
                               distances = "physical",
                               LDmeasure = "r",
                               title = paste0("LD heatmap ", basename(ruta_bfile)),
                               color = colorRampPalette(c("blue", "magenta", "red"))(100))
dev.off()
cat("Heatmap (paquete) guardado: ld_heatmap_paquete.png\n")

# misma logica de reutilizacion que en el script hecho a mano: guardar la
# matriz de r^2 para subconjuntarla despues por cada umbral de poda
saveRDS(list(r2 = resultado_heatmap$LDmatrix, bp = mapa$position, snp_id = mapa$snp.name),
        paste0(ruta_carpeta_salida, "r2_matriz_completa_paquete.rds"))


# ============================================================
# 4. LD decay: r^2 vs distancia fisica, reutilizando la misma matriz
# ============================================================

r2 <- resultado_heatmap$LDmatrix
bp <- mapa$position
indices <- which(upper.tri(r2), arr.ind = TRUE)
df_decay <- data.frame(
  distancia_bp = abs(bp[indices[, 1]] - bp[indices[, 2]]),
  r2 = r2[upper.tri(r2)]
)

p_decay <- ggplot(df_decay, aes(x = distancia_bp, y = r2)) +
  geom_point(alpha = 0.15, size = 0.5, color = "lightgreen") +
  geom_smooth(method = "loess", span = 0.2, color = "red", se = FALSE) +
  labs(title = paste0("LD decay ", basename(ruta_bfile)),
       x = "Distancia fisica (bp)", y = expression(R^2)) +
  theme_minimal(base_size = 12)
ggsave(paste0(ruta_carpeta_salida, "ld_decay_paquete.png"), p_decay, width = 8, height = 6, dpi = 150)
cat("LD decay (paquete) guardado: ld_decay_paquete.png\n")


