# ============================================================
# script de estructura (PCA + DAPC)
# 17 sep 2026
#
# 1ra vez (checkpoint 1, línea base): dejar n_pcs_fijo <- NULL para que xvalDapc decida el número de PCs.
# Checkpoints siguientes (2 a 5): poner en n_pcs_fijo el número que ya decidiste en el checkpoint 1, para correr xvalDapc cada vez.
# ============================================================


# corriendo con todo el dataset para checar estructura antes del filtrado

library(adegenet)

reportar_tiempo <- function(etiqueta, t_referencia) {
  transcurrido <- as.numeric(difftime(Sys.time(), t_referencia, units = "mins"))
  cat(sprintf("[tiempo] %s: %.2f min\n", etiqueta, transcurrido))
  Sys.time()
}
t0 <- Sys.time()
t_inicio_total <- t0

graficar_en_png <- function(ruta, ancho, alto, expr_grafico) {
  png(ruta, width = ancho, height = alto)
  on.exit(dev.off())
  expr_grafico()
}

buscar_K_con_salvaguarda <- function(gl, n_pca, k_max, etiqueta, carpeta_salida) {
  grupos <- find.clusters(gl, n.pca = n_pca, max.n.clust = k_max, choose.n.clust = FALSE)
  
  bic_df <- data.frame(K = seq_along(grupos$Kstat), BIC = as.numeric(grupos$Kstat))
  write.csv(bic_df, paste0(carpeta_salida, "bic_vs_K.csv"), row.names = FALSE)
  graficar_en_png(paste0(carpeta_salida, "bic_vs_K.png"), 800, 500, function() {
    plot(bic_df$K, bic_df$BIC, type = "b", xlab = "K", ylab = "BIC", main = "BIC vs K")
  })
  
  k_elegido <- length(unique(grupos$grp))
  if (k_elegido == k_max) {
    message("[AVISO SERIO] El K elegido automáticamente (", k_elegido, ") choca exactamente ",
            "contra el techo (max.n.clust=", k_max, "). El BIC probablemente no encontró un ",
            "mínimo real dentro del rango probado -- revisa bic_vs_K.png antes de confiar en este K.")
  }
  grupos
}


# ========
# 1. rutas
# ========

ruta_bfile <- "/mnt/data/sur/users/spacheco/data/teosinte/T3604_33929_all"   # sin extensión
ruta_metadata <- "/mnt/data/sur/users/spacheco/data/teosinte/data_teosinte.csv"  # para asignar accession a cada individuo
ruta_carpeta_salida <- "/mnt/data/sur/users/spacheco/results/sep_2026/teo/17_sep/"
 
k_max <- 400   # ajustar según qué tan grande sea el subset (con pocos individuos, bajarlo)


# NULL = correr xvalDapc para decidir el número de PCs (1ra vez).
# Un número = saltarse xvalDapc y usar ese número
n_pcs_fijo <- NULL   

xval_n_pca_max <- 200
xval_n_rep <- 30


dir.create(ruta_carpeta_salida, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(ruta_carpeta_salida)) {
  stop("No se pudo crear/acceder a la carpeta de salida: ", ruta_carpeta_salida,
       " -- revisa mayúsculas/minúsculas y permisos antes de seguir.")
}


# ============================================================
# 2. Convertir .bed/.bim/.fam a .raw con plink2 (tabs -> espacios,porque adegenet::read.PLINK necesita espacios)
# ============================================================

ruta_raw_original <- paste0(ruta_bfile, ".raw")
ruta_raw_convertido <- paste0(ruta_bfile, "_convertido.raw")

if (!file.exists(ruta_raw_convertido)) {
  message("Generando .raw con plink2 --export A...")
  resultado <- system2("plink2", c("--bfile", ruta_bfile, "--export", "A", "--out", ruta_bfile))
  if (resultado != 0 || !file.exists(ruta_raw_original)) {
    stop("plink2 no corrió correctamente. No se puede continuar sin el .raw generado.")
  }
  message("Convirtiendo delimitador (tabs -> espacios)...")
  lineas <- readLines(ruta_raw_original)
  writeLines(gsub("\t", " ", lineas), ruta_raw_convertido)
} else {
  message(ruta_raw_convertido, " ya existe, no se regenera.")
}
t0 <- reportar_tiempo("conversión a .raw", t0)


# ============================================================
# 3. Leer como genlight y asignar población (Accession)
# ============================================================

gl <- read.PLINK(ruta_raw_convertido, quiet = TRUE)
cat("Individuos leídos:", nInd(gl), "| Loci:", nLoc(gl), "\n")

meta <- read.csv(ruta_metadata, stringsAsFactors = FALSE)
orden_gl <- data.frame(IID = indNames(gl))
cruce <- merge(orden_gl, meta[, c("Sample_name", "Accession")],
               by.x = "IID", by.y = "Sample_name", all.x = TRUE, sort = FALSE)
cruce <- cruce[match(indNames(gl), cruce$IID), ]

faltantes <- sum(is.na(cruce$Accession))
if (faltantes > 0) {
  message("[AVISO] ", faltantes, " individuo(s) del genlight no encontraron Accession -- revisar.")
}

pop(gl) <- cruce$Accession
cat("Poblaciones (Accession) distintas asignadas:", length(unique(pop(gl))), "\n")
t0 <- reportar_tiempo("leer genlight + asignar Accession", t0)


# ======
# 4. PCA
# =====

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
pca <- glPca(gl, nf = 10, parallel = (n_cores > 1), n.cores = n_cores)

varianza_pca <- pca$eig / sum(pca$eig) * 100
write.csv(data.frame(eje = paste0("PC", seq_along(varianza_pca)), varianza_pct = varianza_pca)[1:10, ],
          paste0(ruta_carpeta_salida, "pca_varianza.csv"), row.names = FALSE)
write.csv(data.frame(Accession = pop(gl), pca$scores),
          paste0(ruta_carpeta_salida, "pca_scores.csv"), row.names = FALSE)
cat("Varianza explicada, primeros 3 ejes:", round(varianza_pca[1:3], 2), "\n")
t0 <- reportar_tiempo("PCA (glPca)", t0)

library(ggplot2)
n_pops_distintas <- length(unique(pop(gl)))
df_pca <- data.frame(PC1 = pca$scores[, 1], PC2 = pca$scores[, 2], Accession = pop(gl))
p_pca <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Accession)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(x = paste0("PC1 (", round(varianza_pca[1], 1), "%)"),
       y = paste0("PC2 (", round(varianza_pca[2], 1), "%)"),
       title = "PCA -- coloreado por Accession") +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())
if (n_pops_distintas > 30) {
  p_pca <- p_pca + theme(legend.position = "none")
}
ggsave(paste0(ruta_carpeta_salida, "pca_plot.png"), p_pca, width = 9, height = 7, dpi = 150)


# ============================================================
# 5. Número de PCs para el DAPC: fijo, o decidido por xvalDapc
# ============================================================

if (is.null(n_pcs_fijo)) {
  cat("\n--- n_pcs_fijo es NULL -- corriendo xvalDapc para decidir ---\n")
  mat <- as.matrix(gl)
  mat[is.na(mat)] <- 0
  xval <- xvalDapc(mat, pop(gl), n.pca.max = xval_n_pca_max, n.rep = xval_n_rep, xval.plot = FALSE)
  n_pcs <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
  cat("Número de PCs sugerido por xvalDapc:", n_pcs, "\n")
  t0 <- reportar_tiempo(paste0("xvalDapc (hasta ", xval_n_pca_max, " PCs, ", xval_n_rep, " rep.)"), t0)
} else {
  n_pcs <- n_pcs_fijo
  cat("\n--- Usando n_pcs_fijo =", n_pcs, "(sin correr xvalDapc) ---\n")
}

max_pcs_posible <- min(nInd(gl), nLoc(gl)) - 1
if (n_pcs > max_pcs_posible) {
  message("[AVISO] n_pcs (", n_pcs, ") es mayor al máximo posible (", max_pcs_posible, "). Se ajusta.")
  n_pcs <- max_pcs_posible
}


# ========
# 6. DAPC
# =======

grupos <- buscar_K_con_salvaguarda(gl, n_pcs, k_max, "dapc", ruta_carpeta_salida)
t0 <- reportar_tiempo(paste0("find.clusters (", n_pcs, " PCs, K=1 a ", k_max, ")"), t0)

if (length(unique(grupos$grp)) < 2) {
  stop("K encontrado fue 1 -- el DAPC no puede correr con un solo grupo. Revisa bic_vs_K.png.")
}

dapc_obj <- dapc(gl, pop = grupos$grp, n.pca = n_pcs, n.da = length(unique(grupos$grp)) - 1)
saveRDS(dapc_obj, paste0(ruta_carpeta_salida, "dapc_", n_pcs, "pcs.rds"))

k_final <- length(unique(grupos$grp))
cat("\nDAPC con", n_pcs, "PCs -- K encontrado:", k_final, "\n")
cat("Proporción de reasignación correcta:", round(summary(dapc_obj)$assign.prop, 4), "\n")
t0 <- reportar_tiempo(paste0("dapc (", n_pcs, " PCs)"), t0)

paleta <- rainbow(k_final)
graficar_en_png(paste0(ruta_carpeta_salida, "dapc_scatter.png"), 900, 700, function() {
  scatter(dapc_obj, col = paleta, bg = "white", cstar = 0,
          legend = (k_final <= 30), posi.leg = "topright",
          scree.pca = TRUE, posi.pca = "bottomleft",
          main = paste0("DAPC (", n_pcs, " PCs, K=", k_final, ")"))
})
graficar_en_png(paste0(ruta_carpeta_salida, "compoplot.png"), 1200, 600, function() {
  compoplot(dapc_obj, col = paleta, legend = (k_final <= 30),
            main = paste0("Compoplot tipo ADMIXTURE (", n_pcs, " PCs, K=", k_final, ")"))
})

cat(sprintf("\n[tiempo] TOTAL del script: %.2f min\n", as.numeric(difftime(Sys.time(), t_inicio_total, units = "mins"))))
