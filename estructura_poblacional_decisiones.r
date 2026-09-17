# ============================================================
# PCA + DAPC ANTES de cualquier filtro (línea base para comparar contra los
# estructura_poblacional_decisiones.r
# decisiones porque apenas decidiré si vale la pena usar los 3300 componentes con los
# que se quedan ellos o si hay mejores parámteros
# ============================================================

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
  
} # es para que se ejecute sin que se rompa antes de termianr la imagen

# ========
# 1. rutas
# ========
# cluster

# ruta_bfile <- "/mnt/data/sur/users/spacheco/data/teosinte/T3604_33929_all"   # sin extensión
# ruta_metadata <- "/mnt/data/sur/users/spacheco/data/teosinte/data_teosinte.csv"  # para asignar accession a cada individuo
# ruta_carpeta_salida <- "/mnt/data/sur/users/spacheco/results/sep_2026/teo/9_Sep/"

# local cuando corrrí el baby subset para ver las gráficas
ruta_bfile <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_baby_subset"   # sin extensión
ruta_metadata <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"  # para asignar accession a cada individuo
ruta_carpeta_salida <- "/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/teocintle/10_sep_2026/"


# Rango de K a probar en find.clusters (el paper probó 1 a 40)
k_min <- 1
k_max <- 8 #ajustar

# Número de PCs "grande" para replicar el enfoque del paper
# (ahí usaron 3,300 pero no se explica el porqué)
n_pcs_grande <- 3400  # le puse un poco más grande para notar alguna diferencia

# configuración de xvalDapc (validación cruzada para elegir PCs)
xval_n_pca_max <- 200   #rango de PCs a explorar; ajustar si hace falta
xval_n_rep <- 30        # default de adegenet; bajar si tarda demasiado


# ============================================================
# 2. convertir .bed/.bim/.fam a .raw con plink2, y arreglar el delimitador (plink2 usa tabs; adegenet::read.PLINK necesita espacios)
# ============================================================

dir.create(ruta_carpeta_salida, recursive = TRUE, showWarnings = FALSE)
ruta_raw_original <- paste0(ruta_bfile, ".raw")
ruta_raw_convertido <- paste0(ruta_bfile, "_convertido.raw")

if (!file.exists(ruta_raw_convertido)) {
  message("Generando .raw con plink2 --export A...")
  resultado <- system2("plink2", c("--bfile", ruta_bfile, "--export", "A", "--out", ruta_bfile))
  
  if (resultado != 0 || !file.exists(ruta_raw_original)) {
    stop("plink2 no corrió correctamente. No se puede continuar sin el .raw generado.")
  }
  
  message("Convirtiendo delimitador (tabs -> espacios) para que adegenet lo lea...")
  lineas <- readLines(ruta_raw_original)
  lineas <- gsub("\t", " ", lineas)
  writeLines(lineas, ruta_raw_convertido)
} else {
  message(ruta_raw_convertido, " ya existe, no se regenera. Bórralo a mano si cambiaste el .bed/.bim/.fam.")
}
t0 <- reportar_tiempo("conversión a .raw", t0)


# ============================================================
# 3. Leer como genlight y asignar población (Accession)
# ============================================================

gl <- read.PLINK(ruta_raw_convertido, quiet = TRUE)
cat("Individuos leídos:", nInd(gl), "| Loci:", nLoc(gl), "\n")

# Asignar Accession como población, cruzando por Sample_name.
# Esto asume que los IID del .fam coinciden con Sample_name de data_teosinte.csv
meta <- read.csv(ruta_metadata, stringsAsFactors = FALSE)
orden_gl <- data.frame(IID = indNames(gl))
cruce <- merge(orden_gl, meta[, c("Sample_name", "Accession")],
               by.x = "IID", by.y = "Sample_name", all.x = TRUE, sort = FALSE)
# el merge no garantiza el orden, reordenar según indNames(gl)
cruce <- cruce[match(indNames(gl), cruce$IID), ]

faltantes <- sum(is.na(cruce$Accession))
if (faltantes > 0) {
  message("[AVISO] ", faltantes, " individuo(s) del genlight no encontraron ",
          "Accession en ", ruta_metadata, " -- revisar antes de continuar.")
}

pop(gl) <- cruce$Accession
cat("Poblaciones (Accession) distintas asignadas:", length(unique(pop(gl))), "\n")
t0 <- reportar_tiempo("leer genlight + asignar Accession", t0)


# =======
# 4. PCA 
# =======

# glPca() sí soporta varios núcleos (parallel=TRUE, n.cores).
# find.clusters() y xvalDapc() NO tienen esa opción, corren en un solo núcleo. 
# por eso solo esta línea usa varios CPUs; el resto del script no se acelera con más --cpus-per-task.

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat("Usando", n_cores, "núcleo(s) para glPca (detectado de SLURM_CPUS_PER_TASK)\n")
pca <- glPca(gl, nf = 10, parallel = (n_cores > 1), n.cores = n_cores)

varianza_pca <- pca$eig / sum(pca$eig) * 100
write.csv(
  data.frame(eje = paste0("PC", seq_along(varianza_pca)),
             varianza_pct = varianza_pca)[1:10, ],
  paste0(ruta_carpeta_salida, "pca_varianza.csv"), row.names = FALSE
)
write.csv(
  data.frame(Accession = pop(gl), pca$scores),
  paste0(ruta_carpeta_salida, "pca_scores.csv"), row.names = FALSE
)
cat("Varianza explicada, primeros 3 ejes del PCA:", round(varianza_pca[1:3], 2), "\n")
t0 <- reportar_tiempo("PCA (glPca)", t0)

# --- Gráfica de PCA (PC1 vs PC2) ---
library(ggplot2)
n_pops_distintas <- length(unique(pop(gl)))
df_pca <- data.frame(PC1 = pca$scores[, 1], PC2 = pca$scores[, 2], Accession = pop(gl))

p_pca <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Accession)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(x = paste0("PC1 (", round(varianza_pca[1], 1), "%)"),
       y = paste0("PC2 (", round(varianza_pca[2], 1), "%)"),
       title = "PCA -- coloreado por Accession") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

# con muchas poblaciones (piensa en las 276 reales) la leyenda es
# ilegible -- se omite sola si hay demasiadas (ajustable)
if (n_pops_distintas > 30) {
  p_pca <- p_pca + theme(legend.position = "none")
  message("[INFO] ", n_pops_distintas, " poblaciones -- se omitió la leyenda del PCA por ilegible.")
}

ggsave(paste0(ruta_carpeta_salida, "pca_plot.png"), p_pca, width = 9, height = 7, dpi = 150)



# ============================================================
# 5a. DAPC replicando el enfoque del paper (n_pcs_grande fijo)
# ============================================================

buscar_K_con_salvaguarda <- function(gl, n_pca, k_max, etiqueta, carpeta_salida) {
  grupos <- find.clusters(gl, n.pca = n_pca, max.n.clust = k_max, choose.n.clust = FALSE)
  
  # guardar la curva completa de BIC vs K -- el equivalente a la Fig S3
  # del paper, para poder inspeccionarla visualmente, no solo confiar
  # en la elección automática
  bic_df <- data.frame(K = seq_along(grupos$Kstat), BIC = as.numeric(grupos$Kstat))
  write.csv(bic_df, paste0(carpeta_salida, "bic_vs_K_", etiqueta, ".csv"), row.names = FALSE)
  png(paste0(carpeta_salida, "bic_vs_K_", etiqueta, ".png"), width = 800, height = 500)
  plot(bic_df$K, bic_df$BIC, type = "b", xlab = "K", ylab = "BIC",
       main = paste("BIC vs K --", etiqueta))
  dev.off()
  
  k_elegido <- length(unique(grupos$grp))
  if (k_elegido == k_max) {
    message("[AVISO SERIO] Para '", etiqueta, "', el K elegido automáticamente (",
            k_elegido, ") choca exactamente contra el techo (max.n.clust=", k_max,
            "). Esto normalmente significa que el BIC no encontró un mínimo real ",
            "dentro del rango probado -- NO confiar en este K sin revisar el .png ",
            "de la curva a mano. Puede que necesites subir max.n.clust, o elegir K ",
            "manualmente viendo dónde el codo/mínimo de la curva ocurre de verdad.")
  }
  grupos
}

max_pcs_posible <- min(nInd(gl), nLoc(gl)) - 1
if (n_pcs_grande > max_pcs_posible) {
  message("[AVISO] n_pcs_grande (", n_pcs_grande, ") es mayor al máximo posible (",
          max_pcs_posible, "). Se usará el máximo posible en su lugar.")
  n_pcs_grande <- max_pcs_posible
}

grupos <- buscar_K_con_salvaguarda(gl, n_pcs_grande, k_max, "pcs_grande", ruta_carpeta_salida)
t0 <- reportar_tiempo(paste0("find.clusters (", n_pcs_grande, " PCs, K=1 a ", k_max, ")"), t0)

if (length(unique(grupos$grp)) < 2) {
  stop("K encontrado fue 1 (ningún grupo distinto) con n_pcs_grande=", n_pcs_grande,
       " -- el DAPC no puede correr con un solo grupo. Revisa bic_vs_K_pcs_grande.png; ",
       "puede que este número de PCs no esté capturando estructura real, o que de ",
       "verdad no haya estructura distinguible en este punto del pipeline.")
}

dapc_grande <- dapc(gl, pop = grupos$grp, n.pca = n_pcs_grande, n.da = length(unique(grupos$grp)) - 1)

saveRDS(dapc_grande, paste0(ruta_carpeta_salida, "dapc_grande_", n_pcs_grande, "pcs.rds"))
cat("\nDAPC con", n_pcs_grande, "PCs -- K encontrado:", length(unique(grupos$grp)), "\n")
cat("Proporción de reasignación correcta:", round(summary(dapc_grande)$assign.prop, 4), "\n")
t0 <- reportar_tiempo(paste0("dapc (", n_pcs_grande, " PCs)"), t0)

k_grande <- length(unique(grupos$grp))
paleta_grande <- rainbow(k_grande)
graficar_en_png(paste0(ruta_carpeta_salida, "dapc_scatter_pcs_grande.png"), 900, 700, function() {
  scatter(dapc_grande, col = paleta_grande, bg = "white", cstar = 0,
          legend = (k_grande <= 30), posi.leg = "topright",
          scree.pca = TRUE, posi.pca = "bottomleft",
          main = paste0("DAPC (", n_pcs_grande, " PCs, K=", k_grande, ")"))
})

# compoplot: la gráfica de barras tipo ADMIXTURE (probabilidad de pertenencia de cada individuo a cada grupo inferido)
graficar_en_png(paste0(ruta_carpeta_salida, "compoplot_pcs_grande.png"), 1200, 600, function() {
  compoplot(dapc_grande, col = paleta_grande, legend = (k_grande <= 30),
            main = paste0("Compoplot tipo ADMIXTURE (", n_pcs_grande, " PCs, K=", k_grande, ")"))
})
# ============================================================
# 5b. DAPC con el número de PCs sugerido por xvalDapc
# ============================================================

cat("\n--- Corriendo xvalDapc (puede tardar) ---\n")
mat <- as.matrix(gl)
mat[is.na(mat)] <- 0  # xvalDapc no acepta NA -- imputación simple; revisar cuántos NA hay antes

xval <- xvalDapc(mat, pop(gl), n.pca.max = xval_n_pca_max, n.rep = xval_n_rep,
                 xval.plot = FALSE)

n_pcs_xval <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
cat("Número de PCs sugerido por xvalDapc:", n_pcs_xval, "\n")
t0 <- reportar_tiempo(paste0("xvalDapc (hasta ", xval_n_pca_max, " PCs, ", xval_n_rep, " repeticiones)"), t0)

grupos_xval <- buscar_K_con_salvaguarda(gl, n_pcs_xval, k_max, "pcs_xval", ruta_carpeta_salida)

if (length(unique(grupos_xval$grp)) < 2) {
  stop("K encontrado fue 1 (ningún grupo distinto) con n_pcs_xval=", n_pcs_xval,
       " -- el DAPC no puede correr con un solo grupo. Revisa bic_vs_K_pcs_xval.png.")
}

dapc_xval <- dapc(gl, pop = grupos_xval$grp, n.pca = n_pcs_xval,
                  n.da = length(unique(grupos_xval$grp)) - 1)

saveRDS(dapc_xval, paste0(ruta_carpeta_salida, "dapc_xval_", n_pcs_xval, "pcs.rds"))
cat("DAPC con", n_pcs_xval, "PCs (xvalDapc) -- K encontrado:", length(unique(grupos_xval$grp)), "\n")
cat("Proporción de reasignación correcta:", round(summary(dapc_xval)$assign.prop, 4), "\n")
t0 <- reportar_tiempo(paste0("find.clusters + dapc (", n_pcs_xval, " PCs, xvalDapc)"), t0)


k_xval <- length(unique(grupos_xval$grp))
paleta_xval <- rainbow(k_xval)
graficar_en_png(paste0(ruta_carpeta_salida, "dapc_scatter_pcs_xval.png"), 900, 700, function() {
  scatter(dapc_xval, col = paleta_xval, bg = "white", cstar = 0,
          legend = (k_xval <= 30), posi.leg = "topright",
          scree.pca = TRUE, posi.pca = "bottomleft",
          main = paste0("DAPC (", n_pcs_xval, " PCs xvalDapc, K=", k_xval, ")"))
})

graficar_en_png(paste0(ruta_carpeta_salida, "compoplot_pcs_xval.png"), 1200, 600, function() {
  compoplot(dapc_xval, col = paleta_xval, legend = (k_xval <= 30),
            main = paste0("Compoplot tipo ADMIXTURE (", n_pcs_xval, " PCs, K=", k_xval, ")"))
})
cat(sprintf("\n[tiempo] TOTAL del script: %.2f min\n", as.numeric(difftime(Sys.time(), t_inicio_total, units = "mins"))))
# ============================================================
# 6. comparación rápida para decidir qué número de PCs usar en los siguientes checkpoints (paso 3, 4, 5)
# ============================================================

cat("\n=== COMPARACIÓN ===\n")
cat("PCs estilo paper (", n_pcs_grande, "): K =", length(unique(grupos$grp)),
    ", reasignación =", round(summary(dapc_grande)$assign.prop, 4), "\n")
cat("PCs por xvalDapc (", n_pcs_xval, "): K =", length(unique(grupos_xval$grp)),
    ", reasignación =", round(summary(dapc_xval)$assign.prop, 4), "\n")
cat("\nSi la reasignación con PCs grande está muy cerca de 1.0 (100%) y la de xvalDapc es notablemente más baja, 
    es señal de sobreajuste con PCs grande (como en la prueba con datos sintéticos). 
    Decide aquí cuál número de PCs usar de forma FIJA en los checkpoints 2 a 5, para no volver a correr xvalDapc cada vez.\n")


# sí usaré xval