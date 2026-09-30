# ============================================================
# poda_ld.R
# script de poda por desequilibrio de ligamiento (LD pruning) con plink
#
# pasos de plink2:
#   1) --indep-pairwise: recorre los SNPs en ventanas deslizantes y decide cuáles quedan (.prune.in) 
#   y cuáles se podan (.prune.out) por estar correlacionados (r^2) con otro SNP de la misma ventana.
#   2) --extract con el .prune.in: escribe el bfile final, quedándose solo con los SNPs que sobrevivieron la poda.
#
# dirección del umbral
#   r2_umbral bajo  -> poda más agresiva (exige que los SNPs estén casi nada correlacionados entre sí para dejarlos a ambos)
#   r2_umbral alto  -> poda más laxa (tolera SNPs bastante correlacionadosy deja a los dos)
# ============================================================


# ----
# rutas
# ------

# bfile de entrada (tomo el archivo del paso anterior )
archivo_bfile <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/filtrado_inicial"

carpeta_salida <- "/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/teocintle/24_sep/"

# tamaño de la ventana en número de SNPs (no en kb) y el paso con el que se desliza
# el paper de Diana Rivera, se usaron 50 10 0.2 (ventana 50 SNPs, avanza 10 SNPs, umbral 0.2)
tamaño_ventana_snps <- 50
paso_ventana_snps <- 10

# techo de tolerancia de r^2 
# el dataset de Dryad ya viene podado a r2=0.2 por el paper original, así que en teoría, un umbral mas laxo (0.8, 0.9)casi no debería quitar nada
# realmente no sería un filtro real sino confirmación del que ya hicieron, para ver un efecto tendría que bajar el umbral
# pero tampoco es el objetivo
r2_umbral <- 0.8


# -------------------
# funciones auxiliares
# -------------------

dir.create(carpeta_salida, recursive = TRUE, showWarnings = FALSE)

correr_plink2 <- function(args, descripcion) {
  cat("\n>>>", descripcion, "\n")
  resultado <- system2("plink2", args)
  if (resultado != 0) {
    stop("plink2 falló en: ", descripcion, " -- revisar el log de plink2 (.log) para el error real.")
  }
}

contar_lineas <- function(ruta) {
  length(readLines(ruta))
}


# -------------------------------------
#  --indep-pairwise (decide qué se poda)
# -------------------------------------

prefijo_indep <- paste0(carpeta_salida, "indep_pairwise")

correr_plink2(c("--bfile", archivo_bfile,
                "--make-founders",
                "--indep-pairwise", as.character(tamaño_ventana_snps),
                as.character(paso_ventana_snps), as.character(r2_umbral),
                "--out", prefijo_indep),
              paste0("Calculando LD pruning (ventana=", tamaño_ventana_snps,
                     " SNPs, paso=", paso_ventana_snps, ", r2=", r2_umbral, ")"))

# .prune.in  = SNPs que SI se quedan
# .prune.out = SNPs que se podan
n_conservados <- contar_lineas(paste0(prefijo_indep, ".prune.in"))
n_podados <- contar_lineas(paste0(prefijo_indep, ".prune.out"))

cat("\nSNPs a conservar (.prune.in):", n_conservados, "\n")
cat("SNPs a podar (.prune.out):", n_podados, "\n")


# ------------------------------
# aplicar la poda y bfile final
# ----------------------------
archivo_podado <- paste0(carpeta_salida, "podado_ld")

correr_plink2(c("--bfile", archivo_bfile,
                "--extract", paste0(prefijo_indep, ".prune.in"),
                "--make-bed", "--out", archivo_podado),
              "Extrayendo solo los SNPs que sobrevivieron la poda")


# -------
# resumen
# -------

n_snps_original <- contar_lineas(paste0(archivo_bfile, ".bim"))
n_snps_final <- contar_lineas(paste0(archivo_podado, ".bim"))

cat("\n--- resumen ---\n")
cat("SNPs antes de la poda:", n_snps_original, "\n")
cat("SNPs después de la poda:", n_snps_final, "\n")
cat("SNPs eliminados:", n_snps_original - n_snps_final,
    sprintf(" (%.1f%%)\n", 100 * (n_snps_original - n_snps_final) / n_snps_original))
cat("\nArchivo final:", archivo_podado, ".bed/.bim/.fam\n")

