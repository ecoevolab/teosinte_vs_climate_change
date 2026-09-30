# ============================================================
# subset bebé de pruebas bases contiguas
# 25 sept 2026
# 1000 SNPs contiguos de un solo cromosoma, 3 pobs
# es para probar las visualizaciones de LD
# ============================================================

# ----
# rutas
# ------

ruta_bfile <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all"   
ruta_metadata <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"  
ruta_salida_subset <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_baby_subset_contiguo_100"

n_poblaciones <- 3   
n_snps <- 100       
semilla <- 21    #21 para 1000, 100 snps.     

# --------------------------------------------
# elegir 3 poblaciones random con mínimo 8 ind
# ----------------------------------------------

meta <- read.csv(ruta_metadata, stringsAsFactors = FALSE)
conteo_por_pop <- table(meta$Accession)

# solo considerar poblaciones con al menos 8 individuos, para que la prueba se parezca un poco más a una corrida real
candidatas <- names(conteo_por_pop[conteo_por_pop >= 8])

set.seed(semilla)
pops_elegidas <- sort(sample(candidatas, n_poblaciones))
cat("Poblaciones elegidas para la prueba:", paste(pops_elegidas, collapse = ", "), "\n")
cat("Individuos por población:\n")
print(conteo_por_pop[pops_elegidas])

individuos_elegidos <- meta[meta$Accession %in% pops_elegidas, ]


# --------------------------------------------
# escribir el archivo --keep para plink2
# ----------------------------------------------

fam <- read.table(paste0(ruta_bfile, ".fam"), stringsAsFactors = FALSE,
                  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))

fam_filtrado <- fam[fam$IID %in% individuos_elegidos$Sample_name, c("FID", "IID")]
cat("Individuos encontrados en el .fam para esas 3 poblaciones:", nrow(fam_filtrado), "\n")

write.table(fam_filtrado, "keep_individuos_prueba.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# --------------------------------------------
# elegir tramo contiguo de 1000 snps, de 1 solo chr
# ----------------------------------------------

bim <- read.table(paste0(ruta_bfile, ".bim"), stringsAsFactors = FALSE,
                  col.names = c("CHR", "SNP_ID", "CM", "BP", "A1", "A2"))
bim$CHR <- as.character(bim$CHR)

snps_por_cromosoma <- table(bim$CHR)
cromosomas_candidatos <- names(snps_por_cromosoma[snps_por_cromosoma >= n_snps])
if (length(cromosomas_candidatos) == 0) {
  stop("Ningún cromosoma tiene al menos ", n_snps, " SNPs -- bajar n_snps o revisar bim$CHR.")
}

set.seed(semilla)
cromosoma_elegido <- sample(cromosomas_candidatos, 1)

bim_cromosoma <- bim[bim$CHR == cromosoma_elegido, ]
bim_cromosoma <- bim_cromosoma[order(bim_cromosoma$BP), ]

# el tramo que no inicie en la orilla
max_inicio <- nrow(bim_cromosoma) - n_snps + 1
indice_inicio <- sample(100:(max_inicio-100), 1)
bim_tramo <- bim_cromosoma[indice_inicio:(indice_inicio + n_snps - 1), ]

snps_elegidos <- bim_tramo$SNP_ID

cat("Cromosoma elegido:", cromosoma_elegido, "(", nrow(bim_cromosoma), "SNPs totales en ese cromosoma)\n")
cat("Tramo contiguo:", n_snps, "SNPs, posiciones", min(bim_tramo$BP), "a", max(bim_tramo$BP),
    "bp (longitud física:", round((max(bim_tramo$BP) - min(bim_tramo$BP)) / 1000, 1), "kb)\n")

writeLines(snps_elegidos, "extract_snps_prueba.txt")

# --------------------------------------------
# correr plink2 para generar el subset
# ----------------------------------------------

resultado <- system2("plink2", c(
  "--bfile", ruta_bfile,
  "--keep", "keep_individuos_prueba.txt",
  "--extract", "extract_snps_prueba.txt",
  "--make-bed",
  "--out", ruta_salida_subset
))

if (resultado != 0) {
  stop("plink2 no corrió bien armando el subset -- revisar mensaje de arriba.")
}

cat("\nSubset de prueba listo en:", ruta_salida_subset, ".bed/.bim/.fam\n")

