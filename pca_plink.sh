#!/bin/bash
#SBATCH --job-name=pca_plink_0.2
#SBATCH --output=/mnt/data/sur/users/spacheco/results/sep_2026/teo/29_sep/pca_plink_0.2.out
#SBATCH --error=/mnt/data/sur/users/spacheco/results/sep_2026/teo/29_sep/pca_plink_0.2.err
#SBATCH --chdir=/mnt/data/sur/users/spacheco/results/sep_2026/teo/29_sep/
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# este script es porque creo que nunca había subido un ejemplo de job a mi github
# y para explicar cómo voy a estar obteniendo ahora los PCAs porque a través de scripts de R tardan como 1 día y varias horas y con plink directo en terminal, 3 minutos por mucho.

module load plink2

umbral_r2="0.2"
ruta_bfile="/mnt/data/sur/users/spacheco/results/sep_2026/teo/25_sep/r2_${umbral_r2}/podado_ld"
ruta_salida="/mnt/data/sur/users/spacheco/results/sep_2026/teo/29_sep/pca_plink_r2_${umbral_r2}"

# plink2 no  crea la carpeta de salida sola, si no existe, truena 
mkdir -p "${ruta_salida}"

# --make-founders en vez de --bad-freqs: --bad-freqs sola se salta el calculo real de frecuencias alelicas, --make-founders sí las calcula,
# tratando a todos los individuos como fundadores 

plink2 --bfile "${ruta_bfile}" \
       --pca 10 \
       --make-founders \
       --out "${ruta_salida}/pca_r2_${umbral_r2}"

# matriz de relacionalidad cuadrada -- su traza (suma de la diagonal) es igual a la suma de todos los eigenvalores posibles 
# es la forma barata de tener la varianza total real sin pedir muchísimos PCs

plink2 --bfile "${ruta_bfile}" \
       --make-rel square \
       --make-founders \
       --out "${ruta_salida}/rel_r2_${umbral_r2}"

echo "Listo, está en: ${ruta_salida}"
