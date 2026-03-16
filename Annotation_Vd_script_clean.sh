#!/bin/bash
# =============================================================================
# Anotación funcional del genoma ensamblado de Verticillium dahliae (Colon)
# =============================================================================

# --- Parámetros generales -----------------------------------------------------

THREADS=32
ISOLATE="Colon"
SPECIES="Verticillium_dahliae"

# Ensamblado final (producto del pipeline de ensamblado)
GENOME_SCAFFOLD="AliTV/Scaffold/ragtag.scaffold.fasta"

# Directorios de salida
FUNANNOTATE_DIR="06_Funannotate"
EFFECTOR_DIR="06_EffectorP"
TELOMERE_DIR="07_Telomere_tidk"
OCCULTER_DIR="03_OcculterCut"

# =============================================================================
# MÓDULO 1: PRE-PROCESAMIENTO DEL ENSAMBLADO
# =============================================================================

mkdir -p ${FUNANNOTATE_DIR}

# Ordenar contigs por tamaño y normalizar nombres de secuencias
funannotate sort \
    -i ${GENOME_SCAFFOLD} \
    -o ${FUNANNOTATE_DIR}/colon.scaffold.fasta

# Enmascaramiento suave de regiones repetitivas
funannotate mask \
    --cpus ${THREADS} \
    -i ${FUNANNOTATE_DIR}/colon.scaffold.fasta \
    -o ${FUNANNOTATE_DIR}/colon.reptmask.fasta

# =============================================================================
# MÓDULO 2: PREDICCIÓN DE GENES (Funannotate)
# =============================================================================

funannotate predict \
    --cpus ${THREADS} \
    -s ${SPECIES} \
    --isolate ${ISOLATE} \
    --name ${ISOLATE} \
    --busco_seed_species verticillium_longisporum1 \
    -i ${FUNANNOTATE_DIR}/colon.reptmask.fasta \
    -o ${FUNANNOTATE_DIR}/Colon

# =============================================================================
# MÓDULO 3: ANOTACIÓN FUNCIONAL
# =============================================================================

# --- 3.1 InterProScan (análisis de dominios proteicos) ----------------------
interproscan.sh \
    --cpu ${THREADS} \
    -f XML -goterms -pa \
    -i ${FUNANNOTATE_DIR}/Colon/predict_results/${SPECIES}_${ISOLATE}.proteins.fa \
    -d ${FUNANNOTATE_DIR}/Colon/interpro

# Mover resultado al directorio que usa Funannotate
mv ${FUNANNOTATE_DIR}/Colon/interpro/${SPECIES}_${ISOLATE}.proteins.fa.xml \
   ${FUNANNOTATE_DIR}/Colon/annotate_misc/iprscan.xml

# --- 3.2 AntiSMASH (clusters de metabolitos secundarios) --------------------
funannotate remote \
    -i ${FUNANNOTATE_DIR}/Colon \
    -o ${FUNANNOTATE_DIR}/Colon \
    -m antismash

# --- 3.3 Phobius (predicción de péptido señal y dominios transmembrana) -----
phobius -short \
    ${FUNANNOTATE_DIR}/Colon/annotate_misc/genome.proteins.fasta \
    > ${FUNANNOTATE_DIR}/Colon/annotate_misc/phobius.results.txt

# --- 3.4 Anotación final integrando todos los resultados -------------------
funannotate annotate \
    -i ${FUNANNOTATE_DIR}/Colon \
    --cpus ${THREADS}

# =============================================================================
# MÓDULO 4: PREDICCIÓN DE EFECTORES (EffectorP)
# =============================================================================

mkdir -p ${EFFECTOR_DIR}

ANNOTATED_GENOME="${FUNANNOTATE_DIR}/Colon/annotate_results/${SPECIES}_${ISOLATE}.scaffolds.fa"
ANNOTATED_GFF="Verticillium_dahliae_VArg1_sinMT.gff3"
MASKED_GENOME="01_Masked_Genomes/Varg1_sinMT_RM.fasta"

# Extraer proteínas del GFF anotado
gffread ${ANNOTATED_GFF} \
    -g ${MASKED_GENOME} \
    -y ${EFFECTOR_DIR}/Varg1_proteins.fasta

# Filtrar proteínas con señal de secreción (anotadas como SECRETED)
grep "SECRETED" ${ANNOTATED_GFF} \
    | awk '{ for(i=1; i<=NF; i++) if($i ~ /ID=/) { split($i, a, ";"); print a[1] } }' \
    | sed 's/ID=//' \
    > ${EFFECTOR_DIR}/Varg1_secreted_ids.txt

seqkit grep \
    -f ${EFFECTOR_DIR}/Varg1_secreted_ids.txt \
    ${EFFECTOR_DIR}/Varg1_proteins.fasta \
    > ${EFFECTOR_DIR}/Varg1_secreted_proteins.fasta

# Predecir efectores
python tools/EffectorP_3.0.0-beta/EffectorP.py -f \
    -o ${EFFECTOR_DIR}/Varg1_EffectorP_output.tsv \
    -i ${EFFECTOR_DIR}/Varg1_secreted_proteins.fasta

# Generar GFF con las anotaciones de efectores
awk -F "\t" '{print $1}' ${EFFECTOR_DIR}/Varg1_EffectorP_output.tsv \
    > ${EFFECTOR_DIR}/Varg1_EffectorP_ids.txt

gffread --ids ${EFFECTOR_DIR}/Varg1_EffectorP_ids.txt \
    ${ANNOTATED_GFF} \
    > ${EFFECTOR_DIR}/Varg1_EffectorP.gff

# =============================================================================
# MÓDULO 5: PREDICCIÓN DE REGIONES RICAS EN AT (OcculterCut)
# =============================================================================

mkdir -p ${OCCULTER_DIR}
cd ${OCCULTER_DIR}

OcculterCut \
    -f ../01_Masked_Genomes/Varg1_sinMT_RM.fasta \
    -a ../Verticillium_dahliae_VArg1_sinMT.gff3

cd ..

# =============================================================================
# MÓDULO 6: IDENTIFICACIÓN DE TELÓMEROS (tidk)
# =============================================================================

mkdir -p ${TELOMERE_DIR}

TELOMERE_SEQ="TTAGGG"

# Búsqueda del motivo telomérico
tidk search \
    --string ${TELOMERE_SEQ} \
    --extension tsv \
    --dir ${TELOMERE_DIR}/ \
    --output Varg1 \
    01_RepeatMasked/Varg1_sinMT_RM.fasta

tidk search \
    --string ${TELOMERE_SEQ} \
    --extension bedgraph \
    --dir ${TELOMERE_DIR}/ \
    --output Varg1 \
    01_RepeatMasked/Varg1_sinMT_RM.fasta

# Visualización de la distribución telomérica
tidk plot \
    --tsv ${TELOMERE_DIR}/Varg1_telomeric_repeat_windows.tsv \
    --output ${TELOMERE_DIR}/Varg1_telomeric
