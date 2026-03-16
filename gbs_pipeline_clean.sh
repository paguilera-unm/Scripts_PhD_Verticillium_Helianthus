#!/bin/bash
# =============================================================================
# Pipeline de análisis genómico poblacional por GBS (Genotyping-by-Sequencing)
# Verticillium dahliae - aislados locales y comparación con secuencias públicas
# =============================================================================

# --- Parámetros generales ----------------------------------------------------

THREADS=8
REF_GENOME="VdLs17_ref_genome.fna"       # Genoma de referencia VdLs17
POPMAP="popmap.txt"                        # Tabla de poblaciones (formato Stacks)

# Directorios de trabajo
RAW="00_Raw_Data"
MERGED="02_Merged"
TRIMMED="03_Trimmed"
MAPPED="04_Mapped"
VCF_DIR="07_VCF_Fixeds"
MILGROOM_DIR="milgroom_2014"
GENOME_DIR="10_Genome_comparison/ncbi_dataset/data"

# Índices BWA
IDX_LS17="Genome_Index/Refgenomebwaidx"
IDX_LS17_NEW="15_Mapeo_Ls17_new/Index/RefgenomeLs17_new"
IDX_JR2="${GENOME_DIR}/GCA_000400815.2JR2/Refgenomebwaidx"
IDX_85S="${GENOME_DIR}/GCA_004798885.185S/Refgenome85S"

# Función auxiliar: comprimir e indexar un VCF
compress_and_index() {
    local vcf="$1"
    bgzip "$vcf"
    bcftools index -f "${vcf}.gz"
}

# =============================================================================
# MÓDULO 1: PRE-PROCESAMIENTO DE READS CRUDOS
# =============================================================================

# 1.1 - Descomprimir y fusionar reads de múltiples lanes por muestra
gzip -d ${RAW}/*

declare -A SAMPLE_IDS=(
    ["andant"]="sar001-2020"
    ["usa"]="sar002-2020"
    ["colon"]="sar003-2020"
)

for sample in "${SAMPLES[@]}"; do
    raw_id="${SAMPLE_IDS[$sample]}"
    cat ${RAW}/${raw_id}_L00?_R1_001.fastq > ${MERGED}/${sample}-merg_L001_R1_001.fastq
    cat ${RAW}/${raw_id}_L00?_R2_001.fastq > ${MERGED}/${sample}-merg_L001_R2_001.fastq
done

# 1.2 - Control de calidad y demultiplexado con process_radtags
# Enzimas de restricción: PstI (corte raro) y Sau3AI/MboI (corte frecuente)
process_radtags \
    -P \
    -p ./${MERGED} \
    -o ./${TRIMMED} \
    -c -q \
    --renz_1 pstI \
    --renz_2 sau3AI \
    --disable_rad_check


# =============================================================================
# MÓDULO 7: MAPEO Y GENOTIPADO SOBRE EL NUEVO ENSAMBLE VdLs17
# =============================================================================

# 7.1 - Indexar el nuevo genoma de referencia
bwa index -p ${IDX_LS17_NEW} -a bwtsw ${REF_LS17_NEW}

# 7.2 - Mapear las muestras locales
for sample in "${SAMPLES[@]}"; do
    bwa mem -t ${THREADS} ${IDX_LS17_NEW} \
        ${TRIMMED}/${sample}-merg_L001_R1_001.1.fq.gz \
        ${TRIMMED}/${sample}-merg_L001_R2_001.2.fq.gz \
    | samtools view -b \
    | samtools sort --threads 4 > 15_Mapeo_Ls17_new/${sample}.bam
done

# 7.3 - Construir loci y analizar poblaciones
gstacks -I 15_Mapeo_Ls17_new/ -M ${POPMAP} -O 15_Mapeo_Ls17_new/05_Stacks/ -t ${THREADS}

populations \
    -P 15_Mapeo_Ls17_new/05_Stacks/ \
    -M ${POPMAP} \
    -r 0.65 --vcf --genepop --fstats --smooth --hwe \
    -t ${THREADS} \
    -O 15_Mapeo_Ls17_new/06_Populations/

# 7.4 - Corregir alelo de referencia y retener solo FORMAT/GT
bcftools norm -f ${REF_LS17_NEW} -c s \
    15_Mapeo_Ls17_new/06_Populations/populations.snps.vcf \
    > 15_Mapeo_Ls17_new/populations.snps.fix.vcf

bcftools annotate -x INFO,^FORMAT/GT --force \
    15_Mapeo_Ls17_new/populations.snps.fix.vcf \
    -Ov -o 15_Mapeo_Ls17_new/populations.snps.fix.GT.vcf

compress_and_index 15_Mapeo_Ls17_new/populations.snps.fix.GT.vcf

# =============================================================================
# MÓDULO 8: INTEGRACIÓN PROGRESIVA DE TODOS LOS CONJUNTOS DE DATOS
# =============================================================================
# Función: obtener SNPs compartidos entre dos VCFs y fusionarlos
isec_and_merge() {
    local vcf_a="$1"
    local vcf_b="$2"
    local isec_dir="$3"
    local out_vcf="$4"

    mkdir -p "${isec_dir}"
    bcftools isec "${vcf_a}" "${vcf_b}" -p "${isec_dir}"

    compress_and_index "${isec_dir}/0002.vcf"
    compress_and_index "${isec_dir}/0003.vcf"

    bcftools merge -m none -0 \
        "${isec_dir}/0002.vcf.gz" \
        "${vcf_b}" \
        -Ov -o "${out_vcf}"

    compress_and_index "${out_vcf}"
}

LAURA_VCF="15_Mapeo_Ls17_new/Fix_SNPs_lineages_BautistaJalon2020.vcf.gz"
BASE="15_Mapeo_Ls17_new"

# Paso 1: locales + Bautista-Jalón (2020)
isec_and_merge \
    ${BASE}/populations.snps.fix.GT.vcf.gz \
    ${LAURA_VCF} \
    ${BASE}/isec_Locales_vs_Laura \
    ${BASE}/merge_locales_laura.vcf

# Paso 2: + aislado 85S
isec_and_merge \
    ${BASE}/snp_85S_Ls17_fix_norm.vcf.gz \
    ${BASE}/merge_locales_laura.vcf.gz \
    ${BASE}/isec_LocalesLaura_vs_85S \
    ${BASE}/merge_localeslaura_85S.vcf

# Paso 3: + aislado Vd39
isec_and_merge \
    ${BASE}/snp_Vd39_Ls17_fix_norm.vcf.gz \
    ${BASE}/merge_localeslaura_85S.vcf.gz \
    ${BASE}/isec_LocalesLaura_vs_Vd39 \
    ${BASE}/merge_localeslaura85S_Vd39.vcf

# Paso 4: + aislado Vd0991
isec_and_merge \
    ${BASE}/snp_Vd0991_Ls17_fix_norm.vcf.gz \
    ${BASE}/merge_localeslaura85S_Vd39.vcf.gz \
    ${BASE}/isec_LocalesLaura_vs_Vd0991 \
    ${BASE}/merge_localeslaura85SVd39_Vd0991.vcf

# Paso 5: + aislado JR2
isec_and_merge \
    ${BASE}/snp_JR2_Ls17_fix_norm.vcf.gz \
    ${BASE}/merge_localeslaura85SVd39_Vd0991.vcf.gz \
    ${BASE}/isec_LocalesLaura_vs_JR2 \
    ${BASE}/merge_localeslaura85SVd39_Vd0991_JR2.vcf

# Paso 6: + secuencias de origen chino (SRR12769281, SRR12769283)
# (Variantes obtenidas con FreeBayes en Galaxy, filtradas a SNPs)
for srr in SRR12769281 SRR12769283; do
    bcftools view -v snps -Ov \
        -o 15_Mapeo_Ls17_new/${srr}.vcf \
        15_Mapeo_Ls17_new/${srr}_annotated.bcf
    compress_and_index 15_Mapeo_Ls17_new/${srr}.vcf
done

bcftools merge -m none -0 \
    15_Mapeo_Ls17_new/SRR12769281.vcf.gz \
    15_Mapeo_Ls17_new/SRR12769283.vcf.gz \
    -Ov -o 15_Mapeo_Ls17_new/merge_chinos.vcf
compress_and_index 15_Mapeo_Ls17_new/merge_chinos.vcf

isec_and_merge \
    15_Mapeo_Ls17_new/merge_chinos.vcf.gz \
    ${BASE}/merge_localeslaura85SVd39_Vd0991_JR2.vcf.gz \
    15_Mapeo_Ls17_new/isec_Chinos \
    15_Mapeo_Ls17_new/merge_final.vcf

# =============================================================================
# MÓDULO 9: SEGUNDA RONDA DE GBS (9 aislados locales adicionales)
# =============================================================================

mkdir -p 16_GBS2

# 9.1 - Verificar concordancia con el genoma de referencia
bcftools +fixref 16_GBS2/calls.bcf \
    -- -f ${REF_LS17_NEW}

# 9.2 - Integrar los nuevos aislados al conjunto de datos final
isec_and_merge \
    14_S023_S011/merge_Locales_Laura_85Svd39_Chinos_enmbl.vcf.gz \
    16_GBS2/calls.bcf \
    16_GBS2/isec \
    16_GBS2/merge_todos_los_aislados.vcf

# Fin del pipeline
