#!/bin/bash
# =============================================================================
# Ensamblado de novo de genomas de Verticillium dahliae
# Tecnología: PacBio HiFi Sequel II (CCS reads)
# =============================================================================

# --- Parámetros generales -----------------------------------------------------

THREADS=32
GENOME_SIZE="35m"
SAMPLES=("andant" "colon")

# Genomas de referencia
REF_LS17="GCA_000952015.1_ASM95201v1_genomic.fna"
REF_85S="GCA_004798885.1_ASM479888v1_genomic_85S.fna"

# Directorios de trabajo
RAW_DIR="00_Raw_data"
TRIM_DIR="02_Total_reads/02_Trimmed"
MERFIN="tools/merfin/build/bin/merfin"

# =============================================================================
# MÓDULO 1: CONTROL DE CALIDAD DE READS CRUDOS
# =============================================================================

mkdir -p 01_Quality_check

fastqc ${RAW_DIR}/*.fastq.gz -o 01_Quality_check/

# =============================================================================
# MÓDULO 2: TRIMMING DE READS
# =============================================================================
# Se elimina el primer nucleótido de cada read (artefacto de la librería GBS)

mkdir -p ${TRIM_DIR}

for sample in "${SAMPLES[@]}"; do
    trimmomatic SE -phred33 -threads ${THREADS} \
        ${RAW_DIR}/${sample}_hifi_reads.fastq.gz \
        ${TRIM_DIR}/${sample}_HIFI_trimmed.fq.gz \
        HEADCROP:1
done

# =============================================================================
# MÓDULO 3: ENSAMBLADO DE NOVO (cuatro ensambladores comparados)
# =============================================================================

for sample in "${SAMPLES[@]}"; do

    # --- 3.1 Canu ---
    mkdir -p 02_Total_reads/03_Canu_Assembly/${sample^}
    canu \
        -p ${sample} \
        -d 02_Total_reads/03_Canu_Assembly/${sample^} \
        genomeSize=${GENOME_SIZE} useGrid=false \
        -pacbio-hifi ${TRIM_DIR}/${sample}_HIFI_trimmed.fq.gz

    # --- 3.2 Hifiasm ---
    mkdir -p 02_Total_reads/04_Hifiasm_Assembly/${sample^}
    hifiasm \
        -o 02_Total_reads/04_Hifiasm_Assembly/${sample^}/${sample}.asm \
        -t ${THREADS} \
        ${TRIM_DIR}/${sample}_HIFI_trimmed.fq.gz
    # Convertir GFA a FASTA
    awk '/^S/{print ">"$2; print $3}' \
        02_Total_reads/04_Hifiasm_Assembly/${sample^}/${sample}.asm.bp.p_ctg.gfa \
        > 02_Total_reads/04_Hifiasm_Assembly/${sample^}/${sample}.asm.bp.p_ctg.fa

    # --- 3.3 IPA ---
    mkdir -p 02_Total_reads/05_IPA_Assembly/${sample^}
    ipa local \
        --nthreads 20 --njobs 4 \
        --genome-size 35000000 \
        --run-dir 02_Total_reads/05_IPA_Assembly/${sample^} \
        -i ${TRIM_DIR}/${sample}_HIFI_trimmed.fq.gz

    # --- 3.4 Flye ---
    mkdir -p 02_Total_reads/06_Flye_Assembly/${sample^}
    flye \
        --threads ${THREADS} \
        --genome-size ${GENOME_SIZE} \
        --scaffold --iterations 10 \
        --pacbio-hifi ${TRIM_DIR}/${sample}_HIFI_trimmed.fq.gz \
        --out-dir 02_Total_reads/06_Flye_Assembly/${sample^}

done

# =============================================================================
# MÓDULO 4: CONTROL DE CALIDAD DEL ENSAMBLADO
# =============================================================================

# --- 4.1 QUAST: métricas de contiguidad y precisión --------------------------
mkdir -p 02_Total_reads/Quast_results

quast.py \
    02_Total_reads/03_Canu_Assembly/Colon/colon.contigs.fasta \
    02_Total_reads/04_Hifiasm_Assembly/Colon/colon.asm.bp.p_ctg.fa \
    02_Total_reads/05_IPA_Assembly/Colon/19-final/final.p_ctg.fasta \
    02_Total_reads/06_Flye_Assembly/Colon/assembly.fasta \
    ${REF_85S} \
    -r ${REF_LS17} \
    -o 02_Total_reads/Quast_results \
    --fungus

# --- 4.2 BUSCO: evaluación de completitud génica ----------------------------
busco -m geno -c ${THREADS} \
    -i AliTV/input/ \
    -l hypocreales_odb10 \
    -o BUSCO_results \
    --out_path BUSCO_results

# --- 4.3 KAT: análisis de espectro k-mer ------------------------------------
# Se evalúan los cuatro ensamblados del aislado colon
declare -A ASSEMBLIES_COLON=(
    ["Canu"]="02_Total_reads/03_Canu_Assembly/Colon/colon.contigs.fasta"
    ["Hifiasm"]="02_Total_reads/04_Hifiasm_Assembly/Colon/colon.asm.bp.p_ctg.fa"
    ["IPA"]="02_Total_reads/05_IPA_Assembly/Colon/19-final/final.p_ctg.fasta"
    ["Flye"]="02_Total_reads/06_Flye_Assembly/Colon/assembly.fasta"
)

for assembler in "${!ASSEMBLIES_COLON[@]}"; do
    outdir="KAT/Colon/${assembler}"
    mkdir -p ${outdir}
    kat comp -t ${THREADS} \
        -o ${outdir}/colon \
        ${TRIM_DIR}/colon_HIFI_trimmed.fq.gz \
        "${ASSEMBLIES_COLON[$assembler]}"
done

# --- 4.4 Inspector: evaluación de errores de ensamblado ---------------------
for assembler in Canu Hifiasm; do
    if [ "$assembler" == "Canu" ]; then
        asm="02_Total_reads/03_Canu_Assembly/Colon/colon.contigs.fasta"
        outdir="02_Total_reads/03_Canu_Assembly/Colon/inspector"
    else
        asm="02_Total_reads/04_Hifiasm_Assembly/Colon/colon.asm.bp.p_ctg.fa"
        outdir="02_Total_reads/04_Hifiasm_Assembly/Colon/inspector"
    fi
    inspector.py \
        -c ${asm} \
        -r ${TRIM_DIR}/colon_HIFI_trimmed.fq.gz \
        -o ${outdir} \
        --datatype hifi \
        --thread ${THREADS}
done

# --- 4.5 Merfin: evaluación de calidad basada en k-mers (QV*) ---------------
MERYL_DB="MAC/merfin/Colon"
PEAK=134   # Valor de kcov obtenido de GenomeScope

# Construir la base de datos de k-mers de los reads
meryl count k=27 ${TRIM_DIR}/colon_HIFI_trimmed.fq.gz output ${MERYL_DB}/reads.meryl
meryl histogram ${MERYL_DB}/reads.meryl > ${MERYL_DB}/reads.hist
meryl greater-than 1 ${MERYL_DB}/reads.meryl output ${MERYL_DB}/reads.gt1.meryl

# Evaluar QV* de cada ensamblado
declare -A ASSEMBLIES_MERFIN=(
    ["Canu"]="02_Total_reads/03_Canu_Assembly/Colon/colon.contigs.fasta"
    ["Hifiasm"]="02_Total_reads/04_Hifiasm_Assembly/Colon/colon.asm.bp.p_ctg.fa"
    ["IPA"]="MAC/input/colon_IPA.fa"
    ["Flye"]="MAC/input/colon_Flye.fa"
)

for assembler in "${!ASSEMBLIES_MERFIN[@]}"; do
    outdir="MAC/merfin/Colon/${assembler}"
    mkdir -p ${outdir}

    ${MERFIN} -hist \
        -sequence "${ASSEMBLIES_MERFIN[$assembler]}" \
        -readmers ${MERYL_DB}/reads.gt1.meryl \
        -peak ${PEAK} \
        -prob ${MERYL_DB}/lookup_table.txt \
        -output ${outdir}/merfin.hist \
        2> ${outdir}/merfin.hist.log

    ${MERFIN} -completeness \
        -sequence "${ASSEMBLIES_MERFIN[$assembler]}" \
        -readmers ${MERYL_DB}/reads.gt1.meryl \
        -peak ${PEAK} \
        -prob ${MERYL_DB}/lookup_table.txt \
        2> ${outdir}/merfin.completeness
done

# --- 4.6 MUMmer: comparación sinténica con genomas de referencia -------------

MUMMER_OUT="02_Total_reads/06_Flye_Assembly/Mummer"
mkdir -p ${MUMMER_OUT}

# Función: alinear, filtrar y graficar con MUMmer
run_mummer() {
    local ref="$1"
    local query="$2"
    local prefix="$3"

    nucmer --maxgap=500 --minmatch=200 --threads=${THREADS} \
        --prefix="${MUMMER_OUT}/${prefix}" "${ref}" "${query}"
    delta-filter -1 "${MUMMER_OUT}/${prefix}.delta" \
        > "${MUMMER_OUT}/${prefix}.filter"
    mummerplot --png --fat --layout \
        -p "${MUMMER_OUT}/${prefix}_plot" \
        -R "${ref}" \
        "${MUMMER_OUT}/${prefix}.filter"
    gnuplot "${MUMMER_OUT}/${prefix}_plot.gp"
}

run_mummer ${REF_LS17} \
    02_Total_reads/06_Flye_Assembly/Colon/assembly.fasta \
    "Ls17_colon"

run_mummer ${REF_85S} \
    02_Total_reads/06_Flye_Assembly/Colon/assembly.fasta \
    "85S_colon"

# =============================================================================
# MÓDULO 5: SCAFFOLDING CON RAGTAG
# =============================================================================
# El ensamblado seleccionado (Hifiasm) se scaffoldea usando IPA como referencia

mkdir -p AliTV/Scaffold

ragtag.py scaffold \
    -t ${THREADS} -w \
    -o AliTV/Scaffold \
    -e AliTV/exclude_colon_IPA.txt \
    --aligner minimap2 \
    AliTV/colon_HIFIasm_recortado.fa \
    AliTV/input/colon_IPA.fa

# --- Evaluación del scaffold final -------------------------------------------
busco -m geno -c ${THREADS} \
    -i AliTV/Scaffold/ragtag.scaffold.fasta \
    -l hypocreales_odb10 \
    -o BUSCO_scaffold \
    --out_path BUSCO_scaffold

kat comp -t ${THREADS} \
    -o KAT/colon_scaffold \
    ${TRIM_DIR}/colon_HIFI_trimmed.fq.gz \
    AliTV/Scaffold/ragtag.scaffold.fasta

${MERFIN} -hist \
    -sequence AliTV/Scaffold/ragtag.scaffold.fasta \
    -readmers ${MERYL_DB}/reads.gt1.meryl \
    -peak ${PEAK} \
    -prob ${MERYL_DB}/lookup_table.txt \
    -output AliTV/Scaffold/ragtag.scaffold.hist \
    2> AliTV/Scaffold/ragtag.scaffold.hist.log

# Detección de telómeros con Tapestry
weave \
    -a AliTV/Scaffold/ragtag.scaffold.fasta \
    -r ${TRIM_DIR}/colon_HIFI_trimmed.fq.gz \
    -t TTAGGG -c ${THREADS} -l 4000 -f \
    -o AliTV/Scaffold/Colon_Scaffold

# =============================================================================
# Ensamblado final seleccionado: AliTV/Scaffold/ragtag.scaffold.fasta
# =============================================================================
