#!/bin/bash
# =============================================================================
# Pipeline de análisis de RNA-seq dual (Helianthus annuus / Verticillium dahliae)
# Estrategia: mapeo secuencial; los reads no mapeados contra girasol se mapean
#             contra el genoma de V. dahliae para cuantificar la expresión del patógeno
# =============================================================================

# --- Parámetros generales -----------------------------------------------------

THREADS=64
WD="RNA-Seq_Verticillium"

# Genomas y anotaciones de referencia
REF_SUNFLOWER="HanXRQr2.0-SUNRISE-2.1.genome.fasta"
GTF_SUNFLOWER="HanXRQr2.0-SUNRISE-2.1.funcAnnot.gtf"
REF_VD="Verticillium_dahliae_Colon.scaffolds.fa"
GTF_VD="Verticillium_dahliae_Colon.gtf"

# Índices STAR
IDX_SUNFLOWER="STAR_index/HanXRQr21_STAR_index"
IDX_VD="STAR_index/Colon_STAR_index"

# Directorios
RAW="${WD}/00_Raw_data"
QC="${WD}/01_QC"
MAP_SUNF="${WD}/04_Mapping_Sunflower"
MAP_VD="${WD}/07_Mapping_Vd_unmap"
QUANT_SUNF="${WD}/06_Quant_Sunf"
QUANT_VD="${WD}/08_Quant_Vd"
RSEQC="${WD}/09_RSeQC_Sunf"
BED_SUNFLOWER="HanXRQr2.0-SUNRISE-2.1.funcAnnot.bed"

# =============================================================================
# MÓDULO 1: FUSIÓN DE READS POR MUESTRA (múltiples lanes y runs)
# =============================================================================

mkdir -p ${RAW}

for id in $(ls ${RAW_UPSTREAM}/*_R1_001.fastq.gz \
            | awk -F "_" '{print $1"_"$2}' | sort | uniq); do
    cat ${RAW_UPSTREAM}/${id}_L00?_R1_001.fastq.gz > ${RAW}/${id}_L001_R1_001.fastq.gz
done

# =============================================================================
# MÓDULO 2: CONTROL DE CALIDAD
# =============================================================================

mkdir -p ${QC}
fastqc -o ${QC} ${RAW}/*.fastq.gz

# =============================================================================
# MÓDULO 3: INDEXADO DE GENOMAS DE REFERENCIA
# =============================================================================

# Índice para Helianthus annuus (HanXRQr2.1)
mkdir -p ${IDX_SUNFLOWER}
STAR --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${IDX_SUNFLOWER} \
    --genomeFastaFiles ${REF_SUNFLOWER} \
    --sjdbGTFfile ${GTF_SUNFLOWER} \
    --sjdbOverhang 74

# Índice para V. dahliae (genoma ensamblado VArg1)
mkdir -p ${IDX_VD}
STAR --runThreadN ${THREADS} --genomeSAindexNbases 11 \
    --runMode genomeGenerate \
    --genomeDir ${IDX_VD} \
    --genomeFastaFiles ${REF_VD} \
    --sjdbGTFfile ${GTF_VD} \
    --sjdbOverhang 100

# =============================================================================
# MÓDULO 4: MAPEO CONTRA GIRASOL
# =============================================================================

mkdir -p ${MAP_SUNF}

for id in $(ls ${RAW}/*.fastq.gz \
            | awk -F "_" '{print $1"_"$2}' | sort | uniq); do

    STAR --genomeDir ${IDX_SUNFLOWER} \
        --runThreadN ${THREADS} \
        --readFilesIn ${RAW}/${id}_L001_R1_001.fastq.gz \
                      ${RAW}/${id}_L001_R2_001.fastq.gz \
        --outFileNamePrefix ${MAP_SUNF}/${id} \
        --readFilesCommand zcat \
        --quantMode GeneCounts \
        --outSAMtype BAM Unsorted SortedByCoordinate \
        --outSAMattributes Standard \
        --outReadsUnmapped Fastx
done

# Reporte de calidad del mapeo
mkdir -p ${QC}/Mapeo_Sunflower
multiqc -o ${QC}/Mapeo_Sunflower/ ${MAP_SUNF}/*Log.final.out

# =============================================================================
# MÓDULO 5: MAPEO DE READS NO MAPEADOS CONTRA V. DAHLIAE
# =============================================================================
# Los reads que no mapearon contra el transcriptoma del hospedador
# se reasignan al genoma del patógeno

mkdir -p ${MAP_VD}

for id in $(ls ${RAW}/*.fastq.gz \
            | awk -F "_" '{print $1"_"$2}' | sort | uniq); do

    STAR --genomeDir ${IDX_VD} \
        --runThreadN ${THREADS} \
        --readFilesIn ${MAP_SUNF}/${id}Unmapped.out.mate1 \
                      ${MAP_SUNF}/${id}Unmapped.out.mate2 \
        --outFileNamePrefix ${MAP_VD}/${id} \
        --quantMode GeneCounts \
        --outSAMtype BAM Unsorted SortedByCoordinate \
        --outSAMattributes Standard \
        --outReadsUnmapped Fastx
done

# Reporte de calidad del mapeo
mkdir -p ${QC}/Mapeo_Vd_unmap
multiqc -o ${QC}/Mapeo_Vd_unmap ${MAP_VD}/*Log.final.out

# Comprimir reads no mapeados para liberar espacio
for id in $(ls ${RAW}/*.fastq.gz \
            | awk -F "_" '{print $1"_"$2}' | sort | uniq); do
    gzip -9 ${MAP_SUNF}/${id}Unmapped.out.mate1
    gzip -9 ${MAP_SUNF}/${id}Unmapped.out.mate2
done

# =============================================================================
# MÓDULO 6: CUANTIFICACIÓN DE EXPRESIÓN (mmquant)
# =============================================================================

# Girasol
mkdir -p ${QUANT_SUNF}
nombres_sf=$(ls ${MAP_SUNF}/*sorted* | awk -F "/" '{print $2}' | awk -F "_" '{print $1}')

mmquant -t ${THREADS} \
    -a ${GTF_SUNFLOWER} \
    -r ${MAP_SUNF}/*Aligned.sortedByCoord.out.bam \
    -n ${nombres_sf} \
    -o ${QUANT_SUNF}/mmquant_count.tsv

# V. dahliae
mkdir -p ${QUANT_VD}
nombres_vd=$(ls ${MAP_VD}/*sorted* | awk -F "/" '{print $2}' | awk -F "_" '{print $1}')

mmquant -t ${THREADS} \
    -a ${GTF_VD} \
    -r ${MAP_VD}/*Aligned.sortedByCoord.out.bam \
    -n ${nombres_vd} \
    -o ${QUANT_VD}/mmquant_count.tsv

# =============================================================================
# MÓDULO 7: CONTROL DE CALIDAD POST-MAPEO (RSeQC)
# =============================================================================
# Se convierte la anotación de GTF a formato BED requerido por RSeQC

gtfToGenePred ${GTF_SUNFLOWER} HanXRQ.genepred
genePredToBed HanXRQ.genepred ${BED_SUNFLOWER}

samtools faidx ${REF_SUNFLOWER}
awk -F "\t" '{print $1"\t"$2}' ${REF_SUNFLOWER}.fai > HanXRQ.chromsize

mkdir -p ${RSEQC}/{stat,inner_distance,junction_annotation,junction_saturation,tin}

for bam in ${MAP_SUNF}/*Aligned.sortedByCoord.out.bam; do
    id=$(basename ${bam} Aligned.sortedByCoord.out.bam)

    bam_stat.py -i ${bam} \
        > ${RSEQC}/stat/${id}.stat

    inner_distance.py -i ${bam} \
        -o ${RSEQC}/inner_distance/${id} \
        -r ${BED_SUNFLOWER}

    junction_annotation.py -i ${bam} \
        -o ${RSEQC}/junction_annotation/${id} \
        -r ${BED_SUNFLOWER}
done

# Integridad de los transcriptos (TIN) sobre todos los BAMs
tin.py -i mapeos.txt -r ${BED_SUNFLOWER}

# Consolidar resumen de TIN
head -n 1 ${RSEQC}/tin/*.summary.txt > ${RSEQC}/tin/tin_summary.txt
tail -n 1 ${RSEQC}/tin/*.summary.txt | grep ".bam" >> ${RSEQC}/tin/tin_summary.txt

# =============================================================================
# MÓDULO 8: LLAMADO DE VARIANTES EN RNA-seq (girasol)
# =============================================================================

mkdir -p 10_Variant_calling

bcftools mpileup --threads ${THREADS} \
    -f ${REF_SUNFLOWER} \
    -b mapeos.txt \
    | bcftools call --threads ${THREADS} \
        -mv -Ob --ploidy 2 \
        -o 10_Variant_calling/calls.bcf
