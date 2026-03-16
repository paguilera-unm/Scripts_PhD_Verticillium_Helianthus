# Bioinformatic pipelines — *Verticillium dahliae* genomics and transcriptomics

**Doctoral thesis — Pablo Aguilera**  
Instituto Nacional de Tecnología Agropecuaria (INTA)

---

## Overview

This repository contains the bioinformatic scripts used in the doctoral thesis. The work integrates four analytical axes:

1. **De novo genome assembly** of an Argentine isolate of *V. dahliae* using PacBio HiFi long reads.
2. **Structural and functional genome annotation**, including gene prediction, effector identification, and repeat characterization.
3. **Population genomics** by Genotyping-by-Sequencing (GBS), integrating local isolates with publicly available datasets from multiple countries.
4. **Dual RNA-seq** (*Helianthus annuus* host / *V. dahliae* pathogen) to characterize differentially expressed genes during infection at multiple time points post-inoculation.

---

## Repository structure

```
.
├── Assembly_Vd_script_clean.sh       # De novo genome assembly (PacBio HiFi)
├── Annotation_Vd_script_clean.sh     # Structural and functional annotation
├── gbs_pipeline_clean.sh             # GBS population genomics pipeline
├── population_genomic_clean.R        # Population structure, phylogeny and diversity
├── RNA-seq_pipeline_clean.sh         # Dual RNA-seq mapping and quantification
├── EBseq_script_clean.R              # Differential expression analysis (EBSeq)
├── DEG_analisis_clean.R              # DEG filtering and Venn diagrams
└── TIN_analisis_clean.R              # RNA-seq quality control (TIN, fragment size)
```

---

## Scripts

### 1. `Assembly_Vd_script_clean.sh` — Genome assembly

De novo assembly of the *V. dahliae* VArg1 genome from PacBio HiFi CCS reads. Four assemblers were compared and the best-performing assembly was selected for scaffolding.

**Workflow:**

```
Raw HiFi reads
    ├── FastQC              Quality control of raw reads
    └── Trimmomatic         HEADCROP:1 (remove first nucleotide artifact)
            ├── Canu
            ├── Hifiasm     ← Selected for final assembly
            ├── IPA
            └── Flye
                    ├── QUAST       Contiguity metrics (ref: VdLs17)
                    ├── BUSCO       Gene completeness (hypocreales_odb10)
                    ├── KAT         K-mer spectra analysis
                    ├── Inspector   Error evaluation
                    ├── Merfin      QV* and completeness (k-mer based)
                    └── MUMmer      Synteny plots vs. VdLs17 and 85S
                            └── RagTag scaffold
                                    └── BUSCO / KAT / Merfin / Tapestry (re-evaluation)
```

**Key parameters:**

| Parameter | Value |
|---|---|
| Genome size estimate | 35 Mb |
| Threads | 32 |
| BUSCO database | `hypocreales_odb10` |
| Telomere motif | `TTAGGG` |
| Merfin peak (kcov) | 134 (from GenomeScope) |

---

### 2. `Annotation_Vd_script_clean.sh` — Genome annotation

Structural and functional annotation of the scaffolded VArg1 genome using Funannotate as the main framework, complemented by specialized tools for effectors, transposable elements, AT-rich regions, and telomeres.

**Workflow:**

```
ragtag.scaffold.fasta
    ├── funannotate sort / mask    Contig sorting and repeat soft-masking
    └── funannotate predict        Ab initio gene prediction
            ├── InterProScan       Protein domain annotation
            ├── AntiSMASH          Secondary metabolite clusters
            ├── Phobius            Signal peptides and TM domains
            └── funannotate annotate   Final integrated annotation
                    ├── EffectorP 3.0   Effector prediction (secreted proteins)
                    ├── EDTA            Transposable element annotation (Galaxy)
                    ├── OcculterCut     AT-rich region prediction (accessory genome)
                    └── tidk            Telomere identification and visualization
```

**Key parameters:**

| Parameter | Value |
|---|---|
| BUSCO seed species | `verticillium_longisporum1` |
| Isolate name | VArg1 |
| Threads | 32 |

---

### 3. `gbs_pipeline_clean.sh` — GBS population genomics pipeline

Processing of Genotyping-by-Sequencing (GBS) data from local *V. dahliae* isolates. Restriction enzymes: PstI (rare cutter) and Sau3AI/MboI (frequent cutter). The final dataset integrates local isolates with published datasets from Milgroom (2014) and Bautista-Jalón (2020), plus reference genome sequences.

**Workflow:**

```
Raw reads (paired-end, multiple lanes)
    └── process_radtags      Demultiplexing and QC (Stacks)
            └── BWA-MEM      Mapping to VdLs17 reference genome
                    └── gstacks     Locus assembly
                            └── populations     Population statistics (VCF output)
                                    └── bcftools norm / fixref    Correct REF allele
                                            ├── SNP / Indel filtering (bcftools)
                                            └── bcftools isec + merge
                                                    ├── Milgroom 2014 dataset
                                                    ├── Bautista-Jalón 2020 dataset
                                                    ├── Reference genomes (85S, Vd39, Vd0991, JR2, S011)
                                                    │       (SNPs obtained via DNAdiff in Galaxy)
                                                    ├── Chinese isolates (SRR12769281, SRR12769283)
                                                    └── Second GBS batch (9 additional local isolates)
```

**Key parameters:**

| Parameter | Value |
|---|---|
| Restriction enzymes | PstI + Sau3AI |
| Reference genome | VdLs17 (GCA_000952015.1) |
| Minimum representation (populations) | r = 0.65 |
| Threads | 8 |

---

### 4. `population_genomic_clean.R` — Population structure and diversity

Population genomic analysis on the merged SNP dataset. Includes phylogenetic reconstruction, dimensionality reduction, network analysis, linkage disequilibrium, and admixture estimation.

**Analyses:**

| Module | Method | Output |
|---|---|---|
| Genotype heatmap | gplots::heatmap.2 | Visualization of SNP matrix |
| Phylogenetic tree | NJ with 1000 bootstrap replicates (poppr::aboot) | `.nwk` file + annotated figure |
| PCA | glPca (adegenet) | PC1–PC6 scatter plots |
| Minimum Spanning Network | poppr::poppr.msn | Haplotype network |
| Linkage disequilibrium | ia() index of association | IA and rD values |
| Karyotype plot | CMplot | Genomic distribution of markers |
| Admixture / structure | sNMF (LEA) K = 1–15, 10 reps | Cross-entropy + structure barplot |
| Tree comparison | Cophenetic correlation, cophylo | Validation of reduced SNP set |

**Key packages:** `vcfR`, `poppr`, `adegenet`, `LEA`, `dartR`, `ggtree`, `ape`

---

### 5. `RNA-seq_pipeline_clean.sh` — Dual RNA-seq pipeline

Mapping and quantification pipeline for paired-end RNA-seq from sunflower plants inoculated with *V. dahliae*. A sequential mapping strategy is used: reads are first mapped against the host genome; unaligned reads are then mapped against the *V. dahliae* genome to capture pathogen transcripts.

**Workflow:**

```
Raw reads (paired-end, multiple runs)
    └── FastQC              Quality control
            ├── STAR index  HanXRQr2.1 (Helianthus annuus)
            └── STAR index  VArg1 (V. dahliae, this study)
                    └── STAR map → H. annuus (HanXRQr2.1)
                            ├── MultiQC         Mapping report
                            ├── mmquant         Read quantification (sunflower)
                            └── Unmapped reads
                                    └── STAR map → V. dahliae (VArg1)
                                            ├── MultiQC     Mapping report
                                            └── mmquant     Read quantification (Vd)
                                    └── RSeQC
                                            ├── bam_stat.py
                                            ├── inner_distance.py
                                            ├── junction_annotation.py
                                            └── tin.py
                                    └── bcftools mpileup    Variant calling
```

**Key parameters:**

| Parameter | Value |
|---|---|
| Sunflower genome | HanXRQr2.0-SUNRISE-2.1 |
| Vd genome | VArg1 (this study, annotated) |
| sjdbOverhang (sunflower) | 74 |
| sjdbOverhang (Vd) | 100 |
| Threads | 64 |

---

### 6. `EBseq_script_clean.R` — Differential expression analysis

Bayesian differential expression analysis using EBSeq. The analysis is performed independently for each time point (DPI 0, 6, 10, 17). The experimental design includes two sunflower genotypes (RHA266 susceptible, HA89 resistant) × two treatments (Testigo / Inoculado).

**Key concepts:**

- **Patterns**: EBSeq assigns each gene to a pattern based on posterior probabilities. Pattern 4 corresponds to genes differentially expressed in the Genotype × Treatment interaction.
- **Output per DPI**: normalized count tables, log2 fold-change tables per pattern, QQ-plots, density histograms, heatmaps, and Gene Enrichment Analysis (GEA) results.

**Workflow:**

```
mmquant count table + metadata
    └── GetPatterns()      Define expression patterns from experimental design
            └── EBMultiTest()   Bayesian multi-condition test (per DPI)
                    ├── GetMultiPP()    Posterior probabilities per pattern
                    ├── GetMultiFC()    Log2 fold-change per pattern
                    ├── QQP() / DenNHist()   Model fit diagnostics
                    ├── heatmap.2()     Normalized expression heatmaps
                    └── auxfun()        Gene Enrichment Analysis
```

**Key packages:** `EBSeq`, `tidyverse`, `factoextra`

---

### 7. `DEG_analisis_clean.R` — DEG filtering and Venn diagrams

Downstream analysis of EBSeq results. Filters differentially expressed genes (Pattern 4, Genotype × Treatment interaction) and intersects them with a list of candidate genes. Produces Venn diagrams comparing the sets of DEGs across the four time points.

**Input:** `.RData` result objects from `EBseq_script_clean.R`  
**Output:** Venn diagrams (.png) for candidate genes and all DEGs across DPI 0, 6, 10, 17.

**Key packages:** `tidyverse`, `ggVennDiagram`, `patchwork`

---

### 8. `TIN_analisis_clean.R` — RNA-seq quality control

Post-alignment quality control using the Transcript Integrity Number (TIN) computed by RSeQC's `tin.py`. TIN assesses RNA degradation at the transcript level; values close to 100 indicate high RNA integrity. Fragment size distribution is also evaluated.

**Modules:**

| Module | Description | Output |
|---|---|---|
| TIN per sample | Mean TIN ± SD per sample, faceted by DPI | `TIN_plot.png` |
| TIN distribution | Density plot of per-transcript TIN values | `TIN_distribution.png` |
| Fragment size | Count distribution, size density, and size per contig | `fragment_size.png` |

**Key packages:** `tidyverse`, `patchwork`

---

## Dependencies

### Bash / command-line tools

| Tool | Version tested | Usage |
|---|---|---|
| FastQC | 0.11+ | Raw read QC |
| Trimmomatic | 0.39 | Read trimming |
| Canu | 2.2 | Genome assembler |
| Hifiasm | 0.16+ | Genome assembler |
| IPA | 1.8+ | Genome assembler |
| Flye | 2.9+ | Genome assembler |
| QUAST | 5.2+ | Assembly quality |
| BUSCO | 5.4+ | Gene completeness |
| KAT | 2.4+ | K-mer analysis |
| Inspector | 1.2 | Assembly error evaluation |
| Meryl / Merfin | — | K-mer QV* |
| MUMmer (nucmer, mummerplot) | 4.0+ | Synteny alignment |
| RagTag | 2.1+ | Scaffolding |
| Tapestry (weave) | — | Telomere detection |
| Funannotate | 1.8+ | Genome annotation |
| InterProScan | 5.60-92.0 | Protein domain annotation |
| EffectorP | 3.0 | Effector prediction |
| OcculterCut | 1.1 | AT-rich region prediction |
| tidk | 0.2+ | Telomere identification |
| process_radtags | Stacks 2.6+ | GBS demultiplexing |
| BWA | 0.7.17 | Read mapping |
| SAMtools | 1.16+ | BAM processing |
| gstacks / populations | Stacks 2.6+ | GBS locus assembly |
| bcftools | 1.10+ | VCF manipulation |
| bgzip / tabix | htslib 1.10+ | VCF compression and indexing |
| TASSEL | 5 | HapMap to VCF conversion |
| STAR | 2.7.10+ | RNA-seq mapping |
| mmquant | — | Multi-mapping quantification |
| MultiQC | 1.14+ | Mapping QC reports |
| RSeQC | 4.0+ | RNA-seq quality metrics |

### R packages

| Package | Usage |
|---|---|
| `EBSeq` | Differential expression |
| `vcfR` | VCF reading and manipulation |
| `poppr` | Population genetics |
| `adegenet` | Multivariate population genetics |
| `LEA` | Admixture / sNMF |
| `dartR` | Genlight utilities |
| `ggtree` / `treeio` | Phylogenetic tree visualization |
| `ape` | Phylogenetics |
| `phangorn` | Tree comparison |
| `CMplot` | Karyotype / Manhattan plots |
| `ggVennDiagram` | Venn diagrams |
| `ggpubr` | Plot arrangement |
| `patchwork` | Plot composition |
| `tidyverse` | Data manipulation and visualization |

---

## Reference genomes

| Genome | Accession | Usage |
|---|---|---|
| *V. dahliae* VdLs17 | GCA_000952015.1 | Primary reference (GBS, RNA-seq) |
| *V. dahliae* 85S | GCA_004798885.1 | Comparative genomics |
| *V. dahliae* JR2 | GCA_000400815.2 | Comparative genomics |
| *H. annuus* HanXRQr2.1 | — | RNA-seq host reference |
| *V. dahliae* Colon | This study | RNA-seq pathogen reference |

---

## Notes

- Scripts are written to be run sequentially within each analytical axis. Output directories of earlier steps are input for later ones.
- GBS analysis was run on a local Linux workstation. Assembly and RNA-seq pipelines were run on a shared HPC cluster.
- Repeat annotation with EDTA and some InterProScan jobs were executed through the Galaxy Europe platform (`use.galaxy.eu`).
- The `popmap.txt` file (required by Stacks) contains two tab-separated columns: sample name and population ID.

---

## Contact

Pablo Aguilera  
INTA — Instituto Nacional de Tecnología Agropecuaria  
aguilera.pablo@inta.gob.ar
