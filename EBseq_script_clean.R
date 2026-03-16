################################################################################
##  Análisis de expresión diferencial con EBSeq
##  Diseño: comparación Genotipo x Tratamiento a cuatro tiempos post-infección
##  (DPI 0, 6, 10 y 17)
################################################################################

# --- Librerías ----------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(EBSeq)
library(factoextra)

# --- Datos de entrada ---------------------------------------------------------

# Tabla de conteos crudos (genes x muestras)
mmquant <- read_delim("mmquant_unicos.tsv", delim = "\t", trim_ws = TRUE) %>%
  data.frame() %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Metadatos: información de cada muestra
metadata <- read_delim("metadata.tsv", delim = "\t", trim_ws = TRUE) %>%
  mutate(
    Genotipo      = fct_relevel(Genotipo,    "RHA266"),
    Tratamiento   = fct_relevel(Tratamiento, "Testigo"),
    DPI           = as.factor(DPI),
    Grupos        = fct_relevel(paste(Genotipo, Tratamiento, DPI, sep = "_"),
                                "RHA266_Testigo_0"),
    Tratamiento_Genotipo = paste(Tratamiento, Genotipo, sep = "_")
  ) %>%
  column_to_rownames("Sample")

# --- Diseño experimental ------------------------------------------------------
# Se generan todos los patrones de expresión posibles dado el número de grupos
PosParti <- GetPatterns(paste(metadata$Genotipo, metadata$Tratamiento, sep = "_"))

# El diseño se exporta, se edita manualmente para seleccionar patrones biológicamente
# relevantes, y se vuelve a importar
write.table(PosParti, file = "EBseq_diseno.tsv",
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

EBseq_diseno <- read.delim("EBseq_diseno.tsv", row.names = 1) %>% as.matrix()

PlotPattern(EBseq_diseno)

# --- Análisis EBSeq por tiempo post-infección ---------------------------------
# El análisis se realiza de forma independiente para cada DPI

dias_dpi <- as.numeric(levels(metadata$DPI))  # 0, 6, 10, 17

for (dia in dias_dpi) {

  # Subconjunto de datos correspondiente al DPI actual
  table   <- mmquant %>% as.data.frame() %>%
               select(row.names(filter(metadata, DPI == dia))) %>%
               as.matrix()
  diseno  <- filter(metadata, DPI == dia)

  # Normalización y test multivariado
  MultiSize <- MedianNorm(table)
  MultiOut  <- EBMultiTest(table,
                            Conditions   = diseno$Tratamiento_Genotipo,
                            sizeFactors  = MultiSize,
                            AllParti     = EBseq_diseno)

  # Probabilidades posteriores y fold-change
  MultiPP <- GetMultiPP(MultiOut)
  MultiFC <- GetMultiFC(MultiOut)

  # Guardar objeto completo
  result <- list(MultiOut  = MultiOut,
                 MultiPP   = MultiPP,
                 MultiFC   = MultiFC,
                 MultiSize = MultiSize)
  save(result, file = paste0("EBSeq/DPI", dia, ".RData"))

  # --- Diagnóstico del ajuste del modelo -------------------------------------
  outdir <- paste0("EBSeq/DPI", dia)
  dir.create(outdir, showWarnings = FALSE)

  png(filename = paste0(outdir, "/QQplot_DPI", dia, ".png"), width = 512, height = 512)
    par(mfrow = c(2, 2))
    QQP(EBOut = result$MultiOut)
  dev.off()

  png(filename = paste0(outdir, "/DenNHist_DPI", dia, ".png"), width = 512, height = 512)
    par(mfrow = c(2, 2))
    DenNHist(result$MultiOut)
  dev.off()

  # --- Exportar tablas de log fold-change para los patrones de interés ------
  for (patron in c("pattern2", "pattern4")) {

    genes_patron <- names(result$MultiPP$MAP[result$MultiPP$MAP == patron])

    logFC <- data.frame(result$MultiFC$Log2PostFCMat) %>%
      filter(row.names(.) %in% genes_patron) %>%
      mutate(GeneID = gsub("gene:", "", genes_patron)) %>%
      dplyr::select(GeneID, everything())

    write.table(logFC,
                file  = paste0(outdir, "/logFC_DPI_", dia, "_", patron, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }

  # --- Heatmaps de expresión normalizada -------------------------------------
  metadata_filt <- filter(metadata, DPI == dia)

  for (patron in c("pattern1", "pattern2", "pattern3", "pattern4")) {

    png(filename = paste0(outdir, "/heatmap_", patron, "_DPI", dia, ".png"),
        width = 700, height = 550)

    data.frame(result$MultiOut$DataNorm) %>%
      rename_at(vars(colnames(result$MultiOut$DataNorm)),
                ~ paste(metadata_filt$Grupos, metadata_filt$Pool, sep = "_")) %>%
      filter(row.names(.) %in%
               names(result$MultiPP$MAP[result$MultiPP$MAP == patron])) %>%
      as.matrix() %>% t() %>%
      heatmap.2(scale = "column", trace = "none", Colv = TRUE,
                margins = c(10, 15), xlab = "",
                main = paste("Expresión:", patron, "– DPI", dia))

    dev.off()
  }

  # --- Análisis de enriquecimiento funcional (Gene Enrichment Analysis) ------
  for (patron in c("pattern2", "pattern4")) {

    GEA <- auxfun(
      myInterestingGenes = gsub("gene:", "",
                                names(result$MultiPP$MAP[result$MultiPP$MAP == patron])),
      geneID2GO = geneID2GO
    )

    write.table(GEA,
                file  = paste0("GEA_result/EBseq/DPI_", dia, "_", patron, ".tsv"),
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  }
}

# --- Exportar expresión normalizada de todos los DPI -------------------------
# Tabla con cuentas normalizadas y asignación de patrón para cada gen

metadata <- metadata %>%
  mutate(Sample = paste(Grupos, Pool, sep = "_"))

for (dpi_id in c("DPI0", "DPI6", "DPI10", "DPI17")) {

  load(paste0("EBSeq/", dpi_id, ".RData"))

  nDPI <- as.numeric(str_split(dpi_id, "I")[[1]][2])
  metadata_filt <- filter(metadata, DPI == nDPI)

  data.frame(result$MultiOut$DataNorm) %>%
    rename_at(vars(colnames(result$MultiOut$DataNorm)),
              ~ metadata_filt$Sample) %>%
    mutate(GeneID = row.names(.)) %>%
    left_join(data.frame(GeneID  = names(result$MultiPP$MAP),
                         Pattern = result$MultiPP$MAP),
              by = "GeneID") %>%
    select(GeneID, Pattern, everything()) %>%
    write.table(file  = paste0("EBSeq/", dpi_id, "_norm.tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
