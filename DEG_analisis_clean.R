################################################################################
##  Análisis de genes diferencialmente expresados (DEGs)
##  Filtrado a genes candidatos y visualización con diagramas de Venn
################################################################################

# --- Librerías ----------------------------------------------------------------
library(tidyverse)
library(ggVennDiagram)
library(patchwork)

# --- Datos de entrada ---------------------------------------------------------

# Lista de genes candidatos obtenida de la comparación de anotaciones
# entre los genomas HanXRQr2.1 y HanXRQv2 mediante Heliagen
nombres_HanXRQ_v2 <- read_delim("Genes_candidatos/nombres_HanXRQ.v2.txt",
                                 delim = "\t", trim_ws = TRUE)
genes_candidatos <- nombres_HanXRQ_v2$objectA

# --- Extracción de DEGs con patrón 4 por tiempo ------------------------------
# Patrón 4: genes con expresión diferencial en la interacción Genotipo:Tratamiento
# Los resultados de EBSeq están guardados como objetos .RData por DPI

dpis <- c(0, 6, 10, 17)

get_pattern4 <- function(dia, genes_candidatos) {
  load(paste0("EBSeq/DPI", dia, ".RData"))
  genes_p4 <- result$MultiPP$MAP[result$MultiPP$MAP == "pattern4"]
  candidatos <- genes_p4[names(genes_p4) %in% paste0("gene:", genes_candidatos)]
  list(todos      = genes_p4,
       candidatos = candidatos)
}

resultados <- lapply(setNames(dpis, paste0("DPI", dpis)), get_pattern4,
                     genes_candidatos = genes_candidatos)

# --- Diagramas de Venn -------------------------------------------------------

# Función auxiliar para construir el Venn y guardarlo
plot_venn <- function(lista_genes, titulo, subtitulo, archivo) {
  p <- ggVennDiagram(lista_genes) +
    guides(fill = "none") +
    scale_fill_viridis_c() +
    labs(title = titulo, subtitle = subtitulo)

  ggsave(filename = archivo, plot = p,
         device = "png", width = 15, height = 15, dpi = 300, units = "cm")
  p
}

# Venn de genes candidatos con expresión diferencial (patrón 4)
plot_venn(
  lista_genes = lapply(resultados, function(x) names(x$candidatos)),
  titulo      = "Differential expression",
  subtitulo   = "Candidate genes – Pattern 4",
  archivo     = "Genes_candidatos/Venn_pattern4_candidatos.png"
)

# Venn de todos los genes con expresión diferencial (patrón 4)
plot_venn(
  lista_genes = lapply(resultados, function(x) names(x$todos)),
  titulo      = "N° de genes con expresión diferencial",
  subtitulo   = "Interacción Genotipo:Tratamiento – Pattern 4",
  archivo     = "EBSeq/Venn_pattern4_todos.png"
)
