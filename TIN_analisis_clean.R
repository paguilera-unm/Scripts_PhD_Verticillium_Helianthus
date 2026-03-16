################################################################################
##  Control de calidad de RNA-seq: Transcript Integrity Number (TIN)
##  y distribución del tamaño de fragmento (RSeQC)
################################################################################

# --- Librerías ----------------------------------------------------------------
library(tidyverse)
library(patchwork)

# --- Módulo 1: TIN por muestra ------------------------------------------------
# El TIN mide la integridad de los transcriptos a nivel de muestra.
# Valores cercanos a 100 indican RNA de alta integridad.

tin_summary <- read_delim("RSeQC/tin_summary.txt",
                           delim = "\t", trim_ws = TRUE)

metadata <- read_delim("metadata.tsv", delim = "\t", trim_ws = TRUE)

grafico_tin <- tin_summary %>%
  mutate(
    Bam_file = str_remove(Bam_file, "_.*"),
    Bam_file = str_replace_all(Bam_file, "-", ".")
  ) %>%
  left_join(metadata, by = c("Bam_file" = "Sample")) %>%
  ggplot() +
  aes(x      = Bam_file,
      y      = `TIN(mean)`,
      color  = paste(Genotipo, Tratamiento, sep = "_")) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = `TIN(mean)` - `TIN(stdev)`,
                    ymax = `TIN(mean)` + `TIN(stdev)`),
                width = 0.6, linewidth = 2) +
  facet_grid(~ DPI, scales = "free", space = "free") +
  scale_color_brewer(palette = "Set2") +
  ylim(0, 100) +
  labs(x      = "Muestra",
       y      = "TIN",
       title  = "Transcript Integrity Number (TIN)",
       color  = "Genotipo_Tratamiento") +
  guides(x = "none") +
  theme_light()

ggsave(filename = "TIN_plot.png", plot = grafico_tin,
       device = "png", width = 25, height = 10, dpi = 300, units = "cm")

# --- Módulo 2: Distribución del TIN por transcripto --------------------------
# Se evalúa la distribución del TIN a nivel de transcripto individual
# en la muestra con menor profundidad

transcriptos <- read_delim("RSeQC/sar1178-2023-69_S46Aligned.sortedByCoord.out.tin.xls",
                            delim = "\t", trim_ws = TRUE)

TIN_distribution <- transcriptos %>%
  filter(TIN > 0) %>%
  ggplot() +
  aes(x = TIN) +
  geom_density(linewidth = 2) +
  theme_classic() +
  labs(title    = "Distribución del TIN por transcripto",
       subtitle = "Muestra con menor profundidad de secuenciación")

ggsave(filename = "TIN_distribution.png", plot = TIN_distribution,
       device = "png", width = 15, height = 10, dpi = 300, units = "cm")

# --- Módulo 3: Tamaño de fragmento -------------------------------------------

frag <- read_delim("RSeQC/sar1178-2023-69_S4_frag_size.txt",
                   delim = "\t", trim_ws = TRUE)

# Distribución del conteo de anotaciones (escala log10)
conteos <- frag %>%
  filter(frag_count > 0) %>%
  ggplot() +
  aes(x = log10(frag_count)) +
  geom_density(linewidth = 2) +
  theme_light() +
  labs(title = "Conteo de anotaciones (log10)")

# Distribución del tamaño de fragmento (línea de referencia: 180 pb)
dispersion <- frag %>%
  filter(frag_count > 0) %>%
  ggplot() +
  aes(x = frag_mean) +
  geom_density(linewidth = 2) +
  geom_vline(xintercept = 180, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "Distribución del tamaño de fragmento",
       x     = "Tamaño de fragmento (pb)")

# Tamaño de fragmento por cromosoma/contig
cromosomas <- frag %>%
  filter(frag_count > 0) %>%
  ggplot() +
  aes(x = chrom, y = frag_mean) +
  geom_boxplot() +
  geom_hline(yintercept = 180, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_light() +
  labs(title = "Tamaño de fragmento por cromosoma/contig",
       y     = "Tamaño de fragmento (pb)",
       x     = "Cromosoma")

grilla <- conteos / dispersion | cromosomas

ggsave(filename = "fragment_size.png", plot = grilla,
       device = "png", width = 35, height = 25, dpi = 300, units = "cm")
