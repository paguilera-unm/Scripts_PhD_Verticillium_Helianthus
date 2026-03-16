################################################################################
##  Análisis de genómica poblacional de Verticillium dahliae
##  Métodos: árbol filogenético NJ, PCA, red de haplotipos MSN, estructura
##           poblacional (sNMF) y diversidad genética
################################################################################

# --- Librerías ----------------------------------------------------------------
library(readxl)
library(tidyverse)
library(vcfR)
library(poppr)
library(ape)
library(adegenet)
library(LEA)
library(dartR)
library(phylogram)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(ggtree)
library(treeio)
library(phangorn)
library(CMplot)
source("PlotK.R")

# Funciones internas de tidytree requeridas por ggtree
nodeid.tbl_tree         <- utils::getFromNamespace("nodeid.tbl_tree",         "tidytree")
rootnode.tbl_tree       <- utils::getFromNamespace("rootnode.tbl_tree",       "tidytree")
offspring.tbl_tree      <- utils::getFromNamespace("offspring.tbl_tree",      "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item","tidytree")
child.tbl_tree          <- utils::getFromNamespace("child.tbl_tree",          "tidytree")
parent.tbl_tree         <- utils::getFromNamespace("parent.tbl_tree",         "tidytree")

# --- Parámetros ---------------------------------------------------------------

COL_LINAJES <- c("#4DAF4A","#984EA3","#A65628","#FF7F00","#A4D3EE",
                 "#E41A1C","#377EB8","#F781BF","#FFFF33")
OUTDIR      <- "19_Resultados"
N_BOOTSTRAP <- 1000
SET_SEED    <- 9

# =============================================================================
# MÓDULO 1: CARGA Y PREPARACIÓN DE DATOS
# =============================================================================

vcf      <- read.vcfR("18_Denis_10_2022/Locales_laura_filt_fix.vcf.gz")
metadata <- read.delim("metadata.tsv")

# Convertir a genlight (haploides)
gl.rubi  <- vcfR2genlight(vcf)
ploidy(gl.rubi) <- 1

# Alinear metadata con las muestras presentes en el VCF
metadata <- metadata %>%
  slice(match(colnames(vcf@gt), Isolate))

other(gl.rubi) <- metadata$Lineageh

# =============================================================================
# MÓDULO 2: HEATMAP DE GENOTIPOS
# =============================================================================

tabla <- as.data.frame(vcf@gt)

# Codificar genotipos como 0/1/NA
tabla_num <- tabla
for (col in names(tabla)) {
  tabla_num[[col]] <- ifelse(tabla[[col]] == "0/0", 0,
                      ifelse(tabla[[col]] == "1/1", 1, NA))
}

# Subconjunto balanceado: se excluyen muestras del linaje 4B (sobrerrepresentado)
n_4B_excl <- 87
heatmap_df <- tabla_num[-1] %>%
  select(-sample(metadata$Isolate[metadata$Lineageh == "4B" &
                                  metadata$Dataset == "Baustista-Jalon et al. 2021"],
                 n_4B_excl)) %>%
  as.matrix()

metadata_hm <- metadata %>% filter(Isolate %in% colnames(heatmap_df))

png(filename = file.path(OUTDIR, "Heatmap_LNJ_2024.png"),
    width = 60, height = 19, units = "cm", res = 300)

gplots::heatmap.2(heatmap_df,
  dendrogram   = "column",
  col          = c("steelblue", "darkgreen"),
  breaks       = c(0, 0.5, 1),
  key          = FALSE,
  labRow       = FALSE,
  margins      = c(12, 1),
  cexCol       = 1.5,
  trace        = "none",
  ColSideColors = COL_LINAJES[as.factor(metadata_hm$Lineageh)],
  lwid         = c(0.01, 5)
)

dev.off()

# =============================================================================
# MÓDULO 3: ÁRBOL FILOGENÉTICO (NJ con bootstrap)
# =============================================================================

set.seed(SET_SEED)

# Clustering para definir grupos
gl.clust <- find.clusters.genlight(gl.rubi, stat = "BIC",
                                   max.n.clust = 10, n.clust = 9, n.pca = 12)
pop(gl.rubi) <- gl.clust$grp

# Árbol NJ
tree <- aboot(gl.rubi, tree = "nj", distance = bitwise.dist,
              sample = N_BOOTSTRAP, showtree = FALSE, cutoff = 0, quiet = TRUE)
write.tree(tree, file = file.path(OUTDIR, "Todos_nj.nwk"), digits = 10)

arbol <- read.newick(file.path(OUTDIR, "Todos_nj.nwk"))

ggarbol_todos <- root(arbol, node = 366, edgelabel = TRUE) %>%
  drop.tip(sample(metadata$Isolate[metadata$Lineageh == "4B" &
                                   metadata$Dataset == "Baustista-Jalon et al. 2021"],
                  n_4B_excl)) %>%
  groupOTU(., list(
    Baustista_Jalon_2021 = metadata$Isolate[metadata$Dataset == "Baustista-Jalon et al. 2021"],
    NCBI_genomes         = c(metadata$Isolate[metadata$Dataset == "NCBI genomes"], "Getta_Getta"),
    This_study           = metadata$Isolate[metadata$Dataset == "This study"]
  )) %>%
  ggtree(linewidth = 1.2) +
  geom_tiplab(aes(color = group), size = 3.88) +
  geom_label2(aes(subset = node > 135 & (label > 70 | label == 100), label = label),
              fill = "white", label.padding = unit(0.01, "lines"),
              label.size = 0, label.r = unit(0, "lines"), size = 2, alpha = 0.7) +
  geom_cladelab(node  = c(163, 140, 187, 186, 209, 210, 230, 256, 251),
                label = c("2A","4B","2B824","2BR1","6","4A","Local","1A/B","2B334"),
                offset = 0.18, offset.text = 0.01, align = FALSE,
                fontsize = 6, barsize = 0.8) +
  geom_cladelab(node  = c(250, 139, 183),
                label = c("Clade I","Subclade II-1","Subclade II-2"),
                offset = 0.33, offset.text = 0.05, align = TRUE,
                fontsize = 8, barsize = 0.9, angle = 270) +
  geom_rootedge(rootedge = 0.1, linewidth = 1.2) +
  xlim(0, 1.2) +
  scale_color_manual(values = c(Baustista_Jalon_2021 = "grey40",
                                NCBI_genomes         = "#005CAB",
                                This_study           = "#E31B23")) +
  theme_tree2() +
  labs(color = NULL) +
  guides(color = "none")

ggsave(filename = file.path(OUTDIR, "Todos_nj_fix.png"),
       device = "png", width = 25, height = 45, units = "cm", dpi = 300)

# =============================================================================
# MÓDULO 4: ANÁLISIS DE COMPONENTES PRINCIPALES (PCA)
# =============================================================================

set.seed(SET_SEED)
rubi.pca  <- glPca(gl.rubi, nf = 8)
exp_var   <- 100 * rubi.pca$eig[1:8] / sum(rubi.pca$eig[1:8])

pca_scores <- as.data.frame(rubi.pca$scores) %>%
  mutate(linage = other(gl.rubi)[[1]])

# Función para graficar un par de componentes
plot_pca <- function(scores, pc_x, pc_y, exp) {
  ggplot(scores, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y),
                             colour = "linage")) +
    ggforce::geom_mark_ellipse(aes(label = linage, group = linage, fill = linage),
                               alpha = 0.2, linewidth = 1, expand = 0.01) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2.5) +
    guides(fill = "none") +
    theme_bw() +
    labs(x      = paste0("PC", pc_x, " (", round(exp[pc_x], 2), "%)"),
         y      = paste0("PC", pc_y, " (", round(exp[pc_y], 2), "%)"),
         colour = "Lineage",
         title  = "Principal component analysis")
}

pca_plot1_2 <- plot_pca(pca_scores, 1, 2, exp_var)
pca_plot3_4 <- plot_pca(pca_scores, 3, 4, exp_var)
pca_plot5_6 <- plot_pca(pca_scores, 5, 6, exp_var)

var_exp <- ggplot(mapping = aes(x = 1:8, y = exp_var)) +
  geom_col(fill = "#3498DB") +
  theme_classic() +
  scale_x_continuous(breaks = 1:8) +
  labs(title = "PCA Eigenvalues",
       y     = "Percent of variance\nexplained",
       x     = "Eigenvalues")

pca_arreglo <- ggarrange(pca_plot1_2, pca_plot3_4, pca_plot5_6, var_exp,
                         common.legend = TRUE, labels = c("A","B","C","D"),
                         ncol = 2, nrow = 2)

ggsave(plot = pca_arreglo,
       filename = file.path(OUTDIR, "Todos_PCA_fix.png"),
       device = "png", width = 35, height = 28, units = "cm", dpi = 300, bg = "white")

# =============================================================================
# MÓDULO 5: RED DE HAPLOTIPOS (Minimum Spanning Network)
# =============================================================================

SNP_genind     <- vcfR2genind(vcf)
pop(SNP_genind) <- gl.clust$grp

SNP_dist   <- bitwise.dist(SNP_genind, percent = TRUE, missing_match = TRUE)
min_span   <- poppr.msn(SNP_genind, SNP_dist, showplot = FALSE, include.ties = TRUE)

set.seed(69)
png(filename = file.path(OUTDIR, "MSN_Todos.png"),
    width = 25, height = 25, units = "cm", res = 300)
plot_poppr_msn(SNP_genind, min_span,
               inds       = "ALL",
               mlg        = FALSE,
               gadj       = 3,
               nodescale  = 75,
               palette    = COL_LINAJES,
               cutoff     = NULL,
               quantiles  = FALSE,
               beforecut  = TRUE,
               layfun     = igraph::layout_with_lgl)
dev.off()

# =============================================================================
# MÓDULO 6: DESEQUILIBRIO DE LIGAMIENTO
# =============================================================================

pop(SNP_genind) <- metadata$Lineageh

div <- poppr(SNP_genind)

SNP_genind %>%
  popsub(., c("1","2","4","5")) %>%
  clonecorrect() %>%
  ia(sample = 999)

# =============================================================================
# MÓDULO 7: DISTRIBUCIÓN GENÓMICA DE LOS MARCADORES (karyotype plot)
# =============================================================================

stacks <- read.table("18_Denis_10_2022/Locales_laura_filt_fix.vcf.gz",
                     stringsAsFactors = FALSE)
pig <- data.frame(SNP        = 1:nrow(stacks),
                  Chromosome = stacks$V1,
                  Position   = stacks$V2,
                  trait1     = rep(0.5, nrow(stacks)))

CMplot(pig, type = "p", plot.type = "d", bin.size = 1e6,
       chr.den.col = c("darkgreen","yellow","red"),
       file = "tiff", dpi = 300, file.output = TRUE, verbose = TRUE,
       width = 9, height = 6)

# =============================================================================
# MÓDULO 8: ESTRUCTURA POBLACIONAL (sNMF – LEA)
# =============================================================================

gl <- gl.compliance.check(x = gl.rubi) %>% gl.filter.monomorphs()
dartR::gl2geno(x = gl, outfile = "Todos", outpath = OUTDIR)

set.seed(65)
project <- snmf(input.file = file.path(OUTDIR, "Todos.geno"),
                K = 1:15, project = "new", entropy = TRUE, repetitions = 10)

# Gráfico de cross-entropy para selección de K
cross_plot <- PlotK(project)
ggsave(plot = cross_plot,
       filename = file.path(OUTDIR, "cross_entropy_Todos.png"),
       device = "png", width = 20, height = 20, units = "cm", dpi = 300, bg = "white")

# Coeficientes de ancestría para K óptimo (K = 11, basado en cross-entropy)
K_opt <- 11
Qmatrix <- Q(project, K = K_opt, run = which.min(cross.entropy(project, K = K_opt)))
row.names(Qmatrix) <- gl@ind.names

structure_plot <- Qmatrix %>%
  t() %>%
  reshape2::melt() %>%
  mutate(linage = rep(metadata$Lineageh, each = K_opt)) %>%
  filter(!Var2 %in% sample(metadata$Isolate[metadata$Lineageh == "4B" &
                           metadata$Dataset == "Baustista-Jalon et al. 2021"], 87)) %>%
  ggplot() +
  aes(x = Var2, y = value, fill = Var1) +
  geom_col() +
  facet_grid(~ linage, scales = "free", space = "free") +
  scale_fill_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                               "#FFFF33","#A65628","#F781BF","#999999","#66CD00","#A4D3EE")) +
  labs(y = "Admixture coefficient", x = NULL, fill = "Genotype") +
  theme_pubr() +
  guides(fill = "none") +
  theme(axis.text.y = element_text(size = 5))

ggsave(filename = file.path(OUTDIR, "Todos_snmf.png"),
       plot = structure_plot,
       device = "png", width = 20, height = 50, units = "cm", dpi = 300, bg = "white")

# Grilla final combinando árbol, PCA y estructura
grilla_max <- ggarrange(ggarbol_todos, pca_arreglo, structure_plot,
                        ncol = 1, labels = c("A","B","C"))

ggsave(filename = file.path(OUTDIR, "Locales_pop.png"),
       plot = grilla_max,
       device = "png", width = 40, height = 55, units = "cm", dpi = 300, bg = "white")

# =============================================================================
# MÓDULO 9: COMPARACIÓN DE ÁRBOLES (subconjunto de SNPs vs. set completo)
# =============================================================================

tree1 <- read.tree(file.path(OUTDIR, "Laura_nj.nwk"))
tree2 <- read.tree(file.path(OUTDIR, "Laura_filt_nj.nwk"))

# Correlación cofenética
comp1 <- cophenetic.phylo(tree1)
comp2 <- cophenetic.phylo(tree2)
cor.test(comp1, comp2)

# Visualización
png(filename = file.path(OUTDIR, "Laura_vs_Laura_filt.png"),
    width = 26, height = 50, units = "cm", res = 100)
plot(cophylo(multi2di(tree1), multi2di(tree2), rotate.multi = FALSE))
dev.off()
