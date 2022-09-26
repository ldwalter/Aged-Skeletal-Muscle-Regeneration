########################################################
# Cell Type Identification Supplemental Figure
# Lauren Walter - September 2022
########################################################

# Load libraries ####
library(cowplot)
library(ggplot2)
library(Seurat)

# Source file ####
setwd("/FILE/PATH/TO/SEURAT_OBJECT")

seur.obj <- readRDS(file = "SEURAT_OBJECT.rds", refhook = NULL)

# Dot plots with markers used to identify cell types ####
# Define order of cell types displayed on the y-axis
positions_list <- list(c("Erythrocytes", "B cells", "NK cells", "T cells (Cd4+)", "T cells (Cycling; Cd3e+)",  "T cells (Non-cycling; Cd3e+)",
                         "Neutrophils", "Dendritic cells (Cd209a+)", "Dendritic cells (Cd72+)", "Dendritic cells (Xcr1+)", "Dendritic cells (Fscn1+)",
                         "M2 Macrophages (Cx3cr1+)", "M1/M2 Macrophages (Mrc1+)","M1 Macrophages (Ccr2+)", "Monocytes/Macrophages (Patrolling; Ctsa+)",
                         "Monocytes/Macrophages (Cxcl10+)", "Monocytes (Cycling; Cdk1+)"),
                       c("Schwann and Neural/Glial cells", "Tenocytes", "FAPs (Adipogenic)", "FAPs (Pro-remodeling)", "FAPs (Stem)"),
                       c("Endothelial and Myeloid cells", "Endothelial cells (Vein)", "Endothelial cells (Capillary)", "Endothelial cells (Artery)",
                         "Pericytes and Smooth muscle cells", "Myonuclei", "MuSCs and progenitors"))

# Define cell type markers to display
markers_list <- list(c("Smc4", "Cdk1", "Hmgb2", "Cxcl10", "Ctsa", "Ctsb", "Ctsl", "Ctsz", "Fabp5", "Ly6c2", "Ccr2", "Ccl9", "Ccl6",
                       "Tgfbi", "Ccl2", "Cd68", "Ccl8", "Lyz2", "Itgam", "Cd163", "Mrc1", "C1qa", "Cx3cr1", "Mapk14", "Fcgr3",
                       "Csf1r", "Cxcl16", "Ctss", "Cd200r1", "Cd74", "Cd83", "Fscn1", "Flt3", "Xcr1", "Cd72", "Cd209a",
                       "Cxcr4", "Csf1", "S100a8", "S100a9", "Sell", "Cd14", "Mmp9", "Ptprc", "Nkg7", "Ccl5", "Cd8a", "Cd8b1", "Cd3e",
                       "Cd4", "Foxp3", "Cxcr3", "Cd27", "Tbx21", "Gata3", "Gzma", "Klra4", "Klra7", "Klra8", "Klre1", "Klrd1", "Ccr7",
                       "Igkc", "Ms4a1", "Cd19", "Cd22", "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"),
                     c("Igfbp5", "Dpp4", "Thy1", "Tnfaip6", "Cd164", "Egfr", "Cd34", "Il33", "Hmgb2", "Cdk1", "Smc4", "Mcm5", "Tyms",
                       "Rrm2", "Dio2", "Rspo3", "Adam12", "Bmp5", "Myoc", "Col3a1", "Bmp1", "Gsn", "Pdgfra", "Bgn", "Hdlbp", "Mmp14",
                       "Ctsk", "Col1a1", "Dcn", "Mmp2", "Apod", "Tnmd", "Scx", "Ptn", "Mpz"),
                     c("Pax7", "Sdc1", "Acta1", "Myh1", "Myh4", "Acta2", "Apold1", "Ednrb", "Rgs5", "Myl9", "Myh11", "Cspg4",
                       "Pdgfrb", "Mcam", "Alpl", "Hey1", "Cdh5", "Pecam1", "Cxcl12", "Ly6a", "Kdr", "Lpl", "Cd34", "Vwf", "Hif1a",
                       "Icam1", "Lrg1", "Aplnr", "Apln", "S100a8", "S100a9", "Csf1", "Cxcr4", "Itgam", "Smc4", "Cdk1", "Hmgb2"))

# Dot plots 
plot_list <- list()

for (i in 1:3) {
  plot_list[[i]] <- DotPlot(object = seur.obj,
                            assay = "RNA",
                            features = markers_list[[i]],
                            cols = c("grey", "blue"),
                            col.min = -2.5,
                            col.max = 2.5,
                            dot.min = 0,
                            dot.scale = 2,
                            group.by = "Specific_cell_types",
                            split.by = NULL,
                            scale = TRUE,
                            scale.by = "radius",
                            scale.min = NA,
                            scale.max = NA
  ) +
    RotatedAxis() +
    labs(x = NULL, y = NULL) +
    scale_y_discrete(limits = positions_list[[i]]) +
    theme(axis.line = element_line(color = "black", size = 0.25), 
          axis.ticks = element_line(color = "black", size = 0.25),
          axis.text.x = element_text(color = "black", size = 6, face = "italic", angle = 90, vjust = 0.5),
          axis.text.y = element_text(color = "black",  size = 6),
          legend.text = element_text(color = "black", size = 6),
          legend.title = element_text(color = "black", size = 6, hjust = 0.5),
          legend.key.size = unit(0.20, "cm"))
  }

# Organize the plots together
DP <- plot_grid(plot_list[[1]] + theme(legend.position = "none"), 
                plot_list[[2]] + theme(legend.position = "none"),
                plot_list[[3]] + theme(legend.position = "none"),
                labels = c("(A)", "(B)", "(C)"), label_size = 8, nrow = 3, align = "hv", rel_heights = c(3.4, 1, 1.4))

# Add legend
legend <- get_legend(plot_list[[1]])

plot_grid(DP, legend, ncol = 2, rel_widths = c(1,	0.10))

# Save plot as PDF
ggsave(filename = "CT_ID_SuppFig.pdf", plot = last_plot(), width = 8.5, height = 7, units = "in")
