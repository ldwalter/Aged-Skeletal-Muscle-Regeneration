########################################################
# Myonuclei Supplemental Figure
# Lauren Walter - September 2022
########################################################

# Load libraries ####
library(cowplot)
library(ggplot2)
library(Seurat)

# Source file ####
setwd("/FILE/PATH/TO/SEURAT_OBJECT")

myo.seurat <- readRDS(file = "SEURAT_OBJECT.rds", refhook = NULL)

# UMAP with general CT IDs ####
plot1 <- DimPlot(object = myo.seurat,
                 reduction = "umap",
                 group.by = "MyoSub_Gen_CT",
                 shuffle = TRUE,
                 label = FALSE, 
                 repel = TRUE,
                 pt.size = 0.1,
                 ncol = 1) +
  NoLegend() +
  labs(title = NULL, x = "UMAP harmony 1", y = "UMAP harmony 2") +
  scale_color_manual(limits = c("Progenitors_1", "Progenitors_2", "Progenitors_3", "Progenitors_4","Myonuclei (Type IIx/IIb)",
                                "Myonuclei (Type IIb)", "Myonuclei (Type IIx)", "Doublets_1", "Doublets_2"), 
                     values = c("#FA78FA", "#05755D", "#3897E0", "#463F66", "#C85A00", "#5772AD",
                                "#36BF96", "#FAA000", "#9F7EF2")) +
  ylim(-6, 6) +
  theme(plot.title = element_text(color = "black", size = 16),
        plot.subtitle = element_text(color = "black", size = 12, hjust = 0.5),
        axis.line = element_line(color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.title = element_text(color = "black", size = 12, hjust = 0.5, vjust = 1),
        axis.text = element_text(color = "black", size = 12))
  
# Feature plot
genes <- c("Myh1", "Myh4")

FP_list <- list()

for (i in 1:2){
  FP_list[[i]] <- FeaturePlot(object = myo.seurat,
                              features = genes[i],
                              cols = c("grey", "blue"),
                              pt.size =  0.1,
                              reduction = "umap",
                              slot = "data", 
                              label = FALSE,
                              repel = TRUE,
                              combine = TRUE
  ) + 
    labs(title = genes[i], x = "UMAP harmony 1", y = "UMAP harmony 2") +
    ylim(-6, 6) +
    theme(plot.title = element_text(color = "black", size = 16, face = "italic"),
          plot.subtitle = element_text(color = "black", size = 12, hjust = 0.5),
          axis.line = element_line(color = "black", size = 0.5),
          axis.ticks = element_line(color = "black", size = 0.5),
          axis.title = element_text(color = "black", size = 12, hjust = 0.5, vjust = 1),
          axis.text = element_text(color = "black", size = 12),
          legend.position = "right", 
          legend.text = element_text(color = "black", size = 12), 
          legend.key.size = unit(0.40, "cm"))
}

# Organize the UMAP and feature plots together
plot_grid(plot1, FP_list[[1]], FP_list[[2]], labels = c("(A)", "(B)", ""), label_size = 16, ncol = 3, align = "hv")

# Save plot as PDF
ggsave(filename = "Myonuclei_SuppFig_AB.pdf", plot = last_plot(), width = 15, height = 5, units = "in")

# Scatter plots of Myh1 and Myh4 expression in each myonuclei cell type ####
# Define myonuclei cell type IDs
CT <- c("Myonuclei (Type IIx)", "Myonuclei (Type IIb)", "Myonuclei (Type IIx/IIb)")

# Define cell type names
CT_Title <- c("IIx", "IIb", "IIx/IIb")

# Define box annotation for each cell type
Box_Ann <- list(geom_rect(aes(xmin = 0.20, xmax = 6.2, ymin = -0.18, ymax = 0.18), alpha = 0, color = "black", linetype = 2, size = 0.2),
                geom_rect(aes(xmin = -0.18, xmax = 0.18, ymin = 0.35, ymax = 6.9), alpha = 0, color = "black", linetype = 2, size = 0.2),
                geom_rect(aes(xmin = 0.05, xmax = 4.1, ymin = 0.04, ymax = 5.3), alpha = 0, color = "black", linetype = 2, size = 0.2))

# Define percent annotation for each cell type
Perc_Ann <- list(annotate(geom = "text", x = 7.0, y = 0, label = "52%", size = (6*0.36)), 
                 annotate(geom = "text", x = 0.08, y = 7.4, label = "66%", size = (6*0.36)),
                 annotate(geom = "text", x = 4.8, y = 5.0, label = "39%", size = (6*0.36)))

# Scatter plot
SP_list <- list()

Idents(myo.seurat) <- myo.seurat$MyoSub_Gen_CT

for (i in 1:3){
  SP_list[[i]] <- FeatureScatter(object = subset(myo.seurat, idents = CT[i], invert = FALSE),
                                 feature1 = "Myh1",
                                 feature2 = "Myh4",
                                 group.by = "MyoSub_Gen_CT",
                                 slot = "data",
                                 pt.size = 0.1) + 
    NoLegend() +
    labs(title = CT_Title[i]) +
    Box_Ann[[i]] +
    Perc_Ann[[i]] +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8), limits = c(-0.5, 8)) + 
    scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8), limits = c(-0.5, 8)) +
    scale_color_manual(limits = c("Myonuclei (Type IIx/IIb)", "Myonuclei (Type IIb)", "Myonuclei (Type IIx)"),
                       values = c("#C85A00", "#5772AD","#36BF96")) +
    theme(plot.title = element_text(color = "black", size = 6, face = "bold"),
          axis.text.x = element_text(color = "black", size = 6),
          axis.title = element_text(color = "black", size = 6, face = "italic"),
          axis.text.y = element_text(color = "black", size = 6),
          axis.line = element_line(color = "black", size = 0.25), 
          axis.ticks = element_line(color = "black", size = 0.25),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
}

# Plot all three scatter plots together
plot3 <- plot_grid(plotlist = SP_list, labels = c("(C)", "", ""), label_size = 8, ncol = 3, align = "hv")

# Stacked bar plot of Myh1 and Myh4 expression ####
# Define empty lists
Myh1 <- list()
Myh4 <- list()

# For each myonuclei cell type, extract the normalized counts for Myh1 and Myh4
Idents(myo.seurat) <- myo.seurat$MyoSub_Gen_CT

for (i in 1:3){
  Myh1[[i]] <- GetAssayData(object = subset(myo.seurat, idents = CT[i], invert = FALSE), slot = "data")["Myh1",]
  Myh4[[i]] <- GetAssayData(object = subset(myo.seurat, idents = CT[i], invert = FALSE), slot = "data")["Myh4",]
}

# For each myonuclei cell type, generate a dataframe with the normalized counts for Myh1 and Myh4 
df <- list()

for (i in 1:3){
  df[[i]] <- data.frame("Myh1" = Myh1[[i]],
                        "Myh4" = Myh4[[i]])
}

# For each myonuclei cell type, generate a dataframe with the number of cells that either individually express Myh1 or Myh4 or co-express Myh1 and Myh4
results <- list()

for (i in 1:3){
  results[[i]] <- data.frame("Myh4_pos" = c(0,0),
                             "Myh4_neg" = c(0,0))
  
  row.names(results[[i]]) <- c("Myh1_pos", "Myh1_neg")
}

for (i in 1:3){
  for (j in 1:nrow(df[[i]])){
    if (((df[[i]][j,1] == 0) & (df[[i]][j,2] == 0)) == TRUE) {
      results[[i]][2,2] <- (results[[i]][2,2]+1)
    }
    if (((df[[i]][j,1] == 0) & (df[[i]][j,2] > 0)) == TRUE) {
      results[[i]][2,1] <- (results[[i]][2,1]+1)
    }
    if (((df[[i]][j,1] > 0) & (df[[i]][j,2] > 0)) == TRUE) {
      results[[i]][1,1] <- (results[[i]][1,1]+1)
    }
    if (((df[[i]][j,1] > 0) & (df[[i]][j,2] == 0)) == TRUE) {
      results[[i]][1,2] <- (results[[i]][1,2]+1)
    }
  }
}

# Generate a dataframe combining the number of cells that either individually express Myh1 or Myh4 or co-express Myh1 and Myh4 across the myonuclei cell types
Myh1_Myh4_df <- data.frame("Cell_types" = c("Myonuclei (Type IIx)", "Myonuclei (Type IIx)", "Myonuclei (Type IIx)", "Myonuclei (Type IIx)", 
                                            "Myonuclei (Type IIb)", "Myonuclei (Type IIb)", "Myonuclei (Type IIb)", "Myonuclei (Type IIb)",
                                            "Myonuclei (Type IIx/IIb)", "Myonuclei (Type IIx/IIb)", "Myonuclei (Type IIx/IIb)", "Myonuclei (Type IIx/IIb)"),
                           "Exp" = c("Myh1+ Myh4+", "Myh1+ Myh4-", "Myh1- Myh4+", "Myh1- Myh4-",
                                     "Myh1+ Myh4+", "Myh1+ Myh4-", "Myh1- Myh4+", "Myh1- Myh4-",
                                     "Myh1+ Myh4+", "Myh1+ Myh4-", "Myh1- Myh4+", "Myh1- Myh4-"),
                           "Freq" = NA,
                           "Perc" = NA)

Myh1_Myh4_df[1,3] <- results[[1]][1,1]
Myh1_Myh4_df[2,3] <- results[[1]][1,2]
Myh1_Myh4_df[3,3] <- results[[1]][2,1]
Myh1_Myh4_df[4,3] <- results[[1]][2,2]
Myh1_Myh4_df[5,3] <- results[[2]][1,1]
Myh1_Myh4_df[6,3] <- results[[2]][1,2]
Myh1_Myh4_df[7,3] <- results[[2]][2,1]
Myh1_Myh4_df[8,3] <- results[[2]][2,2]
Myh1_Myh4_df[9,3] <- results[[3]][1,1]
Myh1_Myh4_df[10,3] <- results[[3]][1,2]
Myh1_Myh4_df[11,3] <- results[[3]][2,1]
Myh1_Myh4_df[12,3] <- results[[3]][2,2]

# Calculate the fraction of cells that either individually express Myh1 or Myh4 or co-express Myh1 and Myh4 by each myonuclei cell type
for (i in 1:4){
  Myh1_Myh4_df[i,4] <- ((Myh1_Myh4_df[i,3])/sum(Myh1_Myh4_df[1:4,3]))
}

for (i in 5:8){
  Myh1_Myh4_df[i,4] <- ((Myh1_Myh4_df[i,3])/sum(Myh1_Myh4_df[5:8,3]))
}

for (i in 9:12){
  Myh1_Myh4_df[i,4] <- ((Myh1_Myh4_df[i,3])/sum(Myh1_Myh4_df[9:12,3]))
}

# Staked bar plot
plot4 <- ggplot(data = Myh1_Myh4_df, 
                aes(x = Cell_types, 
                    y = Perc, 
                    fill = factor(Exp, levels = c("Myh1- Myh4-", "Myh1- Myh4+", "Myh1+ Myh4-", "Myh1+ Myh4+")))) +
  theme_minimal() +
  geom_bar(position = "fill", stat = "identity") + 
  geom_text(aes(label = paste0(round(Perc*100, digits = 0), "%")), size = 2.16, col = "white", position = position_stack(vjust = .5)) +
  scale_fill_manual(NULL,
                    limits = c("Myh1+ Myh4+", "Myh1+ Myh4-", "Myh1- Myh4+", "Myh1- Myh4-"),
                    values = c("#7a2048", "#408ec6", "#1e2761", "#808080")) + 
  scale_x_discrete(labels = c("Myonuclei (Type IIx)" = "IIx", "Myonuclei (Type IIb)" = "IIb", "Myonuclei (Type IIx/IIb)" = "IIx/IIb"),
                   limits = c("Myonuclei (Type IIx)", "Myonuclei (Type IIb)", "Myonuclei (Type IIx/IIb)")) +
  labs(title = NULL, x = NULL, y = "Fraction of cells") +              
  theme(axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 6, hjust = 1), 
        axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.key.size = unit(0.2, "cm"))

# Organize plots 3 and 4 together
plot_grid(plot3, plot4, labels = c("", "(D)"), label_size = 8, ncol = 2, rel_widths = c(2,1))

# Save plot as PDF
ggsave(filename = "Myonuclei_SuppFig_CD.pdf", plot = last_plot(), width = 7.5, height = 1.88, units = "in")

# Violin plots of select markers in each myonuclei cell type ####
VP_list <- list()

features <- c("Tnnc2", "Tnni2", "Mb", "Cox6a2", "Cox6c", "Atp5e", "Atp5g1", "Ckm", "Myh2", "Myh1", "Myh4", "Chrne", "Col22a1", "percent.mt")

for (i in 1:length(features)){
  VP_list[[i]] <- VlnPlot(object = myo.seurat,
                          features = features[i],
                          cols = NULL,
                          pt.size = 0,
                          idents = NULL,
                          sort = FALSE,
                          assay = "RNA",
                          group.by = "MyoSub_Gen_CT",
                          split.by = NULL,
                          adjust = 1,
                          y.max = NULL,
                          same.y.lims = FALSE,
                          log = FALSE,
                          ncol = 1,
                          slot = "data",
                          split.plot = FALSE,
                          combine = TRUE
  ) +
    NoLegend() +
    scale_fill_manual(limits = c("Myonuclei (Type IIx/IIb)", "Myonuclei (Type IIb)", "Myonuclei (Type IIx)"), 
                      values = c("#C85A00", "#5772AD", "#36BF96")) +
    scale_x_discrete(labels = c("Myonuclei (Type IIx)" = "IIx", "Myonuclei (Type IIb)" = "IIb", "Myonuclei (Type IIx/IIb)" = "IIx/IIb"), 
                     limits = c("Myonuclei (Type IIx)", "Myonuclei (Type IIb)", "Myonuclei (Type IIx/IIb)")) +
    labs(title = features[i], x = NULL) +
    theme(axis.line = element_line(color = "black", size = 0.25),
          axis.ticks = element_line(color = "black", size = 0.25),
          plot.title = element_text(color = "black", size = 6, face = "italic", hjust = 0),
          axis.text.x = element_text(color = "black", size = 6, angle = 0, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 6, hjust = 1),
          axis.title.y = element_text(color = "black", size = 6),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  
}

# Define line size of the violin outline
for (i in 1:14){
  VP_list[[i]]$layers[[1]]$aes_params$size = 0.25
}

# Plot all of the violin plots in a grid
plot_grid(VP_list[[1]], 
          VP_list[[2]] + labs(y = NULL), 
          VP_list[[3]] + labs(y = NULL), 
          VP_list[[4]] + labs(y = NULL), 
          VP_list[[5]] + labs(y = NULL),
          VP_list[[6]] + labs(y = NULL),
          VP_list[[7]] + labs(y = NULL),
          VP_list[[8]],
          VP_list[[9]] + labs(y = NULL),
          VP_list[[10]] + labs(y = NULL),
          VP_list[[11]] + labs(y = NULL),
          VP_list[[12]] + labs(y = NULL),
          VP_list[[13]] + labs(y = NULL),
          VP_list[[14]] + labs(title = "% Mito Reads", y = NULL) + theme(plot.title = element_text(face = "plain")),
          ncol = 7, nrow = 2)

# Save plot as PDF
ggsave(filename = "Myonuclei_SuppFig_EF.pdf", plot = last_plot(), width = 7.5, height = 2, units = "in")
