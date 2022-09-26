########################################################
# Myogenic Doublets Supplemental Figure
# Lauren Walter - September 2022
########################################################

# Load libraries ####
library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)
library(Seurat)

# Source file ####
setwd("/FILE/PATH/TO/SEURAT_OBJECT")

myo.seurat <- readRDS(file = "SEURAT_OBJECT.rds", refhook = NULL)

# Scatter plot of select markers with density curves ####
# From the Seurat object, extract the normalized counts of select markers for each myogenic cell type
Pax7 <- list()
Pecam1 <- list()

Idents(myo.seurat) <- myo.seurat$MyoSub_Gen_CT

CT <- data.frame("CT" = unique(myo.seurat$MyoSub_Gen_CT))

for (i in 1:9){
  Pax7[[i]] <- GetAssayData(object = subset(myo.seurat, idents = CT[i,1], invert = FALSE), slot = "data")["Pax7",]
  Pecam1[[i]] <- GetAssayData(object = subset(myo.seurat, idents = CT[i,1], invert = FALSE), slot = "data")["Pecam1",]
}

# Generate a dataframe with the normalized counts of select markers for each myogenic cell type by cell barcode
df <- list()

for (i in 1:9){
  df[[i]] <- data.frame("Pax7" = Pax7[[i]],
                        "Pecam1" = Pecam1[[i]],
                        "CT" = as.character(i))
}

# Scatter plot of normalized counts of select markers for each myogenic cell type 
plot_list <- list()

for (i in 1:9){
  plot_list[[i]] <- ggplot(data = df[[i]],
                           aes(x = Pax7,
                               y = Pecam1,
                               color = CT)) +
    theme_minimal() +
    geom_point(size = 0.1, alpha = 0.25) +
    labs(title = CT[i,],
         subtitle = paste0(nrow(df[[i]]), " cells"),
         x = NULL, y = NULL) +
    scale_color_manual(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                       values = c("black", "black", "black", "black", "black", "black", "black", "black", "black")) +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(-0.5, 4)) +
    scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(-0.5, 4)) +
    theme(plot.title = element_text(color = "black", size = 6, face = "bold", vjust = 0),
          plot.subtitle = element_text(color = "black", size = 6, hjust = 0),
          legend.position = "none",
          axis.line = element_line(color = "black", size = 0.25),
          axis.ticks = element_line(color = "black", size = 0.25),
          axis.text.x = element_text(color = "black", size = 6),
          axis.text.y = element_text(color = "black", size = 6)) 
}

# Add custom theme settings to select plots
plot_list[[1]] <- plot_list[[1]] + theme(axis.text.x = element_blank())
plot_list[[6]] <- plot_list[[6]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list[[4]] <- plot_list[[4]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list[[7]] <- plot_list[[7]] + theme(axis.text.x = element_blank())
plot_list[[5]] <- plot_list[[5]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list[[2]] <- plot_list[[2]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list[[8]] <- plot_list[[8]] + theme(axis.text.y = element_blank())
plot_list[[9]] <- plot_list[[9]] + theme(axis.text.y = element_blank())

# Add density curve to each plot
for (i in 1:9){
  plot_list[[i]] <- ggMarginal(plot_list[[i]], type = "density", groupFill = TRUE)
}

# Organize the plots together
plot1 <- annotate_figure(p = ggarrange(plot_list[[1]], plot_list[[6]], plot_list[[4]], plot_list[[7]], plot_list[[5]], plot_list[[2]],
                                       plot_list[[3]], plot_list[[8]], plot_list[[9]], ncol = 3, nrow = 3),
                         bottom = text_grob("Pax7", size = 8, hjust = 0.5, face = "italic"),
                         left = text_grob("Pecam1", size = 8, face = "italic", rot = 90))

# From the Seurat object, extract the normalized counts of select markers for each myogenic cell type
Acta1 <- list()
C1qa <- list()

Idents(myo.seurat) <- myo.seurat$MyoSub_Gen_CT

CT <- data.frame("CT" = unique(myo.seurat$MyoSub_Gen_CT))

for (i in 1:9){
  Acta1[[i]] <- GetAssayData(object = subset(myo.seurat, idents = CT[i,1], invert = FALSE), slot = "data")["Acta1",]
  C1qa[[i]] <- GetAssayData(object = subset(myo.seurat, idents = CT[i,1], invert = FALSE), slot = "data")["C1qa",]
}

# Generate a dataframe with the normalized counts of select markers for each myogenic cell type by cell barcode
df2 <- list()

for (i in 1:9){
  df2[[i]] <- data.frame("Acta1" = Acta1[[i]],
                         "C1qa" = C1qa[[i]],
                         "CT" = as.character(i))
}

# Scatter plot of normalized counts of select markers for each myogenic cell type
plot_list2 <- list()

for (i in 1:9){
  plot_list2[[i]] <- ggplot(data = df2[[i]],
                            aes(x = Acta1,
                                y = C1qa,
                                color = CT)) +
    theme_minimal() +
    geom_point(size = 0.1, alpha = 0.25) +
    labs(title = CT[i,],
         subtitle = paste0(nrow(df[[i]]), " cells"),
         x = NULL, y = NULL) +
    scale_color_manual(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                       values = c("black", "black", "black", "black", "black", "black", "black", "black", "black")) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(-0.5, 8)) +
    scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(-0.5, 8)) +
    theme(plot.title = element_text(color = "black", size = 6, face = "bold", vjust = 0),
          plot.subtitle = element_text(color = "black", size = 6, hjust = 0),
          legend.position = "none",
          axis.line = element_line(color = "black", size = 0.25),
          axis.ticks = element_line(color = "black", size = 0.25),
          axis.text.x = element_text(color = "black", size = 6),
          axis.text.y = element_text(color = "black", size = 6)) 
}

# Add custom theme settings to select plots
plot_list2[[1]] <- plot_list2[[1]] + theme(axis.text.x = element_blank())
plot_list2[[6]] <- plot_list2[[6]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list2[[4]] <- plot_list2[[4]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list2[[7]] <- plot_list2[[7]] + theme(axis.text.x = element_blank())
plot_list2[[5]] <- plot_list2[[5]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list2[[2]] <- plot_list2[[2]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list2[[8]] <- plot_list2[[8]] + theme(axis.text.y = element_blank())
plot_list2[[9]] <- plot_list2[[9]] + theme(axis.text.y = element_blank())

# Add density curve to each plot
for (i in 1:9){
  plot_list2[[i]] <- ggMarginal(plot_list2[[i]], type = "density", groupFill = TRUE)
}

# Organize the plots together
plot2 <- annotate_figure(p = ggarrange(plot_list2[[1]], plot_list2[[6]], plot_list2[[4]], plot_list2[[7]], plot_list2[[5]], plot_list2[[2]],
                                       plot_list2[[3]], plot_list2[[8]], plot_list2[[9]], ncol = 3, nrow = 3),
                         bottom = text_grob("Acta1", size = 8, hjust = 0.5, face = "italic"),
                         left = text_grob("C1qa", size = 8, face = "italic", rot = 90))

# Organize plot1 and plot2 together
plot_grid(plot1, plot2, labels = c("(A)", "(B)"), label_size = 8, ncol = 2, nrow = 1, align = "hv")

# Save plot as PDF
ggsave(filename = "Myogenic_Doublets_SuppFig.pdf", plot = last_plot(), width = 7.5, height = 5, units = "in")
