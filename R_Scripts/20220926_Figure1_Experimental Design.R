########################################################
# Figure 1 - Experimental Design
# Lauren Walter - September 2022
########################################################

# Load libraries ####
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Seurat)

# Source file ####
setwd("/FILE/PATH/TO/SEURAT_OBJECT")

seur.obj <- readRDS(file = "SEURAT_OBJECT.rds", refhook = NULL)

# Pie chart of the fraction of cells within each age group ####
# Generate dataframe with the number of cells in each age group
df <- data.frame(table(seur.obj$Age.Word))

# Calculate the fraction of cells in each age group
df$Perc <- df$Freq/sum(df$Freq)

# Define position of labels
df <- df %>%
  arrange(desc(Var1)) %>%
  mutate(lab.ypos = cumsum(Perc) - 0.5*Perc)

# Plot pie chart 
plot1 <- ggplot(data = df,
                aes(x = "",
                    y = Perc,
                    fill = Var1)) +
  theme_void() +
  geom_bar(position = 'fill', stat = 'identity', width = 1, color = "white") +
  geom_text(aes(y = lab.ypos, label = paste0(round(Perc*100, digits = 1), "%")), size = 2, col = "white", fontface = "bold") +
  scale_fill_manual("Age", limits = c("Young", "Old", "Geriatric"), values = c('#F08205', '#4E94C7', '#4C026E')) +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  coord_polar("y", start=0)

# Stacked bar plot of the fraction of cells in each time point by age group ####
# Generate dataframe with the number of cells in each age group and time point
df2 <- data.frame(table(seur.obj$Age.Word, 
                        seur.obj$Time.Point))

# Stacked bar plot
plot2 <- ggplot(data = df2,
                aes(x = Var1,
                    y = Freq,
                    fill = factor(Var2, levels = c("0", "1", "2", "3.5", "5", "7")))) +
  theme_minimal() +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual("DPI", limits = c('0', '1', '2', '3.5', '5', '7'), 
                    values = c('#700324', '#e3104f', '#f584a6', '#024769', '#097cb5', '#60b6e0')) +
  labs(title = NULL, x = NULL, y = 'Fraction of cells by age') +
  theme(axis.text.x = element_text(color = "black", size = 6, angle = 0, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 6, hjust = 1),
        axis.title.y = element_text(color = "black", size = 6),
        legend.position = "none") +
  scale_x_discrete(limits = c("Young", "Old", "Geriatric")) +
  coord_fixed(3)

# UMAP colored by specific cell type IDs ####
plot3 <- DimPlot(object = seur.obj,
                 reduction = 'umap',
                 group.by = "Specific_cell_types",
                 shuffle = TRUE,
                 label = FALSE,
                 repel = TRUE,
                 pt.size = 0.1,
                 ncol = 1) +
  NoLegend() +
  labs(title = NULL, x = "UMAP harmony 1", y = "UMAP harmony 2") +
  scale_color_manual(limits = c('B cells', 'Dendritic cells (Cd209a+)', 'Dendritic cells (Cd72+)', 'Dendritic cells (Fscn1+)', 'Dendritic cells (Xcr1+)', 
                                'Endothelial and Myeloid cells', 'Endothelial cells (Artery)', 'Endothelial cells (Capillary)', 
                                'Endothelial cells (Vein)', 'Erythrocytes', 'FAPs (Adipogenic)', 'FAPs (Pro-remodeling)', 'FAPs (Stem)', 
                                'M1 Macrophages (Ccr2+)', 'M1/M2 Macrophages (Mrc1+)', 'M2 Macrophages (Cx3cr1+)', 'Monocytes (Cycling; Cdk1+)',  
                                'Monocytes/Macrophages (Cxcl10+)', 'Monocytes/Macrophages (Patrolling; Ctsa+)', 'MuSCs and progenitors',  
                                'Myonuclei', 'Neutrophils', 'NK cells', 'Pericytes and Smooth muscle cells', 'Schwann and Neural/Glial cells',  
                                'T cells (Cd4+)', 'T cells (Cycling; Cd3e+)', 'T cells (Non-cycling; Cd3e+)', 'Tenocytes'), 
                     values = c('#F20A53', '#4D1F82', '#8F4ECC', '#3D78E0', '#D6475F', '#000000', '#940A1D', '#FAA000', '#447173', '#3C0AAA', 
                                '#5AB40A', '#362354', '#065B66', '#0AB4A6', '#5C946A', '#82FAA0', '#278A4B', '#B84B85', '#354852', '#FA78FA', 
                                '#C85A00', '#709EB5', '#076594', '#826E00', '#5078FA', '#D15115', '#871E8F', '#83CCBD', '#793EA8')) +
  coord_fixed() +
  theme(plot.title = element_text(color = "black", size = 8),
        plot.subtitle = element_text(color = "black", size = 6, hjust = 0.5),
        axis.line = element_line(color = "black", size = 0.25),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.title = element_text(color = "black", size = 6, hjust = 0.5, vjust = 1),
        axis.text = element_text(color = "black", size = 6))

# UMAP colored by age group ####
plot_list <- list()
Age <- c("Young", "Old", "Geriatric")
color <- list(scale_color_manual(limits = c("Young", "Old", "Geriatric"), values = c('#F08205', '#D3D3D3', '#D3D3D3')),
              scale_color_manual(limits = c("Young", "Old", "Geriatric"), values = c('#D3D3D3', '#4E94C7', '#D3D3D3')),
              scale_color_manual(limits = c("Young", "Old", "Geriatric"), values = c('#D3D3D3', '#D3D3D3', '#4C026E')))

for (i in 1:3) {
  plot_list[[i]] <- DimPlot(object = seur.obj,
                            reduction = 'umap',
                            group.by = "Age.Word",
                            shuffle = TRUE,
                            label = FALSE,
                            repel = TRUE,
                            pt.size = 0.1,
                            ncol = 1) +  
    NoLegend() +
    labs(title = Age[i], x = "UMAP harmony 1", y = "UMAP harmony 2") +
    color[i] +
    coord_fixed() +
    theme(plot.title = element_text(color = "black", size = 8),
          plot.subtitle = element_text(color = "black", size = 6, hjust = 0.5),
          axis.line = element_line(color = "black", size = 0.25),
          axis.ticks = element_line(color = "black", size = 0.25),
          axis.title = element_text(color = "black", size = 6, hjust = 0.5, vjust = 1),
          axis.text = element_text(color = "black", size = 6))
}

plot4 <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 1, align = "hv")

# Organize the plots together
left <- plot_grid(plot1, plot2, labels = c('(D)', '(E)'), label_size = 8, align = 'hv', nrow = 2)
top <- plot_grid(left, plot3, labels = c('', '(F)'), label_size = 8, ncol = 2, rel_widths = c(1,2))

plot_grid(top, plot4, labels = c('', '(G)'), label_size = 8, ncol = 1, rel_heights = c(2,1))

# Save plot as PDF
ggsave(filename = "Figure1_Experimental_Design.pdf", plot = last_plot(), width = 7, height = 7, units = "in")
