########################################################
# Cell Cycle Scoring
# Lauren Walter - June 2022
########################################################

# Load libraries ####
library(biomaRt) 
library(Seurat)
library(ggplot2)
library(ggpubr)
library(hash)
library(scales)
library(plyr)
library(dplyr)

# Source file ####
setwd("/FILE/PATH/TO/SEURAT_OBJECT")

seur.obj <- readRDS(file = "SEURAT_OBJECT.rds", refhook = NULL)

# Cell Cycle Scoring following Seurat's workflow ####
# https://satijalab.org/seurat/archive/v3.0/cell_cycle_vignette.html

# S phase and G2/M phase markers from Tirosh et al, 2015
S.humGenes <- cc.genes$s.genes
G2M.humGenes <- cc.genes$g2m.genes

# Convert human genes to mouse genes
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

s.genes <- convertHumanGeneList(S.humGenes)
g2m.genes <- convertHumanGeneList(G2M.humGenes)

# Assign cell cycle scores
seur.obj <- CellCycleScoring(seur.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Converting Cartesian coordinates to polar coordinates ####

df  <- data.frame("CT_ID" = colnames(seur.obj),
                  "S.Score" = seur.obj$S.Score, 
                  "G2M.Score" = seur.obj$G2M.Score, 
                  "Phase" = seur.obj$Phase, 
                  "r" = 0,
                  "theta" = 0,
                  "Cdk1" = GetAssayData(object = seur.obj, slot = "data")['Cdk1',],
                  "Cdk4" = GetAssayData(object = seur.obj, slot = "data")['Cdk4',],
                  "Age" = seur.obj$Age.Word,
                  "Time.Point" = seur.obj$Time.Point)

for (i in 1:nrow(df)){
  df[i,5] <- sqrt((df[i,2]^2) + (df[i,3]^2)) # calculate r
  df[i,6] <- atan(df[i,2]/df[i,3]) # calculate theta
}

polar_df <- mutate(df, theta = ifelse(df[,2] < 0, atan(df[,3] / df[,2]) + pi,
                                      ifelse(df[,3] < 0 , atan(df[,3] / df[,2]) + 2*pi, atan(df[,3] / df[,2])))) %>%
  arrange(theta)

ggplot(data = polar_df,
       aes(x = theta, 
           y = r,
           color = Phase)) +
  theme_grey() +
  geom_point(size = 0.1, alpha = 0.5) +
  scale_color_manual("Cell Cycle\nPhase", limits = c('G1', 'S', 'G2M'), values = c('#6CCAD9', '#369C5F', '#C93A95')) +
  labs(x = "Theta", y = "r") +
  scale_x_continuous(breaks = c(0, pi/2, pi, (3*pi)/2, 2*pi), labels = label_number(accuracy = 0.001)) +
  theme(text = element_text(color = "black", size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color = "black", size = 6, hjust = 0),
        plot.subtitle = element_text(color = "black", size = 6, hjust = 0),
        axis.title = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(color = "black", size = 0.25), 
        axis.ticks = element_line(color = "black", size = 0.25),
        legend.position = "none") 

# Define the quadrant each cell is in
polar_df$Quadrant <- 0

for (i in 1:nrow(polar_df)){
  if (((polar_df[i,6] >= 0) & (polar_df[i,6] < (pi/2))) == TRUE){
    polar_df[i,11] <- 1
  }
  if (((polar_df[i,6] >= (pi/2)) & (polar_df[i,6] < pi)) == TRUE){
    polar_df[i,11] <- 2
  }
  if (((polar_df[i,6] >= pi) & (polar_df[i,6] < ((3*pi)/2))) == TRUE){
    polar_df[i,11] <- 3
  }
  if (((polar_df[i,6] >= ((3*pi)/2)) & (polar_df[i,6] < (2*pi))) == TRUE){
    polar_df[i,11] <- 4
  }
}

# Rescale within each quadrant. This is the Normalized theta value.

quadrants <- split(polar_df, polar_df$Quadrant)

quad1 <- nrow(quadrants[[1]])
quad2 <- nrow(quadrants[[2]])
quad3 <- nrow(quadrants[[3]])
quad4 <- nrow(quadrants[[4]])

quadrants[[3]]$Rescale <- rescale(quadrants[[3]][,6], to = c(0, 0.25))
quadrants[[4]]$Rescale <- rescale(quadrants[[4]][,6], to = c(0.25, 0.5))
quadrants[[1]]$Rescale <- rescale(quadrants[[1]][,6], to = c(0.5, 0.75))
quadrants[[2]]$Rescale <- rescale(quadrants[[2]][,6], to = c(0.75, 1))

tmp <- rbind(quadrants[[1]], quadrants[[2]], quadrants[[3]], quadrants[[4]])

# Define non-G1 cells (Normalized theta value > 0.375)

tmp$G1_Status <- TRUE

for (i in 1:nrow(tmp)){
  if ((tmp[i,12] > 0.375) == TRUE) {
    tmp[i,13] <- FALSE
  }
}

# Add G1 status to Seurat object
df2 <- data.frame("CT_IDs" = colnames(seur.obj),
                  "G1_Status" = NA)

df2$CT_IDs <- gsub("-1", "", df2$CT_IDs)

tmp$CT_ID <- gsub("-1", "", tmp$CT_ID)

h <- hash( keys = tmp[,1], values = tmp[,13] )

for (i in 1:nrow(df2)){
  df2[i,2] <- h[[ df2[i,1] ]]
}

df2$CT_IDs <- paste0(df2$CT_IDs, "-1") 

seur.obj$G1_Status <- df2[,2]
