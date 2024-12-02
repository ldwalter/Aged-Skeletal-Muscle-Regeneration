########################################################
# Senescence Scoring
# Lauren Walter - September 2022
########################################################

# Load packages ####
library(escape)
library(Seurat)

# Source file ####
setwd("/FILE/PATH/TO/SEURAT_OBJECT")

# Seurat object
seur.obj <- readRDS(file = "SEURAT_OBJECT.rds", refhook = NULL)

# Two-way Senescence Score ####

# Define gene set of interest. As an example, this is the FBR gene set.
FBR <- read.csv(file = "Supplemental_Table_2_FBR.csv", header = TRUE)

# Specifically for the FBR gene list, remove genes that have a FDR higher than 0.05
FBR <- subset(FBR, FDR <= 0.05)

# Z-Score all genes within a given cell type cluster
# Define the active identity
Idents(seur.obj) <- seur.obj$Specific_cell_types

# Make a dataframe with the cell type names
CT_IDs <- data.frame(unique(seur.obj$Specific_cell_types))

# Define the genes to z-score
all.genes <- rownames(seur.obj)

# Subset the full Seurat object by cell type and z-score all genes within the subset
seur_list <- list()

for (i in 1:29){
  seur_list[[i]] <- subset(seur.obj, idents = CT_IDs[i,1], invert = FALSE)
  seur_list[[i]] <- ScaleData(object = seur_list[[i]], assay = "RNA", features = all.genes)
}

# Two-way Senescence Score function was developed by Chris Cherry in the Elisseeff Lab
# https://github.com/Chris-Cherry/sctools/blob/master/R/signed_set_scoring.R

#' Scores cell expression of genes in set taking into account sign of fold change
#'
#' Takes in a gene set and uses normalized scaled gene expression data to calculate a score describing extent to which
#' cells are expressing genes in the set. Seurat object must contain scale_data assay. Gene set must be organized
#' such that "genes" are genes and "FC" is log(foldchange). All genes will be assumed to be significant.
#'
#' @param geneset   A list of genes and their fold changes. 
#' @param ser       A Seurat object to be scored by the gene set. Must contain scaled data
#' @param from_gene Original gene type
#' @param to_gene   Gene type for conversion
#' @param scaled    Boolean to determine whether to scale score by the number of genes in the set
#' @importFrom      stats na.omit
#' @return          Outputs a vector of score values named as cell names
#' @export

signed_set_scoring <- function(ser, geneset, from_gene = "MGI", to_gene = "MGI", scaled = TRUE){
  
  # Gene conversion if required
  if(to_gene != from_gene){
    gene_conv = convert_genes(geneset$genes, from_gene, to_gene)
    mid = match(geneset$genes, gene_conv[,1])
    geneset$genes = gene_conv[mid, 2]
    geneset = na.omit(geneset)
  }
  
  # Calculate summed z-scores, taking into account directionality of fold change
  pos = geneset$genes[which(geneset$FC>0)]
  neg = geneset$genes[which(geneset$FC<0)]
  scale_dat = GetAssayData(object = ser, slot = "scale.data")
  pos_ind = match(pos, rownames(scale_dat))
  if (!identical(which(is.na(pos_ind)), integer(0))){
    pos_ind = pos_ind[-which(is.na(pos_ind))]}
  pos_gene_subset = scale_dat[pos_ind,]
  neg_ind = match(neg, rownames(scale_dat))
  if (!identical(which(is.na(neg_ind)), integer(0))){
    neg_ind = neg_ind[-which(is.na(neg_ind))]}
  neg_gene_subset = scale_dat[neg_ind,]
  if (is.null(dim(neg_gene_subset))){
    neg_score = neg_gene_subset * -1
  } else {
    neg_score = colSums(neg_gene_subset*-1)
  }
  if (is.null(dim(pos_gene_subset))){
    pos_score = pos_gene_subset
  } else {
    pos_score = colSums(pos_gene_subset)
  }
  scores = neg_score + pos_score
  names(scores) = colnames(ser)
  
  # Scale by number of genes
  if (scaled) {scores = scores/length(geneset$genes)}
  
  return(scores)
}

# Define the gene list to be used
genes <- data.frame("genes" = FBR$Gene,  # Mouse gene IDs
                    "FC" = FBR$logFC)    # Fold-change of the gene

# Calculate the Two-way Senescence Score for each cell in every subset (based on cell type ID)
score_list <- list()

for (i in 1:29){
  score_list[[i]] <- signed_set_scoring(ser = seur_list[[i]],
                                       geneset = genes,
                                       from_gene = "MGI",
                                       to_gene = "MGI",
                                       scaled = TRUE)
}

# Add the Two-Way Senescence Score to the Seurat object
for (i in 1:29){
   seur_list[[i]] <- AddMetaData(object = seur_list[[i]],
                                metadata = score_list[[i]],
                                col.name="TwoWay_FBR_SenScore")
}

seur.obj$TwoWay_FBR_SenScore <- 0

for (i in 1:29){
  seur.obj$TwoWay_FBR_SenScore[colnames(seur_list[[i]])] <- paste(seur_list[[i]]$TwoWay_FBR_SenScore)
}

seur.obj$TwoWay_FBR_SenScore <- as.numeric(seur.obj$TwoWay_FBR_SenScore)

# One-way Senescence Score ####

# Define gene set of interest. As an example, this is the SenMayo gene set.
gene.sets <- list(Senescence = c("Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3", "Bmp2", "Bmp6", "C3", "Ccl1", "Ccl2", "Ccl20", "Ccl24", "Ccl26", 
                                 "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cd55", "Cd9", "Csf1", "Csf2", "Csf2rb", "Cst10", "Ctnnb1", "Ctsb", "Cxcl1", "Cxcl10", 
                                 "Cxcl12", "Cxcl16", "Cxcl2", "Cxcl3", "Cxcr2", "Dkk1", "Edn1", "Egf", "Egfr", "Ereg", "Esm1", "Ets2", "Fas", "Fgf1", "Fgf2", 
                                 "Fgf7", "Gdf15", "Gem", "Gmfg", "Hgf", "Hmgb1", "Icam1", "Icam5", "Igf1", "Igfbp1", "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", 
                                 "Igfbp6", "Igfbp7", "Il10", "Il13", "Il15", "Il18", "Il1a", "Il1b", "Il2", "Il6", "Il6st", "Il7", "Inha", "Iqgap2", "Itga2", 
                                 "Itpka", "Jun", "Kitl", "Lcp1", "Mif", "Mmp13", "Mmp10", "Mmp12", "Mmp13", "Mmp14", "Mmp2", "Mmp3", "Mmp9", "Nap1l4", "Nrg1", 
                                 "Pappa", "Pecam1", "Pgf", "Pigf", "Plat", "Plau", "Plaur", "Ptbp1", "Ptger2", "Ptges", "Rps6ka5", "Scamp4", "Selplg", "Sema3f", 
                                 "Serpinb3a", "Serpine1", "Serpine2", "Spp1", "Spx", "Timp2", "Tnf", "Tnfrsf11b", "Tnfrsf1a", "Tnfrsf1b", "Tubgcp2", "Vegfa", 
                                 "Vegfc", "Vgf", "Wnt16", "Wnt2"))

# Calculate enrichment score
ES <- enrichIt(obj = seur.obj, 
               gene.sets = gene.sets, 
               method = "ssGSEA",
               groups = 1000, 
               cores = 15, 
               min.size = 5, 
               ssGSEA.norm = FALSE)

# Add the One-way Senescence Score to the Seurat object
seur.obj <- AddMetaData(object = seur.obj, metadata = ES, col.name = "SenMayo")
