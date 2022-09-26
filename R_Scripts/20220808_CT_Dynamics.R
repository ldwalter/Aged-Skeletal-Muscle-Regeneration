########################################################
# Cell Type Dynamics: Plotting and Statistical Analysis
# Lauren Walter - August 2022
########################################################

# Load libraries ####
library(ggplot2)
library(ggpubr)
library(nlme)
library(plyr)
library(rlist)
library(scales)
library(Seurat)
library(stringr)

# Source file ####
setwd("/FILE/PATH/TO/SEURAT_OBJECT")

seur.obj <- readRDS(file = "SEURAT_OBJECT.rds", refhook = NULL)

# For each sample, calculate the percent of cells in each cell type ####
Idents(seur.obj) <- seur.obj$Specific_cell_types

df <- data.frame(table(seur.obj$Sample.ID, seur.obj$Specific_cell_types))

df_list <- split(df,df$Var2)

for (i in 1:29){
  colnames(df_list[[i]]) <- c("Sample.ID", "Cell.type", as.character(df_list[[i]][1,2]))
  
  df_list[[i]][,2] <- NULL
}

cell_count <- list.cbind(df_list)

cell_count <- cell_count[c(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58)]

colnames(cell_count) <- c("Sample.ID", as.character(unique(df$Var2)))

# Remove the erythrocytes from the data frame
cell_count$Erythrocytes <- NULL

# Count the number of cells in each sample
cell_count$Total <- 0

for (i in 1:nrow(cell_count)){
  cell_count[i,30] <- sum(cell_count[i, 2:29])
}

# Divide each value in a given row by the total cell count for that row
CT_dynamics <- data.frame(matrix(NA, nrow = 65, ncol = 28))

for (i in 1:65){
  for (j in 2:29){
    CT_dynamics[i,j-1] <- cell_count[i,j]/cell_count[i,30]
  }
}

# Make the column labels the cell type annotation
colnames(CT_dynamics) <- colnames(cell_count)[2:29]

# Add sample ID, age, and time point columns
CT_dynamics$Sample.ID <- cell_count$Sample.ID
CT_dynamics$Age <- c("Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", 
                     "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Ger", "Old", "Old", "Old", "Old", "Old", "Old", "Old", "Old", "Old", "Old", 
                     "Old", "Old", "Old", "Old", "Old", "Old", "Old", "Old", "Old", "Old", "Old", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", 
                     "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng", "Yng")
CT_dynamics$Time.Point <- c(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3.5, 3.5, 3.5, 3.5, 5, 5, 5, 5, 7, 7, 7, 7, 
                            0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3.5, 3.5, 3.5, 3.5, 5, 5, 5, 7, 7, 7, 
                            0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3.5, 3.5, 3.5, 3.5, 5,5, 5, 7, 7, 7)
CT_dynamics <- CT_dynamics[,c(29:31,1:28)]

# Remove punctuation and spaces from column names
colnames(CT_dynamics) <- str_replace_all(colnames(CT_dynamics), "[[:punct:]]", ".")
colnames(CT_dynamics) <- str_replace_all(colnames(CT_dynamics), " ", ".")
colnames(CT_dynamics) <- str_replace_all(colnames(CT_dynamics), "[-|=|\\+]", ".")

CT_names <- colnames(cell_count)[2:29]

write.csv(x = CT_dynamics, file = "CT_dynamics.csv", row.names = FALSE)
write.csv(x = CT_names, file = "CT_names.csv", row.names = FALSE)

# Statistical analysis using Non-linear modeling ####

# Plot theme
line.theme <- theme(text = element_text(color = "black", size = 6),
                    plot.title = element_text(color = "black", size = 6, hjust = 0),
                    plot.subtitle = element_text(color = "black", size = 6, hjust = 0),
                    axis.title = element_text(color = "black", size = 6),
                    axis.text.x = element_text(color = "black"), 
                    axis.text.y = element_text(color = "black"),
                    axis.line = element_line(color = "black", size = 0.25), 
                    axis.ticks = element_line(color = "black", size = 0.25),
                    legend.key.size = unit(0.25, "cm")) 

# Determine which non-linear equation best fits the data
# For every cell type, group data by age and calculate the coefficients for the quaternary, cubic, and quadratic equations

df_list <- list()
grp_list <- list()
quart_results <- list()
cubic_results <- list()
quad_results <- list()

for (i in 1:28){
  
  # Split the cell type dynamics data frame by cell type
  df_list[[i]] <- CT_dynamics[1:3]
  df_list[[i]][,4] <- CT_dynamics[i+3]
  colnames(df_list[[i]])[4] <- "Cell.Type"
  
  # For every cell type, group the data by age
  grp_list[[i]] <- groupedData(Cell.Type ~ Time.Point | Age, data = df_list[[i]])
  
  # For every cell type, calculate the coefficients for the non-linear equation that best fits the data
  quart_results[[i]] <- nlsList(Cell.Type ~ a*Time.Point^4 + b*Time.Point^3 + c*Time.Point^2 + d*Time.Point + e, data = grp_list[[i]],
                                start = list(a = 0, b = 0, c = 0, d = 0, e = 0))
  
  cubic_results[[i]] <- nlsList(Cell.Type ~ a*Time.Point^3 + b*Time.Point^2 + c*Time.Point + d, data = grp_list[[i]],
                                start = list(a = 0, b = 0, c = 0, d = 0))
  
  quad_results[[i]] <- nlsList(Cell.Type ~ a*Time.Point^2 + b*Time.Point + c, data = grp_list[[i]],
                               start = list(a = 0, b = 0, c = 0))
}

# The type of equation used for each cell type was selected based on the confidence interval and significance (p<0.05) of the leading coefficient
# No modeling equation went below the second degree

# View the coefficients (change number in double brackets to designate which cell type you want to see)
summary(quart_results[[1]])
summary(cubic_results[[1]])
summary(quad_results[[1]])

# View the confidence intervals for each coefficient (change number in double brackets to designate which cell type you want to see)  
plot(intervals(quart_results[[1]]))
plot(intervals(cubic_results[[1]]))
plot(intervals(quad_results[[1]]))

# If the leading coefficient was significantly different than zero and the confidence interval did not include zero for two out of the three age groups, it was concluded that the leading coefficient was needed  
# Otherwise the leading coefficient was not needed and the degree of the equation went down one
# The non-linear equation used for each cell type is detailed below

# Determine if the quadratic formula is different by age group
# Followed the method described here: https://stats.stackexchange.com/questions/26611/how-to-test-the-effect-of-a-grouping-variable-with-a-non-linear-model

quad_dyn <- subset(CT_dynamics, select = c("Sample.ID", "Age", "Time.Point", "B.cells", "Dendritic.cells..Cd209a..", "Dendritic.cells..Cd72..", 
                                           "Dendritic.cells..Xcr1..", "Endothelial.and.Myeloid.cells", "Endothelial.cells..Capillary.",
                                           "FAPs..Adipogenic.", "FAPs..Stem.", "Monocytes.Macrophages..Cxcl10..", "Myonuclei", 
                                           "Pericytes.and.Smooth.muscle.cells", "Schwann.and.Neural.Glial.cells", "T.cells..Cd4..", 
                                           "T.cells..Cycling..Cd3e..", "T.cells..Non.cycling..Cd3e..", "Tenocytes"))

quad_names <- data.frame("Cell_types" = c("B cells", "Dendritic cells (Cd209a+)", "Dendritic cells (Cd72+)", "Dendritic cells (Xcr1+)", 
                                          "Endothelial and Myeloid cells", "Endothelial cells (Capillary)", "FAPs (Adipogenic)", "FAPs (Stem)", 
                                          "Monocytes/Macrophages (Cxcl10+)", "Myonuclei", "Pericytes and Smooth muscle cells", 
                                          "Schwann and Neural/Glial cells", "T cells (Cd4+)", "T cells (Cycling; Cd3e+)", "T cells (Non-cycling; Cd3e+)",
                                          "Tenocytes"))

# Define null hypothesis: Age groups have the same coefficients
nls.null <- list()

for (i in 4:ncol(quad_dyn)){
  nls.null[[i-3]] <- nls(formula = quad_dyn[[i]] ~ a*Time.Point^2 + b*Time.Point + c,
                         data = quad_dyn,
                         start=list(a = 0, b = 0, c = 0))
  
}

# Convert the null hypothesis list into a data frame
nls.null.df <- data.frame("Cell_types" = quad_names$Cell_types,
                          "a" = NA, "b" = NA, "c" = NA)

for (i in 1:length(nls.null)){
  nls.null.df[i,2] <- coef(nls.null[[i]])[[1]]
  nls.null.df[i,3] <- coef(nls.null[[i]])[[2]]
  nls.null.df[i,4] <- coef(nls.null[[i]])[[3]]
}

# Convert the quad_dyn data frame into a list of data frames (each data frame is one cell type)
CTdyn_list <- list()

for (i in 4:ncol(quad_dyn)){
  CTdyn_list[[i-3]] <- data.frame(quad_dyn[1],
                                  quad_dyn[2],
                                  quad_dyn[3],
                                  quad_dyn[i])
}

# Change the column names
for (i in 1:length(CTdyn_list)){
  colnames(CTdyn_list[[i]])[4] <- "cell_type"
}

# Plots for every cell type with the null hypothesis
null_plots <- list()

for (i in 1:length(CTdyn_list)){
  null_plots[[i]] <- ggplot(data = CTdyn_list[[i]],
                            aes(x = as.numeric(Time.Point),
                                y = cell_type,
                                color = Age)) +
    theme_classic() +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(aes(fill = Age), alpha = 0.2, linetype = 0) +
    scale_color_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_fill_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_x_continuous(breaks = seq(0,7,1)) +
    line.theme +
    labs(title = quad_names[i,1],
         subtitle = "Null Hypothesis",
         x = "Days post-injury",
         y = "Fraction of total cells") 
}

# Combine multiple plots into one ggplot and add the null hypothesis equation
ggarrange(null_plots[[1]] + stat_function(fun=function(x) nls.null.df[1,2]*x^2 + nls.null.df[1,3]*x + nls.null.df[1,4], color="black"),
          null_plots[[2]] + stat_function(fun=function(x) nls.null.df[2,2]*x^2 + nls.null.df[2,3]*x + nls.null.df[2,4], color="black"),
          null_plots[[3]] + stat_function(fun=function(x) nls.null.df[3,2]*x^2 + nls.null.df[3,3]*x + nls.null.df[3,4], color="black"),
          null_plots[[4]] + stat_function(fun=function(x) nls.null.df[4,2]*x^2 + nls.null.df[4,3]*x + nls.null.df[4,4], color="black"),
          null_plots[[5]] + stat_function(fun=function(x) nls.null.df[5,2]*x^2 + nls.null.df[5,3]*x + nls.null.df[5,4], color="black"),
          null_plots[[6]] + stat_function(fun=function(x) nls.null.df[6,2]*x^2 + nls.null.df[6,3]*x + nls.null.df[6,4], color="black"),
          null_plots[[7]] + stat_function(fun=function(x) nls.null.df[7,2]*x^2 + nls.null.df[7,3]*x + nls.null.df[7,4], color="black"),
          null_plots[[8]] + stat_function(fun=function(x) nls.null.df[8,2]*x^2 + nls.null.df[8,3]*x + nls.null.df[8,4], color="black"),
          null_plots[[9]] + stat_function(fun=function(x) nls.null.df[9,2]*x^2 + nls.null.df[9,3]*x + nls.null.df[9,4], color="black"),
          null_plots[[10]] + stat_function(fun=function(x) nls.null.df[10,2]*x^2 + nls.null.df[10,3]*x + nls.null.df[10,4], color="black"),
          null_plots[[11]] + stat_function(fun=function(x) nls.null.df[11,2]*x^2 + nls.null.df[11,3]*x + nls.null.df[11,4], color="black"),
          null_plots[[12]] + stat_function(fun=function(x) nls.null.df[12,2]*x^2 + nls.null.df[12,3]*x + nls.null.df[12,4], color="black"),
          null_plots[[13]] + stat_function(fun=function(x) nls.null.df[13,2]*x^2 + nls.null.df[13,3]*x + nls.null.df[13,4], color="black"),
          null_plots[[14]] + stat_function(fun=function(x) nls.null.df[14,2]*x^2 + nls.null.df[14,3]*x + nls.null.df[14,4], color="black"),
          null_plots[[15]] + stat_function(fun=function(x) nls.null.df[15,2]*x^2 + nls.null.df[15,3]*x + nls.null.df[15,4], color="black"),
          null_plots[[16]] + stat_function(fun=function(x) nls.null.df[16,2]*x^2 + nls.null.df[16,3]*x + nls.null.df[16,4], color="black"),
          ncol = 6,
          nrow = 3,
          common.legend = TRUE,
          legend = "right")

#  Define alternative hypothesis: Age groups have different coefficients
nls.alt <- list()

for (i in 4:ncol(quad_dyn)){
  nls.alt[[i-3]] <- nls(formula = quad_dyn[[i]] ~ 
                          as.numeric(Age=="Ger")* aGer * Time.Point^2 + as.numeric(Age=="Ger")* bGer *Time.Point + as.numeric(Age=="Ger")* cGer +
                          as.numeric(Age=="Old")* aOld * Time.Point^2 + as.numeric(Age=="Old")* bOld *Time.Point + as.numeric(Age=="Old")* cOld +
                          as.numeric(Age=="Yng")* aYng * Time.Point^2 + as.numeric(Age=="Yng")* bYng *Time.Point + as.numeric(Age=="Yng")* cYng ,
                        data = quad_dyn, 
                        start = list(aGer = 0, bGer = 0, cGer = 0, 
                                     aOld = 0, bOld = 0, cOld = 0,  
                                     aYng = 0, bYng = 0, cYng = 0))
}

# Convert the alternative hypothesis list into a data frame
nls.alt.df <- data.frame("Cell_types" = quad_names$Cell_types,
                         "aGer" = NA, "bGer" = NA, "cGer" = NA,
                         "aOld" = NA, "bOld" = NA, "cOld" = NA,
                         "aYng" = NA, "bYng" = NA, "cYng" = NA)

for (i in 1:length(nls.alt)){
  nls.alt.df[i,2] <- coef(nls.alt[[i]])[[1]]
  nls.alt.df[i,3] <- coef(nls.alt[[i]])[[2]]
  nls.alt.df[i,4] <- coef(nls.alt[[i]])[[3]]
  nls.alt.df[i,5] <- coef(nls.alt[[i]])[[4]]
  nls.alt.df[i,6] <- coef(nls.alt[[i]])[[5]]
  nls.alt.df[i,7] <- coef(nls.alt[[i]])[[6]]
  nls.alt.df[i,8] <- coef(nls.alt[[i]])[[7]]
  nls.alt.df[i,9] <- coef(nls.alt[[i]])[[8]]
  nls.alt.df[i,10] <- coef(nls.alt[[i]])[[9]]
}

# Plots for every cell type with the alternative hypothesis
alt_plots <- list()

for (i in 1:length(CTdyn_list)){
  alt_plots[[i]] <- ggplot(data = CTdyn_list[[i]],
                           aes(x = as.numeric(Time.Point),
                               y = cell_type,
                               color = Age)) +
    theme_classic() +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(aes(fill = Age), alpha = 0.2, linetype = 0) +
    scale_color_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_fill_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_x_continuous(breaks = seq(0,7,1)) +
    line.theme +
    labs(title = quad_names[i,1],
         subtitle = "Alternative Hypothesis",
         x = "Days post-injury",
         y = "Fraction of total cells") 
}

# Combine multiple plots into one ggplot and add the alternative hypothesis equations
null_alt_plots <- list()

for (i in 1:16){
  null_alt_plots[[i]] <- ggarrange(null_plots[[i]] + stat_function(fun=function(x) nls.null.df[i,2]*x^2 + nls.null.df[i,3]*x + nls.null.df[i,4], color="black"),
                                   alt_plots[[i]] + stat_function(fun=function(x) nls.alt.df[i,2]*x^2 + nls.alt.df[i,3]*x + nls.alt.df[i,4], color="#4c026e") +
                                     stat_function(fun=function(x) nls.alt.df[i,5]*x^2 + nls.alt.df[i,6]*x + nls.alt.df[i,7], color="#4e94c7") +
                                     stat_function(fun=function(x) nls.alt.df[i,8]*x^2 + nls.alt.df[i,9]*x + nls.alt.df[i,10], color="#f08205"),
                                   ncol = 2,
                                   nrow = 1,
                                   common.legend = TRUE,
                                   legend = "right")
}

ggarrange(null_alt_plots[[2]], null_alt_plots[[3]], null_alt_plots[[4]], null_alt_plots[[9]], ncol = 2, nrow = 2)
ggarrange(null_alt_plots[[1]], null_alt_plots[[13]], null_alt_plots[[14]], null_alt_plots[[15]],  ncol = 2, nrow = 2)
ggarrange(null_alt_plots[[7]], null_alt_plots[[8]], null_alt_plots[[12]], null_alt_plots[[16]],  ncol = 2, nrow = 2)
ggarrange(null_alt_plots[[5]], null_alt_plots[[6]], null_alt_plots[[10]], null_alt_plots[[11]],  ncol = 2, nrow = 2)

# Likelihood Ratio Test to see if the alternative hypothesis fits the data significantly better than the null hypothesis
anova_summary <- list()

for (i in 1:16){
  anova_summary[[i]] <- anova(nls.null[[i]], nls.alt[[i]])
}

# Organize results from ANOVA into a data frame
quadratic_stats <- data.frame("Cell_types" = quad_names$Cell_types,
                              "pval" = NA,
                              "FDR_pval" = NA,
                              "Significant" = NA)

for (i in 1:16){
  quadratic_stats[i,2] <- anova_summary[[i]][6][2,1]
}

for (i in 1:16){
  quadratic_stats[i,3] <- p.adjust(p = quadratic_stats[i,2], method = "fdr", n = 28)
}

for (i in 1:16){
  if ((quadratic_stats[i,3] <= 0.05) == TRUE){
    quadratic_stats[i,4] <- "T"
  } else
    quadratic_stats[i,4] <- "F"
}

# Determine if the cubic formula is different by age group

cub_dyn <- subset(CT_dynamics, select = c("Sample.ID", "Age", "Time.Point", "Dendritic.cells..Fscn1..", "Endothelial.cells..Artery.", 
                                          "Endothelial.cells..Vein.", "FAPs..Pro.remodeling.", "M1.M2.Macrophages..Mrc1..", "M2.Macrophages..Cx3cr1..",
                                          "MuSCs.and.progenitors", "NK.cells"))

cub_names <- data.frame("Cell_types" = c("Dendritic cells (Fscn1+)", "Endothelial cells (Artery)", "Endothelial cells (Vein)", "FAPs (Pro-remodeling)",
                                         "M1/M2 Macrophages (Mrc1+)", "M2 Macrophages (Cx3cr1+)", "MuSCs and progenitors", "NK cells"))

# Define null hypothesis: Age groups have the same coefficients
nls.null <- list()

for (i in 4:ncol(cub_dyn)){
  nls.null[[i-3]] <- nls(formula = cub_dyn[[i]] ~ a*Time.Point^3 + b*Time.Point^2 + c*Time.Point + d,
                         data = cub_dyn,
                         start=list(a = 0, b = 0, c = 0, d = 0))
  
}

# Convert the null hypothesis list into a data frame
nls.null.df <- data.frame("Cell_types" = cub_names$Cell_types,
                          "a" = NA, "b" = NA, "c" = NA, "d" = NA)

for (i in 1:length(nls.null)){
  nls.null.df[i,2] <- coef(nls.null[[i]])[[1]]
  nls.null.df[i,3] <- coef(nls.null[[i]])[[2]]
  nls.null.df[i,4] <- coef(nls.null[[i]])[[3]]
  nls.null.df[i,5] <- coef(nls.null[[i]])[[4]]
}

# Convert the cub_dyn data frame into a list of data frames (each data frame is one cell type)
CTdyn_list <- list()

for (i in 4:ncol(cub_dyn)){
  CTdyn_list[[i-3]] <- data.frame(cub_dyn[1],
                                  cub_dyn[2],
                                  cub_dyn[3],
                                  cub_dyn[i])
}

# Change the column names
for (i in 1:length(CTdyn_list)){
  colnames(CTdyn_list[[i]])[4] <- "cell_type"
}

# Plots for every cell type with the null hypothesis
null_plots <- list()

for (i in 1:length(CTdyn_list)){
  null_plots[[i]] <- ggplot(data = CTdyn_list[[i]],
                            aes(x = as.numeric(Time.Point),
                                y = cell_type,
                                color = Age)) +
    theme_classic() +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(aes(fill = Age), alpha = 0.2, linetype = 0) +
    scale_color_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_fill_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_x_continuous(breaks = seq(0,7,1)) +
    line.theme +
    labs(title = cub_names[i,1],
         subtitle = "Null Hypothesis",
         x = "Days post-injury",
         y = "Fraction of total cells") 
}

# Combine multiple plots into one ggplot and add the null hypothesis equation
ggarrange(null_plots[[1]] + stat_function(fun=function(x) nls.null.df[1,2]*x^3 + nls.null.df[1,3]*x^2 + nls.null.df[1,4]*x + nls.null.df[1,5], color="black"),
          null_plots[[2]] + stat_function(fun=function(x) nls.null.df[2,2]*x^3 + nls.null.df[2,3]*x^2 + nls.null.df[2,4]*x + nls.null.df[2,5], color="black"),
          null_plots[[3]] + stat_function(fun=function(x) nls.null.df[3,2]*x^3 + nls.null.df[3,3]*x^2 + nls.null.df[3,4]*x + nls.null.df[3,5], color="black"),
          null_plots[[4]] + stat_function(fun=function(x) nls.null.df[4,2]*x^3 + nls.null.df[4,3]*x^2 + nls.null.df[4,4]*x + nls.null.df[4,5], color="black"),
          null_plots[[5]] + stat_function(fun=function(x) nls.null.df[5,2]*x^3 + nls.null.df[5,3]*x^2 + nls.null.df[5,4]*x + nls.null.df[5,5], color="black"),
          null_plots[[6]] + stat_function(fun=function(x) nls.null.df[6,2]*x^3 + nls.null.df[6,3]*x^2 + nls.null.df[6,4]*x + nls.null.df[6,5], color="black"),
          ncol = 3,
          nrow = 2,
          common.legend = TRUE,
          legend = "right")

#  Define alternative hypothesis: Age groups have different coefficients
nls.alt <- list()

for (i in 4:ncol(cub_dyn)){
  nls.alt[[i-3]] <- nls(formula = cub_dyn[[i]] ~ 
                          as.numeric(Age=="Ger")* aGer * Time.Point^3 + as.numeric(Age=="Ger")* bGer *Time.Point^2 + as.numeric(Age=="Ger")* cGer *Time.Point + as.numeric(Age=="Ger")* dGer +
                          as.numeric(Age=="Old")* aOld * Time.Point^3 + as.numeric(Age=="Old")* bOld *Time.Point^2 + as.numeric(Age=="Old")* cOld *Time.Point + as.numeric(Age=="Old")* dOld +
                          as.numeric(Age=="Yng")* aYng * Time.Point^3 + as.numeric(Age=="Yng")* bYng *Time.Point^2 + as.numeric(Age=="Yng")* cYng *Time.Point + as.numeric(Age=="Yng")* dYng,
                        data = cub_dyn, 
                        start = list(aGer = 0, bGer = 0, cGer = 0, dGer = 0, 
                                     aOld = 0, bOld = 0, cOld = 0, dOld = 0,   
                                     aYng = 0, bYng = 0, cYng = 0, dYng = 0)) 
}

# Convert the alternative hypothesis list into a data frame
nls.alt.df <- data.frame("Cell_types" = cub_names$Cell_types,
                         "aGer" = NA, "bGer" = NA, "cGer" = NA, "dGer" = NA,
                         "aOld" = NA, "bOld" = NA, "cOld" = NA, "dOld" = NA,
                         "aYng" = NA, "bYng" = NA, "cYng" = NA, "dYng" = NA)

for (i in 1:length(nls.alt)){
  nls.alt.df[i,2] <- coef(nls.alt[[i]])[[1]]
  nls.alt.df[i,3] <- coef(nls.alt[[i]])[[2]]
  nls.alt.df[i,4] <- coef(nls.alt[[i]])[[3]]
  nls.alt.df[i,5] <- coef(nls.alt[[i]])[[4]]
  nls.alt.df[i,6] <- coef(nls.alt[[i]])[[5]]
  nls.alt.df[i,7] <- coef(nls.alt[[i]])[[6]]
  nls.alt.df[i,8] <- coef(nls.alt[[i]])[[7]]
  nls.alt.df[i,9] <- coef(nls.alt[[i]])[[8]]
  nls.alt.df[i,10] <- coef(nls.alt[[i]])[[9]]
  nls.alt.df[i,11] <- coef(nls.alt[[i]])[[10]]
  nls.alt.df[i,12] <- coef(nls.alt[[i]])[[11]]
  nls.alt.df[i,13] <- coef(nls.alt[[i]])[[12]]
}

# Plots for every cell type with the alternative hypothesis
alt_plots <- list()

for (i in 1:length(CTdyn_list)){
  alt_plots[[i]] <- ggplot(data = CTdyn_list[[i]],
                           aes(x = as.numeric(Time.Point),
                               y = cell_type,
                               color = Age)) +
    theme_classic() +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(aes(fill = Age), alpha = 0.2, linetype = 0) +
    scale_color_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_fill_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_x_continuous(breaks = seq(0,7,1)) +
    line.theme +
    labs(title = cub_names[i,1],
         subtitle = "Alternative Hypothesis",
         x = "Days post-injury",
         y = "Fraction of total cells") 
}

# Combine multiple plots into one ggplot and add the alternative hypothesis equations
null_alt_plots <- list()

for (i in 1:8){
  null_alt_plots[[i]] <- ggarrange(null_plots[[i]] + stat_function(fun=function(x) nls.null.df[i,2]*x^3 + nls.null.df[i,3]*x^2 + nls.null.df[i,4]*x + nls.null.df[i,5],  color="black"),
                                   alt_plots[[i]] + stat_function(fun=function(x) nls.alt.df[i,2]*x^3 + nls.alt.df[i,3]*x^2 + nls.alt.df[i,4]*x + nls.alt.df[i,5], color="#4c026e") +
                                     stat_function(fun=function(x) nls.alt.df[i,6]*x^3 + nls.alt.df[i,7]*x^2 + nls.alt.df[i,8]*x + nls.alt.df[i,9], color="#4e94c7") +
                                     stat_function(fun=function(x) nls.alt.df[i,10]*x^3 + nls.alt.df[i,11]*x^2 + nls.alt.df[i,12]*x + nls.alt.df[i,13], color="#f08205"),
                                   ncol = 2,
                                   nrow = 1,
                                   common.legend = TRUE,
                                   legend = "right")
}

ggarrange(null_alt_plots[[1]], null_alt_plots[[5]], null_alt_plots[[6]], null_alt_plots[[8]], ncol = 2, nrow = 2)
ggarrange(null_alt_plots[[2]], null_alt_plots[[3]], null_alt_plots[[4]], null_alt_plots[[7]], ncol = 2, nrow = 2)

# Likelihood Ratio Test to see if the alternative hypothesis fits the data significantly better than the null hypothesis
anova_summary <- list()

for (i in 1:8){
  anova_summary[[i]] <- anova(nls.null[[i]], nls.alt[[i]])
}

# Organize results from ANOVA into a data frame
cubic_stats <- data.frame("Cell_types" = cub_names$Cell_types,
                          "pval" = NA,
                          "FDR_pval" = NA,
                          "Significant" = NA)

for (i in 1:8){
  cubic_stats[i,2] <- anova_summary[[i]][6][2,1]
}

for (i in 1:8){
  cubic_stats[i,3] <- p.adjust(p = cubic_stats[i,2], method = "fdr", n = 28)
}

for (i in 1:8){
  if ((cubic_stats[i,3] <= 0.05) == TRUE){
    cubic_stats[i,4] <- "T"
  } else
    cubic_stats[i,4] <- "F"
}

# Determine if the quaternary formula is different by age group

quat_dyn <- subset(CT_dynamics, select = c("Sample.ID", "Age", "Time.Point", "M1.Macrophages..Ccr2..", "Monocytes..Cycling..Cdk1..",
                                           "Monocytes.Macrophages..Patrolling..Ctsa..", "Neutrophils"))

quat_names <- data.frame("Cell_types" = c("M1 Macrophages (Ccr2+)", "Monocytes (Cycling; Cdk1+)", 
                                          "Monocytes/Macrophages (Patrolling; Ctsa+)", "Neutrophils"))

# Define null hypothesis: Age groups have the same coefficients
nls.null <- list()

for (i in 4:ncol(quat_dyn)){
  nls.null[[i-3]] <- nls(formula = quat_dyn[[i]] ~ a*Time.Point^4 + b*Time.Point^3 + c*Time.Point^2 + d*Time.Point + e,
                         data = quat_dyn,
                         start=list(a = 0, b = 0, c = 0, d = 0, e = 0))
  
}

# Convert the null hypothesis list into a data frame
nls.null.df <- data.frame("Cell_types" = quat_names$Cell_types,
                          "a" = NA, "b" = NA, "c" = NA, "d" = NA, "e" = NA)

for (i in 1:length(nls.null)){
  nls.null.df[i,2] <- coef(nls.null[[i]])[[1]]
  nls.null.df[i,3] <- coef(nls.null[[i]])[[2]]
  nls.null.df[i,4] <- coef(nls.null[[i]])[[3]]
  nls.null.df[i,5] <- coef(nls.null[[i]])[[4]]
  nls.null.df[i,6] <- coef(nls.null[[i]])[[5]]
}

# Convert the quat_dyn data frame into a list of data frames  (each data frame is one cell type)
CTdyn_list <- list()

for (i in 4:ncol(quat_dyn)){
  CTdyn_list[[i-3]] <- data.frame(quat_dyn[1],
                                  quat_dyn[2],
                                  quat_dyn[3],
                                  quat_dyn[i])
}

# Change the column names
for (i in 1:length(CTdyn_list)){
  colnames(CTdyn_list[[i]])[4] <- "cell_type"
}

# Plots for every cell type with the null hypothesis
null_plots <- list()

for (i in 1:length(CTdyn_list)){
  null_plots[[i]] <- ggplot(data = CTdyn_list[[i]],
                            aes(x = as.numeric(Time.Point),
                                y = cell_type,
                                color = Age)) +
    theme_classic() +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(aes(fill = Age), alpha = 0.2, linetype = 0) +
    scale_color_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_fill_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_x_continuous(breaks = seq(0,7,1)) +
    line.theme +
    labs(title = quat_names[i,1],
         subtitle = "Null Hypothesis",
         x = "Days post-injury",
         y = "Fraction of total cells") 
}

# Combine multiple plots into one ggplot and add the null hypothesis equation
ggarrange(null_plots[[1]] + stat_function(fun=function(x) nls.null.df[1,2]*x^4 + nls.null.df[1,3]*x^3 + nls.null.df[1,4]*x^2 + nls.null.df[1,5]*x + nls.null.df[1,6], color="black"),
          null_plots[[2]] + stat_function(fun=function(x) nls.null.df[2,2]*x^4 + nls.null.df[2,3]*x^3 + nls.null.df[2,4]*x^2 + nls.null.df[2,5]*x + nls.null.df[2,6], color="black"),
          null_plots[[3]] + stat_function(fun=function(x) nls.null.df[3,2]*x^4 + nls.null.df[3,3]*x^3 + nls.null.df[3,4]*x^2 + nls.null.df[3,5]*x + nls.null.df[3,6], color="black"),
          null_plots[[4]] + stat_function(fun=function(x) nls.null.df[4,2]*x^4 + nls.null.df[4,3]*x^3 + nls.null.df[4,4]*x^2 + nls.null.df[4,5]*x + nls.null.df[4,6], color="black"),
          ncol = 2,
          nrow = 2,
          common.legend = TRUE,
          legend = "right")

#  Define alternative hypothesis: Age groups have different coefficients
nls.alt <- list()

for (i in 4:ncol(quat_dyn)){
  nls.alt[[i-3]] <- nls(formula = quat_dyn[[i]] ~ 
                          as.numeric(Age=="Ger")* aGer * Time.Point^4 + as.numeric(Age=="Ger")* bGer *Time.Point^3 + as.numeric(Age=="Ger")* cGer *Time.Point^2 + as.numeric(Age=="Ger")* dGer *Time.Point + as.numeric(Age=="Ger")* eGer + 
                          as.numeric(Age=="Old")* aOld * Time.Point^4 + as.numeric(Age=="Old")* bOld *Time.Point^3 + as.numeric(Age=="Old")* cOld *Time.Point^2 + as.numeric(Age=="Old")* dOld *Time.Point + as.numeric(Age=="Old")* eOld +
                          as.numeric(Age=="Yng")* aYng * Time.Point^4 + as.numeric(Age=="Yng")* bYng *Time.Point^3 + as.numeric(Age=="Yng")* cYng *Time.Point^2 + as.numeric(Age=="Yng")* dYng *Time.Point + as.numeric(Age=="Yng")* eYng, 
                        data = quat_dyn, 
                        start = list(aGer = 0, bGer = 0, cGer = 0, dGer = 0, eGer = 0, 
                                     aOld = 0, bOld = 0, cOld = 0, dOld = 0, eOld = 0,   
                                     aYng = 0, bYng = 0, cYng = 0, dYng = 0, eYng = 0)) 
}

# Convert the alternative hypothesis list into a data frame
nls.alt.df <- data.frame("Cell_types" = quat_names$Cell_types,
                         "aGer" = NA, "bGer" = NA, "cGer" = NA, "dGer" = NA, "eGer" = NA,
                         "aOld" = NA, "bOld" = NA, "cOld" = NA, "dOld" = NA, "eOld" = NA,
                         "aYng" = NA, "bYng" = NA, "cYng" = NA, "dYng" = NA, "eYng" = NA)

for (i in 1:length(nls.alt)){
  nls.alt.df[i,2] <- coef(nls.alt[[i]])[[1]]
  nls.alt.df[i,3] <- coef(nls.alt[[i]])[[2]]
  nls.alt.df[i,4] <- coef(nls.alt[[i]])[[3]]
  nls.alt.df[i,5] <- coef(nls.alt[[i]])[[4]]
  nls.alt.df[i,6] <- coef(nls.alt[[i]])[[5]]
  nls.alt.df[i,7] <- coef(nls.alt[[i]])[[6]]
  nls.alt.df[i,8] <- coef(nls.alt[[i]])[[7]]
  nls.alt.df[i,9] <- coef(nls.alt[[i]])[[8]]
  nls.alt.df[i,10] <- coef(nls.alt[[i]])[[9]]
  nls.alt.df[i,11] <- coef(nls.alt[[i]])[[10]]
  nls.alt.df[i,12] <- coef(nls.alt[[i]])[[11]]
  nls.alt.df[i,13] <- coef(nls.alt[[i]])[[12]]
  nls.alt.df[i,14] <- coef(nls.alt[[i]])[[13]]
  nls.alt.df[i,15] <- coef(nls.alt[[i]])[[14]]
  nls.alt.df[i,16] <- coef(nls.alt[[i]])[[15]]
}

# Plots for every cell type with the alternative hypothesis
alt_plots <- list()

for (i in 1:length(CTdyn_list)){
  alt_plots[[i]] <- ggplot(data = CTdyn_list[[i]],
                           aes(x = as.numeric(Time.Point),
                               y = cell_type,
                               color = Age)) +
    theme_classic() +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(aes(fill = Age), alpha = 0.2, linetype = 0) +
    scale_color_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_fill_manual("Age", limits = c("Yng", "Old", "Ger"), values = c("#f08205", "#4e94c7", "#4c026e")) +
    scale_x_continuous(breaks = seq(0,7,1)) +
    line.theme +
    labs(title = quat_names[i,1],
         subtitle = "Alternative Hypothesis",
         x = "Days post-injury",
         y = "Fraction of total cells") 
}

# Combine multiple plots into one ggplot and add the alternative hypothesis equations
null_alt_plots <- list()

for (i in 1:4){
  null_alt_plots[[i]] <- ggarrange(null_plots[[i]] + stat_function(fun=function(x) nls.null.df[i,2]*x^4 + nls.null.df[i,3]*x^3 + nls.null.df[i,4]*x^2 + nls.null.df[i,5]*x + nls.null.df[i,6], color= "black"),
                                   alt_plots[[i]] + stat_function(fun=function(x) nls.alt.df[i,2]*x^4 + nls.alt.df[i,3]*x^3 + nls.alt.df[i,4]*x^2 + nls.alt.df[i,5]*x + nls.alt.df[i,6], color= "#4c026e") +
                                     stat_function(fun=function(x) nls.alt.df[i,7]*x^4 + nls.alt.df[i,8]*x^3 + nls.alt.df[i,9]*x^2 + nls.alt.df[i,10]*x + nls.alt.df[i,11], color="#4e94c7") +
                                     stat_function(fun=function(x) nls.alt.df[i,12]*x^4 + nls.alt.df[i,13]*x^3 + nls.alt.df[i,14]*x^2 + nls.alt.df[i,15]*x + nls.alt.df[i,16], color= "#f08205"),
                                   ncol = 2,
                                   nrow = 1,
                                   common.legend = TRUE,
                                   legend = "right")
}

ggarrange(null_alt_plots[[1]], null_alt_plots[[2]], null_alt_plots[[3]], null_alt_plots[[4]], ncol = 2, nrow = 2)

# Likelihood Ratio Test to see if the null hypothesis fits the data significantly better than the alternative hypothesis
anova_summary <- list()

for (i in 1:4){
  anova_summary[[i]] <- anova(nls.null[[i]], nls.alt[[i]])
}

# Organize results from ANOVA into a data frame
quaternary_stats <- data.frame("Cell_types" = quat_names$Cell_types,
                               "pval" = NA,
                               "FDR_pval" = NA,
                               "Significant" = NA)

for (i in 1:4){
  quaternary_stats[i,2] <- anova_summary[[i]][6][2,1]
}

for (i in 1:4){
  quaternary_stats[i,3] <- p.adjust(p = quaternary_stats[i,2], method = "fdr", n = 28)
}

for (i in 1:4){
  if ((quaternary_stats[i,3] <= 0.05) == TRUE){
    quaternary_stats[i,4] <- "T"
  } else
    quaternary_stats[i,4] <- "F"
}

# Combine statistics data frames for all three methods
stats_sum <- rbind(quadratic_stats, cubic_stats, quaternary_stats)
stats_sum$Cell_types <- as.character(stats_sum$Cell_types)
stats_sum <- stats_sum[order(stats_sum$Cell_types),]

# Calculate the mean and standard deviation for each group ####

# Helper function 
# http://sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func,
                    varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# Group the data by age and time point
table_list <- list()

for (i in 4:ncol(CT_dynamics)){
  table_list[[i-3]] <- data_summary(CT_dynamics,
                                    varname = colnames(CT_dynamics[i]),
                                    groupnames = c("Age", "Time.Point"))
}

# Calculate the ymin and ymax values for each time point in each age group
for (i in 1:length(table_list)){
  table_list[[i]][5] <- table_list[[i]][3]-table_list[[i]][4]  # Calculate ymin
  table_list[[i]][6] <- table_list[[i]][3]+table_list[[i]][4]  # Calculate ymax
  
  # Change the column names
  colnames(table_list[[i]])[3] <- "cell_type"
  colnames(table_list[[i]])[5] <- "y_min"
  colnames(table_list[[i]])[6] <- "y_max"
}

# Add column with age and replicate ID
CT_dynamics$Age_Rep <- c("Ger_A", "Ger_B", "Ger_C", "Ger_D", "Ger_A", "Ger_B", "Ger_C", "Ger_D", "Ger_A", "Ger_B", "Ger_C", "Ger_D", "Ger_A", 
                         "Ger_B", "Ger_C", "Ger_D", "Ger_A", "Ger_B", "Ger_C", "Ger_D", "Ger_A", "Ger_B", "Ger_C", "Ger_D", "Old_A", "Old_B", 
                         "Old_C", "Old_D", "Old_A", "Old_B", "Old_C", "Old_D", "Old_A", "Old_B", "Old_C", "Old_A", "Old_B", "Old_C", "Old_D", 
                         "Old_A", "Old_B", "Old_C", "Old_A", "Old_B", "Old_C", "Yng_A", "Yng_B", "Yng_C", "Yng_A", "Yng_B", "Yng_C", "Yng_D", 
                         "Yng_A", "Yng_B", "Yng_C", "Yng_A", "Yng_B", "Yng_C", "Yng_D", "Yng_A", "Yng_B", "Yng_C", "Yng_A", "Yng_B", "Yng_C")

CT_dynamics <- CT_dynamics[,c(1,32,2:31)]

# Group the data by age and replicate ID
table_list2 <- list()

for (i in 5:ncol(CT_dynamics)){
  table_list2[[i-4]] <- data_summary(CT_dynamics,
                                     varname = colnames(CT_dynamics[i]),
                                     groupnames = c("Age_Rep", "Age", "Time.Point"))
}

# Change the column names
for (i in 1:length(table_list2)){
  colnames(table_list2[[i]])[4] <- "cell_type"
}

# Line plots for each cell type ####
plot_list <- list()

for (i in 1:length(table_list)){
  plot_list[[i]] <- ggplot(data = table_list[[i]],
                           aes(x = as.numeric(Time.Point),
                               y = cell_type,
                               group = Age)) +
    theme_classic() +
    geom_ribbon(aes(ymin = y_min, ymax = y_max, fill = Age), alpha = 0.2) +
    geom_line(aes(color = Age), size = 0.5) +
    geom_point(data = table_list2[[i]], 
               aes(x = as.numeric(Time.Point),
                   y = cell_type,
                   group = Age_Rep,
                   color = Age), 
               size = 0.5, alpha = 0.5) +
    scale_color_manual("Age", labels = c("Yng" = "Young", "Old" = "Old", "Ger" = "Geriatric"),
                       limits = c("Yng", "Old", "Ger"), 
                       values = c("#f08205", "#4e94c7", "#4c026e")) +   
    scale_fill_manual("Age", labels = c("Yng" = "Young", "Old" = "Old", "Ger" = "Geriatric"),
                      limits = c("Yng", "Old", "Ger"), 
                      values = c("#f08205", "#4e94c7", "#4c026e")) +    
    scale_x_continuous(breaks = seq(0,7,1)) +
    line.theme +
    labs(title = CT_names[i], 
         subtitle = paste0("FDR adj. p-value = ", scientific(stats_sum[i,3], digits = 3)),
         x = "Days post-injury", 
         y = "Fraction of total cells")
}

# Plot the cell types that have a statistically significant difference in their cell type dynamics
ggarrange(plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]], 
          plot_list[[18]] + labs(title = "Monocytes/Macrophages \n(Patrolling; Ctsa+)"), 
          plot_list[[22]], plot_list[[27]], plot_list[[1]], plot_list[[19]],  
          common.legend = T, legend = "right", ncol = 3, nrow =3, align = "hv")

# Plot the cell types that do not have a statistically significant difference in their cell type dynamics
ggarrange(plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], 
          plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]], 
          plot_list[[12]], plot_list[[17]], plot_list[[20]], plot_list[[21]], plot_list[[23]], 
          plot_list[[24]], plot_list[[25]], plot_list[[26]], plot_list[[28]],
          common.legend = T, legend = "right", ncol = 4, nrow = 5, align = "hv")

# Plot summarizing the p-values for all cell types ####

ggplot(data = stats_sum,
       aes(x = -log10(FDR_pval),
           y = reorder(Cell_types, -FDR_pval),
           width = 0.75)) +
  theme_classic() +
  geom_bar(stat = "identity", color = "black", fill = "black", size = 0.25) +
  geom_vline(xintercept = -log10(0.05), color = "red") +
  labs(x = "-log10(FDR adj. p-value)", y = NULL) +
  theme(axis.ticks = element_line(size = 0.25), 
        axis.line = element_line(size = 0.25),
        axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(color = "black", size = 8))

# Stacked bar plots of the fraction of cells per cell type ####

# Plot theme
SBP_theme <- theme(text = element_text(color = "black", size = 6),
                   plot.title = element_text(color = "black", size = 6, hjust = 0),
                   plot.subtitle = element_text(color = "black", size = 6, hjust = 0),
                   axis.title = element_text(color = "black", size = 6),
                   axis.text.x = element_text(color = "black"), 
                   axis.text.y = element_text(color = "black"),
                   axis.line = element_line(color = "black", size = 0.25), 
                   axis.ticks = element_line(color = "black", size = 0.25),
                   legend.key.size = unit(0.15, "cm")) 

# Count the number of cells in each age group, time point, and specific cell type ID
df <- data.frame(table(seur.obj$Age.Word, seur.obj$Time.Point, seur.obj$Specific_cell_types))

# Remove cells with the "Erythrocytes" ID
df2 <- df %>%
  filter(!row_number() %in% c(163:180))

# Assign a broader ID to all immune cells
df2$Immune_Broad <- as.character(df2$Var3)

CT_immune <- c("B cells", "Dendritic cells (Cd209a+)", "Dendritic cells (Cd72+)", "Dendritic cells (Fscn1+)", "Dendritic cells (Xcr1+)", "Endothelial and Myeloid cells", "M1 Macrophages (Ccr2+)",
               "M1/M2 Macrophages (Mrc1+)", "M2 Macrophages (Cx3cr1+)", "Monocytes (Cycling; Cdk1+)", "Monocytes/Macrophages (Cxcl10+)", 
               "Monocytes/Macrophages (Patrolling; Ctsa+)", "Neutrophils", "NK cells", "T cells (Cd4+)", "T cells (Cycling; Cd3e+)", "T cells (Non-cycling; Cd3e+)")

for (i in 1:504){
  for (j in 1:17){
    if ((df2[i,5] == CT_immune[j]) == TRUE){
      df2[i,5] <- "Immune cells"
    }
  }
}

# Assign a broader ID to all endothelial cells
df2$Endo_Broad <- df2$Immune_Broad

CT_endo <- c("Endothelial cells (Artery)", "Endothelial cells (Capillary)", "Endothelial cells (Vein)")

for (i in 1:504){
  for (j in 1:3){
    if ((df2[i,6] == CT_endo[j]) == TRUE){
      df2[i,6] <- "Endothelial cells"
    }
  }
}

# Assign a broader ID to all FAPs
df2$FAPs_Broad <- df2$Endo_Broad

CT_FAPs <- c("FAPs (Adipogenic)", "FAPs (Pro-remodeling)", "FAPs (Stem)")

for (i in 1:504){
  for (j in 1:3){
    if ((df2[i,7] == CT_FAPs[j]) == TRUE){
      df2[i,7] <- "FAPs"
    }
  }
}

# Plot the fraction of cells for each cell type with broad immune cell, endothelial cell, and FAPs IDs
df3 <- df2 %>%
  group_by(Var2, FAPs_Broad) %>%
  summarise(Count = sum(Freq))

ggplot(data = df3, 
       aes(x = as.character(Var2),
           y = Count,
           fill = FAPs_Broad)) +
  theme_classic() +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days post-injury", y = "% of total cells", fill = "Broad CT IDs") +
  scale_fill_manual(limits = c("Endothelial cells", "FAPs", "Immune cells", "MuSCs and progenitors", "Myonuclei",  "Pericytes and Smooth muscle cells", "Schwann and Neural/Glial cells",  "Tenocytes"), 
                    values = c("#2E3B65", "#8BA690", "#B9DDFF", "#FA78FA", "#C85A00", "#826E00", "#5078FA", "#793EA8")) +
  SBP_theme

ggsave(filename = "Broad_CT_Dynamics_SBP.pdf", plot = last_plot(), width = 4, height = 2, units = "in")

# Plot the fraction of cells for each cell type with broad immune cell IDs
df4 <- df2 %>%
  group_by(Var2, Immune_Broad) %>%
  summarise(Count = sum(Freq))

ggplot(data = df4, 
       aes(x = as.character(Var2),
           y = Count,
           fill = Immune_Broad)) +
  theme_classic() +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days post-injury", y = "% of total cells", fill = "Broad CT IDs") +
  scale_fill_manual(limits = c("Endothelial cells (Artery)", "Endothelial cells (Capillary)", "Endothelial cells (Vein)", 
                               "FAPs (Adipogenic)", "FAPs (Pro-remodeling)", "FAPs (Stem)", "Immune cells",
                               "MuSCs and progenitors",  "Myonuclei", "Pericytes and Smooth muscle cells", "Schwann and Neural/Glial cells",  
                               "Tenocytes"), 
                    values = c("#940A1D", "#FAA000", "#447173", "#5AB40A", "#362354", "#065B66", "#B9DDFF", "#FA78FA", 
                               "#C85A00", "#826E00", "#5078FA", "#793EA8")) +
  SBP_theme

ggsave(filename = "ImmuneBroad_CT_Dynamics_SBP.pdf", plot = last_plot(), width = 4, height = 2, units = "in")

# Plot the fraction of cells for ONLY the immune cells
df5 <- df2 %>%
  group_by(Var2, Var3) %>%
  summarise(Count = sum(Freq))

df6 <- subset(df5, !(Var3 %in% c("Endothelial cells (Artery)", "Endothelial cells (Capillary)", "Endothelial cells (Vein)", 
                                 "FAPs (Adipogenic)", "FAPs (Pro-remodeling)", "FAPs (Stem)",  "MuSCs and progenitors",  "Myonuclei", 
                                 "Pericytes and Smooth muscle cells", "Schwann and Neural/Glial cells", "Tenocytes")))

ggplot(data = df6, 
       aes(x = as.character(Var2),
           y = Count,
           fill = Var3)) +
  theme_classic() +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days post-injury", y = "% of immune cells", fill = "Immune CT IDs") +
  scale_fill_manual(limits = c("B cells", "Dendritic cells (Cd209a+)", "Dendritic cells (Cd72+)", "Dendritic cells (Fscn1+)", "Dendritic cells (Xcr1+)", "Endothelial and Myeloid cells",
                               "M1 Macrophages (Ccr2+)", "M1/M2 Macrophages (Mrc1+)", "M2 Macrophages (Cx3cr1+)", 
                               "Monocytes (Cycling; Cdk1+)", "Monocytes/Macrophages (Cxcl10+)", "Monocytes/Macrophages (Patrolling; Ctsa+)", "Neutrophils", "NK cells",  
                               "T cells (Cd4+)", "T cells (Cycling; Cd3e+)", "T cells (Non-cycling; Cd3e+)"), 
                    values = c("#F20A53", "#4D1F82", "#8F4ECC", "#3D78E0", "#D6475F", "#000000", "#0AB4A6", "#5C946A",  "#82FAA0", "#278A4B", 
                               "#B84B85", "#354852", "#709EB5", "#076594", "#D15115", "#871E8F",
                               "#83CCBD")) +
  SBP_theme

ggsave(filename = "ImmuneOnly_CT_Dynamics_SBP.pdf", plot = last_plot(), width = 4, height = 2, units = "in")

# Assign a broader ID to all dendritic cells
df6$Immune_Broad <- as.character(df6$Var3)

CT_Dendritic <- c("Dendritic cells (Cd209a+)", "Dendritic cells (Cd72+)", "Dendritic cells (Fscn1+)", "Dendritic cells (Xcr1+)")

for (i in 1:102){
  for (j in 1:4){
    if ((df6[i,4] == CT_Dendritic[j]) == TRUE){
      df6[i,4] <- "Dendritic cells"
    }
  }
}

# Assign a broader ID to all T cells
CT_TCells <- c("T cells (Cd4+)", "T cells (Cycling; Cd3e+)", "T cells (Non-cycling; Cd3e+)")

for (i in 1:102){
  for (j in 1:3){
    if ((df6[i,4] == CT_TCells[j]) == TRUE){
      df6[i,4] <- "T cells"
    }
  }
}

# Plot the fraction of cells for ONLY the immune cells with broad dendritic and T cell IDs
df7 <- df6 %>%
  group_by(Var2, Immune_Broad) %>%
  summarise(Count = sum(Count))

ggplot(data = df7, 
       aes(x = as.character(Var2),
           y = Count,
           fill = Immune_Broad)) +
  theme_classic() +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days post-injury", y = "% of immune cells", fill = "Immune CT IDs") +
  scale_fill_manual(limits = c("B cells", "Dendritic cells", "Endothelial and Myeloid cells",
                               "M1 Macrophages (Ccr2+)", "M1/M2 Macrophages (Mrc1+)", "M2 Macrophages (Cx3cr1+)", 
                               "Monocytes (Cycling; Cdk1+)", "Monocytes/Macrophages (Cxcl10+)", "Monocytes/Macrophages (Patrolling; Ctsa+)", "Neutrophils", "NK cells",  
                               "T cells"), 
                    values = c("#F20A53", "#9C6DA5", "#000000", "#0AB4A6", "#5C946A",  "#82FAA0", "#278A4B", 
                               "#B84B85", "#354852", "#709EB5", "#076594", "#AFE4DE")) +
  SBP_theme

ggsave(filename = "ImmuneOnly_CT_Dynamics_SBP2.pdf", plot = last_plot(), width = 4, height = 2, units = "in")
