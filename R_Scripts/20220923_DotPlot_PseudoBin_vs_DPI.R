########################################################
# Dotplot with pseudotime bins on the x-axis and the DPI on the y-axis
# Lauren Walter - September 2022
########################################################

# Load libraries ####
library(dplyr)
library(ggplot2)
library(Seurat)

# Source file ####
setwd("/FILE/PATH/TO/SEURAT_OBJECT")

myo.seurat <- readRDS(file = "SEURAT_OBJECT.rds", refhook = NULL)

# Dotplot with pseudotime bins on the x-axis and the DPI on the y-axis ####
# Generate dataframe with the number of cells at each age, time point, and pseudotime bin
df <- data.frame(table(myo.seurat$Age.Word, 
                       myo.seurat$Time.Point, 
                       myo.seurat$Pseudo.bins2)) 

# Define column names                 
colnames(df) <- c("Age", "Time.Point", "Pseudo.bins2", "Freq")

# Add up the total number of cells by time point
df <- df %>%
  group_by(Time.Point) %>%
  mutate(Sum = sum(Freq))

# Calculate the fraction of cells per age, time point, and pseudotime bin by the total number of cells at the corresponding time point
df$Fraction <- df$Freq/df$Sum

# Generate dataframe to define how the dots should be colored
df2 <- data.frame("Age" = myo.seurat$Age.Word,
                  "Time.Point" = myo.seurat$Time.Point, 
                  "Pseudo.bins2" = myo.seurat$Pseudo.bins2,
                  "Color_by" =  GetAssayData(object = myo.seurat)["Isg15",]) # Can change this to another gene or metadata column

# Calculate the average color by value per age, time point, and pseudotime bin
df2 <- df2 %>%
  group_by(Age, Time.Point, Pseudo.bins2) %>%
  mutate(Mean = mean(Color_by))

# Remove column from dataframe
df2$Color_by <- NULL

# Remove all rows that are not unique (we only need one row for every age, time point, and pseudotime bin because we calculated the average value)
df2 <- unique(df2)

# Add a column to the original dataframe of the variable to color by
df$Color_by <- NA

for (i in 1:nrow(df)){
  for (j in 1:nrow(df2)){
      if (((df[i,1] == df2[j,1]) & (df[i,2] == df2[j,2]) & (df[i,3] == df2[j,3])) == TRUE)
        df[i,7] <- df2[j,4]
  }
}

# Remove rows that are incomplete (ie. contain NAs)
data_complete <- df[complete.cases(df), ]

# Dotplot with pseudotime bins on the x-axis and the DPI on the y-axis
ggplot(data = data_complete, 
       aes(x = Pseudo.bins2,
           y = Time.Point,
           size = Fraction*100)) +
  theme_minimal() +
  geom_point(color = "black", shape = 21, stroke = 1) +
  geom_point(aes(color = (as.numeric(data_complete$Color_by)))) + 
  scale_color_gradient(low = "white", high = "#5f1654", na.value = NA) + 
  labs(title = "Colored by the average expression of Isg15",
       x = "Pseudotime Bins with equal # of cells", 
       y = "Days Post-Injury", 
       size = "% of cells in a\ngiven bin by DPI", 
       color = "Avg. Isg15\nexpression") +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")) +
  scale_y_discrete(limits = c("7", "5", "3.5", "2", "1", "0")) +
  theme(panel.border = element_rect(color = "black", size = 0.5, linetype = "solid", fill = NA),
        axis.ticks = element_line(color = "black", size = 0.25),
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.key.size = unit(0.20, "cm"),
        legend.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 6), 
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(size = 6, color = "black"), 
        axis.title.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, color = "black", hjust = 1),
        axis.title.y = element_text(size = 6),
        strip.text = element_text(size = 8)) +
  facet_grid(~ factor(data_complete$Age, levels = c("Young", "Old", "Geriatric")), scales = "fixed") +
  guides(size = guide_legend(override.aes = list(shape = 21))) +
  scale_size_area(breaks = c(0, 5, 10, 15), limits = c(0, NA), max_size = 4)

# Save plot as PDF
ggsave(filename = "PseudoBins_DPI_Isg15_DP.pdf", plot = last_plot(), width = 8.5, height = 2, units = "in")

# Alternate option to color by the fraction of non-G1 cells (or some other binary variable) ####
# Generate dataframe with the number of cells at each age, time point,pseudotime bin, and G1 status
df <- data.frame(table(myo.seurat$Age.Word, 
                       myo.seurat$Time.Point, 
                       myo.seurat$Pseudo.bins2, 
                       myo.seurat$G1_Status))  # Can change this to another metadata column

# Define column names                 
colnames(df) <- c("Age", "Time.Point", "Pseudo.bins2", "G1.Status", "Freq")

# Add up the total number of cells by time point
df <- df %>%
  group_by(Time.Point) %>%
  mutate(Sum = sum(Freq))

# Calculate the fraction of cells per age, time point, and pseudotime bin by the total number of cells at the corresponding time point
df$Fraction <- df$Freq/df$Sum

# Add up the total number of cells that are in G1 and are not in G1 per age, time point, and pseudotime bin
df <- df %>%
  group_by(Age, Time.Point, Pseudo.bins2) %>%
  mutate(Phase.Sum = sum(Freq))

# Calculate the fraction of cells that are in G1 and are not in G1 per age, time point, and pseudotime bin
df$NonG1 <- (df$Freq)/(df$Phase.Sum)

# Focus on the cells not in G1
df <- subset(df, G1.Status == FALSE)

# Remove rows that are incomplete (ie. contain NAs)
data_complete <- df[complete.cases(df), ]

# Dotplot with pseudotime bins on the x-axis and the DPI on the y-axis
ggplot(data = data_complete, 
       aes(x = Pseudo.bins2,
           y = Time.Point,
           size = Fraction*100)) +
  theme_minimal() +
  geom_point(color = "black", shape = 21, stroke = 1) +
  geom_point(aes(color = (as.numeric(data_complete$NonG1*100)))) + 
  scale_color_gradient(low = "white", high = "#5f1654", na.value = NA) + 
  labs(title = "Colored by the fraction of cells in S, G2, or M",
       x = "Pseudotime Bins with equal # of cells", 
       y = "Days Post-Injury", 
       size = "% of cells in a\ngiven bin by DPI", 
       color = "Fraction of\ncells in S, G2, or M") +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")) +
  scale_y_discrete(limits = c("7", "5", "3.5", "2", "1", "0")) +
  theme(panel.border = element_rect(color = "black", size = 0.5, linetype = "solid", fill = NA),
        axis.ticks = element_line(color = "black", size = 0.25),
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.key.size = unit(0.20, "cm"),
        legend.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 6), 
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(size = 6, color = "black"), 
        axis.title.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, color = "black", hjust = 1),
        axis.title.y = element_text(size = 6),
        strip.text = element_text(size = 8)) +
  facet_grid(~ factor(data_complete$Age, levels = c("Young", "Old", "Geriatric")), scales = "fixed") +
  guides(size = guide_legend(override.aes = list(shape = 21))) +
  scale_size_area(breaks = c(0, 5, 10, 15), limits = c(0, NA), max_size = 4)

# Save plot as PDF
ggsave(filename = "PseudoBins_DPI_NonG1_DP.pdf", plot = last_plot(), width = 8.5, height = 2, units = "in")

# Dotplot with pseudotime bins on the x-axis and the DPI on the y-axis by replicate ####
# Generate dataframe with the number of cells at each age, time point, pseudotime bin, and sample ID
df <- data.frame(table(myo.seurat$Age.Word, 
                       myo.seurat$Time.Point, 
                       myo.seurat$Pseudo.bins2,
                       myo.seurat$Sample.ID))

# Define column names                 
colnames(df) <- c("Age", "Time.Point", "Pseudo.bins2", "Sample.ID", "Freq")

# Add up the total number of cells by time point
df <- df %>%
  group_by(Time.Point) %>%
  mutate(Sum = sum(Freq))

# Calculate the fraction of cells per age, time point, pseudotime bin, and sample ID by the total number of cells at the corresponding time point
df$Fraction <- df$Freq/df$Sum

# Remove and replace text within the sample IDs
df$Sample.ID <- gsub("26mo_", "", df$Sample.ID)
df$Sample.ID <- gsub("Old_", "", df$Sample.ID)
df$Sample.ID <- gsub("Yng_", "", df$Sample.ID)
df$Sample.ID <- gsub("D0_", "", df$Sample.ID)
df$Sample.ID <- gsub("D1_", "", df$Sample.ID)
df$Sample.ID <- gsub("D2_", "", df$Sample.ID)
df$Sample.ID <- gsub("D3.5_", "", df$Sample.ID)
df$Sample.ID <- gsub("D5_", "", df$Sample.ID)
df$Sample.ID <- gsub("D7_", "", df$Sample.ID)

# Add a column defining the age and time point
df$Age_TP <- paste0(df$Age, " DPI ", df$Time.Point)

# Dotplot with pseudotime bins on the x-axis and the DPI on the y-axis by replicate 
ggplot(data = df, 
       aes(x = Pseudo.bins2,
           y = Sample.ID,
           size = Fraction*100,
           stroke = 0,
           color = Time.Point)) + # Can change this to "Sex" if you want to color by that 
  theme_minimal() +
  geom_point() +
  scale_color_manual("DPI", limits = c("0", "1", "2", "3.5", "5", "7"), values = c("#700324", "#e3104f", "#f584a6", "#024769", "#097cb5", "#60b6e0")) +
  labs(x = "Pseudotime Bins with equal # of cells", 
       y = "Days Post-Injury", 
       size = "% of cells in a\ngiven bin by DPI", 
       color = "Time point") +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")) +
  scale_y_discrete(limits = c("D", "C", "B", "A")) +
  theme(panel.border = element_rect(color = "#D9D9D9", size = 0.5, linetype = "solid", fill = NA),
        axis.ticks = element_line(color = "#D9D9D9", size = 0.25),
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.key.size = unit(0.20, "cm"),
        legend.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 6), 
        axis.text.x = element_text(size = 6, color = "black"), 
        axis.title.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, color = "black", hjust = 1),
        axis.title.y = element_text(size = 6),
        strip.text = element_text(size = 8)) +
  guides(size = guide_legend(override.aes = list(shape = 21))) +
  scale_size_area(breaks = c(0, 4, 8), limits = c(0, NA), max_size = 4) +
  facet_wrap(~ factor(Age_TP, levels = c("Young DPI 0", "Old DPI 0", "Geriatric DPI 0",
                                         "Young DPI 1", "Old DPI 1", "Geriatric DPI 1", 
                                         "Young DPI 2", "Old DPI 2", "Geriatric DPI 2", 
                                         "Young DPI 3.5", "Old DPI 3.5", "Geriatric DPI 3.5", 
                                         "Young DPI 5", "Old DPI 5", "Geriatric DPI 5", 
                                         "Young DPI 7", "Old DPI 7", "Geriatric DPI 7")), scales = "fixed", ncol = 3, nrow = 6)

# Save plot as PDF
ggsave(filename = "PseudoBins_DPI_Rep_TP_DP.pdf", plot = last_plot(), width = 8.5, height = 5, units = "in")
