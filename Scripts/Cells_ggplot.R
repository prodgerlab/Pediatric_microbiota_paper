## Pathological Phimosis is Associated with Foreskin Immune Cell Infiltration but Not Microbiota Composition
## Microbial diversity analysis - alpha diversiy
## Prepared by: Rachel Penney
## Reviewed by: Jorge Rojas-Vargas


# Load libraries
library(readxl)
library(scales)
library(ggplot2)
library(dplyr)

# Load data
metadata <- read_excel("Complete_PDFS_cells.xlsx")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]

#Get medians for abstract

# Replace 'your_data' with the actual name of your data frame
EL_subset <- metadata %>%
  dplyr::filter(Cohort == "EL")
EL_subset_numeric <- EL_subset %>%
  dplyr::select(-1)
EL_subset_numeric <- EL_subset_numeric %>%
  dplyr::mutate(across(everything(), as.numeric))

median_values <- EL_subset_numeric %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))
print(median_values)*1000000

PA_subset <- metadata %>%
  dplyr::filter(Cohort == "PA")
PA_subset_numeric <- PA_subset %>%
  dplyr::select(-1)
PA_subset_numeric <- PA_subset_numeric %>%
  dplyr::mutate(across(everything(), as.numeric))

median_values <- PA_subset_numeric %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))
print(median_values)*1000000

# Ensure variables are correctly formatted
metadata$Cohort <- as.factor(metadata$Cohort)

# Define custom colors
box_colors <- c("EL" = "steelblue1", "NE" = "darkorange", "PA" = "firebrick")

# Custom function to format scientific notation as "3 × 10^-4"
scientific_labels <- function(x) {
  ifelse(x == 0, "0", parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x))))
}


y_breaks <- function(limits) {
  all_breaks <- scales::extended_breaks()(limits)  # Generate standard breaks
  all_breaks[seq(1, length(all_breaks), by = 2)]   # Keep every other one
}

#This is the code used for the boxplots in my paper
# Create the boxplot
# Filter the data (same logic as your ggplot call)
metadata$CD11c_epi_IFS <- as.numeric(metadata$CD11c_epi_IFS) 

filtered_data <- na.omit(metadata[, c("Cohort", "CD11c_epi_IFS")])

# Create the boxplot with counts shown
p1 <- ggplot(filtered_data, aes(x = Cohort, y = CD11c_epi_IFS, fill = Cohort)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.2, color = "black", size = 2, alpha = 0.7) +  
  scale_fill_manual(values = box_colors) +  
  labs(y = NULL, x = NULL) +  
  scale_y_continuous(breaks = y_breaks, labels = scientific_labels) +  #switch c(0,0.1) back to y_breaks to automatic y-axis instead of manual
  theme_classic() +  
  theme(
    axis.title.y = element_blank(),  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.y = element_text(size = 22, color = "black", angle = 90, hjust = 0.5),  
    legend.position = "none",  
    panel.border = element_blank(),  
    axis.line = element_line(linewidth = 1),  
    plot.margin = margin(t = 50, r = 20, b = 50, l = 20)  
  )

p1

citation("Hmisc")




# Now import a new data frame to remove NE and multiply all data by 1M to make it um^2 instead of mm^2
metadata <- read_excel("Cells_NoNE_um^2.xlsx")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]

# Ensure variables are correctly formatted
metadata$Cohort <- as.factor(metadata$Cohort)

# Define custom colors
box_colors <- c("EL" = "steelblue1", "PA" = "firebrick")

y_breaks <- function(limits) {
  all_breaks <- scales::extended_breaks()(limits)  # Generate standard breaks
  all_breaks[seq(1, length(all_breaks), by = 2)]   # Keep every other one
}

#This is the code used for the boxplots in my paper
# Create the boxplot
# Filter the data (same logic as your ggplot call)
metadata$CD68_derm_IFS <- as.numeric(metadata$CD68_derm_IFS) 

filtered_data <- na.omit(metadata[, c("Cohort", "CD68_derm_IFS")])

# Create the boxplot with counts shown
p1 <- ggplot(filtered_data, aes(x = Cohort, y = CD68_derm_IFS, fill = Cohort)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.2, color = "black", size = 2, alpha = 0.7) +  
  scale_fill_manual(values = box_colors) +  
  labs(y = NULL, x = NULL) +  
  scale_y_continuous(breaks = y_breaks) +  
  theme_classic() +  
  theme(
    axis.title.y = element_blank(),  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.y = element_text(size = 22, color = "black", angle = 90, hjust = 0.5),  
    legend.position = "none",  
    panel.border = element_blank(),  
    axis.line = element_line(linewidth = 1),  
    plot.margin = margin(t = 50, r = 20, b = 50, l = 20)  
  )

p1




