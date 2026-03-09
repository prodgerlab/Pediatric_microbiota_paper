## Pathological Phimosis is Associated with Foreskin Immune Cell Infiltration but Not Microbiota Composition
## Abundance correlations and comparisons
## Prepared by: Rachel Penney
## Reviewed by: Jorge Rojas-Vargas

# Load libraries
library(pheatmap)
library(corrplot)
library(Hmisc)
library(dplyr)
library(tidyr)
library(tibble)
library(pald)

# Load data
load("abund_PEDS_genus.RData")
load("abund_PC_genus.RData")

process_genus <- function(df, top = 30) {
  
  # eliminar Unassigned
  df <- df[rownames(df) != "Unassigned", , drop = FALSE]
  
  # ordenar por abundancia total
  row_sums <- rowSums(df[, 2:ncol(df)])
  df <- df[order(row_sums, decreasing = TRUE), ]
  
  # top N
  df <- head(df, top)
  
  # transponer
  df <- as.data.frame(t(df))
  
  
}


## -------------------------------------------------
##
## Spearman's correlation heatmap
##
## -------------------------------------------------

## ---------------
## PED heatmap
## ---------------

abd_ped_30 <- process_genus(abund_PEDS_genus)
cor_matrix_ped <- cor(abd_ped_30, method = "spearman", use = "pairwise.complete.obs")
abs_ped_30 <- abs(cor_matrix_ped[lower.tri(cor_matrix_ped)])

# Perform hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix_ped), method = "complete")
plot(hc, main = "Horizontal Dendrogram", horiz = TRUE)

# Determine the number of clusters dynamically based on a height cutoff (e.g., 0.5)
clusters <- cutree(hc, h = 1.3)
clusters

#h = 1.3, clusters = 2
#h = 1.2, clusters = 5
#h = 1.1, clusters = 6
#h = 1, clusters = 9

# Count the number of unique clusters
num_clusters <- length(unique(clusters))

#Now try adding p-values
cor_test <- cor.mtest(cor_matrix_ped, conf.level = 0.95)

# Apply FDR correction to the p-value matrix
adjusted_p <- matrix(p.adjust(cor_test$p, method = "fdr"), 
                     nrow = nrow(cor_test$p), 
                     ncol = ncol(cor_test$p), 
                     byrow = FALSE)

# Verify dimensions
if (!all(dim(cor_matrix_ped) == dim(adjusted_p))) {
  stop("Dimensions of cor_matrix and adjusted_p do not match.")
}

rownames(adjusted_p) <- rownames(cor_matrix_ped)
colnames(adjusted_p) <- colnames(cor_matrix_ped)


# Add significance stars to the correlation plot using FDR-adjusted p-values
# Generate the corrplot
corrplot(cor_matrix_ped, 
         method = "color", 
         order = "hclust", 
         hclust.method = "complete",
         tl.col = "black",  
         tl.cex = 0.8,
         col = colorRampPalette(c("dodgerblue4", "white", "firebrick"))(200),
         number.cex = 0.4,          # Set text size for significance stars
         p.mat = adjusted_p,        # FDR-adjusted p-values
         sig.level = 0.05,          # Threshold for significance
         insig = "label_sig",       # Add significance stars
         pch.cex = 0.8)             # Adjust size of stars



## ---------------
## ADULT heatmap
## ---------------

abd_adult_30 <- process_genus(abund_PC_genus)
cor_matrix_adult <- cor(abd_adult_30, method = "spearman", use = "pairwise.complete.obs")
abs_adult_30 <-abs(cor_matrix_adult[lower.tri(cor_matrix_adult)])

# Perform hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix_adult), method = "complete")
plot(hc, main = "Horizontal Dendrogram", horiz = TRUE)

# Determine the number of clusters dynamically based on a height cutoff (e.g., 0.5)
clusters <- cutree(hc, h = 1.3)
clusters

# Count the number of unique clusters
num_clusters <- length(unique(clusters))

#Now try adding p-values
cor_test <- cor.mtest(cor_matrix_adult, conf.level = 0.95)

# Apply FDR correction to the p-value matrix
adjusted_p <- matrix(p.adjust(cor_test$p, method = "fdr"), 
                     nrow = nrow(cor_test$p), 
                     ncol = ncol(cor_test$p), 
                     byrow = FALSE)

# Verify dimensions
if (!all(dim(cor_matrix_adult) == dim(adjusted_p))) {
  stop("Dimensions of cor_matrix and adjusted_p do not match.")
}

rownames(adjusted_p) <- rownames(cor_matrix_adult)
colnames(adjusted_p) <- colnames(cor_matrix_adult)


# Add significance stars to the correlation plot using FDR-adjusted p-values
# Generate the corrplot
corrplot(cor_matrix_adult, 
         method = "color", 
         order = "hclust", 
         hclust.method = "complete",
         tl.col = "black",  
         tl.cex = 0.8,
         col = colorRampPalette(c("dodgerblue4", "white", "firebrick"))(200),
         number.cex = 0.4,          # Set text size for significance stars
         p.mat = adjusted_p,        # FDR-adjusted p-values
         sig.level = 0.05,          # Threshold for significance
         insig = "label_sig",       # Add significance stars
         pch.cex = 0.8)             # Adjust size of stars



## ------------------------------------------------------
##
## Comparing Spearman's values between peds and adults
##
## ------------------------------------------------------


## -------------------------------------------------
## Basic statistics & Wilcoxon test of the top 30
## -------------------------------------------------

stats_30 <- tibble(
  cohort      = c("Paediatric", "Adult"),
  median_abs_r  = c(median(abs_ped_30),   median(abs_adult_30)),
  min_abs_r= c(min(abs_ped_30), min(abs_adult_30)),
  max_abs_r= c(max(abs_ped_30), max(abs_adult_30))
)
print(stats_30)
# cohort     median_abs_r min_abs_r max_abs_r
# Paediatric      0.114  0.000187     0.518
# Adult           0.236  0.000992     0.755

# Wilcoxon test
wilcox.test(abs_ped_30, abs_adult_30)
## Results: p-value <2e-16


## -------------------------------------------------
## Identify genera shared by both datasets
## -------------------------------------------------
shared_genera <- intersect(colnames(abd_ped_30), colnames(abd_adult_30))
length(shared_genera)   # should be 22

## -------------------------------------------------
## Correlation matrices for the shared genera
## -------------------------------------------------
# Paediatric
cor_ped_shared <- cor(abd_ped_30[ , shared_genera], method = "spearman", use    = "pairwise.complete.obs")
cor_ped_test <- cor.mtest(cor_ped_shared, conf.level = 0.95)
adjusted_p_ped <- matrix(p.adjust(cor_ped_test$p, method = "fdr"), 
                         nrow = nrow(cor_ped_test$p), 
                         ncol = ncol(cor_ped_test$p), 
                         byrow = FALSE)

# Adult
cor_adult_shared <- cor(abd_adult_30[ , shared_genera], method = "spearman", use    = "pairwise.complete.obs")
cor_adult_test <- cor.mtest(cor_adult_shared, conf.level = 0.95)
adjusted_p_adult <- matrix(p.adjust(cor_adult_test$p, method = "fdr"), 
                           nrow = nrow(cor_adult_test$p), 
                           ncol = ncol(cor_adult_test$p), 
                           byrow = FALSE)

# Absolute correlations (lower‑triangle, no diagonal)
abs_ped   <- abs(cor_ped_shared[lower.tri(cor_ped_shared)])
abs_adult <- abs(cor_adult_shared[lower.tri(cor_adult_shared)])

# Real correlations (lower‑triangle, no diagonal)
real_ped   <- cor_ped_shared[lower.tri(cor_ped_shared)]
real_adult <- cor_adult_shared[lower.tri(cor_adult_shared)]


## -------------------------------------------------
## Basic statistics & Wilcoxon test of shared taxa
## -------------------------------------------------
stats <- tibble(
  cohort      = c("Paediatric", "Adult"),
  median_abs_r  = c(median(abs_ped),   median(abs_adult)),
  min_abs_r= c(min(abs_ped), min(abs_adult)),
  max_abs_r= c(max(abs_ped), max(abs_adult))
)
print(stats)
# cohort     median_abs_r min_abs_r max_abs_r
# Paediatric      0.117  0.00141      0.518
# Adult           0.269  0.000992     0.755

wilcox_res <- wilcox.test(abs_ped, abs_adult)
print(wilcox_res)
## Results: p-value <2e-16


## -------------------------------------------------
## Which genus‑pairs change the most?
## -------------------------------------------------
# Vector of pair names for clarity
get_pair_names <- function(nms) {
  # returns "GenusA|GenusB" for each lower‑tri pair
  idx <- which(lower.tri(matrix(0, length(nms), length(nms))), arr.ind = TRUE)
  paste(nms[idx[,1]], nms[idx[,2]], sep = "|")
}

pair_names <- get_pair_names(shared_genera)

delta_df <- tibble(
  Pair        = pair_names,
  abs_r_ped   = abs_ped,
  abs_r_adult = abs_adult,
  Diff        = abs_r_adult - abs_r_ped         # positive = stronger in adults
) %>%
  arrange(desc(abs(Diff)))                      # largest shifts first

# Show top 10 shifts
print(head(delta_df, 10))
# Pair                            abs_r_ped abs_r_adult  Diff
# 1 Propionimicrobium|Campylobacter   0.0611        0.636 0.574
# 2 Dialister|Staphylococcus          0.0107        0.533 0.522
# 3 Actinotignum|Staphylococcus       0.0677        0.575 0.508
# 4 Winkia|Finegoldia                 0.260         0.755 0.495
# 5 Campylobacter|Hoylesella          0.00893       0.490 0.482
# 6 Propionimicrobium|Varibaculum     0.0171        0.496 0.479
# 7 Murdochiella|Ezakiella            0.00141       0.476 0.475
# 8 Actinotignum|Anaerococcus         0.00177       0.472 0.471
# 9 Staphylococcus|Hoylesella         0.0649        0.528 0.463
# 10 Mobiluncus|Finegoldia            0.0562        0.515 0.459


delta_df_real <- tibble(
  Pair        = pair_names,
  real_r_ped   = real_ped,
  real_r_adult = real_adult,
  Diff        = real_r_adult - real_r_ped         # positive = stronger in adults
) %>%
  arrange(desc(abs(Diff)))                      # largest shifts first

# Show top 10 shifts
print(head(delta_df_real, 20))
#    Pair                             real_r_ped real_r_adult   Diff
# 1 Propionimicrobium|Porphyromonas     -0.212         0.423  0.635
# 2 Mobiluncus|Ezakiella                -0.149         0.475  0.624
# 3 Staphylococcus|Hoylesella            0.0649       -0.528 -0.593
# 4 Propionimicrobium|Staphylococcus     0.129        -0.462 -0.591
# 5 Propionimicrobium|Campylobacter      0.0611        0.636  0.574
# 6 Mobiluncus|Finegoldia                0.0562       -0.515 -0.571
# 7 Actinotignum|Hoylesella             -0.0528        0.505  0.557
# 8 Dialister|Staphylococcus            -0.0107       -0.533 -0.522
# 9 Actinotignum|Staphylococcus         -0.0677       -0.575 -0.508
# 10 Ezakiella|Hoylesella               -0.124         0.375  0.499
# 11 Fusobacterium|Campylobacter        -0.144          0.353  0.497
# 12 Winkia|Finegoldia                   0.260          0.755  0.495
# 13 Campylobacter|Hoylesella            0.00893        0.490  0.482
# 14 Negativicoccus|Prevotella           0.209         -0.271 -0.480
# 15 Actinotignum|Varibaculum           -0.0670         0.412  0.479
# 16 Propionimicrobium|Varibaculum       0.0171         0.496  0.479
# 17 Murdochiella|Ezakiella             -0.00141        0.476  0.478
# 18 Mobiluncus|Corynebacterium          0.176         -0.300 -0.476
# 19 Streptococcus|Prevotella           -0.212          0.260  0.472
# 20 Actinotignum|Anaerococcus          -0.00177       -0.472 -0.471





## ------------------------------------------------------
##
## Cohesion analysis for peds and adults
##
## ------------------------------------------------------


#------------------------------------------------------------
# Function to compute cohesion according to Herren & McMahon 2017
#------------------------------------------------------------
compute_cohesion <- function(rel.d, iter = 500) {
  # rel.d: matrix (samples × taxa) of relative abundances
  # iter: number of permutations for the null model
  
  taxa <- colnames(rel.d)
  nt <- ncol(rel.d)
  
  # Observed correlation
  cor.true <- cor(rel.d, method = "spearman", use = "pairwise.complete.obs")
  
  # Matrix to store the median of null correlations
  med.null <- matrix(NA, nrow = nt, ncol = nt,
                     dimnames = list(taxa, taxa))
  
  # Loop over focal taxa
  for (k in seq_len(nt)) {
    # matrix (nt × iter) storing all null correlations for taxon k
    perm.cor.mat <- matrix(NA, nrow = nt, ncol = iter)
    
    for (i in seq_len(iter)) {
      # copy rel.d and permute each row except column k
      perm.d <- rel.d
      for (j in seq_len(nrow(rel.d))) {
        pos <- which(rel.d[j, ] > 0)
        nf  <- setdiff(pos, k)
        if (length(nf) > 1) {
          perm.d[j, nf] <- sample(rel.d[j, nf])
        }
      }
      # null correlation
      perm.cor.mat[, i] <- cor(perm.d, method = "spearman")[, k]
    }
    
    # median of null correlations for taxon k
    med.null[, k] <- apply(perm.cor.mat, 1, median, na.rm = TRUE)
  }
  
  # Observed – expected matrix
  obs.exp <- cor.true - med.null
  diag(obs.exp) <- 0
  
  # Connectedness functions
  pos.mean <- function(x) if(any(x > 0)) mean(x[x > 0], na.rm = TRUE) else 0
  neg.mean <- function(x) if(any(x < 0)) mean(x[x < 0], na.rm = TRUE) else 0
  
  connected.pos <- apply(obs.exp, 2, pos.mean)
  connected.neg <- apply(obs.exp, 2, neg.mean)
  
  # Sample-level cohesion
  pos.cohesion <- as.vector(rel.d %*% connected.pos)
  neg.cohesion <- as.vector(rel.d %*% connected.neg)
  
  tibble::tibble(
    sample        = rownames(rel.d),
    pos.cohesion  = pos.cohesion,
    neg.cohesion  = neg.cohesion
  )
}

compute_connectedness <- function(rel.d, iter = 500) {
  taxa <- colnames(rel.d)
  nt <- ncol(rel.d)
  
  # Observed correlation
  cor.true <- cor(rel.d, method = "spearman", use = "pairwise.complete.obs")
  
  # Matrix to store the median of null correlations
  med.null <- matrix(NA, nrow = nt, ncol = nt,
                     dimnames = list(taxa, taxa))
  
  # Loop over focal taxa
  for (k in seq_len(nt)) {
    # matrix (nt × iter) storing all null correlations for taxon k
    perm.cor.mat <- matrix(NA, nrow = nt, ncol = iter)
    
    for (i in seq_len(iter)) {
      # copy rel.d and permute each row except column k
      perm.d <- rel.d
      for (j in seq_len(nrow(rel.d))) {
        pos <- which(rel.d[j, ] > 0)
        nf  <- setdiff(pos, k)
        if (length(nf) > 1) {
          perm.d[j, nf] <- sample(rel.d[j, nf])
        }
      }
      # null correlation
      perm.cor.mat[, i] <- cor(perm.d, method = "spearman")[, k]
    }
    
    # median of null correlations for taxon k
    med.null[, k] <- apply(perm.cor.mat, 1, median, na.rm = TRUE)
  }
  
  # Observed – expected matrix
  obs.exp <- cor.true - med.null
  diag(obs.exp) <- 0
  
  # Connectedness functions
  pos.mean <- function(x) if(any(x>0)) mean(x[x>0]) else 0
  neg.mean <- function(x) if(any(x<0)) mean(x[x<0]) else 0
  
  connected.pos <- apply(obs.exp, 2, pos.mean)
  connected.neg <- apply(obs.exp, 2, neg.mean)
  
  tibble(
    Taxon             = colnames(rel.d),
    pos_connectedness = connected.pos,
    neg_connectedness = connected.neg
  )
}


#------------------------------------------------------------
# Apply to pediatric and adult cohorts
#------------------------------------------------------------

# Prepare relative abundance tables
ped.rel   <- as.matrix(abd_ped_30)
adult.rel <- as.matrix(abd_adult_30)

# Compute cohesion per sample
coh_ped   <- compute_cohesion(ped.rel,   iter = 1000) %>% mutate(cohort = "Pediatric")
coh_adult <- compute_cohesion(adult.rel, iter = 1000) %>% mutate(cohort = "Adult")

coh_all <- bind_rows(coh_ped, coh_adult)

# Compute genus-level connectedness vectors
conn_ped   <- compute_connectedness(ped.rel,   iter = 1000) %>% mutate(cohort = "Pediatric")
conn_adult <- compute_connectedness(adult.rel, iter = 1000) %>% mutate(cohort = "Adult")

conn_all <- bind_rows(conn_ped, conn_adult)


#------------------------------------------------------------
# Compare cohesion distributions
#------------------------------------------------------------
wilcox_pos <- wilcox.test(pos.cohesion ~ cohort, data = coh_all)
wilcox_neg <- wilcox.test(neg.cohesion ~ cohort, data = coh_all)

print(wilcox_pos)
# W = 2543, p-value = 0.04
print(wilcox_neg)
# W = 715, p-value = 1e-10

# Cohort-level summary
coh_all %>% 
  group_by(cohort) %>% 
  summarise(
    median_pos = median(pos.cohesion),
    median_neg = median(neg.cohesion)
  ) %>% print()

# cohort     median_pos median_neg
# 1 Adult         0.152     -0.134
# 2 Pediatric     0.147     -0.112


#------------------------------------------------------------
# Compare connectedness distributions
#------------------------------------------------------------
wilcox_pos <- wilcox.test(pos_connectedness ~ cohort, data = conn_all)
wilcox_neg <- wilcox.test(neg_connectedness ~ cohort, data = conn_all)

print(wilcox_pos)
# W = 592, p-value = 0.04
print(wilcox_neg)
# W = 237, p-value = 0.001

# Cohort-level summary
conn_all %>% 
  group_by(cohort) %>% 
  summarise(
    median_pos = median(pos_connectedness),
    median_neg = median(neg_connectedness)
  ) %>% print()

# cohort    median_pos median_neg
# 1 Adult          0.169     -0.148
# 2 Pediatric      0.143     -0.104



