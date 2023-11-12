###### Amy Willis, August 21 2023
###### lyPMA

## Steps
## 1. set up data
## 2. rerun tinyvamp
## 3. save

######## Step 1

library(dplyr)
library(tibble)
library(magrittr)
library(phyloseq)
library(tinyvamp)
ps <- readRDS("../data/PHY_20230521_MGK_host_tidy.rds")
ps_lympa_genus <- ps$phyloseq_count %>% 
  tax_glom("Genus") %>%
  phyloseq::subset_samples(treatment == "lyPMA")  %>%
  phyloseq::prune_taxa(taxa_sums(.) > 0, .)
zero_counts <- which(ps_lympa_genus %>% sample_sums() == 0) %>% names
ps_lympa_genus <- ps_lympa_genus %>%
  phyloseq::subset_samples(!(baylor_other_id %in% zero_counts)) %>% # has zero reads
  phyloseq::prune_taxa(taxa_sums(.) > 0, .) 
W <- ps_lympa_genus %>% otu_table %>% as.matrix %>% t
jj <- ncol(W) 
nn <- nrow(W) 
sample_ids <- ps_lympa_genus %>% sample_data %>% as_tibble %>% pull(sample_id)
subject_ids <- ps_lympa_genus %>% sample_data %>% as_tibble %>% pull(subject_id)
treatment_ids <- ps_lympa_genus %>% sample_data %>% as_tibble %>% pull(treatment)
specimen_ids <- paste(subject_ids, treatment_ids, sep = "_")
sample_specimen_map <- tibble(sample_ids, specimen_ids)
KK <- specimen_ids %>% unique %>% length
Z <- matrix(0, nrow = nn, ncol = KK)
rownames(Z) <- sample_ids
colnames(Z) <- unique(specimen_ids) 
for (i in 1:nn) {
  Z[i, sample_specimen_map[i, 2] %>% unlist] <- 1
}

stopifnot(all(apply(Z, 1, sum) == 1))


ps_lympa_genus %>% 
  subset_samples(subject_id == "Mock") %>%
  taxa_sums() %>% sort %>% tail(10)


P <- matrix(0, nrow = KK, ncol = jj)
rownames(P) <- unique(specimen_ids)
colnames(P) <- ps_lympa_genus %>% otu_table %>% rownames
P <- ((t(Z) %*% W ) / (t(Z) %*% W %>% apply(1, sum)))
mock_samples <- which((rownames(P) %>% substr(1, 4)) == "Mock")
neg_samples <- which((rownames(P) %>% substr(1, 4)) == "Neg.")
P[mock_samples, ] <- 0 
P[mock_samples, "Listeria_monocytogenes"] <- 0.12
P[mock_samples, "Pseudomonas_aeruginosa_group"] <- 0.12
P[mock_samples, "Bacillus_intestinalis"] <- 0.12 ## note misclassification
P[mock_samples, "Escherichia_coli"] <- 0.12 
P[mock_samples, "Salmonella_enterica"] <- 0.12 
P[mock_samples, "Lactobacillus_fermentum"] <- 0.12 
P[mock_samples, "Enterococcus_faecalis"] <- 0.12 
P[mock_samples, "Staphylococcus_aureus"] <- 0.12 
P[mock_samples, "Saccharomyces_cerevisiae"] <- 0.02 
P[mock_samples, "Cryptococcus_neoformans"] <- 0.02 
P[neg_samples, ] <- 0

P[mock_samples, ] %>% sum ### checked

P_fixed_indices <- matrix(FALSE, nrow = KK, ncol = jj)
rownames(P_fixed_indices) <- unique(specimen_ids)
colnames(P_fixed_indices) <- colnames(W)
P_fixed_indices[mock_samples, ] <- TRUE
P_fixed_indices[neg_samples, ] <- TRUE

old_bhats <- c(Staphylococcus_aureus = -0.365, Escherichia_coli = 0.478, 
               Pseudomonas_aeruginosa_group = -0.765, 
               Lactobacillus_fermentum = 0.136, Bacillus_intestinalis = 0.3, 
               Listeria_monocytogenes = -0.006, Salmonella_enterica = 0.388, 
               Saccharomyces_cerevisiae = 0.091, Cryptococcus_neoformans = -0.269)
tmts <- treatment_ids %>% unique
pp <- tmts %>% length # 1
B <- matrix(0, nrow = pp, ncol = jj)
colnames(B) <- colnames(W)
rownames(B) <- tmts
B[, names(old_bhats)] <- matrix(old_bhats, nrow = pp, ncol = length(old_bhats),
                                byrow = T)

B_fixed_indices <- matrix(TRUE, nrow = pp, ncol = jj) ## in general, don't estimate efficiencies...
colnames(B_fixed_indices) <- colnames(W)
rownames(B_fixed_indices) <- tmts
B_fixed_indices[, c("Listeria_monocytogenes", "Pseudomonas_aeruginosa_group", 
                    "Bacillus_intestinalis", "Escherichia_coli", 
                    "Salmonella_enterica", "Lactobacillus_fermentum", 
                    "Enterococcus_faecalis", "Staphylococcus_aureus", 
                    "Saccharomyces_cerevisiae", "Cryptococcus_neoformans")] <- FALSE
B_fixed_indices[, "Enterococcus_faecalis"] <- TRUE

X <- matrix(1, nrow = nn, ncol = pp)
gammas_fixed_indices <- rep(FALSE, nn)
gammas <- log(apply(W, 1, sum)) + rnorm(nn)

KK_tilde <- pp
P_tilde <- matrix(0, ncol = jj, nrow = KK_tilde)
P_tilde_fixed_indices <- matrix(FALSE, ncol = jj, nrow = KK_tilde)
neg_profile <- (t(Z) %*% W)[neg_samples, ] 
P_tilde <- matrix(neg_profile / sum(neg_profile), ncol = jj, nrow = KK_tilde)


DNA_bac_ng_uLs <- ps_lympa_genus %>% sample_data %>% as_tibble %>% pull(DNA_bac_ng_uL)
Zis <- 1/DNA_bac_ng_uLs
Z_tilde <- matrix(0, nrow = nn, ncol = KK_tilde)
rownames(Z_tilde) <- rownames(Z)
colnames(Z_tilde) <- tmts
specimen_treatment_map <- tibble(specimen_ids, treatment_ids) %>% 
  mutate(treatment_ids = as.character(treatment_ids))
for (i in 1:nn) {
  Z_tilde[i, specimen_treatment_map[i, "treatment_ids"] %>% unlist] <- Zis[i]/exp(mean(log(Zis)))
}

Z_tilde_gamma_cols <- 1

gamma_tilde <- matrix(0, nrow = KK_tilde, ncol = 1)
gamma_tilde_fixed_indices <- matrix(FALSE, nrow = KK_tilde, ncol = 1)

X_tilde <- matrix(0, nrow = KK_tilde, ncol = pp) 

######## Step 2
lypma_fix_p <- estimate_parameters(W = W,
                                   X = X,
                                   Z = Z,
                                   Z_tilde = Z_tilde,
                                   Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                   gammas = gammas,
                                   gammas_fixed_indices = gammas_fixed_indices,
                                   P = P,
                                   P_fixed_indices = P_fixed_indices,
                                   B = B,
                                   B_fixed_indices = B_fixed_indices,
                                   X_tilde = X_tilde,
                                   P_tilde = P_tilde,
                                   P_tilde_fixed_indices = P_tilde_fixed_indices,
                                   gamma_tilde = gamma_tilde,
                                   gamma_tilde_fixed_indices = gamma_tilde_fixed_indices, 
                                   profile_P = TRUE, 
                                   verbose=TRUE, 
                                   max_barrier=1e50, 
                                   profiling_maxit = 200,
                                   criterion = "Poisson")
system("say done")
lypma_fix_p$optimization_status


########################################################
######## Step 3
########################################################

lypma_phats <- lypma_fix_p$P
colnames(lypma_phats) <- colnames(lypma_fix_p$W)
rownames(lypma_phats) <- colnames(lypma_fix_p$Z)
lypma_phats <- lypma_phats[!grepl("Neg.", rownames(lypma_phats)) & 
                             !grepl("Mock", rownames(lypma_phats)), ]
lypma_phats %>% apply(1, sum)
(lypma_phats > 0) %>% apply(1, sum) %>% summary

saveRDS(lypma_phats, "../output/lypma_p_hats_all_v3.RDS")
