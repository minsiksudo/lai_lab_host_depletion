###### Amy Willis, August 21 2023
###### untreated

## Steps
## 1. set up data
## 2. run tinyvamp
## 3. save


######## Step 1

library(dplyr)
library(tibble)
library(magrittr)
library(phyloseq)
library(tinyvamp)
ps <- readRDS("../data/PHY_20230521_MGK_host_tidy.rds")
ps_untreated_genus <- ps$phyloseq_count %>% 
  tax_glom("Genus") %>%
  phyloseq::subset_samples(treatment == "Untreated")  %>%
  phyloseq::prune_taxa(taxa_sums(.) > 0, .)
zero_counts <- which(ps_untreated_genus %>% sample_sums() == 0) %>% names
ps_untreated_genus <- ps_untreated_genus %>%
  phyloseq::subset_samples(!(baylor_other_id %in% zero_counts)) %>% # has zero reads
  phyloseq::prune_taxa(taxa_sums(.) > 0, .) 
W <- ps_untreated_genus %>% otu_table %>% as.matrix %>% t
jj <- ncol(W) 
nn <- nrow(W) 
sample_ids <- ps_untreated_genus %>% sample_data %>% as_tibble %>% pull(sample_id)
subject_ids <- ps_untreated_genus %>% sample_data %>% as_tibble %>% pull(subject_id)
treatment_ids <- ps_untreated_genus %>% sample_data %>% as_tibble %>% pull(treatment)
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


ps_untreated_genus %>% 
  subset_samples(subject_id == "Mock") %>%
  taxa_sums() %>% sort %>% tail(10)


P <- matrix(0, nrow = KK, ncol = jj)
rownames(P) <- unique(specimen_ids)
colnames(P) <- ps_untreated_genus %>% otu_table %>% rownames
P <- ((t(Z) %*% W ) / (t(Z) %*% W %>% apply(1, sum)))
mock_samples <- which((rownames(P) %>% substr(1, 4)) == "Mock")
neg_samples <- which((rownames(P) %>% substr(1, 4)) == "Neg.")
##########################################
##### this is the line that I needed to fix
##########################################
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


DNA_bac_ng_uLs <- ps_untreated_genus %>% sample_data %>% as_tibble %>% pull(DNA_bac_ng_uL)
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

### these initializations arose from previous runs and are provided here to maximize reproducibility
gammas_start <- structure(c(13.1053910105608, 15.7110184530917, 15.39606379788, 
                            15.9428631259283, 13.0547245150994, 15.1054575356105, 12.668631455152, 
                            15.3275004832138, 14.1264568678292, 13.5872940185279, 12.0230613920716, 
                            11.3688212940528, 10.3499719606036, 10.6292355000468, 11.6706685125541, 
                            10.0524627358166, 9.26955120853209, 7.71027089615233, 11.4712636244606, 
                            17.3599255355884, 14.0356837133997, 16.9267422030901, 17.4230870798191, 
                            10.148485406341, 9.81229417582491, 10.5071230125039, 17.4800829720119, 
                            17.4621742830996, 17.5107639365228, 9.65007364648706, 11.9454665314861), dim = c(31L, 1L))
bs_start <- structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.140474346842468, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.67963414074045, 
                        0, 0.336199227408064, 0, 0, 0, 0, 0, 0.0255087306003175, 0, 0, 
                        0, 0.855503936522122, -0.0501454909323506, 0, 0, 0, 0, 1.43524449346565, 
                        0, 0, 0.00834133122850654, 0.15363704997865, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0), dim = c(1L, 67L))
gamma_tilde_start <- structure(-7.41026570596103, dim = c(1L, 1L))
P_start <- structure(c(0, 0, 0, 1.72740721959321e-05, 0, 0, 0, 0, 0, 0, 
                       0.00349429671820065, 0, 0.029816760499058, 0.00425862262970908, 
                       0.0910413326053007, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.34828869169366e-05, 
                       0, 5.24608969764844e-05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0.75279194764397, 0.0656785980409625, 0.2729249044004, 0.96444260185339, 
                       0.175768784990338, 0.922020701787871, 0.0986434571418946, 0.830425588151898, 
                       0.301857142788117, 0.347517167438421, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0178497322115479, 0, 
                       0.0387425933099556, 0.0435057016377144, 0.0121872403047347, 0, 
                       0, 0, 0, 0, 0, 0.0256765746974149, 0.829256200311256, 0.592127548884878, 
                       0.0193889812222199, 0.152872943873033, 0.0254130021090739, 0.40757163760958, 
                       0.0432406996865129, 0.553883240288182, 0.0942512650372467, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0.0139294946994098, 0, 0, 0, 0.00176019089494582, 0, 
                       0, 0, 1.78762497678273e-05, 8.43929355750148e-06, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.81573331690281e-06, 
                       0, 2.14662352747685e-05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0.000829968744320491, 0.000670615089860729, 0, 0, 
                       0, 0.00145767304418648, 0, 0.000740654846884929, 0, 0.00228149957628716, 
                       0.0062136705592197, 0, 0.013404982596443, 0.0107866044181167, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0.00874445768166974, 0.00904138323331914, 0.00893714975266351, 
                       0.0018957404930357, 0.0297328021189437, 0.00828588077906704, 
                       0.332846262603287, 0.0093702845996482, 0.00510100457667748, 0.0290230543049182, 
                       0, 0.558650456584772, 0, 0.000974461621651069, 0.0420069027548727, 
                       0.81879605081101, 0, 0, 0.252677125043202, 0.12, 0, 0.19690876636498, 
                       0, 0, 0.00844003529121944, 0.624879657810459, 0.0313128197570985, 
                       0, 0.100697362137923, 0, 0.519226142231866, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0.0132711980386289, 0.12, 0, 0, 0, 0, 0, 0, 0, 0.00229960489324769, 
                       0, 0, 0, 0.0217583012204586, 0.00360346694191939, 0.117507469222118, 
                       0.160424250196425, 0.296311354146692, 0, 1.00000000005883, 0, 
                       0.0830206010259554, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0.00852338179625593, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0.00323683537173552, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0107268211786847, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0.00494728209673458, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0.0328240303363271, 0.0215834737825634, 0.0616279790020976, 
                       0, 0, 0, 0, 0, 0, 0, 0.00361610997648596, 0.00140874452883069, 
                       0, 0, 0, 0.00955341321218189, 0, 0.00404612929981549, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00998066221515121, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0160500246606518, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0880838954752431, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0.00466877458306849, 0.00334447603511549, 0.00120271256525783, 
                       0.00555971088158454, 0.016745811295818, 0.0126612889577958, 0.129930550360778, 
                       0.0157198549616428, 0.00288251614167806, 0.00998237104776778, 
                       0, 0, 0, 0, 0, 0.089500169710473, 0, 1.00000000008578, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00505532274484008, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00293846161060117, 0, 
                       0, 0, 0, 0, 0, 0, 0.0580387086483091, 0, 0, 0, 0, 0.12, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00232352720004347, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.852783237942414, 
                       0.0632511547882078, 0.75596325214641, 0.678394808540622, 0.195673670758288, 
                       0, 0, 0, 0, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0.0112094791217526, 0.0882153874626672, 
                       0.12271988555083, 0.00024684053624019, 0, 0.000241357544792675, 
                       0.00477827738731308, 0.000493749620180179, 0.131489312124886, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0.29659737253541, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0.00637018697446453, 0, 0.0251458945841694, 
                       0.0393680626028057, 0.177701407854813, 0, 0, 0, 0.325010394129189, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0917037795533182, 
                       0, 0, 0, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00996671999681602, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00824154652628728, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.368281251199463, 
                       0, 0.0327623581861485, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.12, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.12, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00532327828720942, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0.0066658116588676, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), dim = c(21L, 67L))

######## Step 2

untreated_fix_p <- estimate_parameters(W = W,
                                       X = X,
                                       Z = Z,
                                       Z_tilde = Z_tilde,
                                       Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                       gammas = gammas_start,
                                       gammas_fixed_indices = gammas_fixed_indices,
                                       P = P_start,
                                       P_fixed_indices = P_fixed_indices,
                                       B = bs_start,
                                       B_fixed_indices = B_fixed_indices,
                                       X_tilde = X_tilde,
                                       P_tilde = P_tilde,
                                       P_tilde_fixed_indices = P_tilde_fixed_indices,
                                       gamma_tilde = gamma_tilde_start,
                                       gamma_tilde_fixed_indices = gamma_tilde_fixed_indices, 
                                       max_barrier=1e50, 
                                       barrier_maxit=Inf, 
                                       criterion = "Poisson", 
                                       profile_P = TRUE, 
                                       verbose=TRUE)
system("say done")
untreated_fix_p$optimization_status

########################################################
######## Step 3
########################################################

untreated_phats <- untreated_fix_p$P
colnames(untreated_phats) <- colnames(untreated_fix_p$W)
rownames(untreated_phats) <- colnames(untreated_fix_p$Z)
untreated_phats <- untreated_phats[!grepl("Neg.", rownames(untreated_phats)) & 
                                     !grepl("Mock", rownames(untreated_phats)), ]
untreated_phats %>% apply(1, sum)
(untreated_phats > 0) %>% apply(1, sum) %>% summary

# saveRDS(untreated_phats, "../output/untreated_p_hats_v2.RDS")
