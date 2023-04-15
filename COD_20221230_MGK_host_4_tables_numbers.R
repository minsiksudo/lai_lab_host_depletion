#Required packages
library(readxl)
library(ggplot2)
library(phyloseq)
library(vegan)
library(tidyverse)
library(microbiome)
library(viridis)
library(ggpubr)
library(kableExtra)
library(ggtext)
library(lme4)
library(lmerTest)

rm(list = ls())

# Loading files -----------------------------------------------------------
#loading tidy phyloseq object
phyloseq <- read_rds("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221129_MGK_host_tidy_tax.rds")
phyloseq_path <- read_rds("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221229_MGK_host_tidy_path.rds")
#sample data loading


#Formattings
sample_data <- sample_data(phyloseq$phyloseq_count) %>% data.frame(check.names = F)
sample_data$treatment

#phyloseq object



# Statistics of samples for main text -------------------------------------


#FIgure 1 numbers
#qPCR result by sample type and control (1/0)

sample_data %>% data.frame() %>% 
        dplyr::filter(sample_type %in% c("Sputum", "nasal_swab", "BAL")) %>% 
        group_by (sample_type, control) %>%
        arrange(proportion) %>%
        summarise(N = n(),
                  `Total DNA<br>(median [IQR])<br>[&mu;g mL<sup>-1</sup>]` =
                          paste(format(round(median(DNA_host_nondil + DNA_bac_nondil),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(DNA_host_nondil + DNA_bac_nondil, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(DNA_host_nondil + DNA_bac_nondil, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
#                  `Total DNA (mean ± SD)` = paste(format(round(mean(DNA_host_nondil + DNA_bac_nondil),2), nsmall = 2, big.mark = ","), " ± ", format(round(sd(DNA_host_nondil + DNA_bac_nondil, 0.25),2), nsmall = 2, big.mark = ","), sep = ""),
                  `Host DNA<br>(median [IQR])<br>[&mu;g mL<sup>-1</sup>]` =
        paste(format(round(median(DNA_host_nondil),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(DNA_host_nondil, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(DNA_host_nondil, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
 #                 `Host DNA (mean ± SD)` = paste(format(round(mean(DNA_host_nondil),2), nsmall = 2, big.mark = ","), " ± ", format(round(sd(DNA_host_nondil, 0.25),2), nsmall = 2, big.mark = ","), sep = ""),
                  `Bacterial DNA<br>(median [IQR])<br>[&mu;g mL<sup>-1</sup>]` =
        paste(format(round(median(DNA_bac_nondil),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(DNA_bac_nondil, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(DNA_bac_nondil, 0.75),), nsmall = 2, big.mark = ","), "]", sep = ""),
  #                `Bacterial DNA (mean ± SD)` = paste(format(round(mean(DNA_bac_nondil),2), nsmall = 2, big.mark = ","), " ± ", format(round(sd(DNA_bac_nondil, 0.25),2), nsmall = 2, big.mark = ","), sep = ""),
                  `DNA proportion<br>(median [IQR])<br>[%]`=
        paste(format(round(median(proportion),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(proportion, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(proportion, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
   #               `DNA proportion (mean ± SD)` = paste(format(round(mean(proportion),2), nsmall = 2, big.mark = ","), " ± ", format(round(sd(proportion, 0.25),2), nsmall = 2, big.mark = ","), sep = ""),
        ) %>% data.frame(check.names = F) %>%
        arrange(sample_type, control) %>% mutate_all(linebreak) %>% kbl(format = "html", escape = F) %>% kable_styling(full_width = 0)


# #qPCR result by treatment method ----------------------------------------

sample_data %>% data.frame() %>% 
        dplyr::filter(sample_type %in% c("Sputum", "nasal_swab", "BAL")) %>% 
        mutate(host_perc_seq = Host_mapped/Raw_reads * 100) %>% 
        mutate(`Sample type` = case_when(sample_type == "BAL" ~ 'BAL',
                                         sample_type == "nasal_swab" ~ 'Nasal swab',
                                         sample_type == "Sputum" ~ 'Sputum',)) %>%
        mutate(Treatment = case_when(treatment == "control" ~ 'Control',
                                     treatment == "lypma" ~ 'lyPMA',
                                     treatment == "benzonase" ~ 'Benzonase',
                                     treatment == "host_zero" ~ 'Host-zero',
                                     treatment == "molysis" ~ 'Molysis',
                                     treatment == "qiaamp" ~ 'QIAamp',)) %>%
        group_by (`Sample type`, `Treatment`) %>%
        arrange(proportion) %>%
        summarise(N = n(),
                  `Total DNA<br>(median [IQR])<br>[&mu;g mL<sup>-1</sup>]` = paste(format(round(median(DNA_host_nondil + DNA_bac_nondil),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(DNA_host_nondil + DNA_bac_nondil, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(DNA_host_nondil + DNA_bac_nondil, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  #`Total DNA (mean ± SD)` = paste(format(round(mean(DNA_host_nondil + DNA_bac_nondil),2), nsmall = 2, big.mark = ","), " ± ", format(round(sd(DNA_host_nondil + DNA_bac_nondil, 0.25),2), nsmall = 2, big.mark = ","), sep = ""),
                  `Host DNA<br>(median [IQR])<br>[&mu;g mL<sup>-1</sup>]` = paste(format(round(median(DNA_host_nondil),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(DNA_host_nondil, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(DNA_host_nondil, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  #`Host DNA (mean ± SD)` = paste(format(round(mean(DNA_host_nondil),2), nsmall = 2, big.mark = ","), " ± ", format(round(sd(DNA_host_nondil, 0.25),2), nsmall = 2, big.mark = ","), sep = ""),
                  `Bacterial DNA<br>(median [IQR])<br>[&mu;g mL<sup>-1</sup>]` = paste(format(round(median(DNA_bac_nondil),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(DNA_bac_nondil, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(DNA_bac_nondil, 0.75),), nsmall = 2, big.mark = ","), "]", sep = ""),
                  #`Bacterial DNA (mean ± SD)` = paste(format(round(mean(DNA_bac_nondil),2), nsmall = 2, big.mark = ","), " ± ", format(round(sd(DNA_bac_nondil, 0.25),2), nsmall = 2, big.mark = ","), sep = ""),
                  `DNA proportion<br>(median [IQR])<br>[%]` = paste(format(round(median(proportion * 100),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(proportion * 100, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(proportion * 100, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  #`DNA proportion (mean ± SD)` = paste(format(round(mean(proportion),2), nsmall = 2, big.mark = ","), " ± ", format(round(sd(proportion, 0.25),2), nsmall = 2, big.mark = ","), sep = ""),
        ) %>% data.frame(check.names = F) %>% 
        mutate(Treatment = factor(Treatment, levels = c("Control", "lyPMA", "Benzonase", "Host-zero", "Molysis", "QIAamp"))) %>% 
        arrange(`Sample type`, `Treatment`) %>%  mutate_all(linebreak) %>% kbl(format = "html", escape = F) %>% kable_styling(full_width = 0)




# #Table2 - numbers -------------------------------------------------------

#seqeuncing results of all samples

sample_data %>% data.frame() %>% 
        dplyr::filter(sample_type %in% c("Sputum", "nasal_swab", "BAL")) %>% 
        arrange(proportion) %>% 
        summarise(N = n(),
                  `Raw reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Raw_reads/10000000),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(Raw_reads/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Raw_reads/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  `Host reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Host_mapped/10000000),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(Host_mapped/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Host_mapped/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  `Host reads proportion<br>(median [IQR])<br>[%]` = paste(format(round(median(sequencing_host_prop * 100),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(sequencing_host_prop * 100, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(sequencing_host_prop * 100, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  `Final reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Final_reads/10000000),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(Final_reads/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Final_reads/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
        ) %>% data.frame(check.names = F) %>% mutate_all(linebreak) %>% kbl(format = "html", escape = F) %>% kable_styling(full_width = 0)


#sequencing result by sample type and control (1/0)


sample_data %>% data.frame() %>% 
        dplyr::filter(sample_type %in% c("Sputum", "nasal_swab", "BAL")) %>% 
        mutate(host_perc_seq = Host_mapped/Raw_reads * 100) %>% 
        mutate(`Sample type` = case_when(sample_type == "BAL" ~ 'BAL',
                                         sample_type == "nasal_swab" ~ 'Nasal swab',
                                         sample_type == "Sputum" ~ 'Sputum',)) %>%
        mutate(Treatment = case_when(treatment == "control" ~ 'Control',
                                     treatment == "lypma" ~ 'lyPMA',
                                     treatment == "benzonase" ~ 'Benzonase',
                                     treatment == "host_zero" ~ 'Host-zero',
                                     treatment == "molysis" ~ 'Molysis',
                                     treatment == "qiaamp" ~ 'QIAamp',)) %>%
        group_by (`Sample type`, control) %>%
        arrange(proportion) %>% 
        summarise(N = n(),
                  `Raw reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Raw_reads/10000000),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(Raw_reads/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Raw_reads/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  `Host reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Host_mapped/10000000),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(Host_mapped/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Host_mapped/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  `Host reads proportion<br>(median [IQR])<br>[%]` = paste(format(round(median(sequencing_host_prop * 100),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(sequencing_host_prop * 100, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(sequencing_host_prop * 100, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
                  `Final reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Final_reads/10000000),2), nsmall = 2, big.mark = ","), " [", format(round(quantile(Final_reads/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Final_reads/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
        ) %>% data.frame(check.names = F) %>% mutate_all(linebreak) %>% kbl(format = "html", escape = F) %>% kable_styling(full_width = 0)


#Figure3 or table (Table 1)
sample_data$sample_type
sample_data %>% data.frame() %>% 
        dplyr::filter(sample_type %in% c("Sputum", "Nasal swab", "BAL")) %>% 
        mutate(host_perc_seq = Host_mapped/Raw_reads * 100) %>% 
        mutate(`Sample type` = case_when(sample_type == "BAL" ~ 'BAL',
                                         sample_type == "Nasal swab" ~ 'Nasal swab',
                                         sample_type == "Sputum" ~ 'Sputum',)) %>%
        mutate(Treatment = case_when(treatment == "Control" ~ 'Control',
                                     treatment == "lyPMA" ~ 'lyPMA',
                                     treatment == "Benzonase" ~ 'Benzonase',
                                     treatment == "Host zero" ~ 'Host-zero',
                                     treatment == "Molysis" ~ 'Molysis',
                                     treatment == "QIAamp" ~ 'QIAamp',)) %>%
        group_by (`Sample type`, `Treatment`) %>%
        arrange(host_perc_seq) %>% 
        summarise(N = n(),
              #    `Raw reads (mean ± SD)` = paste(format(round(mean(Raw_reads),0), nsmall = 2, big.mark = ","), format(round(sd(Raw_reads),0), nsmall = 2, big.mark = ","), sep = " ± "),
              #    `Host reads (mean ± SD)` = paste(format(round(mean(Host_mapped),0), nsmall = 2, big.mark = ","), format(round(sd(Host_mapped),0), nsmall = 2, big.mark = ","), sep = " ± "),
              #    `% host reads (mean ± SD)` = paste(round(mean(host_perc_seq), 2), round(sd(host_perc_seq), 2), sep = " ± "),
              #    `Final reads (mean ± SD)` = paste(format(round(mean(Final_reads),0), nsmall = 2, big.mark = ","), format(round(sd(Final_reads),0), nsmall = 2, big.mark = ","), sep = " ± "),
              `Raw reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Raw_reads/10000000),2), nsmall = 2, big.mark = ","), "<br>[", format(round(quantile(Raw_reads/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Raw_reads/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
              `Host reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Host_mapped/10000000),2), nsmall = 2, big.mark = ","), "<br>[", format(round(quantile(Host_mapped/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Host_mapped/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
              `Host reads proportion<br>(median [IQR])<br>[%]` = paste(format(round(median(sequencing_host_prop * 100),2), nsmall = 2, big.mark = ","), "<br>[", format(round(quantile(sequencing_host_prop * 100, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(sequencing_host_prop * 100, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
              `Final reads<br>(median [IQR])<br>[reads x 10<sup>7</sup>]` = paste(format(round(median(Final_reads/10000000),2), nsmall = 2, big.mark = ","), "<br>[", format(round(quantile(Final_reads/10000000, 0.25),2), nsmall = 2, big.mark = ","), ", ", format(round(quantile(Final_reads/10000000, 0.75),2), nsmall = 2, big.mark = ","), "]", sep = ""),
        ) %>% data.frame(check.names = F) %>% 
        mutate(Treatment = factor(Treatment, levels = c("Control", "lyPMA", "Benzonase", "Host-zero", "Molysis", "QIAamp"))) %>% 
        arrange(`Sample type`, `Treatment`) %>%
        mutate_all(linebreak) %>% kbl(format = "html", escape = F) %>% kable_styling(full_width = 0) %>%
        save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/Table1_IQR.html", self_contained = T)



#Supplemenatary tables

# Supplementary tables - model of sequencing result -----------------------


#Model with final reads


#Table S1
lmer(log10(Final_reads) ~ sample_type * treatment + (1|original_sample), data = sample_data) %>%
        summary() 

lmer(log10(Final_reads) ~ sample_type * treatment + (1|original_sample), data = sample_data) %>%
        summary() %>%
        .$coefficients %>%
        kbl %>%
        save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/TableS1.html", self_contained = T)
#Table S2
lmer(sequencing_host_prop ~ sample_type * treatment + (1|original_sample), data = sample_data) %>% summary()

displlmer(sequencing_host_prop ~ sample_type * treatment + (1|original_sample), data = sample_data) %>% summary() %>%
        .$coefficients %>%
        kbl %>%
        save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/TableS2.html", self_contained = T)


# Supplementary tables - model of alpha diversity ----------------------------------


#Table S3
lmer(S.obs ~ sample_type * treatment + log10 (Final_reads) + (1|original_sample), data = sample_data) %>% summary()

lmer(S.obs ~ sample_type * treatment + log10 (Final_reads) + (1|original_sample), data = sample_data) %>% summary() %>%
        .$coefficients %>%
        kbl %>%
        save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/TableS3.html", self_contained = T)

#Table S4
lmer(S.obs ~ treatment + log10 (Final_reads) + (1|original_sample), data = subset(sample_data, sample_data$sample_type == "Nasal swab")) %>% summary() %>% .$coefficients %>%
        kbl %>% save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/TableS4.html", self_contained = T)
#Table S5
lmer(S.obs ~ treatment + log10 (Final_reads) + (1|original_sample), data = subset(sample_data, sample_data$sample_type == "BAL")) %>% summary()%>% .$coefficients %>%
        kbl %>% save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/TableS5.html", self_contained = T)
#Table S6
lmer(S.obs ~ treatment + log10 (Final_reads) + (1|original_sample), data = subset(sample_data, sample_data$sample_type == "Sputum")) %>% summary()%>% .$coefficients %>%
        kbl %>% save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/TableS6.html", self_contained = T)


#Table S7
lmer(data_shannon ~ sample_type * treatment + log10 (Final_reads) + (1|original_sample), data = sample_data) %>% summary()

lmer(data_shannon ~ sample_type * treatment + log10 (Final_reads) + (1|original_sample), data = sample_data) %>% summary() %>%
        .$coefficients %>%
        kbl %>%
        save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/TableS7.html", self_contained = T)

#Table S8
lmer(dbp ~ sample_type * treatment + log10 (Final_reads) + (1|original_sample), data = sample_data) %>% summary()

lmer(dbp ~ sample_type * treatment + log10 (Final_reads) + (1|original_sample), data = sample_data) %>% summary() %>%
        .$coefficients %>%
        kbl %>%
        save_kable(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/TableS8.html", self_contained = T)





# Table 2 - beta diversity permanova result -------------------------------


# Table for beta diversity


bray_perm_ <- vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + treatment + log10(Final_reads) + original_sample,
                             data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 

bray_perm_inter <- vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type * treatment + log10(Final_reads) + original_sample,
                                  data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 

bray_perm_ns <- vegan::adonis2(distance(subset_samples(phyloseq_rel_nz, sample_type == "nasal_swab"), method="bray") ~ lypma + benzonase + host_zero + molysis + qiaamp + log10(Final_reads) + original_sample,
                               data = subset_samples(phyloseq_rel_nz, sample_type == "nasal_swab") %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 

bray_perm_bal  <- vegan::adonis2(distance(subset_samples(phyloseq_rel_nz, sample_type == "BAL"), method="bray") ~  lypma + benzonase + host_zero + molysis + qiaamp + log10(Final_reads) + original_sample,
                                 data = subset_samples(phyloseq_rel_nz, sample_type == "BAL") %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 

bray_perm_spt <- vegan::adonis2(distance(subset_samples(phyloseq_rel_nz, sample_type == "Sputum"), method="bray") ~ lypma + benzonase + host_zero + molysis + qiaamp + log10(Final_reads) + original_sample,
                                data = subset_samples(phyloseq_rel_nz, sample_type == "Sputum") %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 


bray_perm <- vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + log10(Final_reads) + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                            data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 


bray_perm %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
        mutate(row.names = case_when(row.names == "sample_type" ~ 'Sample type',
                                     row.names == "lypma" ~ 'lyPMA',
                                     row.names == "benzonase" ~ 'Benzonase',
                                     row.names == "host_zero" ~ 'Host zero',
                                     row.names == "molysis" ~ 'Molysis',
                                     row.names == "qiaamp" ~ 'QIAamp',
                                     row.names == "original_sample" ~ 'Subject',
                                     row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                     row.names == "Residual" ~ 'Residual',
                                     row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
        round(4) %>% select(c(`R2`, `Pr(>F)`))

bray_perm_ %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
        mutate(row.names = case_when(row.names == "sample_type" ~ 'Sample type',
                                     row.names == "treatment" ~ 'Treatment',
                                     row.names == "original_sample" ~ 'Subject',
                                     row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                     row.names == "Residual" ~ 'Residual',
                                     row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
        round(4) %>% select(c(`R2`, `Pr(>F)`))
bray_perm_inter
tabe_inter <- bray_perm_inter %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
        mutate(row.names = case_when(row.names == "sample_type" ~ 'Sample type',
                                     row.names == "treatment" ~ 'Treatment',
                                     row.names == "original_sample" ~ 'Subject',
                                     row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                     row.names == "sample_type:treatment" ~ 'Sample type X treatment',
                                     row.names == "Residual" ~ 'Residual',
                                     row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
        round(4) %>% select(c(`R2`, `Pr(>F)`))


bray_perm_ns

table_ns <- bray_perm_ns %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
        mutate(row.names = case_when(row.names == "lypma" ~ 'lyPMA',
                                     row.names == "benzonase" ~ 'Benzonase',
                                     row.names == "host_zero" ~ 'Host zero',
                                     row.names == "molysis" ~ 'Molysis',
                                     row.names == "qiaamp" ~ 'QIAamp',
                                     row.names == "original_sample" ~ 'Subject id',
                                     row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                     row.names == "Residual" ~ 'Residual',
                                     row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
        round(4) %>% select(c(`R2`, `Pr(>F)`))

table_bal <- bray_perm_bal %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
        mutate(row.names = case_when(row.names == "lypma" ~ 'lyPMA',
                                     row.names == "benzonase" ~ 'Benzonase',
                                     row.names == "host_zero" ~ 'Host zero',
                                     row.names == "molysis" ~ 'Molysis',
                                     row.names == "qiaamp" ~ 'QIAamp',
                                     row.names == "original_sample" ~ 'Subject id',
                                     row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                     row.names == "Residual" ~ 'Residual',
                                     row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
        round(3) %>% select(c(`R2`, `Pr(>F)`))

table_spt <- bray_perm_spt %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
        mutate(row.names = case_when(row.names == "lypma" ~ 'lyPMA',
                                     row.names == "benzonase" ~ 'Benzonase',
                                     row.names == "host_zero" ~ 'Host zero',
                                     row.names == "molysis" ~ 'Molysis',
                                     row.names == "qiaamp" ~ 'QIAamp',
                                     row.names == "original_sample" ~ 'Subject id',
                                     row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                     row.names == "Residual" ~ 'Residual',
                                     row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
        round(3) %>% select(c(`R2`, `Pr(>F)`))

table_all <- bray_perm %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
        mutate(row.names = case_when(row.names == "sample_type" ~ 'Sample type',
                                     row.names == "lypma" ~ 'lyPMA',
                                     row.names == "benzonase" ~ 'Benzonase',
                                     row.names == "host_zero" ~ 'Host zero',
                                     row.names == "molysis" ~ 'Molysis',
                                     row.names == "qiaamp" ~ 'QIAamp',
                                     row.names == "original_sample" ~ 'Subject id',
                                     row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                     row.names == "Residual" ~ 'Residual',
                                     row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
        round(4) %>% select(c(`R2`, `Pr(>F)`))

#Figure 6 models

#Observed species

subset(sample_data(phyloseq_path), sample_data(phyloseq_path)$sample_type %in% c("Sputum", "nasal_swab", "BAL")) %>%
        sample_data %>% data.frame() %>% lmer(S.obs ~ sample_type + treatment + (1|original_sample), data = .) %>% summary
#Shannon
subset(sample_data(phyloseq_path), sample_data(phyloseq_path)$sample_type %in% c("Sputum", "nasal_swab", "BAL")) %>%
        sample_data %>% data.frame() %>% lmer(data_shannon ~ sample_type + treatment + (1|original_sample), data = .) %>% summary
#Invsimpson
subset(sample_data(phyloseq_path), sample_data(phyloseq_path)$sample_type %in% c("Sputum", "nasal_swab", "BAL")) %>% sample_data %>%
        data.frame() %>% lmer(data_invsimpson ~ sample_type + treatment + (1|original_sample), data = .) %>% summary
#Berger-Parker
subset(sample_data(phyloseq_path), sample_data(phyloseq_path)$sample_type %in% c("Sputum", "nasal_swab", "BAL")) %>% sample_data %>%
        data.frame() %>% lmer(dbp ~ sample_type + treatment + (1|original_sample), data = .) %>% summary
#Beta - perm


                        bray_perm_path_ <- vegan::adonis2(distance(phyloseq_path_rel_nz, method="bray") ~ sample_type + treatment + log10(Final_reads) + original_sample,
                                                     data = phyloseq_path_rel_nz %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 
                        
                        bray_perm_path <- vegan::adonis2(distance(phyloseq_path_rel_nz, method="bray") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp + log10(Final_reads) + original_sample,
                                                            data = phyloseq_path_rel_nz %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 
                        
                        
                        bray_perm_path_inter <- vegan::adonis2(distance(phyloseq_path_rel_nz, method="bray") ~ sample_type * treatment + log10(Final_reads) + original_sample,
                                                          data = phyloseq_path_rel_nz %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 
                        
                        bray_perm_path_ns <- vegan::adonis2(distance(subset_samples(phyloseq_path_rel_nz, sample_type == "nasal_swab"), method="bray") ~ lypma + benzonase + host_zero + molysis + qiaamp + log10(Final_reads) + original_sample,
                                                       data = subset_samples(phyloseq_path_rel_nz, sample_type == "nasal_swab") %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 
                        
                        bray_perm_path_bal  <- vegan::adonis2(distance(subset_samples(phyloseq_path_rel_nz, sample_type == "BAL"), method="bray") ~  lypma + benzonase + host_zero + molysis + qiaamp + log10(Final_reads) + original_sample,
                                                         data = subset_samples(phyloseq_path_rel_nz, sample_type == "BAL") %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 
                        
                        bray_perm_path_spt <- vegan::adonis2(distance(subset_samples(phyloseq_rel_nz, sample_type == "Sputum"), method="bray") ~ lypma + benzonase + host_zero + molysis + qiaamp + log10(Final_reads) + original_sample,
                                                        data = subset_samples(phyloseq_path_rel_nz, sample_type == "Sputum") %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 
                        
                        
                        bray_perm <- vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + log10(Final_reads) + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                    data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), permutations = 10000) 
                        
                        
                        
                        bray_perm_path %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
                                mutate(row.names = case_when(row.names == "sample_type" ~ 'Sample type',
                                                             row.names == "lypma" ~ 'lyPMA',
                                                             row.names == "benzonase" ~ 'Benzonase',
                                                             row.names == "host_zero" ~ 'Host zero',
                                                             row.names == "molysis" ~ 'Molysis',
                                                             row.names == "qiaamp" ~ 'QIAamp',
                                                             row.names == "original_sample" ~ 'Subject',
                                                             row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                                             row.names == "Residual" ~ 'Residual',
                                                             row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
                                round(4) %>% select(c(`R2`, `Pr(>F)`))
                        
                        bray_perm_path %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
                                mutate(row.names = case_when(row.names == "sample_type" ~ 'Sample type',
                                                             row.names == "treatment" ~ 'Treatment',
                                                             row.names == "original_sample" ~ 'Subject',
                                                             row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                                             row.names == "Residual" ~ 'Residual',
                                                             row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
                                round(4) %>% select(c(`R2`, `Pr(>F)`))
                        bray_perm_path_inter
                        tabe_path_inter <- bray_perm_path_inter %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
                                mutate(row.names = case_when(row.names == "sample_type" ~ 'Sample type',
                                                             row.names == "treatment" ~ 'Treatment',
                                                             row.names == "original_sample" ~ 'Subject',
                                                             row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                                             row.names == "sample_type:treatment" ~ 'Sample type X treatment',
                                                             row.names == "Residual" ~ 'Residual',
                                                             row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
                                round(4) %>% select(c(`R2`, `Pr(>F)`))
                        
                        
                        table_path_ns <- bray_perm_path_ns %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
                                mutate(row.names = case_when(row.names == "lypma" ~ 'lyPMA',
                                                             row.names == "benzonase" ~ 'Benzonase',
                                                             row.names == "host_zero" ~ 'Host zero',
                                                             row.names == "molysis" ~ 'Molysis',
                                                             row.names == "qiaamp" ~ 'QIAamp',
                                                             row.names == "original_sample" ~ 'Subject id',
                                                             row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                                             row.names == "Residual" ~ 'Residual',
                                                             row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
                                round(4) %>% select(c(`R2`, `Pr(>F)`))
                        
                        table_path_bal <- bray_perm_path_bal %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
                                mutate(row.names = case_when(row.names == "lypma" ~ 'lyPMA',
                                                             row.names == "benzonase" ~ 'Benzonase',
                                                             row.names == "host_zero" ~ 'Host zero',
                                                             row.names == "molysis" ~ 'Molysis',
                                                             row.names == "qiaamp" ~ 'QIAamp',
                                                             row.names == "original_sample" ~ 'Subject id',
                                                             row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                                             row.names == "Residual" ~ 'Residual',
                                                             row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
                                round(3) %>% select(c(`R2`, `Pr(>F)`))
                        
                        table_path_spt <- bray_perm_path_spt %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
                                mutate(row.names = case_when(row.names == "lypma" ~ 'lyPMA',
                                                             row.names == "benzonase" ~ 'Benzonase',
                                                             row.names == "host_zero" ~ 'Host zero',
                                                             row.names == "molysis" ~ 'Molysis',
                                                             row.names == "qiaamp" ~ 'QIAamp',
                                                             row.names == "original_sample" ~ 'Subject id',
                                                             row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                                             row.names == "Residual" ~ 'Residual',
                                                             row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
                                round(3) %>% select(c(`R2`, `Pr(>F)`))
                        
                        table_path_all <- bray_perm %>% data.frame(check.names = F) %>% rownames_to_column('row.names') %>% 
                                mutate(row.names = case_when(row.names == "sample_type" ~ 'Sample type',
                                                             row.names == "lypma" ~ 'lyPMA',
                                                             row.names == "benzonase" ~ 'Benzonase',
                                                             row.names == "host_zero" ~ 'Host zero',
                                                             row.names == "molysis" ~ 'Molysis',
                                                             row.names == "qiaamp" ~ 'QIAamp',
                                                             row.names == "original_sample" ~ 'Subject id',
                                                             row.names == "log10(Final_reads)" ~ 'log10(Final reads)',
                                                             row.names == "Residual" ~ 'Residual',
                                                             row.names == "Total" ~ 'Total')) %>% column_to_rownames('row.names') %>% 
                                round(4) %>% select(c(`R2`, `Pr(>F)`))
                        
                        
                        
#Beta-distance by group
#distances of betadiversity - boxplots
                                                        library(harrietr)
                                                        
                                                        bray_dist_long_path <- distance(phyloseq_path_rel_nz, method="bray") %>% as.matrix() %>% melt_dist() #making long data of distance matrices
                                                        #Adding sample type and treatment name. 
                                                        #this can be also done by merging metadata into the `bray_dist_long`
                                                        names <- data.frame(str_split_fixed(bray_dist_long_path$iso1, "_", 3))
                                                        names2 <- data.frame(str_split_fixed(bray_dist_long_path$iso2, "_", 3))
                                                        bray_dist_long_path$sample_id_1 <- paste(names$X1, names$X2, sep = "_")
                                                        bray_dist_long_path$method_1 <- ifelse(grepl("control", bray_dist_long_path$iso1),"control", 
                                                                                               ifelse(grepl("lyPMA", bray_dist_long_path$iso1),"lypma", 
                                                                                                      ifelse(grepl("benzonase", bray_dist_long_path$iso1),"benzonase", 
                                                                                                             ifelse(grepl("host", bray_dist_long_path$iso1),"host_zero", 
                                                                                                                    ifelse(grepl("qia", bray_dist_long_path$iso1),"qiaamp", 
                                                                                                                           ifelse(grepl("moly", bray_dist_long_path$iso1),"molysis", 
                                                                                                                                  NA))))))
                                                        #Adding data for iso 2 also should be done
                                                        bray_dist_long_path$sample_id_2 <- paste(names2$X1, names2$X2, sep = "_")
                                                        bray_dist_long_path$method_2 <-ifelse(grepl("control", bray_dist_long_path$iso2),"control", 
                                                                                              ifelse(grepl("lyPMA", bray_dist_long_path$iso2),"lypma", 
                                                                                                     ifelse(grepl("benzonase", bray_dist_long_path$iso2),"benzonase", 
                                                                                                            ifelse(grepl("host", bray_dist_long_path$iso2),"host_zero", 
                                                                                                                   ifelse(grepl("qia", bray_dist_long_path$iso2),"qiaamp", 
                                                                                                                          ifelse(grepl("moly", bray_dist_long_path$iso2),"molysis", 
                                                                                                                                 NA))))))
                                                        #subsetting distances of my interest
                                                        path_bray_dist_long_within_sampleid <- subset(bray_dist_long_path, bray_dist_long_path$sample_id_1 == bray_dist_long_path$sample_id_2)
                                                        path_bray_dist_long_within_sampleid_from_control <- subset(path_bray_dist_long_within_sampleid, path_bray_dist_long_within_sampleid$method_1 == "control" | path_bray_dist_long_within_sampleid$method_2 == "control" )
                                                        path_bray_dist_long_within_sampleid_from_control$treatment <- path_bray_dist_long_within_sampleid_from_control$method_1
                                                        path_bray_dist_long_within_sampleid_from_control$treatment <- ifelse(path_bray_dist_long_within_sampleid_from_control$treatment == "control", path_bray_dist_long_within_sampleid_from_control$method_2, path_bray_dist_long_within_sampleid_from_control$treatment)
                                                        path_bray_dist_long_within_sampleid_from_control$sample_type <- ifelse(grepl("NS", path_bray_dist_long_within_sampleid_from_control$iso1), "nasal_swab",
                                                                                                                               ifelse(grepl("CFB", path_bray_dist_long_within_sampleid_from_control$iso1), "Sputum",
                                                                                                                                      ifelse(grepl("BAL", path_bray_dist_long_within_sampleid_from_control$iso1), "BAL", NA)))

path_bray_dist_long_within_sampleid_from_control %>% lm(dist ~ sample_type + treatment, data = .) %>% summary



# Summary table -----------------------------------------------------------



## Final results summary {.tabset}

### Sequencing results

matrix(nrow=3,ncol=5) %>% data.frame() %>% rename(lyPMA = X1, Benzonase = X2, `Host zero` = X3, Molysis = X4, QIAamp = X5) %>%
        rownames_to_column("x") %>% mutate(x = c("BAL", "Nasal swab", "Sputum"),
                                           lyPMA = c("No increase in final reads<br>Taxa beta changed",
                                                     "No increase in final reads",
                                                     "No increase in final reads"),
                                           Benzonase = c("No decrease in host %",
                                                         "No decrease in host %",
                                                         "No decrease in host %"),
                                           `Host zero` = c(NA,
                                                           NA,
                                                           NA),
                                           Molysis = c("No decrease in host %",
                                                       "High cahnge of failure in library pep",
                                                       NA),
                                           QIAamp = c("No decrease in host %",
                                                      NA,
                                                      "No decrease in host %")) %>% column_to_rownames("x") %>%
        kbl(format = "html", caption = "Table of issues of each treatment method") %>%
        kable_styling(full_width = 0, html_font = "serif")


matrix(nrow=3,ncol=5) %>% data.frame() %>% rename(lyPMA = X1, Benzonase = X2, `Host zero` = X3, Molysis = X4, QIAamp = X5) %>%
        rownames_to_column("x") %>% mutate(x = c("BAL", "Nasal swab", "Sputum"),
                                           lyPMA = c(NA,
                                                     "Shannon +"),
                                           Benzonase = c(NA,
                                                         NA,
                                                         "Richness + Shannon + InvSimp +"),
                                           `Host zero` = c(NA,
                                                           "Richness + Shannon + InvSimp + BPI -",
                                                           NA),
                                           Molysis = c(NA,
                                                       "Richness + Shannon + InvSimp + BPI -",
                                                       "Beta changed"),
                                           QIAamp = c("Beta changed",
                                                      NA,
                                                      "Beta  changed")) %>% column_to_rownames("x") %>%
        kbl(format = "html", caption = "Table of community changes induced by each treatment method") %>%
        kable_styling(full_width = 0, html_font = "serif")


matrix(nrow=3,ncol=5) %>% data.frame() %>% rename(lyPMA = X1, Benzonase = X2, `Host zero` = X3, Molysis = X4, QIAamp = X5) %>%
        rownames_to_column("x") %>% mutate(x = c("BAL", "Nasal swab", "Sputum"),
                                           lyPMA = c(NA,
                                                     NA,
                                                     "Shannon +"),
                                           Benzonase = c(NA,
                                                         NA,
                                                         "Shannon +"),
                                           `Host zero` = c(NA,
                                                           "Richness +",
                                                           "Shannon +"),
                                           Molysis = c(NA,
                                                       "Richness + InvSimp + BPI +",
                                                       "Shannon +"),
                                           QIAamp = c(NA,
                                                      "Richness + Shannon +",
                                                      "Shannon +")) %>% column_to_rownames("x") %>%
        kbl(format = "html", caption = "Table of functional diversity changes induced by each treatment method") %>%
        kable_styling(full_width = 0, html_font = "serif")


matrix(nrow=3,ncol=5) %>% data.frame() %>% rename(lyPMA = X1, Benzonase = X2, `Host zero` = X3, Molysis = X4, QIAamp = X5) %>%
        rownames_to_column("x") %>% mutate(x = c("BAL", "Nasal swab", "Sputum"),
                                           lyPMA = c("Listeria",
                                                     "Listeria",
                                                     "Listeria, Candida, Corynebacterium"),
                                           Benzonase = c("Listeria",
                                                         "Listeria",
                                                         "Listeria, Candida, Corynebacterium"),
                                           `Host zero` = c("Listeria",
                                                           "Listeria",
                                                           "Listeria, Candida, Corynebacterium"),
                                           Molysis = c("Streptococcaceae, Listeria",
                                                       "Streptococcaceae, Listeria",
                                                       "Streptococcaceae, Listeria, Candida, Corynebacterium"),
                                           QIAamp = c("Listeria",
                                                      "Listeria",
                                                      "Listeria, Candida, Corynebacterium")) %>% column_to_rownames("x") %>%
        kbl(format = "html", caption = "Table of potential contaminants identified by decontam and DA analysis") %>%
        kable_styling(full_width = 0, html_font = "serif") %>%
        column_spec(2:6, italic = T) #%>%

