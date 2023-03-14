#Required packages
library(readxl)
library(ggplot2)
library(phyloseq)
library(vegan)
library(tidyverse)
library(microbiome)
library(ggpubr)
library(viridis)
library(decontam)
library(lme4)
library(lmerTest)
library(writexl)
library(harrietr)
library(Maaslin2)

rm(list = ls())

# Loading files -----------------------------------------------------------
        #loading tidy phyloseq object
        phyloseq <- read_rds("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221129_MGK_host_tidy_tax.rds")
        #sample data loading
        sample_data <- sample_data(phyloseq)
        
        
# Initail QC --------------------------------------------------------------
        #Quesetions - QC
        
                                
                                hist(sample_data$Final_reads)
                                #overview - FInal reads
                                
                                hist(sample_data$total_reads)
                                #sample_data$treatment <- factor(sample_data$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                #sample_data$sample_type <- factor(sample_data$sample_type, levels = c("nasal_swab", "sputum", "bal", "pos_control", "neg_control"))
                                
        #Q0. How many samples failed in sequencing
                                
                                sample_data$total_reads
                                hist((log10(sample_data$total_reads + 1)),100) # one somaple showed 0 reads
                                sample_data %>% data.frame %>% filter(total_reads == 0) # two BAL sampels showed 0 total reads
                                
                                #how were the samples failed in library prep?
                                
                                #BAL
                                lib_fail_list <- c("NS_6_b_lyPMA", "NS_8_b_lyPMA", "NS_17_b_lyPMA", "NS_21_b_host_zero", "NS_22_b_host_zero", "NS_19_d_molysis", "NS_21_d_molysis", "NS_22_d_molysis", "NS_23_d_molysis", "NS_26_b_lyPMA", "BAL_073_b_host_zero", "BAL_073_c_molysis", "BAL_078_d_lyPMA")
                                sample_data(phyloseq)$lib_failed <- row.names(sample_data(phyloseq)) %in% lib_fail_list

                                #I have no idea how Baylor decisded those sample as "library prep failed samples
        
        # Q0. How was the positive control and negative control? were there any contaminants? 
                                sample_data
                                
                                phyloseq_control <- subset_samples(phyloseq, grepl("2022", Row.names))
                                phyloseq_sample <- subset_samples(phyloseq, !grepl("2022", Row.names))
                                        sample_data(phyloseq_control)$is.neg <- 1
                                        sample_data(phyloseq_sample)$is.neg <- 0
                                        
                                
                                phyloseq_control <- prune_taxa(taxa_sums(phyloseq_control) != 0, phyloseq_control) 
                                        #Mock community data
                                        
                                        zymo_mock<-read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_sicas2/data_raw/DAR_20210929_zymo_mock_data.xlsx")
                                        zymo_mock<-data.frame(zymo_mock, row.names = T)
                                        names(zymo_mock) <- "zymo_mock"
                                        physeq_mock <- phyloseq(otu_table(zymo_mock, taxa_are_rows = T), tax_table(phyloseq_control))
                                        mock_sample_data <- rbind(c("zymo_mock", "mock")) %>% data.frame(.) 
                                        names(mock_sample_data) <- c("Row.names", "sample_type")
                                        row.names(mock_sample_data) <- mock_sample_data$Row.names
                                        physeq_mock <- merge_phyloseq(sample_data(mock_sample_data), physeq_mock)
                                        #Mock community & controls phyloseq
                                        controls_phyloseq <- merge_phyloseq(physeq_mock, phyloseq_control)
                                        #Mock community - relative abundance phyloseq
                                        phyloseq_control_rel <- transform_sample_counts(controls_phyloseq, function(x){x / sum(x)})
                                #Species richness of each control groups
                                S.obs <- rowSums(t(otu_table(phyloseq_control)) != 0)
                                ps.gen <- phyloseq::tax_glom(phyloseq_control_rel, "Genus", NArm = TRUE)
                                
                                sample_data(phyloseq_control_rel)$detail <- c("Mock, 10 spp.", "Pos., 41 spp.", "Neg., 3 spp.")
                                sample_data(phyloseq_control_rel)$detail <- factor(sample_data(phyloseq_control_rel)$detail, levels = c("Mock, 10 spp.", "Pos., 41 spp.", "Neg., 3 spp."))
                                
                                plot_bar(phyloseq_control_rel, fill="Species") + 
                                        facet_wrap (~ detail, scales= "free_x", nrow=1) +
                                        ylab("Relative abudnacne") +
                                        title("Species richenss") +
                                        theme_classic(base_size = 11, base_family = "serif")
                                sample_data(ps.gen)$detail <- c("Mock, 10 Genus", "Pos., 28 Genus", "Neg., 3 Genus")
                                sample_data(ps.gen)$detail <- factor(sample_data(ps.gen)$detail, levels = c("Mock, 10 Genus", "Pos., 28 Genus", "Neg., 3 Genus"))
                                plot_bar(ps.gen, fill="Genus") + 
                                        facet_wrap (~ detail, scales= "free_x", nrow=1) +
                                        ylab("Relative abudnacne") +
                                        title("Species richenss") +
                                        theme_classic(base_size = 11, base_family = "serif")
                                #there could be opportunistic pathogens...
                                
                                #Mo-mock control data
                                        ps.gen_nomock <- subset_taxa(ps.gen, otu_table(ps.gen)[,1] == 0)
                                        phyloseq_control_nomock <- subset_taxa(phyloseq_control, (tax_table(phyloseq_control) %>% data.frame() %>% .$Genus) %in% (tax_table(ps.gen_nomock) %>% data.frame() %>% .$Genus))
                                        sample_data(phyloseq_control_nomock)
                                        phyloseq_control_nomock <- subset_samples(phyloseq_control_nomock, sample_type != "mock")
                                        phyloseq_control_nomock

# Decontam package --------------------------------------------------------

                
                        #running decontam
                                        
                                phyloseq_nomock <- merge_phyloseq(phyloseq_sample, phyloseq_control_nomock)
                                
                                        #species level decontama
                                        #Control phylose object
                                        #decontam - frequency method
                                        
                                        sample_data(phyloseq_nomock)$total_dna_ngul <- sample_data(phyloseq_nomock)$DNA_bac_nondil + sample_data(phyloseq_nomock)$DNA_host_nondil
                                        phy_decontam <- subset_samples(phyloseq_nomock, total_reads != 0)
                                        sample_data(phy_decontam)$DNA_bac_nondil <- ifelse(sample_data(phy_decontam)$DNA_bac_nondil %>% is.na(),
                                                                                           0.01, 
                                                                                           sample_data(phy_decontam)$DNA_bac_nondil)
                                        sample_data(phy_decontam)$total_dna_ngul <- ifelse(sample_data(phy_decontam)$total_dna_ngul %>% is.na(),
                                                                                           0.01, 
                                                                                           sample_data(phy_decontam)$total_dna_ngul)
                                        phy_decontam_freq <- subset_samples(phy_decontam, !is.na(DNA_bac_nondil))
                                        contamdf.freq <- isContaminant(phy_decontam, method="frequency", conc="DNA_bac_nondil")
                                        contamdf.freq$contaminant %>% sum()
                                        which(contamdf.freq$contaminant)#list of contam
                                        
                                        plot_frequency(phy_decontam, taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))], conc="DNA_bac_nondil") + 
                                                xlab("DNA concentration (qPCR; ) [ng / uL]")        
                                        
                                        #decontam - prevalence method
                                        sample_data(phy_decontam)$is.neg <- sample_data(phy_decontam)$is.neg == 1
                                        
                                        contamdf.prev <- isContaminant(phy_decontam, method="prevalence", neg="is.neg", threshold = 0.5)
                                        
                                        table(contamdf.prev$contaminant)
                                        ps.pa <- transform_sample_counts(phy_decontam, function(abund) 1*(abund>0))
                                        ps.pa.neg <- prune_samples(sample_data(ps.pa)$is.neg == 1, ps.pa)
                                        ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.neg == 0, ps.pa)
                                        # Make data.frame of prevalence in positive and negative samples
                                        df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                                                            contaminant=contamdf.prev$contaminant)
                                        ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
                                                xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
                                        
                                        #decontam - combined method
                                        phy_decontam_comb <- subset_samples(phy_decontam, total_reads != 0)
                                        phy_decontam_comb <- subset_samples(phy_decontam_comb, !is.na(DNA_bac_nondil))
                                        contamdf.comb <- isContaminant(phy_decontam_comb, method="combined", neg="is.neg", conc = "DNA_bac_nondil")
                                        table(contamdf.comb$contaminant)
                                        #using combined method, 6 potential candidates were made.
                                        
                                        otu_table <- otu_table(phyloseq)
                                        otu_table_decontam <- subset(otu_table, !(row.names(otu_table) %in% row.names(subset(contamdf.comb, contamdf.comb$contaminant))))
                                        
                                        phyloseq_decontam <- merge_phyloseq(sample_data(phyloseq), tax_table(phyloseq), otu_table(otu_table_decontam, taxa_are_rows = T))
                                        phyloseq_decontam
                                        #Getting info of both methods
                                        contamdf.freq_true <- subset(contamdf.freq, contamdf.freq$contaminant)
                                        contamdf.prev_true <- subset(contamdf.prev, contamdf.prev$contaminant)
                                        contamdf.comb_true <- subset(contamdf.comb, contamdf.comb$contaminant)
                                        
                                        otu_table_decontam2 <- subset(otu_table, !(row.names(otu_table) %in% row.names(contamdf.prev_true)) & !(row.names(otu_table) %in% row.names(contamdf.freq_true)))
                                        phyloseq_decontam2 <- merge_phyloseq(sample_data(phyloseq), tax_table(phyloseq), otu_table(otu_table_decontam, taxa_are_rows = T))

# Decontam - sensitivity analysis -----------------------------------------

                                        
                                        #decontam-sensitivity analysis
                                        
                                        
                                                        phy_decontam <- subset_samples(phy_decontam, !grepl("control", Row.names))
                                                        contamdf.freqs <- isContaminant(phy_decontam, method="frequency", conc="DNA_bac_nondil")
                                                        contamdf.freqs$contaminant %>% sum()
                                                        which(contamdf.freqs$contaminant)#list of contam
                                                        plot_frequency(phy_decontam, taxa_names(phy_decontam_freq)[c(which(contamdf.freqs$contaminant))], conc="DNA_bac_nondil") + 
                                                                xlab("DNA concentration (pigogreen fluorescent intensity) [ng / uL]")        
                                                        
                                                        #decontam - prevalence method
                                                        sample_data(phy_decontam)$is.neg <- sample_data(phy_decontam)$is.neg == 1
                                                        
                                                        contamdf.prev_s <- isContaminant(phy_decontam, method="prevalence", neg="is.neg", threshold = 0.5)
                                                        
                                                        table(contamdf.prev_s$contaminant)
                                                        ps.pa <- transform_sample_counts(phy_decontam, function(abund) 1*(abund>0))
                                                        ps.pa.neg <- prune_samples(sample_data(ps.pa)$is.neg == 1, ps.pa)
                                                        ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.neg == 0, ps.pa)
                                                        # Make data.frame of prevalence in positive and negative samples
                                                        df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                                                                            contaminant=contamdf.prev$contaminant)
                                                        ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
                                                                xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
                                                        
                                                        #decontam - combined method
                                                        phy_decontam_comb_s <- subset_samples(phy_decontam, total_reads != 0)
                                                        phy_decontam_comb_s <- subset_samples(phy_decontam_comb_s, !is.na(DNA_bac_nondil))
                                                        contamdf.comb <- isContaminant(phy_decontam_comb, method="combined", neg="is.neg", conc = "DNA_bac_nondil")
                                                        table(contamdf.comb$contaminant)
                                                        #using combined method, 6 potential candidates were made.
                                                        
                                                        otu_table_decontam <- subset(otu_table, !(row.names(otu_table) %in% row.names(subset(contamdf.comb, contamdf.comb$contaminant))))
                                                        phyloseq_decontam <- merge_phyloseq(sample_data(phyloseq), tax_table(phyloseq), otu_table(otu_table_decontam, taxa_are_rows = T))
                                                        
                                                        #Getting info of both methods
                                                        
                                                        contamdf.freq_true_s <- subset(contamdf.freqs, contamdf.freqs$contaminant)
                                                        contamdf.prev_true_s <- subset(contamdf.prev_s, contamdf.prev_s$contaminant)
                                                        contamdf.comb_true_s <- subset(contamdf.comb, contamdf.comb$contaminant)
                                                        

# Decontam - for each sample type -----------------------------------------

                                                        
                                #Running decontam for differnet sample types
                                                        
                                                        phy_decontam %>% sample_data
                                              
                                                        #species level decontama
                                                        #decontam - frequency method
                                                        
                                                        sample_data(phyloseq_nomock)$total_dna_ngul <- sample_data(phyloseq_nomock)$DNA_bac_nondil + sample_data(phyloseq_nomock)$DNA_host_nondil
                                                        phy_decontam <- subset_samples(phy_decontam, total_reads != 0)
                                                        phy_decontam_ns <- subset_samples(phy_decontam, sample_type == "nasal_swab")
                                                        phy_decontam_ns_prev <- subset_samples(phy_decontam, sample_type != "BAL" & sample_type != "Sputum" )
                                                        phy_decontam_bal <- subset_samples(phy_decontam, sample_type == "BAL")
                                                        phy_decontam_bal_prev <- subset_samples(phy_decontam, sample_type != "nasal_swab" & sample_type != "Sputum")
                                                        phy_decontam_spt <- subset_samples(phy_decontam, sample_type == "Sputum")
                                                        phy_decontam_spt_prev <- subset_samples(phy_decontam, sample_type != "BAL" & sample_type != "nasal_swab")
                                                        
                                                        contamdf.freq$contaminant %>% sum()
                                                        which(contamdf.freq$contaminant)#list of contam
                                                        
                                                        #negative control
                                                        phyloseq_control %>% otu_table()
                                                        phyloseq %>% otu_table() %>% subset(.,grepl("Sutter", row.names(.)))
                                                        phyloseq %>% otu_table() %>% subset(.,grepl("Sutter", row.names(.)))
                                        #frequency method
                                                #All sample type
                                                        plot_frequency(phy_decontam, taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))], conc="DNA_bac_nondil") + 
                                                                xlab("Bacteiral DNA concentration (qPCR) [ng / uL]") 

                                                        
                                                #Nasal swab                        
                                                        #Bacterial DNA
                                                        contamdf.freq <- isContaminant(phy_decontam_ns, method="frequency", conc="DNA_bac_nondil")
                                                        plot_frequency(phy_decontam_ns, taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))], conc="DNA_bac_nondil") + 
                                                                xlab("Bacteiral DNA concentration (qPCR) [ng / uL]")        
                                                        #Total DNA
                                                        contamdf.freq <- isContaminant(phy_decontam_ns, method="frequency", conc="total_dna_ngul")
                                                        plot_frequency(phy_decontam_ns, taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))], conc="DNA_bac_nondil") + 
                                                                xlab("Total DNA concentration (qPCR) [ng / uL]")        
                                                        taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))]
                                                        
                                                        
                                                #BAL
                                                        #Bacterial DNA
                                                        contamdf.freq <- isContaminant(phy_decontam_bal, method="frequency", conc="DNA_bac_nondil")
                                                        plot_frequency(phy_decontam_bal, taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))], conc="DNA_bac_nondil") + 
                                                                xlab("Bacteiral DNA concentration (qPCR) [ng / uL]")        
                                                        #Total DNA
                                                        contamdf.freq <- isContaminant(phy_decontam_bal, method="frequency", conc="total_dna_ngul")
                                                        plot_frequency(phy_decontam_bal, taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))], conc="DNA_bac_nondil") + 
                                                                xlab("Total DNA concentration (qPCR) [ng / uL]")        
                                                        taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))]
                                                        
                                                #Sputum
                                                        #Bacterial DNA
                                                        contamdf.freq <- isContaminant(phy_decontam_spt, method="frequency", conc="DNA_bac_nondil")
                                                        plot_frequency(phy_decontam_spt, taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))], conc="DNA_bac_nondil") + 
                                                                xlab("Bacteiral DNA concentration (qPCR) [ng / uL]")        
                                                        #Total DNA
                                                        
                                                        contamdf.freq <- isContaminant(phy_decontam_spt, method="frequency", conc="total_dna_ngul")
                                                        plot_frequency(phy_decontam_spt, taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))], conc="DNA_bac_nondil") + 
                                                                xlab("Total DNA concentration (qPCR) [ng / uL]")        
                                                        taxa_names(phy_decontam_freq)[c(which(contamdf.freq$contaminant))]
                                                         
                                        #decontam - prevalence method
                                                        decontam_prev <- function (data) {
                                                                sample_data(data)$is.neg <- sample_data(data)$is.neg == 1
                                                                contamdf.prev <- isContaminant(data, method="prevalence", neg="is.neg", threshold = 0.5)
                                                                table(contamdf.prev$contaminant)
                                                                ps.pa <- transform_sample_counts(data, function(abund) 1*(abund>0))
                                                                ps.pa.neg <- prune_samples(sample_data(ps.pa)$is.neg == 1, ps.pa)
                                                                ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.neg == 0, ps.pa)
                                                                # Make data.frame of prevalence in positive and negative samples
                                                                df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                                                                                    contaminant=contamdf.prev$contaminant)
                                                                ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
                                                                        xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
                                                         #       subset(df.pa, df.pa$contaminant)
                                                        }
                                                        decontam_prev(phy_decontam_bal_prev)
                                                        decontam_prev(phy_decontam_ns_prev)
                                                        decontam_prev(phy_decontam_spt_prev)
                                                #decontam - combined method
                                                        phy_decontam_comb <- subset_samples(phy_decontam, total_reads != 0)
                                                        phy_decontam_comb <- subset_samples(phy_decontam_comb, !is.na(DNA_bac_nondil))
                                                        contamdf.comb <- isContaminant(phy_decontam_comb, method="combined", neg="is.neg", conc = "DNA_bac_nondil")
                                                        table(contamdf.comb$contaminant)
                                                        
                                                        

# Aanalysis 1. qPCR result vs seqeucning result ---------------------------


        #Q1. is qPCR result consistent with sequencing result
                                
                                ggplot(sample_data, aes(x = proportion, y = sequencing_host_prop)) +
                                        geom_point () + 
                                        xlab("Host DNA proportion by qPCR") +
                                        ylab("(Host DNA / mapped bacterial DNA) by seqeuncing")
                                
                                #is the total reads the same with final reads?
                                ggplot(sample_data, aes(x = log10(Final_reads), y = log10(total_reads))) +
                                        geom_point () + 
                                        xlab("Final reads in Kneaddata") +
                                        ylab("Sum of reads of OTU table")
                                
                                
                                lm_eqn <- function(df, a, b){
                                        m <- lm(a~ b, df);
                                        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                                         list(a = format(unname(coef(m)[1]), digits = 2),
                                                              b = format(unname(coef(m)[2]), digits = 2),
                                                              r2 = format(summary(m)$r.squared, digits = 3)))
                                        as.character(as.expression(eq));
                                }
                                sample_data <- data.frame(sample_data)
                                lm_eqn(sample_data, sample_data$Host_mapped / (sample_data$Host_mapped + sample_data$Final_reads), sample_data$proportion)
                                
                                ggplot(sample_data, aes(x = proportion, y = Host_mapped / (Host_mapped + Final_reads))) +
                                        geom_point () + 
                                        xlab("Host DNA proportion by qPCR") +
                                        ylab("(Host reads / host + final reads) by seqeuncing")
                                #......? what is the meaning of the metaphan mapped reads...?
                                ggplot(sample_data, aes(x = sample_type, y = proportion)) +
                                        geom_boxplot ()
                                        
                                #......? what is the meaning of the metaphan mapped reads...?
                                
                                

# Host DNA proportion change ----------------------------------------------

                                
        # Q2. Did Host DNA proportion changed by treatment at each group, in sequecing data?
                                
                                ggplot(sample_data, aes(x = treatment, y = sequencing_host_prop)) +
                                        xlab("Sample type") +
                                        ylab("(Host DNA / final reads) by seqeuncing") +
                                        geom_boxplot() +
                                        theme_bw() +
                                        facet_wrap(~ sample_type, scales= "free_x", ncol = 1)
                                
                                #facetted plot
                                ggplot(sample_data, aes(x = sample_type, y = sequencing_host_prop, fill = treatment)) +
                                        xlab("Sample type") +
                                        ylab("(Host DNA / final reads) by seqeuncing") +
                                        geom_boxplot() 
                                
                                ggplot(sample_data, aes(x = sequencing_host_prop, y = log10(Final_reads), fill = treatment)) +
                                        xlab("sequencing_host_prop") +
                                        ylab("Final reads") +
                                        geom_point() 
                                
                                

# Sequencing depth change in each treatment group -------------------------

                                
        # Q3. Did seqeuncing depth (total reads) chagned by treatment in each group?
                                
                                #It seems that qPCR captured a lot of bacterial DNA, however, sequencing pipeline was unable to map all the reads for taxonomic annotation.
                                #1. species richness and sequencing depth?
                                #calculation of alpha diversity indices
                                
                                alpha_diversity <- function(data) {
                                        otu_table <- otu_table(data)
                                        S.obs <- rowSums(t(otu_table) != 0)
                                        sample_data <- sample_data(data)
                                        data_evenness <- vegan::diversity(t(otu_table)) / log(vegan::specnumber(t(otu_table)))                # calculate evenness index using vegan package
                                        data_shannon <- vegan::diversity(t(otu_table), index = "shannon")                           # calculate Shannon index using vegan package
                                        data_hill <- exp(data_shannon)                           # calculate Hills index
                                        data_dominance <- microbiome::dominance(otu_table, index = "all", rank = 1, relative = A, aggregate = TRUE) # dominance (Berger-Parker index), etc.
                                        data_invsimpson <- vegan::diversity(t(otu_table), index = "invsimpson")                          # calculate Shannon index using vegan package
                                        alpha_diversity <- cbind(S.obs, data_shannon, data_hill, data_invsimpson, data_evenness,data_dominance) # combine all indices in one data table
                                        alpha_diversity$Row.names <- row.names(alpha_diversity)
                                        sample_data <- merge(data.frame(sample_data), alpha_diversity, by = "Row.names")
                                        sample_data$sample_type <- as.character(sample_data$sample_type)
                                        row.names(sample_data) <- sample_data$Row.names
                                        sample_data
                                }
                                
                                        sample_data(phyloseq_sample) <- sample_data(alpha_diversity(phyloseq_sample)) 
                                        sample_data(phyloseq_decontam) <- sample_data(alpha_diversity(phyloseq_decontam))
                                        sample_data(phyloseq_decontam2) <- sample_data(alpha_diversity(phyloseq_decontam2))
                                        
                                #sample_data$sample_type <- ifelse(grepl("BAL", sample_data$sample_id), "bal",
                                #                                  ifelse(grepl("NS", sample_data$sample_id), "nasal_swab",
                                #                                         ifelse(grepl("CFB", sample_data$sample_id), "sputum", sample_data$sample_type)))
                                #sample_data


# Summary figures - facet -------------------------------------------------

                                
        #Summary figures - for the ease of overview
                #Peggy's script - @Project_SICAS2_microbiome/5_Scripts/PSL/COD_20221104_PSL_Host_Depletion_Sequencing.R
                        ## figures -----------------------------------------------------------------
                                
                                sample_data %>% filter(sample_type == "nasal_swab") %>% filter(treatment == "qiaamp") %>% pull(Final_reads) %>% sort /1000
                                        sample_data %>% filter(sample_type == "nasal_swab") %>% filter(treatment == "host_zero") %>% pull(Final_reads) %>% sort /1000
                                        
                                facet_figure <- function(data, target) {
                                        sample_data <- sample_data(phyloseq_sample) %>% data.frame(check.names = F)
                                        sample_data$treatment <- factor(sample_data$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                        sample_data$sample_type <- tolower(sample_data$sample_type)
                                        names(sample_data) <- ifelse(names(sample_data) == "S.obs", "species_richness", names(sample_data))
                                        names(sample_data) <- ifelse(names(sample_data) == "dbp", "berger_parker", names(sample_data))
                                        sample_data$sample_type <- factor(sample_data$sample_type, levels = c("bal", "nasal_swab", "sputum"))
                                        file1_long <- sample_data %>%
                                                filter(sample_type %in% c("bal", "nasal_swab", "sputum")) %>%
                                                select(sample_type, treatment, Raw_reads, Host_mapped, Final_reads, "log10(Final_reads)",Metaphlan_mapped, sequencing_host_prop, total_reads, species_richness, data_invsimpson, data_shannon, berger_parker) %>%
                                                pivot_longer(cols = c("Raw_reads", "Host_mapped", "Final_reads", "log10(Final_reads)", "total_reads", "Metaphlan_mapped", "sequencing_host_prop", "species_richness", "data_shannon", "data_invsimpson", "berger_parker"), names_to = "feature", values_to = "value") %>% 
                                                mutate(feature = factor(feature, levels = c("Raw_reads", "Host_mapped", "Final_reads","log10(Final_reads)", "total_reads", "Metaphlan_mapped", "sequencing_host_prop", "species_richness", "data_shannon", "data_invsimpson", "berger_parker"))) 
                                        p2 <- file1_long %>% 
                                                filter(feature == target) %>%
                                                ggboxplot(x = "treatment", y = "value", fill= "treatment", palette = viridis(6)) +
                                                guides(fill=guide_legend(title="Treatment"))
                                        
                                        facet(p2 + theme_classic(), facet.by = c("sample_type", "feature"), scales = "free_y") +
                                                xlab("Treatment") +
                                                #guides(x =  guide_axis(angle = 90)) +
                                                aes(fill = viridis(6))
                                                # + guides(x =  guide_axis(angle = 90))                                                
                                }
                                sample_data(phyloseq_sample)$`log10(Final_reads)` <- sample_data(phyloseq_sample)$Final_reads %>% log10()
                                
                                facet_figure(phyloseq_sample, "Final_reads")
                                facet_figure(phyloseq_sample, "log10(Final_reads)")
                                facet_figure(phyloseq_sample, "Metaphlan_mapped")
                                facet_figure(phyloseq_sample, "sequencing_host_prop")
                                facet_figure(phyloseq_sample, "species_richness")
                                facet_figure(phyloseq_sample, "data_shannon")
                                
                                facet_figure(phyloseq_sample, "data_invsimpson")
                                facet_figure(phyloseq_sample, "berger_parker")
                                #data %>% select(sample_type, treatment, Raw_reads, Host_mapped, Final_reads, Metaphlan_mapped, sequencing_host_prop, total_reads, species_richness, data_invsimpson, data_shannon, berger_parker) %>%
                                #                mutate_all( - mean(a)) / sd()
                                
                                facet_figure(phyloseq_decontam, "species_richness")
                                facet_figure(phyloseq_decontam, "data_shannon")
                                facet_figure(phyloseq_decontam, "data_invsimpson")
                                facet_figure(phyloseq_decontam, "berger_parker")
                                

# Summary figures - facet and z-score -------------------------------------

                                
                                
                                
                                
                                
                                facet_figure_z <- function(data, target) {
                                        sample_data <- sample_data(data) %>% data.frame()
                                        sample_data <- do.call(data.frame,                      # Replace Inf in data by NA
                                                           lapply(sample_data,
                                                                  function(x) replace(x, is.infinite(x), NA)))
                                        zscore <- sample_data
                                        zscore$`log(Raw_reads)` <- log(zscore$Raw_reads)
                                        zscore$`log(Host_mapped)` <- log(zscore$Host_mapped)
                                        zscore$`log(Final_reads)` <- log(zscore$Final_reads)
                                        zscore$`log(Metaphlan_mapped)` <- log(zscore$Metaphlan_mapped)
                                        zscore <- zscore %>% select(`log(Raw_reads)`, `log(Host_mapped)`, `log(Final_reads)`, `log(Metaphlan_mapped)`, sequencing_host_prop, total_reads, S.obs, data_invsimpson, data_shannon, dbp) %>% scale %>% data.frame(check.names = F)
                                        sample_data <- sample_data %>% select(-c(Raw_reads, Host_mapped, Final_reads, Metaphlan_mapped, sequencing_host_prop, total_reads, S.obs, data_invsimpson, data_shannon, dbp))
                                        sample_data <- cbind (zscore, sample_data)
                                        sample_data$treatment <- factor(sample_data$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                        sample_data$sample_type <- tolower(sample_data$sample_type)
                                        names(sample_data) <- ifelse(names(sample_data) == "S.obs", "species_richness", names(sample_data))
                                        names(sample_data) <- ifelse(names(sample_data) == "dbp", "berger_parker", names(sample_data))
                                        sample_data$sample_type <- factor(sample_data$sample_type, levels = c("bal", "nasal_swab", "sputum"))
                                        file1_long <- sample_data %>%
                                                filter(sample_type %in% c("bal", "nasal_swab", "sputum")) %>%
                                                select(sample_type, treatment, `log(Raw_reads)`, `log(Host_mapped)`, `log(Final_reads)`, `log(Metaphlan_mapped)`, sequencing_host_prop, total_reads, species_richness, data_invsimpson, data_shannon, berger_parker) %>%
                                                pivot_longer(cols = c("log(Raw_reads)", "log(Host_mapped)", "log(Final_reads)", "log(Metaphlan_mapped)", "sequencing_host_prop", "species_richness", "data_shannon", "data_invsimpson", "berger_parker"), names_to = "feature", values_to = "value") %>%
                                                mutate(feature = factor(feature, levels = c("log(Raw_reads)", "log(Host_mapped)", "log(Final_reads)", "log(Metaphlan_mapped)", "sequencing_host_prop", "species_richness", "data_shannon", "data_invsimpson", "berger_parker"))) 
                                        p2 <- file1_long %>% 
                                                filter(feature %in% target) %>%
                                                ggboxplot(x = "treatment", y = "value", fill= "treatment", palette = viridis(6)) +
                                                guides(fill=guide_legend(title="Treatment"))
                                        
                                        facet(p2 + theme_classic(), facet.by = c("sample_type", "feature")) +
                                                ylab("z score") +
                                                xlab("Treatment") +
                                                guides(x =  guide_axis(angle = 90)) +
                                                aes(fill = viridis(6))
                                                
                                        }
                                facet_figure_z(phyloseq_sample, c("log(Raw_reads)", "log(Host_mapped)", "log(Final_reads)", "log(Metaphlan_mapped)", "sequencing_host_prop", "species_richness", "data_shannon", "data_invsimpson", "berger_parker"))
                                facet_figure_z(phyloseq_sample, c("log(Raw_reads)", "log(Host_mapped)", "log(Final_reads)", "log(Metaphlan_mapped)", "sequencing_host_prop"))
                                
                                
                                
                                
                                
                                ggplot(data.frame(sample_data(phyloseq_sample), check.names = F), aes(x = log10(Final_reads), y = S.obs, col = sample_type)) +
                                        geom_point() +
                                        scale_color_manual(values = viridis(3)) +
                                        theme_bw()
                                
                                ggplot(data.frame(sample_data(phyloseq_sample), check.names = F), aes(x = Final_reads, y = S.obs, col = treatment)) +
                                        geom_point() +
                                        #scale_color_manual(values = viridis(6)) +
                                        theme_bw() +
                                        ylab("Species richness") +
                                        facet_wrap(~sample_type)

                                ggplot(data.frame(sample_data(phyloseq_sample), check.names = F), aes(x = Final_reads, y = S.obs, col = sample_type)) +
                                        geom_point() +
                                        #scale_color_manual(values = viridis(6)) +
                                        theme_bw() +
                                        ylab("Species richness") +
                                        facet_wrap(~treatment)

# Modeling - trimming data ------------------------------------------------

                                
                        #modeling final reads
                        #using all the data.
                                sample_data <- sample_data(phyloseq_sample) %>% data.frame()
                                sample_data$treatment <- factor(sample_data$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                
                                sample_data$sample_type <- tolower(sample_data$sample_type)
                                names(sample_data) <- ifelse(names(sample_data) == "S.obs", "species_richness", names(sample_data))
                                names(sample_data) <- ifelse(names(sample_data) == "dbp", "berger_parker", names(sample_data))
                                
                                sample_data$sample_type <- factor(sample_data$sample_type, levels = c("bal", "nasal_swab", "sputum"))


# Modeling - Final reads ~ input ------------------------------------------

                                                
                        #final reads
                                
                                
                                lmer(log10(Final_reads) ~ treated + (1|original_sample), data = sample_data) %>% summary
                                #control
                                lmer(log10(Final_reads) ~ treatment + (1|original_sample), data = sample_data) %>% summary()
                                #treatment method
                                lmer(log10(Final_reads) ~ sample_type + treatment + (1|original_sample), data = sample_data) %>% summary()
                                lmer(log10(Final_reads) ~ sample_type * treatment + (1|original_sample), data = sample_data) %>% summary()
                        

# Modeling - Alpha diversity ~ input --------------------------------------

                                
                        #alpha diversity
                        #Species richness
                                #control
                                lmer(species_richness ~ treated + (1|original_sample), data = sample_data) %>% summary()
                                #treatment method
                                lmer(species_richness ~ sample_type + treatment + (1|original_sample), data = sample_data) %>% summary()
                                lmer(species_richness ~ sample_type * treatment + (1|original_sample), data = sample_data) %>% summary()
                                lmer(species_richness ~ Final_reads + sample_type * treatment + (1|original_sample), data = sample_data) %>% summary()
                        #Shannon
                                #control
                                lmer(data_shannon ~ treated + (1|original_sample), data = sample_data) %>% summary()
                                #treatment method
                                lmer(data_shannon ~ sample_type + treatment + (1|original_sample), data = sample_data) %>% summary()
                                lmer(data_shannon ~ sample_type * treatment + (1|original_sample), data = sample_data) %>% summary()
                                
                                

# Modeling - using decontam data ------------------------------------------

                                        
                                        
                                #using decontam data.
                                
                                sample_data_decontam <- sample_data(phyloseq_decontam) %>% data.frame()
                                        sample_data_decontam$treatment <- factor(sample_data_decontam$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                        sample_data_decontam$sample_type <- tolower(sample_data_decontam$sample_type)
                                        names(sample_data_decontam) <- ifelse(names(sample_data_decontam) == "S.obs", "species_richness", names(sample_data_decontam))
                                        names(sample_data_decontam) <- ifelse(names(sample_data_decontam) == "dbp", "berger_parker", names(sample_data_decontam))
                                        sample_data_decontam$sample_type <- factor(sample_data_decontam$sample_type, levels = c("bal", "nasal_swab", "sputum"))
                                        
                                        #final reads
                                        #control
                                        lmer(Final_reads ~ sample_type + treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        #treatment method
                                        lmer(Final_reads ~ sample_type + treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        lmer(Final_reads ~ sample_type * treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        
                                        #alpha diversity
                                        #Species richness
                                        #control
                                        lmer(species_richness ~ sample_type + treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        
                                        lmer(species_richness ~ log10(Final_reads) + sample_type + treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        lmer(species_richness ~ log10(Final_reads) + sample_type + treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        lmer(species_richness ~ log10(Final_reads) + sample_type * treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        #treatment method
                                        lmer(species_richness ~ sample_type + treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        lmer(species_richness ~ sample_type * treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        #Shannon
                                        #control
                                        lmer(data_shannon ~ sample_type + treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        #treatment method
                                        lmer(data_shannon ~ sample_type * treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        
                                #using decontam data 2 (harsher condition.
                                        sample_data_decontam <- sample_data(phyloseq_decontam2) %>% data.frame()
                                        sample_data_decontam$treatment <- factor(sample_data_decontam$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                        sample_data_decontam$sample_type <- tolower(sample_data_decontam$sample_type)
                                        names(sample_data_decontam) <- ifelse(names(sample_data_decontam) == "S.obs", "species_richness", names(sample_data_decontam))
                                        names(sample_data_decontam) <- ifelse(names(sample_data_decontam) == "dbp", "berger_parker", names(sample_data_decontam))
                                        sample_data_decontam$sample_type <- factor(sample_data_decontam$sample_type, levels = c("bal", "nasal_swab", "sputum"))
                                        
                                        #final reads
                                        #control
                                        lmer(Final_reads ~ sample_type + treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        #treatment method
                                        lmer(Final_reads ~ sample_type + treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        lmer(Final_reads ~ sample_type * treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        
                                        #alpha diversity
                                        #Species richness
                                        #control
                                        lmer(species_richness ~ sample_type + treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        #treatment method
                                        lmer(species_richness ~ sample_type + treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        lmer(species_richness ~ sample_type * treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        #Shannon
                                        #control
                                        lmer(data_shannon ~ sample_type + treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        #treatment method
                                        lmer(data_shannon ~ sample_type * treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                        
                                #decontam for nasal swab samples
                                        
                                        sample_data_decontam <- subset(sample_data_decontam, sample_data_decontam$sample_type == "nasal_swab")
                                                #final reads
                                                #control
                                                lmer(Final_reads ~ treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                                #treatment method
                                                lmer(Final_reads ~  treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                                
                                                
                                                #alpha diversity
                                                #Species richness
                                                #control
                                                lmer(species_richness ~  treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                                #treatment method
                                                lmer(species_richness ~ treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                                #Shannon
                                                #control
                                                lmer(data_shannon ~ treated + (1|original_sample), data = sample_data_decontam) %>% summary()
                                                #treatment method
                                                lmer(data_shannon ~ treatment + (1|original_sample), data = sample_data_decontam) %>% summary()
                                                
                                        
                                        
                #using Nasal swabs
                                        phyloseq_ns <- subset_samples(phyloseq_sample, sample_type == "nasal_swab")
                                        sample_data <- sample_data(subset_samples(phyloseq_ns)) %>% data.frame()
                                        sample_data$treatment <- factor(sample_data$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                        sample_data$sample_type <- tolower(sample_data$sample_type)
                                        names(sample_data) <- ifelse(names(sample_data) == "S.obs", "species_richness", names(sample_data))
                                        names(sample_data) <- ifelse(names(sample_data) == "dbp", "berger_parker", names(sample_data))
                                        
                                        
                                        #final reads
                                        #control
                                        lmer(Final_reads ~ treated + (1|original_sample), data = sample_data) %>% summary()
                                        #treatment method
                                        lmer(Final_reads ~ treatment + (1|original_sample), data = sample_data) %>% summary()
                                        
                                        
                                        #alpha diversity
                                        #Species richness
                                        lmer(species_richness ~ treatment + (1|original_sample), data = sample_data) %>% summary()
                                        lmer(species_richness ~ log10(Final_reads) + treatment + (1|original_sample), data = sample_data) %>% summary()
                #Nasal swabs decontam
                                        
                                        
                                        lmer(species_richness ~ treatment + (1|original_sample), data = subset(sample_data_decontam, sample_data$sample_type =="nasal_swab")) %>% summary()
                                        
                                        
                                        #Shannon
                                        lmer(data_shannon ~ treatment + (1|original_sample), data = sample_data) %>% summary()
                                        lmer(data_shannon ~ treatment + (1|original_sample), data = subset(sample_data_decontam, sample_data_decontam$sample_type =="nasal_swab")) %>% summary()
                                        
                                        par(mfrow = c(1, 2))
                                        boxplot(sample_data$species_richness ~ sample_data$treatment, xlab = "Treatment", ylab = "Species richness", main = "Nasal swab, before decontam", ylim = c(0,40)) 
                                        boxplot(sample_data_decontam$species_richness ~ sample_data_decontam$treatment, xlab = "Treatment", ylab = "Species richness", main = "Nasal swab, after decontam", ylim = c(0,40))
                                        dev.off()
                                        
                                        par(mfrow = c(1, 2))
                                        boxplot(sample_data$data_shannon ~ sample_data$treatment, xlab = "Treatment", ylab = "Shannon", main = "Nasal swab, before decontam", ylim = c(0,2)) 
                                        boxplot(sample_data_decontam$data_shannon ~ sample_data_decontam$treatment, xlab = "Treatment", ylab = "Shannon", main = "Nasal swab, after decontam", ylim = c(0,2))
                                        dev.off()

# Beta diversity - plotting -----------------------------------------------

                                        
        #beta-diversity
                                sample_data(phyloseq)
                                phyloseq_rel = transform_sample_counts(phyloseq_sample, function(x){x / sum(x)})
                                phyloseq_rel_nz = subset_samples(phyloseq_rel, S.obs != 0)
                                phyloseq_rel_nz = subset_samples(phyloseq_rel_nz, sample_type != "pos_control" & sample_type != "neg_control")
                                sample_data(phyloseq_rel_nz)
                                
                                #pca_bray_rel <- ordinate(phyloseq_rel,  method = "PCoA", distance = "bray")
                                
                                pca_bray <- ordinate(phyloseq_rel_nz,  method = "PCoA", distance = "bray")
                                NMDS_bray <- ordinate(phyloseq_rel_nz,  method = "NMDS", distance = "bray")
                                NMDS_jsd <- ordinate(phyloseq_rel_nz,  method = "NMDS", distance = "jsd")
                                
                                
                                ordinate(phyloseq_rel_nz,  method = "NMDS", distance = "bray") %>% plot_ordination(phyloseq_rel_nz, ., col = "original_sample", shape = "treatment" ) + geom_point(size = 4)

# Beta diversity - PERMANOVA ----------------------------------------------


                                #Jake/Laura comment on 11142022
                                #Look subset of the data for PERMANOVA
                                
                                        #1. interactin term
                                
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type * treatment + original_sample,
                                               data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                               permutations = 10000)
                                
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type * lypma + sample_type * benzonase + sample_type * host_zero + sample_type * molysis + sample_type * qiaamp + original_sample,
                                               data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                               permutations = 10000)
                                
                                
                                #Orignial analyses
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + treatment + original_sample,
                                               data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                               permutations = 10000)
                                
                                
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + treatment,
                                                           strata = sample_data(phyloseq_rel_nz)$original_sample, data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                           permutations = 10000)
                                
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp,
                                                            strata = sample_data(phyloseq_rel_nz)$original_sample, data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                            permutations = 10000)
                                        
                                
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + log10(Final_reads) + lypma + benzonase + host_zero + molysis + qiaamp,
                                               strata = sample_data(phyloseq_rel_nz)$original_sample, data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                               permutations = 10000)
                                
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + log10(Final_reads) + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                               data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), permutations = 10000)
                                
                                
                                
                                
                                
                                
                                
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp,
                                                       strata = sample_data(phyloseq_rel_nz)$original_sample, data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                       permutations = 10000)
                                #sample_id as fixed effect
                                vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                            data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                            permutations = 10000)
                                        #other distances
                                        vegan::adonis2(distance(phyloseq_rel_nz, method="jaccard") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                       data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                       permutations = 10000)
                                        vegan::adonis2(distance(phyloseq_rel_nz, method="jsd") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                       data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                       permutations = 10000)
                                        vegan::adonis2(distance(phyloseq_rel_nz, method="kulczynski") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                       data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                       permutations = 10000)
                                        vegan::adonis2(distance(phyloseq_rel_nz, method="manhattan") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                       data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                       permutations = 10000)
                                        
                                        vegan::adonis2(distance(phyloseq_rel_nz, method="manhattan") ~ Final_reads + sample_type + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                       data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                       permutations = 10000)
                                        
                                        dist_methods <- unlist(distanceMethodList)
                                        print(dist_methods)
                                        
                                        
                                perm_bray

# Beta diversity - boxplots -----------------------------------------------

                                
                        #distances of betadiversity - boxplots
                                bray_dist_long <- distance(phyloseq, method="bray") %>% as.matrix() %>% melt_dist() #making long data of distance matrices
                                #Adding sample type and treatment name. 
                                #this can be also done by merging metadata into the `bray_dist_long`
                                names <- data.frame(str_split_fixed(bray_dist_long$iso1, "_", 3))
                                names2 <- data.frame(str_split_fixed(bray_dist_long$iso2, "_", 3))
                                bray_dist_long$sample_id_1 <- paste(names$X1, names$X2, sep = "_")
                                bray_dist_long$method_1 <- ifelse(grepl("control", bray_dist_long$iso1),"control", 
                                                                  ifelse(grepl("lyPMA", bray_dist_long$iso1),"lyPMA", 
                                                                         ifelse(grepl("benzonase", bray_dist_long$iso1),"benzonase", 
                                                                                ifelse(grepl("host", bray_dist_long$iso1),"host_zero", 
                                                                                       ifelse(grepl("qia", bray_dist_long$iso1),"qiaamp", 
                                                                                              ifelse(grepl("moly", bray_dist_long$iso1),"molysis", 
                                                                                                     NA))))))
                                #Adding data for iso 2 also should be done
                                bray_dist_long$sample_id_2 <- paste(names2$X1, names2$X2, sep = "_")
                                bray_dist_long$method_2 <-ifelse(grepl("control", bray_dist_long$iso2),"control", 
                                                                 ifelse(grepl("lyPMA", bray_dist_long$iso2),"lyPMA", 
                                                                        ifelse(grepl("benzonase", bray_dist_long$iso2),"benzonase", 
                                                                               ifelse(grepl("host", bray_dist_long$iso2),"host_zero", 
                                                                                      ifelse(grepl("qia", bray_dist_long$iso2),"qiaamp", 
                                                                                             ifelse(grepl("moly", bray_dist_long$iso2),"molysis", 
                                                                                                    NA))))))
                                #subsetting distances of my interest
                                bray_dist_long_within_sampleid <- subset(bray_dist_long, bray_dist_long$sample_id_1 == bray_dist_long$sample_id_2)
                                bray_dist_long_within_sampleid_from_control <- subset(bray_dist_long_within_sampleid, bray_dist_long_within_sampleid$method_1 == "control" | bray_dist_long_within_sampleid$method_2 == "control" )
                                bray_dist_long_within_sampleid_from_control$treatment <- bray_dist_long_within_sampleid_from_control$method_1
                                bray_dist_long_within_sampleid_from_control$treatment <- ifelse(bray_dist_long_within_sampleid_from_control$treatment == "control", bray_dist_long_within_sampleid_from_control$method_2, bray_dist_long_within_sampleid_from_control$treatment)
                                bray_dist_long_within_sampleid_from_control$treatment <- factor(bray_dist_long_within_sampleid_from_control$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                bray_dist_long_within_sampleid_from_control$treatment <- factor(bray_dist_long_within_sampleid_from_control$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                
                                bray_dist_long_within_sampleid_from_control$sample_type <- ifelse(grepl("NS", bray_dist_long_within_sampleid_from_control$iso1), "nasal_swab",
                                                                                                        ifelse(grepl("CFB", bray_dist_long_within_sampleid_from_control$iso1), "sputum",
                                                                                                                     ifelse(grepl("BAL", bray_dist_long_within_sampleid_from_control$iso1), "bal", NA)))
                                
                                ggplot(bray_dist_long_within_sampleid_from_control, aes(y = dist, x = treatment)) +
                                        geom_boxplot() +
                                        theme_bw() +
                                        facet_wrap(~sample_type) +
                                        guides(x =  guide_axis(angle = 90))      
                                
                                

w# Beta diversity - for NS subset data -------------------------------------

                                
                                        
                                        
        #beta-diversity - NS decontam
                                        sample_data(phyloseq_decontam2)
                                        phyloseq_rel = transform_sample_counts(phyloseq_decontam2, function(x){x / sum(x)})
                                        phyloseq_rel_nz = subset_samples(phyloseq_rel, S.obs != 0)
                                        phyloseq_rel_nz_decontam = subset_samples(phyloseq_rel_nz, sample_type == "nasal_swab")
                                        phyloseq_rel = transform_sample_counts(phyloseq, function(x){x / sum(x)})
                                        phyloseq_rel_nz = subset_samples(phyloseq_rel, S.obs != 0)
                                        phyloseq_rel_nz = subset_samples(phyloseq_rel_nz, sample_type == "nasal_swab")
                                        #pca_bray_rel <- ordinate(phyloseq_rel,  method = "PCoA", distance = "bray")
                                        
                                        NMDS_bray <- ordinate(phyloseq_rel_nz,  method = "NMDS", distance = "bray")
                                        NMDS_bray_decon <- ordinate(phyloseq_rel_nz_decontam,  method = "NMDS", distance = "bray")
                                        
                                        ordinate(phyloseq_rel_nz,  method = "NMDS", distance = "bray") %>% plot_ordination(phyloseq_rel_nz, ., col = "original_sample", shape = "treatment" ) + geom_point(size = 4)
                                        ordinate(phyloseq_rel_nz_decontam,  method = "NMDS", distance = "bray") %>% plot_ordination(phyloseq_rel_nz_decontam, ., col = "original_sample", shape = "treatment" ) + geom_point(size = 4)
                                        
                                        ordinate(phyloseq_rel_nz,  method = "PCoA", distance = "bray") %>% plot_ordination(phyloseq_rel_nz, ., col = "original_sample", shape = "treatment" ) + geom_point(size = 4)
                                        ordinate(phyloseq_rel_nz_decontam,  method = "PCoA", distance = "bray") %>% plot_ordination(phyloseq_rel_nz_decontam, ., col = "original_sample", shape = "treatment" ) + geom_point(size = 4)
                                        
                                        vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ treatment + original_sample, data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                                    permutations = 10000)
                                        vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ lypma + benzonase + host_zero + molysis + qiaamp,
                                                                    strata = sample_data(phyloseq_rel_nz)$original_sample, data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                                    permutations = 10000)
                                        vegan::adonis2(distance(phyloseq_rel_nz, method="bray") ~ lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                                    data = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F),
                                                                    permutations = 10000)
                                        vegan::adonis2(distance(phyloseq_rel_nz_decontam, method="bray") ~ lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                                       data = phyloseq_rel_nz_decontam %>% sample_data %>% data.frame(check.names = F),
                                                       permutations = 10000)
                                        
                                        (harrietr)
                                        bray_dist_long <- distance(phyloseq_rel_nz_decontam, method="bray") %>% as.matrix() %>% melt_dist()
                                        str_split_fixed(bray_dist_long$iso1, "_", 5)
                                        names <- data.frame(str_split_fixed(bray_dist_long$iso1, "_", 3))
                                        names2 <- data.frame(str_split_fixed(bray_dist_long$iso2, "_", 3))
                                        bray_dist_long$sample_id_1 <- paste(names$X1, names$X2, sep = "_")
                                        bray_dist_long$method_1 <- ifelse(grepl("control", bray_dist_long$iso1),"control", 
                                                                          ifelse(grepl("lyPMA", bray_dist_long$iso1),"lyPMA", 
                                                                                 ifelse(grepl("benzonase", bray_dist_long$iso1),"benzonase", 
                                                                                        ifelse(grepl("host", bray_dist_long$iso1),"host_zero", 
                                                                                               ifelse(grepl("qia", bray_dist_long$iso1),"qiaamp", 
                                                                                                      ifelse(grepl("moly", bray_dist_long$iso1),"molysis", 
                                                                                                             NA))))))
                                        
                                        bray_dist_long$sample_id_2 <- paste(names2$X1, names2$X2, sep = "_")
                                        bray_dist_long$method_2 <-ifelse(grepl("control", bray_dist_long$iso2),"control", 
                                                                         ifelse(grepl("lyPMA", bray_dist_long$iso2),"lyPMA", 
                                                                                ifelse(grepl("benzonase", bray_dist_long$iso2),"benzonase", 
                                                                                       ifelse(grepl("host", bray_dist_long$iso2),"host_zero", 
                                                                                              ifelse(grepl("qia", bray_dist_long$iso2),"qiaamp", 
                                                                                                     ifelse(grepl("moly", bray_dist_long$iso2),"molysis", 
                                                                                                            NA))))))
                                        bray_dist_long_within_sampleid
                                        bray_dist_long_within_sampleid <- subset(bray_dist_long, bray_dist_long$sample_id_1 == bray_dist_long$sample_id_2)
                                        bray_dist_long_within_sampleid_from_control <- subset(bray_dist_long_within_sampleid, bray_dist_long_within_sampleid$method_1 == "control" | bray_dist_long_within_sampleid$method_2 == "control" )
                                        bray_dist_long_within_sampleid_from_control$treatment <- bray_dist_long_within_sampleid_from_control$method_1
                                        bray_dist_long_within_sampleid_from_control$treatment <- ifelse(bray_dist_long_within_sampleid_from_control$treatment == "control", bray_dist_long_within_sampleid_from_control$method_2, bray_dist_long_within_sampleid_from_control$treatment)
                                        bray_dist_long_within_sampleid_from_control$treatment <- factor(bray_dist_long_within_sampleid_from_control$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                        bray_dist_long_within_sampleid_from_control$treatment <- factor(bray_dist_long_within_sampleid_from_control$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                        
                                        bray_dist_long_within_sampleid_from_control$sample_type <- ifelse(grepl("NS", bray_dist_long_within_sampleid_from_control$iso1), "nasal_swab",
                                                                                                          ifelse(grepl("CFB", bray_dist_long_within_sampleid_from_control$iso1), "sputum",
                                                                                                                 ifelse(grepl("BAL", bray_dist_long_within_sampleid_from_control$iso1), "bal", NA)))
                                        bray_dist_long_within_sampleid_from_control
                                        ggplot(bray_dist_long_within_sampleid_from_control, aes(y = dist, x = treatment)) +
                                                geom_boxplot() +
                                                theme_bw() +
                                                facet_wrap(~sample_type) +
                                                guides(x =  guide_axis(angle = 90))      
                                        
                                
                                        decontam_change <- merge(sample_data(phyloseq_sample) %>% data.frame() %>% select(c("Row.names", "sample_type", "original_sample", "treatment", "S.obs")),
                                                sample_data(phyloseq_decontam2) %>% data.frame() %>% select(c("Row.names", "S.obs")),
                                              by = "Row.names")
                                        
                                        decontam_change$change <- decontam_change$S.obs.x - decontam_change$S.obs.y
                                        decontam_change
                                        ggplot(decontam_change, aes(x = treatment, y = change)) +
                                                geom_boxplot() +
                                                facet_wrap(~ sample_type, ncol = 1) +
                                                ylab("Number of contmainants") +
                                                xlab("Treatment methods")
                                                #guides(x =  guide_axis(angle = 90))      
                                        
                                        decontam_change$treatment <- factor(decontam_change$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                        decontam_change$sample_type <- tolower(decontam_change$sample_type)
                                        decontam_change$sample_type <- factor(decontam_change$sample_type, levels = c("bal", "nasal_swab", "sputum"))
                                        lm(change ~ treatment, data = decontam_change) %>% summary()
                                        lmer(change ~ treatment + (1|original_sample), data = decontam_change) %>% summary()
                                        lmer(change ~ treatment * sample_type + (1|original_sample), data = decontam_change) %>% summary()
                                        
                                        lm(change ~ treatment + sample_type, data = decontam_change) %>% summary()
                                        lm(change ~ treatment * sample_type, data = decontam_change) %>% summary()

                                        lm(change ~ treatment, data = subset(decontam_change, decontam_change$sample_type == "nasal_swab")) %>% summary()
                                
                                        
                                #library_fiailed_samples
                                ordinate(phyloseq_rel_nz,  method = "NMDS", distance = "bray") %>% plot_ordination(phyloseq_rel_nz, ., shape = "lib_failed", col = "original_sample" ) + geom_point(size = 4)
                                ordinate(phyloseq_rel_nz,  method = "NMDS", distance = "bray") %>% plot_ordination(phyloseq_rel_nz, ., shape = "lib_failed", col = "original_sample" ) + geom_point(size = 4)
                                
                                
                                
                                lib_fail_list
                                
                                #plot_abundance
                                        

# Additional models - S.obs ~ sequencing depth + others -------------------

                                        
                #Seqeuncing depth - modeling
                                sample_data_decontam <- phyloseq_decontam2 %>% sample_data %>% data.frame()
                                sample_data_decontam$treatment <- factor(sample_data_decontam$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                sample_data_decontam$sample_type <- tolower(sample_data_decontam$sample_type)
                                sample_data_decontam$sample_type <- factor(sample_data_decontam$sample_type, levels = c("bal", "nasal_swab", "sputum", "pos_control", "neg_control"))
                                lmer(S.obs ~ log10(Final_reads) + sample_type + (1|original_sample), data = sample_data_decontam) %>% summary
                                lmer(S.obs ~ log10(Final_reads) + treated + (1|original_sample), data = sample_data_decontam) %>% summary
                                lmer(S.obs ~ log10(Final_reads) + treatment + (1|original_sample), data = sample_data_decontam) %>% summary
                                lmer(S.obs ~ log10(Final_reads) + treatment + sample_type + (1|original_sample), data = sample_data_decontam) %>% summary
                                lmer(S.obs ~ log10(Final_reads) + treatment * sample_type + (1|original_sample), data = sample_data_decontam) %>% summary
                                lmer(S.obs ~ log10(Final_reads) + treatment + sample_type + treatment * sample_type + (1|original_sample), data = sample_data_decontam) %>% summary
                                lmer(S.obs ~ log10(Final_reads) + treatment + (1|original_sample), data = subset(sample_data_decontam, sample_data_decontam$sample_type == "nasal_swab")) %>% summary
                                lmer(S.obs ~ log10(Final_reads) + treatment + (1|original_sample), data = subset(sample_data_decontam, sample_data_decontam$sample_type == "bal")) %>% summary
                                lmer(S.obs ~ log2(Final_reads) + treatment + (1|original_sample), data = subset(sample_data_decontam, sample_data_decontam$sample_type == "sputum")) %>% summary
                                
                                
                                
                        #Extraction date
                                
                                lmer(species_richness ~ extraction_date + (1|original_sample), data = sample_data) %>% summary
                                lmer(species_richness ~ extraction_date + sample_type + (1|original_sample), data = sample_data) %>% summary
                                
                                #There could be differences between extraction date, howver, thier difference was much less than that of sample_type
                        #depletion date
                                lmer(species_richness ~ depletion_date + (1|original_sample), data = sample_data) %>% summary
                                lmer(species_richness ~ depletion_date + sample_type + (1|original_sample), data = sample_data) %>% summary
                                lmer(species_richness ~ depletion_date + sample_type + treatment + (1|original_sample), data = sample_data) %>% summary
                                
                                lmer(species_richness ~ depletion_date + sample_type + log10(Final_reads) + (1|original_sample), data = sample_data) %>% summary

                                
                                

# Contamination check - bacterial DNA amount - by treatment ---------------

        ggplot(sample_data, aes(x = treatment, y = DNA_bac_nondil)) +
                                        geom_boxplot() +
                                        facet_wrap (~ original_sample, scales= "free_x")
                                        
                                sample_data
                                
                        #depletion date
                                lmer(DNA_bac_nondil ~ depletion_date + (1|original_sample), data = sample_data) %>% summary
                                lmer(DNA_bac_nondil ~ depletion_date + treatment + (1|original_sample), data = sample_data) %>% summary
                                lmer(DNA_bac_nondil ~ depletion_date + sample_type + (1|original_sample), data = sample_data) %>% summary
                                lmer(DNA_bac_nondil ~ depletion_date + treatment + sample_type + (1|original_sample), data = sample_data) %>% summary
                                
                                
# DA analsyis - MaAslin ---------------------------------------------------

                                
                                
                                
                                
                                
        #DA analysis - MaAslin
                                
                                ?Maaslin2
                                input_data = system.file(
                                        "extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
                                input_data
                                input_metadata = system.file(
                                        "extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file
                                df_input_data = read.table(file = input_data, header = TRUE, sep = "\t",
                                                           row.names = 1,
                                                           stringsAsFactors = FALSE)
                                
                                sample_data_sample <- sample_data(phyloseq_rel_nz) %>% data.frame()
                                otu_data_sample <- otu_table(phyloseq_rel_nz) %>% t %>% data.frame()
                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                

                #Running MaAslin for all sample without decontam
                                
                                #model 1 (all samples): taxa ~ sample_type + final_reads + lypma + benzonase + host_zero + molysis + qiaamp + (1 | original_sample)
                                #for taxa differentially abundant by host depletion method, look to see which ones overlap with potential contaminant taxa
                                #model 2 (nasal only): taxa ~ sample_type + final_reads + lypma + benzonase + host_zero + molysis + qiaamp + (1 | original_sample)
                                
                                

# Maaslin - # y ~ sample_type + treatment  --------------------------------

                
                                #all samples
                                fit_data = Maaslin2(
                                        input_data = otu_data_sample, 
                                        input_metadata = sample_data_sample, 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("sample_type", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # default
                                        random_effects = c("original_sample"), 
                                        reference = c("sample_type,BAL"))
                                fit_data$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()
                                
                       
                        #NS
                        # y ~ + sample_type + treatment 
                                
                                phyloseq_rel_nz_ns <- subset_samples(phyloseq_rel_nz, sample_type == "nasal_swab")
                                sample_data_sample <- sample_data(phyloseq_rel_nz_ns) %>% data.frame()
                                otu_data_sample <- otu_table(phyloseq_rel_nz_ns) %>% t %>% data.frame()
                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                
                                fit_data_ns = Maaslin2(
                                        input_data = otu_data_sample, 
                                        input_metadata = sample_data_sample, 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # default
                                        random_effects = c("original_sample"))
                                
                                fit_data_ns$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()
                                fit_data_ns$results
                        
                        #bal
                        # y ~ + sample_type + treatment 
                                phyloseq_rel_nz_bal <- subset_samples(phyloseq_rel_nz, sample_type == "BAL")
                                sample_data_sample <- sample_data(phyloseq_rel_nz_bal) %>% data.frame()
                                otu_data_sample <- otu_table(phyloseq_rel_nz_bal) %>% t %>% data.frame()
                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                
                                fit_data_bal = Maaslin2(
                                        input_data = otu_data_sample, 
                                        input_metadata = sample_data_sample, 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # default
                                        random_effects = c("original_sample"))
                                
                                fit_data_bal$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()
                                
                        #sputum
                        # y ~ + sample_type + treatment 
                                phyloseq_rel_nz_spt <- subset_samples(phyloseq_rel_nz, sample_type == "Sputum")
                                sample_data_sample <- sample_data(phyloseq_rel_nz_spt) %>% data.frame()
                                otu_data_sample <- otu_table(phyloseq_rel_nz_spt) %>% t %>% data.frame()
                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                
                                fit_data_spt = Maaslin2(
                                        input_data = otu_data_sample, 
                                        input_metadata = sample_data_sample, 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # default
                                        random_effects = c("original_sample"))
                                
                                fit_data_spt$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()
                                install.packages("writexl")
                                
                                sheets <- list(
                                        "all_sample" = fit_data$results,
                                        "nasal_swab" =fit_data_ns$results,
                                        "bal" = fit_data_bal$results,
                                        "sputum" = fit_data_spt$results)
                                write_xlsx(sheets, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/6_Results/host_depleiton/maaslin/DAT_20221207_MGK_maaslin_treatment.xlsx")
                                

# Maaslin - # # y ~ log(final reads) + sample_type + treatment  -----------

        
                                
                                #all samples
                                sample_data_sample <- sample_data(phyloseq_rel_nz) %>% data.frame()
                                otu_data_sample <- otu_table(phyloseq_rel_nz) %>% t %>% data.frame()
                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                
                                fit_data = Maaslin2(
                                        input_data = otu_data_sample, 
                                        input_metadata = sample_data_sample, 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("sample_type", "log10.Final_reads.", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # default
                                        random_effects = c("original_sample"), 
                                        reference = c("sample_type,BAL"))
                                
                                fit_data_1$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()
                                
                                #NS
                                # # y ~ log(final reads) + sample_type + treatment 
                                
                                phyloseq_rel_nz_ns <- subset_samples(phyloseq_rel_nz, sample_type == "nasal_swab")
                                sample_data_sample <- sample_data(phyloseq_rel_nz_ns) %>% data.frame()
                                otu_data_sample <- otu_table(phyloseq_rel_nz_ns) %>% t %>% data.frame()
                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                
                                fit_data_ns = Maaslin2(
                                        input_data = otu_data_sample, 
                                        input_metadata = sample_data_sample, 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads.","lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # default
                                        random_effects = c("original_sample"))
                                
                                fit_data_ns$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()
                                fit_data_ns$results
                                
                                #bal
                                # # y ~ log(final reads) + sample_type + treatment 
                                phyloseq_rel_nz_bal <- subset_samples(phyloseq_rel_nz, sample_type == "BAL")
                                sample_data_sample <- sample_data(phyloseq_rel_nz_bal) %>% data.frame()
                                otu_data_sample <- otu_table(phyloseq_rel_nz_bal) %>% t %>% data.frame()
                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                
                                fit_data_bal = Maaslin2(
                                        input_data = otu_data_sample, 
                                        input_metadata = sample_data_sample, 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads.","lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # default
                                        random_effects = c("original_sample"))
                                
                                fit_data_bal$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()
                                
                                #sputum
                                # # y ~ log(final reads) + sample_type + treatment 
                                phyloseq_rel_nz_spt <- subset_samples(phyloseq_rel_nz, sample_type == "Sputum")
                                sample_data_sample <- sample_data(phyloseq_rel_nz_spt) %>% data.frame()
                                otu_data_sample <- otu_table(phyloseq_rel_nz_spt) %>% t %>% data.frame()
                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                
                                fit_data_spt = Maaslin2(
                                        input_data = otu_data_sample, 
                                        input_metadata = sample_data_sample, 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads.","lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # default
                                        random_effects = c("original_sample"))
                                
                                fit_data_spt$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()
                                
                                sheets <- list(
                                        "all_sample" = fit_data$results,
                                        "nasal_swab" =fit_data_ns$results,
                                        "bal" = fit_data_bal$results,
                                        "sputum" = fit_data_spt$results)
                                write_xlsx(sheets, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/6_Results/host_depleiton/maaslin/DAT_20221207_MGK_maaslin_log(Final_reads)+sample_type+treatment.xlsx")
                                
                                
                                
                                

                                
                                

# #Overlap between DA anaysis by treatment and potential contaminant --------

                                
                                
                                fit_data_1_result_sig_host <- subset(fit_data_1_result_sig, fit_data_1_result_sig$value == "host_zero"  & fit_data_1_result_sig$qval < 0.1 )
                                fit_data_1_result_sig_qiaamp <- subset(fit_data_1_result_sig, fit_data_1_result_sig$value == "qiaamp" & fit_data_1_result_sig$qval < 0.1 )
                                fit_data_1_result_sig_ben <- subset(fit_data_1_result_sig, fit_data_1_result_sig$value == "benzonase" & fit_data_1_result_sig$qval < 0.1 )
                                fit_data_1_result_sig_lypma <- subset(fit_data_1_result_sig, fit_data_1_result_sig$value == "lypma" & fit_data_1_result_sig$qval < 0.1 )
                                fit_data_1_result_sig_molysis <- subset(fit_data_1_result_sig, fit_data_1_result_sig$value == "molysis" & fit_data_1_result_sig$qval < 0.1 )
                                
                                
                                
                                
                                
                                fit_data_1_result_sig_host$decontam_freq <- fit_data_1_result_sig_host$feature %in% c(row.names(contamdf.freq_true))
                                fit_data_1_result_sig_host$decontam_prev <- fit_data_1_result_sig_host$feature %in% c(row.names(contamdf.prev_true))
                                fit_data_1_result_sig_host$decontam_comb <- fit_data_1_result_sig_host$feature %in% c(row.names(contamdf.comb_true))
                                
                                
                                fit_data_1_result_sig_qiaamp$decontam_freq <- fit_data_1_result_sig_qiaamp$feature %in% c(row.names(contamdf.freq_true))
                                fit_data_1_result_sig_qiaamp$decontam_prev <- fit_data_1_result_sig_qiaamp$feature %in% c(row.names(contamdf.prev_true))
                                fit_data_1_result_sig_qiaamp$decontam_comb <- fit_data_1_result_sig_qiaamp$feature %in% c(row.names(contamdf.comb_true))
                                
                                fit_data_1_result_sig_host %>% select(c("feature", "metadata", "coef", "qval", "decontam_freq", "decontam_prev", "decontam_comb"))
                                fit_data_1_result_sig_qiaamp %>% select(c("feature", "metadata", "coef", "qval", "decontam_freq", "decontam_prev", "decontam_comb"))
                                
                                                               
                                #Sensitivity analysis -
                                        #interaction term
                                
                                
                                
                                
                                        #BAL and sputum
                                        #BAL
                                                phyloseq_rel_bal <- subset_samples(phyloseq_rel_nz, sample_type == "BAL") 
                                                sample_data_sample <- sample_data(phyloseq_rel_bal) %>% data.frame()
                                                otu_data_sample <- otu_table(phyloseq_rel_nz) %>% t %>% data.frame()
                                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                                
                                                fit_data_bal = Maaslin2(
                                                        input_data = otu_data_sample, 
                                                        input_metadata = sample_data_sample, 
                                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                                        fixed_effects = c("Final_reads", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                                        random_effects = c("original_sample"), 
                                                        reference = c("sample_type,BAL"))
                                                ggplot(fit_data_bal$results, aes(y = -log10(qval), x = coef, col = metadata)) +
                                                        geom_point() +
                                                        xlab("MaAslin coefficient") +
                                                        ylab("-log10(q-value)") +
                                                        geom_hline(yintercept = 1, col = "blue") +
                                                        geom_vline(xintercept = 0, col = "blue") +
                                                        geom_text(aes( -1, 2, label = "Line: q = 0.1, fc = 0", vjust = -1), col = "blue", size = 4)
                                                
                                                
                                        #sputum
                                                phyloseq_rel_sputum <- subset_samples(phyloseq_rel_nz, sample_type == "Sputum") 
                                                sample_data_sample <- sample_data(phyloseq_rel_sputum) %>% data.frame()
                                                otu_data_sample <- otu_table(phyloseq_rel_nz) %>% t %>% data.frame()
                                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                                
                                                fit_data_sputum = Maaslin2(
                                                        input_data = otu_data_sample, 
                                                        input_metadata = sample_data_sample, 
                                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                                        fixed_effects = c("Final_reads", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                                        random_effects = c("original_sample"), 
                                                        reference = c("sample_type,BAL"))
                                                
                                                ggplot(fit_data_sputum$results, aes(y = -log10(qval), x = coef, col = metadata)) +
                                                        geom_point() +
                                                        xlab("MaAslin coefficient") +
                                                        ylab("-log10(q-value)") +
                                                        geom_hline(yintercept = 1, col = "blue") +
                                                        geom_vline(xintercept = 0, col = "blue") +
                                                        geom_text(aes( -1, 2, label = "Line: q = 0.1, fc = 0", vjust = -1), col = "blue", size = 4)
                                                
                                        
                                                fit_data_sputum$results %>% filter(qval < 0.1) %>% .$metadata %>% table
                                        
                                        #sputum and BAL
                                                phyloseq_rel_sputum <- subset_samples(phyloseq_rel_nz, sample_type != "nasal_swab") 
                                                sample_data_sample <- sample_data(phyloseq_rel_sputum) %>% data.frame()
                                                otu_data_sample <- otu_table(phyloseq_rel_nz) %>% t %>% data.frame()
                                                sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")
                                                
                                                fit_data_sputum_bal = Maaslin2(
                                                        input_data = otu_data_sample, 
                                                        input_metadata = sample_data_sample, 
                                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                                        fixed_effects = c("sample_type", "Final_reads", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                                        random_effects = c("original_sample"), 
                                                        reference = c("sample_type,BAL"))
                                                
                                                ggplot(fit_data_sputum_bal$results, aes(y = -log10(qval), x = coef, col = metadata)) +
                                                        geom_point() +
                                                        xlab("MaAslin coefficient") +
                                                        ylab("-log10(q-value)") +
                                                        geom_hline(yintercept = 1, col = "blue") +
                                                        geom_vline(xintercept = 0, col = "blue") +
                                                        geom_text(aes( -1, 2, label = "Line: q = 0.1, fc = 0", vjust = -1), col = "blue", size = 4)
                                                
                                                fit_data_sputum_bal$results %>% filter(qval < 0.1) %>% .$metadata %>% table
                                                
                                                
                                                
                                ns_sum_abundance <- phyloseq_ns %>% transform_sample_counts(function(x) {x/sum(x)}) %>% otu_table %>% rowSums()
                                ns_prevalnece <- phyloseq_ns %>% otu_table != 0
                                ns_prevalnece <- ns_prevalnece %>% rowSums()
                                decontam_comb1 <- names(ns_prevalnece) %in% row.names(contamdf.comb_true)
                                
                                bal_sum_abundance <- phyloseq_rel_bal %>% otu_table %>% rowSums()
                                bal_prevalnece <- phyloseq_rel_bal %>% otu_table != 0
                                bal_prevalnece <- bal_prevalnece %>% rowSums()
                                decontam_comb2 <- names(bal_prevalnece) %in% row.names(contamdf.comb_true)
                                
                                sputum_sum_abundance <- phyloseq_rel_sputum %>% otu_table %>% rowSums()
                                sputum_prevalnece <- phyloseq_rel_sputum %>% otu_table != 0
                                sputum_prevalnece <- sputum_prevalnece %>% rowSums()
                                decontam_comb3 <- names(sputum_prevalnece) %in% row.names(contamdf.comb_true)
                                
                                
                                cbind (ns_sum_abundance, ns_prevalnece, decontam_comb1) %>% .[order(ns_sum_abundance, decreasing = T),] %>% head(., 25)
                                cbind (bal_sum_abundance, bal_prevalnece, decontam_comb2) %>% .[order(bal_sum_abundance, decreasing = T),] %>% head(., 25)
                                cbind (sputum_sum_abundance, sputum_prevalnece, decontam_comb3) %>% .[order(sputum_sum_abundance, decreasing = T),] %>% head(., 25)
                                
                                
                                

# Additional analysis - CFB / insilico mock community ---------------------

                                
                                
                                
        #Double-check:
                #1. Diversity indices after removing newly found taxa
                                phyloseq_taxa_in_control <- phyloseq_sample
                                
                                phyloseq_taxa_in_control %>% sample_data
                                #phyloseq_taxa_in_control <- subset_samples(phyloseq_taxa_in_control, S.obs != 0)
                                
                                phyloseq_taxa_in_control %>% sample_data %>% .$original_sample %>% table()
                                
                                for (i in 1:20) {
                                sample_data <- sample_data(phyloseq_taxa_in_control)
                                sample_id_temp <- sample_data$original_sample %>% table() %>% names() %>%.[i]
                                sample_data %>% data.frame %>% subset(., .$original_sample == sample_id_temp & grepl("control", sample_data$Row.names))
                                sample_data %>% data.frame %>% subset(., .$original_sample == sample_id_temp)
                                control_group <- ifelse((sample_data %>% data.frame %>% subset(., .$original_sample == sample_id_temp & grepl("control", sample_data$Row.names)) %>% .$S.obs) == 0,
                                       "qiaamp", "control")
                                phyloseq_dummy <- subset_samples(phyloseq_sample, original_sample == sample_id_temp)
                                phyloseq_taxa_in_control <- subset_samples(phyloseq_taxa_in_control, !(original_sample %in% sample_data(phyloseq_dummy)$original_sample))
                                number_control <- otu_table(phyloseq_dummy) %>% data.frame() %>% names() %>% grep(., pattern = control_group)
                                number_control <- names(otu_table(phyloseq_dummy) %>% data.frame())[number_control]
                                otu_table(phyloseq_dummy)[, number_control] %>% subset(., .[,1] != 0) %>% taxa_names()
                                control_taxa <- otu_table(phyloseq_dummy)[, number_control] %>% subset(., .[,1] != 0) %>% taxa_names()
                                phyloseq_dummy <- subset_taxa(phyloseq_dummy, taxa_names(phyloseq_dummy) %in% control_taxa)
                                phyloseq_taxa_in_control <- merge_phyloseq(phyloseq_dummy, phyloseq_taxa_in_control)
                                }
                                
                                read_sum <- phyloseq_taxa_in_control %>% otu_table() %>% rowSums() %>% data.frame()
                                non_zero_taxa <- subset(read_sum, read_sum$. !=0 ) %>% row.names()
                                phyloseq_taxa_in_control <- subset_taxa(phyloseq_taxa_in_control, taxa_names(phyloseq_taxa_in_control) %in% non_zero_taxa)
                                sample_data(phyloseq_taxa_in_control) <- sample_data(phyloseq)
                                sample_data(phyloseq_taxa_in_control) <- sample_data(alpha_diversity(phyloseq_taxa_in_control)) 
                                
                                
                                
                                facet_figure_z(phyloseq_taxa_in_control, c("species_richness", "data_shannon", "data_invsimpson", "berger_parker"))
                                phyloseq_taxa_in_control = transform_sample_counts(phyloseq_taxa_in_control, function (x) {x/sum(x)})
                                phyloseq_taxa_in_control = subset_samples(phyloseq_taxa_in_control, S.obs != 0)
                                
                                
                                ordinate(phyloseq_taxa_in_control,  method = "NMDS", distance = "bray") %>% plot_ordination(phyloseq_taxa_in_control, ., col = "original_sample", shape = "treatment" ) + geom_point(size = 4)
                                ordinate(phyloseq_taxa_in_control,  method = "NMDS", distance = "bray") %>% plot_ordination(phyloseq_taxa_in_control, ., col = "sample_type", shape = "treatment" ) + geom_point(size = 4)
                                
                                vegan::adonis2(distance(phyloseq_taxa_in_control, method="bray") ~ sample_type + lypma + benzonase + host_zero + molysis + qiaamp + original_sample,
                                               data = phyloseq_taxa_in_control %>% sample_data %>% data.frame(check.names = F),
                                               permutations = 10000)
                                
                                #Almost the same trend is observed.

        #Double-check:
                #2. CF project - reads
                                
                                
                                CFB_phyloseq <- readRDS("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/DAT_20220719_JMV_CFB_Sputum_Microbial_AbsAbund.rds")
                                CFB_phyloseq %>% sample_data
                                sputum_origin <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220407_MGK_sputum_BAL_host_depletion.xlsx", sheet = 2) %>% data.frame()
                                sputum_origin$ID <- gsub("_", ".", sputum_origin$ID)
                                row.names(sputum_origin) <- sputum_origin$ID
                                
                                CFB_a <- subset_samples(CFB_phyloseq, grepl(paste(row.names(subset(sputum_origin, sputum_origin$aliqout_group == "a")), collapse = "|"), sample_names(CFB_phyloseq)))
                                CFB_b <- subset_samples(CFB_phyloseq, grepl(paste(row.names(subset(sputum_origin, sputum_origin$aliqout_group == "b")), collapse = "|"), sample_names(CFB_phyloseq)))
                                CFB_c <- subset_samples(CFB_phyloseq, grepl(paste(row.names(subset(sputum_origin, sputum_origin$aliqout_group == "c")), collapse = "|"), sample_names(CFB_phyloseq)))
                                CFB_d <- subset_samples(CFB_phyloseq, grepl(paste(row.names(subset(sputum_origin, sputum_origin$aliqout_group == "d")), collapse = "|"), sample_names(CFB_phyloseq)))
                                CFB_e <- subset_samples(CFB_phyloseq, grepl(paste(row.names(subset(sputum_origin, sputum_origin$aliqout_group == "e")), collapse = "|"), sample_names(CFB_phyloseq)))
                                
                                CFB_a <- CFB_a %>% otu_table() %>% colSums %>% data.frame
                                CFB_b <- CFB_b %>% otu_table() %>% colSums %>% data.frame
                                CFB_c <- CFB_c %>% otu_table() %>% colSums %>% data.frame
                                CFB_d <- CFB_d %>% otu_table() %>% colSums %>% data.frame
                                CFB_e <- CFB_e %>% otu_table() %>% colSums %>% data.frame
                                CFB_control <- cbind(CFB_a, CFB_b, CFB_c, CFB_d, CFB_e)
                                names(CFB_control) <- c("CFB_a", "CFB_b", "CFB_c", "CFB_d", "CFB_e")
                                
                                CFB_control_surv <- CFB_control > 0
                                dummy <- CFB_control_surv %>% colSums() %>% data.frame()
                                dummy$treatment <- "CFB_sequencing"
                                dummy$original_sample <- row.names(dummy)
                                names(dummy) <- c("S.obs", "treatment", "original_sample")
                                sample_data <- sample_data(phyloseq_sample) %>% data.frame()
                                sample_data <- sample_data %>% subset(., .$sample_type == "Sputum") 
                                sample_data <- select(sample_data, c("treatment", "S.obs", "original_sample"))
                                sample_data <- rbind(dummy, sample_data)
                                
                                sample_data$treatment <- factor(sample_data$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp", "CFB_sequencing"))
                                lmer(S.obs ~ treatment + (1|original_sample), data = sample_data) %>% summary
                                lm(S.obs ~ treatment, data = sample_data) %>% summary
                                
                                sample_data %>% 
                                        ggplot(aes(x = treatment, y = S.obs)) +
                                        geom_boxplot()
                                
                                
                                
                                
                                
                                
                                
                                