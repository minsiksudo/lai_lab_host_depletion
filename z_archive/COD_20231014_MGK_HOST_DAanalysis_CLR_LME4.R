library(phyloseq)
library(Maaslin2)
library(tidyverse)
library(decontam)
library(microbiome)
library(qvalue)
library(lme4)



# Loading data ------------------------------------------------------------


phyloseq <- readRDS("/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20230521_MGK_host_tidy.rds")

# Adding alpha diversity indices ------------------------------------------


alpha_diversity <- function(data) {
        otu_table <- otu_table(data) %>% .[colSums(.) !=0]
        S.obs <- rowSums(t(otu_table) != 0)
        sample_data <- sample_data(data)
        data_evenness <- vegan::diversity(t(otu_table)) / log(vegan::specnumber(t(otu_table))) # calculate evenness index using vegan package
        data_shannon <- vegan::diversity(t(otu_table), index = "shannon") # calculate Shannon index using vegan package
        data_hill <- exp(data_shannon)                           # calculate Hills index
        
        data_dominance <- microbiome::dominance(otu_table, index = "all", rank = 1, aggregate = TRUE) # dominance (Berger-Parker index), etc.
        data_invsimpson <- vegan::diversity(t(otu_table), index = "invsimpson")                          # calculate Shannon index using vegan package
        alpha_diversity <- cbind(S.obs, data_shannon, data_hill, data_invsimpson, data_evenness,data_dominance) # combine all indices in one data table
        sample_data <- merge(data.frame(sample_data), alpha_diversity, by = 0, all = T) %>% column_to_rownames(var = "Row.names")
}

sample_data(phyloseq$phyloseq_rel) <- sample_data(alpha_diversity(phyloseq$phyloseq_count))
sample_data(phyloseq$phyloseq_count) <- sample_data(alpha_diversity(phyloseq$phyloseq_count)) 
sample_data(phyloseq$phyloseq_path_rpk) <- sample_data(alpha_diversity(phyloseq$phyloseq_path_rpk))

#Chaning name of function taxa
#taxa_names(phyloseq$phyloseq_path_rpk) <- tax_table(phyloseq$phyloseq_path_rpk) %>% data.frame %>% .$group

phyloseq_unfiltered <- phyloseq
# Decontam ----------------------------------------------------------------
sample_data(phyloseq_unfiltered$phyloseq_rel)$is.neg <- grepl("Neg", sample_data(phyloseq_unfiltered$phyloseq_rel)$sample_type)

phyloseq_decontam_bal <- phyloseq_unfiltered$phyloseq_rel %>% 
        subset_samples(S.obs != 0) %>%
        subset_samples(sample_type == "Neg." | sample_type == "BAL")
phyloseq_decontam_ns <- phyloseq_unfiltered$phyloseq_rel %>% 
        subset_samples(S.obs != 0) %>%
        subset_samples(sample_type == "Neg." | sample_type == "Nasal")
phyloseq_decontam_spt <- phyloseq_unfiltered$phyloseq_rel %>% 
        subset_samples(S.obs != 0) %>%
        subset_samples(sample_type == "Neg." | sample_type == "Sputum")


contaminant_combined_bal <- 
        data.frame("BAL", fix.empty.names = F, 
                   isContaminant(phyloseq_decontam_bal, method="combined", neg = "is.neg", threshold = 0.1, conc = "DNA_bac_ng_uL") %>% subset(.,.$contaminant) %>% row.names
        )

contaminant_combined_ns <- 
        data.frame("Nasal swab", fix.empty.names = F, 
                   isContaminant(phyloseq_decontam_ns, method="combined", neg = "is.neg", threshold = 0.1, conc = "DNA_bac_ng_uL") %>% subset(.,.$contaminant) %>% row.names
        )


contaminant_combined_spt <- 
        data.frame("Sputum", fix.empty.names = F, 
                   isContaminant(phyloseq_decontam_spt, method="combined", neg = "is.neg", threshold = 0.1, conc = "DNA_bac_ng_uL") %>% subset(.,.$contaminant) %>% row.names
        )



contaminants <- rbind(contaminant_combined_bal, contaminant_combined_ns, contaminant_combined_spt)

names(contaminants) <- c("Sample type", "Taxa")

species_italic <- function(data) {
        names <- gsub("_", " ", data)
        names <- gsub("[]]|[[]", "", names)
        names <- gsub(" sp", " sp.", names)
        names <- gsub(" sp.", "* sp.", names)
        names <- gsub(" group", "* group.", names)
        names <- ifelse(grepl("[*]", names), paste("*", names, sep = ""), paste("*", names, "*", sep = ""))
        names
}


contaminants %>% 
        mutate(Taxa = species_italic(Taxa))


contaminant <- contaminants$Taxa



# Filtering samples failed sequencing -------------------------------------

phyloseq$phyloseq_count <- subset_samples(phyloseq$phyloseq_count, S.obs != 0 & sample_type %in% c("BAL", "Nasal", "Sputum", "Mock"))
phyloseq$phyloseq_rel <- subset_samples(phyloseq$phyloseq_rel, S.obs != 0 & sample_type %in% c("BAL", "Nasal", "Sputum", "Mock"))


# Prevalence filtering  ---------------------------------------------------

phyloseq_unfiltered$phyloseq_rel <- transform_sample_counts(phyloseq_unfiltered$phyloseq_rel,
                                                            function(x){x/sum(x)})
taxa_qc <- data.frame("species" =
                              otu_table(
                                      subset_samples(
                                              phyloseq_unfiltered$phyloseq_rel,S.obs != 0 &
                                                      sample_type %in% c("Mock", "BAL", "Nasal", "Sputum"))) %>%
                              t() %>% colnames(),
                      "prevalence" =
                              ifelse(subset_samples(phyloseq_unfiltered$phyloseq_rel, S.obs != 0 & 
                                                            sample_type %in% c("Mock", "BAL", "Nasal", "Sputum")) %>%
                                             otu_table() > 0, 1, 0)%>% 
                              t() %>%
                              colSums(), #Prevalence of taxa
                      "mean_rel_abd" = 
                              subset_samples(phyloseq_unfiltered$phyloseq_rel,
                                             S.obs != 0 & sample_type %in% c("Mock", "BAL", "Nasal", "Sputum")) %>%
                              otu_table() %>%
                              t() %>%
                              colMeans(na.rm = T) #mean relativ abundacne 
)


function_qc <- data.frame("function" =
                                  otu_table(
                                          subset_samples(
                                                  phyloseq_unfiltered$phyloseq_path_rpk,
                                                  S.obs != 0 & 
                                                          sample_type %in% 
                                                          c("Mock", "BAL", "Nasal", "Sputum")
                                          )
                                  ) %>% 
                                  t() %>% 
                                  colnames(),
                          "prevalence" = 
                                  ifelse(subset_samples(phyloseq_unfiltered$phyloseq_path_rpk,
                                                        S.obs != 0 & 
                                                                sample_type %in% 
                                                                c("Mock", "BAL", "Nasal", "Sputum")
                                  ) %>% 
                                          otu_table() > 0, 
                                  1, 
                                  0
                                  ) %>% 
                                  t() %>% 
                                  colSums(), #Prevalence of taxa
                          "mean_rpk" = 
                                  subset_samples(phyloseq_unfiltered$phyloseq_path_rpk, 
                                                 S.obs != 0 & 
                                                         sample_type %in% 
                                                         c("Mock", "BAL", "Nasal", "Sputum")
                                  ) %>% 
                                  otu_table() %>% 
                                  t() %>% 
                                  colMeans(na.rm = T), #mean relativ abundacne 
                          unidentified = 
                                  ifelse((subset_samples(phyloseq_unfiltered$phyloseq_path_rpk,
                                                         S.obs != 0 & sample_type %in% c("Mock", "BAL", "Nasal", "Sputum")) %>% 
                                                  otu_table() > 0) %>%
                                                 row.names() %in% c("UNMAPPED", "UNINTEGRATED")
                                         , 1, 0)
)


# Making a list of filters ------------------------------------------------


red_flag_taxa <- data.frame(species = taxa_qc$species,
                            prevalence = taxa_qc$prevalence,
                            mean_rel_abd = taxa_qc$mean_rel_abd,
                            red_flag_prev_abd = 
                                    ifelse(taxa_qc$prevalence < 
                                                   otu_table(
                                                           subset_samples(
                                                                   phyloseq_unfiltered$phyloseq_rel,
                                                                   S.obs != 0 & 
                                                                           sample_type %in% 
                                                                           c("Mock", "BAL", "Nasal", "Sputum"))) %>%
                                                   t %>% rownames() %>%
                                                   length * 0.05 &
                                                   #Removing taxa with zero prevalence - taxa from nasal swabs
                                                   taxa_qc$mean_rel_abd <
                                                   taxa_qc %>%
                                                   subset(., .$prevalence != 0) %>%
                                                   .$mean_rel_abd %>%
                                                   quantile(., 0.75), 1,0),
                            red_flag_prev =
                                    ifelse(taxa_qc$prevalence < 
                                                   otu_table(
                                                           subset_samples(
                                                                   phyloseq_unfiltered$phyloseq_rel,
                                                                   S.obs != 0 & 
                                                                           sample_type %in% 
                                                                           c("Mock", "BAL", "Nasal", "Sputum"))) %>%
                                                   t %>% rownames() %>%
                                                   length * 0.05,
                                           1,
                                           0)) %>%
        mutate(red_flag_decontam = species %in% (contaminants$Taxa %>% unique()))

subset(red_flag_taxa, red_flag_taxa$red_flag_prev == 1 & red_flag_taxa$red_flag_prev_abd == 0)

#Unampped function were removed

red_flag_function <- 
        data.frame(function. = function_qc$function.,
                   prevalence = function_qc$prevalence,
                   mean_rel_abd = function_qc$mean_rpk,
                   red_flag_prev_abd = 
                           ifelse(function_qc$prevalence < 
                                          otu_table(
                                                  subset_samples(
                                                          phyloseq_unfiltered$phyloseq_path_rpk,
                                                          S.obs != 0 & 
                                                                  sample_type %in% 
                                                                  c("Mock", "BAL", "Nasal", "Sputum"))) %>%
                                          t %>% 
                                          rownames() %>%
                                          length * 0.05 &
                                          #Removing taxa with zero prevalence - taxa from nasal swabs
                                          function_qc$mean_rpk <
                                          function_qc %>%
                                          subset(., .$prevalence != 0) %>%
                                          .$mean_rpk %>%
                                          quantile(., 0.75), 1,0),
                   red_flag_prev =
                           ifelse(function_qc$prevalence <
                                          otu_table(
                                                  subset_samples(
                                                          phyloseq_unfiltered$phyloseq_path_rpk,
                                                          S.obs != 0 & 
                                                                  sample_type %in% 
                                                                  c("Mock", "BAL", "Nasal", "Sputum"))) %>%
                                          t %>% 
                                          rownames() %>%
                                          length * 0.05,
                                  1,
                                  0)) %>%
        mutate(red_flag_prev_abd = case_when(function. %in% c("UNMAPPED", "UNINTEGRATED") ~ 1,
                                             .default = red_flag_prev_abd))

subset(red_flag_function, red_flag_function$red_flag_prev == 1 & red_flag_function$red_flag_prev_abd == 0)

# Adding variables for MaAsLin --------------------------------------------

sample_data(phyloseq$phyloseq_rel)$log10.Final_reads <- log10(sample_data(phyloseq$phyloseq_rel)$Final_reads)
sample_data(phyloseq$phyloseq_rel)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_rel)$sample_type, sample_data(phyloseq$phyloseq_rel)$treatment, sep = ":")
sample_data(phyloseq$phyloseq_count)$log10.Final_reads <- log10(sample_data(phyloseq$phyloseq_count)$Final_reads)
sample_data(phyloseq$phyloseq_count)$log10.Final_reads <- log10(sample_data(phyloseq$phyloseq_count)$Final_reads)
sample_data(phyloseq$phyloseq_path_rpk)$log10.Final_reads <- log10(sample_data(phyloseq$phyloseq_path_rpk)$Final_reads)
sample_data(phyloseq$phyloseq_path_rpk)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_path_rpk)$sample_type, sample_data(phyloseq$phyloseq_path_rpk)$treatment, sep = ":")


# Making new phyloseq object ----------------------------------------------

phyloseq$phyloseq_count_filtered <- prune_taxa(subset(red_flag_taxa,
                                                      red_flag_taxa$red_flag_prev_abd != 1)$species,
                                               phyloseq$phyloseq_count)

phyloseq$phyloseq_path_rpk_filtered <- prune_taxa(subset(red_flag_taxa, red_flag_function$red_flag_prev_abd != 1)$function., phyloseq$phyloseq_path_rpk)
phyloseq$phyloseq_rel_filtered <- prune_taxa(subset(red_flag_taxa,
                                                    red_flag_taxa$red_flag_prev_abd != 1)$species,
                                             phyloseq$phyloseq_rel) %>% 
        transform_sample_counts(., function (x) {x/sum(x)})


# CLR transform with `microbiome` package ---------------------------------



# Using phyloseq$phyloseq_rel_filtered, do CLR transform 

physeq0_clr_all <- phyloseq$phyloseq_count_filtered %>%
        microbiome::transform(., 'clr')

        #Do I have to add pseudo counts?
                physeq0_clr__pseudocnt_all <- phyloseq$phyloseq_count_filtered %>% {
                        otu_table(phyloseq$phyloseq_count_filtered) <- phyloseq$phyloseq_count_filtered %>% otu_table + 1
                        phyloseq$phyloseq_count_filtered
                } %>%
                        microbiome::transform(., 'clr')
                physeq0_clr__pseudocnt_all %>% otu_table
                
                #It seems like adding pseudo counts does not change the association, as they were `centered`.
                
        
        #How about using relative abundance?
                physeq0_clr_all_rel <- phyloseq$phyloseq_rel_filtered %>%
                microbiome::transform(., 'clr')
        
              
                
# CLR transformation was conducted after stratifying the data by sample type



physeq0_clr_spt <- subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "Sputum") %>% 
        prune_taxa(taxa_sums(.) > 0, .) %>%
        microbiome::transform(., 'clr')

physeq0_clr_bal <- subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "BAL") %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        microbiome::transform(., 'clr')

physeq0_clr_ns <- subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "Nasal") %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        microbiome::transform(., 'clr')

#Sanity check
        # Original data
        #phyloseq$phyloseq_count_filtered %>% otu_table %>% t %>% data.frame() %>% view
        # CLR transformed data, all samples at once 
        #physeq0_clr_all %>% otu_table %>% t %>% data.frame() %>% view
        # CLR transformed data, partial samples
        #physeq0_clr_spt %>% otu_table

#It seems like CLR-transfomred data may result in changes from 0-0 reads.


# LME4 for CLR normalized data --------------------------------------------

#Running LME4 for each taxa

#BAL


lmer_taxa_bal <- data.frame()

a <- physeq0_clr_bal %>% taxa_sums() %>% length()

for(i in 1:a) {
        #Creating a data frame that includes CLR transformed data of i-th bug.
        #making differnt otu tables for each sample type
        otu_table_bal <- physeq0_clr_bal %>% otu_table
        
        #tax table for different sample type
        taxa_data_bal <- otu_table_bal[i] %>% t %>% data.frame()
        
        #Making a merged dataframe having sample data and CLR transformed output
        lme_data_bal <- merge(taxa_data_bal, sample_data(physeq0_clr_bal), by = 0) %>% column_to_rownames("Row.names")
        
        #generating a character of formula.
        #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
        lme4_formula <- paste(names(taxa_data_bal), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
        
        #BAL stratified analysis
        BAL_result <- lme4::lmer(formula = lme4_formula,
                                 data = lme_data_bal) %>% 
                lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                summary %>%
                .$coefficients %>%
                data.frame(check.names = F) %>% 
                mutate(feature = names(taxa_data_bal),
                       sample_type = "BAL") %>% 
                rownames_to_column("metadata")  %>%
                subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        
        #row binding all the associations of i-th taxa to one data frame
        lmer_taxa_bal <- rbind(lmer_taxa_bal,
                           BAL_result) %>% 
                remove_rownames()
        
}

#Nasal

lmer_taxa_ns <- data.frame()

b <- physeq0_clr_ns %>% taxa_sums() %>% length()


for(i in 1:b) {
        #Creating a data frame that includes CLR transformed data of i-th bug.
        #making differnt otu tables for each sample type
        otu_table_ns <- physeq0_clr_ns %>% otu_table
        
        #tax table for different sample type
        taxa_data_ns <- otu_table_ns[i] %>% t %>% data.frame()
        
        #Making a merged dataframe having sample data and CLR transformed output
        lme_data_ns <- merge(taxa_data_ns, sample_data(physeq0_clr_ns), by = 0) %>% column_to_rownames("Row.names")
        
        #generating a character of formula.
        #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
        lme4_formula <- paste(names(taxa_data_ns), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
        
        #BAL stratified analysis
        Nasal_result <- lme4::lmer(formula = lme4_formula,
                                 data = lme_data_ns) %>% 
                lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                summary %>%
                .$coefficients %>%
                data.frame(check.names = F) %>% 
                mutate(feature = names(taxa_data_ns),
                       sample_type = "Nasal") %>% 
                rownames_to_column("metadata")  %>%
                subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        
        #row binding all the associations of i-th taxa to one data frame
        lmer_taxa_ns <- rbind(lmer_taxa_ns,
                           Nasal_result) %>% 
                remove_rownames()
        

}

lmer_taxa_spt <- data.frame()

c <- physeq0_clr_spt %>% taxa_sums() %>% length()


for(i in 1:c) {
        #Creating a data frame that includes CLR transformed data of i-th bug.
        
        #making differnt otu tables for each sample type
        otu_table_spt <- physeq0_clr_spt %>% otu_table
        
                #tax table for different sample type
        taxa_data_spt <- otu_table_spt[i] %>% t %>% data.frame()
        
        #Making a merged dataframe having sample data and CLR transformed output
        lme_data_spt <- merge(taxa_data_spt, sample_data(physeq0_clr_spt), by = 0) %>% column_to_rownames("Row.names")
        
        #generating a character of formula.
        #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
        lme4_formula <- paste(names(taxa_data_spt), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
        
        #stratified analysis
        Sputum_result <- lme4::lmer(formula = lme4_formula,
                                   data = lme_data_spt) %>% 
                lmerTest::as_lmerModLmerTest() %>%
                summary %>%
                .$coefficients %>%
                data.frame(check.names = F) %>% 
                mutate(feature = names(taxa_data_spt),
                       sample_type = "Sputum") %>% 
                rownames_to_column("metadata")  %>%
                subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        
        #row binding all the associations of i-th taxa to one data frame
        lmer_taxa_spt <- rbind(lmer_taxa_spt,
                           Sputum_result) %>% 
                remove_rownames()
        
        
}

#Calculating q-value

lmer_taxa_bal <- lmer_taxa_bal %>% subset(., .$sample_type == "BAL") %>%
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  

lmer_taxa_ns <-  lmer_taxa_ns %>% subset(., .$sample_type == "Nasal") %>%
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  

lmer_taxa_spt <-  lmer_taxa_spt %>% subset(., .$sample_type == "Sputum") %>%
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  





lmer_taxa_bal %>% subset(., .$q_val < 0.1 )

lmer_taxa_ns %>% subset(., .$q_val < 0.1 ) %>% .$feature %>% table

lmer_taxa_spt %>% subset(., .$q_val < 0.1 ) %>% .$feature %>% table
lmer_taxa_spt %>% subset(., .$q_val < 0.1 ) 

lmer_taxa_ns %>% subset(., .$feature == "Staphylococcus_aureus" & .$metadata == "lypma" & .$sample_type == "Nasal") 
maaslin_ns %>% subset(., .$feature == "Staphylococcus_aureus" & .$metadata == "lypma") 

# In this thread, MaAsLin showed different result with a CLR transformed data tested with GLM.
# https://forum.biobakery.org/t/different-results-running-glm-and-maaslin2-using-same-methods-transformations/1327/3



# MaAslin with CLR --------------------------------------------------------



#Running maaslin to compare with for-loop

        #for bal

capture.output(
        maaslin_all2 = 
                Maaslin2(input_data = subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "BAL") %>%
                                 otu_table %>% t %>% data.frame(),
                         input_metadata = subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "BAL") %>% 
                                 sample_data %>% data.frame(check.names = F), 
                         output = "/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw", 
                         fixed_effects = c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                         transform = "NONE", 
                         # CLR does the normalization and transformation at once.
                         #https://forum.biobakery.org/t/questions-about-normalization-method-clr-vs-tss/4352/2
                         normalization = "CLR",
                         random_effects = c("subject_id"), 
                         #reference = c("sample_type,BAL"), 
                         plot_heatmap = F,
                         plot_scatter = F))

maaslin_bal <- read.csv("/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw/all_results.tsv", sep = "\t")
maaslin_bal_transformed_otu <- read.csv("/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw/features/filtered_data_norm_transformed.tsv",
                                       sep = "\t")

#feature ~ lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)
#for nasal

capture.output(
        maaslin_all2 = 
                Maaslin2(input_data = subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "Nasal") %>%
                                 otu_table %>% t %>% data.frame(),
                         input_metadata = subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "Nasal") %>% 
                                 sample_data %>% data.frame(check.names = F), 
                         output = "/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw", 
                         fixed_effects = c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                         transform = "NONE", 
                         # CLR does the normalization and transformation at once.
                         #https://forum.biobakery.org/t/questions-about-normalization-method-clr-vs-tss/4352/2
                         normalization = "CLR",
                         random_effects = c("subject_id"),
                         plot_heatmap = F,
                         plot_scatter = F))

maaslin_ns <- read.csv("/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw/all_results.tsv", sep = "\t")
maaslin_ns_transformed_otu <- read.csv("/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw/features/filtered_data_norm_transformed.tsv",
                                        sep = "\t")



#for sputum
capture.output(
        maaslin_all2 = 
                Maaslin2(input_data = subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "Sputum") %>%
                                 otu_table %>% t %>% data.frame(),
                         input_metadata = subset_samples(phyloseq$phyloseq_count_filtered, sample_type == "Sputum") %>% 
                                 sample_data %>% data.frame(check.names = F), 
                         output = "/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw", 
                         fixed_effects = c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                         transform = "NONE", 
                         # CLR does the normalization and transformation at once.
                         #https://forum.biobakery.org/t/questions-about-normalization-method-clr-vs-tss/4352/2
                         normalization = "CLR",
                         random_effects = c("subject_id"),
                         plot_heatmap = F,
                         plot_scatter = F))

maaslin_sputum <- read.csv("/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw/all_results.tsv", sep = "\t")
maaslin_spt_transformed_otu <- read.csv("/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw/features/filtered_data_norm_transformed.tsv",
                                        sep = "\t")



#Comparing lme4 otuput with maaslin output

        #with BAL
sig_bal %>% subset(., .$q_val < 0.1) %>% .$taxa %>% unique
maaslin_bal %>% subset(., .$qval <0.1) %>% .$feature %>% unique


#Associations are different....Why?
#Maaslin may use different CLR calculation with microbiome::transform().
        #Maaslin2 normalization
        maaslin_spt_transformed_otu <- read.csv("/Users/minsikkim/Dropbox/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_raw/features/filtered_data_norm_transformed.tsv",
                 sep = "\t")
        #microbiome normalization - count
        clr_spt_transformed_otu <- otu_table(physeq0_clr_spt) %>% t %>% data.frame()
        #microbiome normalization - relative abundance
        clr_spt_transformed_otu_rel <- physeq0_clr_all_rel %>% otu_table %>% t %>% data.frame()
        
        #both normalization - for read counts and relative abundances with microbiome package - showed the same output.
        #Maaslin showed slightly different data.


# Masslin normalization for LMER ------------------------------------------

        

# LMER on maaslin transformed data
        lmer_maaslin_bal <- data.frame()
        
        for(i in 1:54) {
                #Creating a data frame that includes CLR transformed data of i-th bug.
                
                #making differnt otu tables for each sample type
                otu_table_bal <- maaslin_bal_transformed_otu %>% column_to_rownames("feature")
                
                #tax table for different sample type
                taxa_data_bal <- otu_table_bal[i] %>% data.frame()
                
                #Making a merged dataframe having sample data and CLR transformed output
                lme_data_bal <- merge(taxa_data_bal, sample_data(physeq0_clr_bal), by = 0) %>% column_to_rownames("Row.names")
                
                #generating a character of formula.
                #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
                lme4_formula <- paste(names(taxa_data_bal), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
                
                #BAL stratified analysis
                BAL_result <- lme4::lmer(formula = lme4_formula,
                                         data = lme_data_bal) %>% 
                        lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                        summary %>%
                        .$coefficients %>%
                        data.frame(check.names = F) %>% 
                        mutate(taxa = names(taxa_data_bal),
                               sample_type = "BAL") %>% 
                        rownames_to_column("metadata")  %>%
                        subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
                
                
                #row binding all the associations of i-th taxa to one data frame
                lmer_maaslin_bal <- rbind(lmer_maaslin_bal,
                                   BAL_result %>% remove_rownames()
                )
                
        }
        
        lmer_maaslin_ns <- data.frame()
        
        for(i in 1:29) {
                #Creating a data frame that includes CLR transformed data of i-th bug.
                
                #making differnt otu tables for each sample type
                otu_table_ns <- maaslin_ns_transformed_otu %>% column_to_rownames("feature")
                
                #tax table for different sample type
                taxa_data_ns <- otu_table_ns[i] %>% data.frame()
                
                #Making a merged dataframe having sample data and CLR transformed output
                lme_data_ns <- merge(taxa_data_ns, sample_data(physeq0_clr_ns), by = 0) %>% column_to_rownames("Row.names")
                
                #generating a character of formula.
                #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
                lme4_formula <- paste(names(taxa_data_ns), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
                
                #BAL stratified analysis
                NS_result <- lme4::lmer(formula = lme4_formula,
                                         data = lme_data_ns) %>% 
                        lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                        summary %>%
                        .$coefficients %>%
                        data.frame(check.names = F) %>% 
                        mutate(taxa = names(taxa_data_ns),
                               sample_type = "Nasal") %>% 
                        rownames_to_column("metadata")  %>%
                        subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
                
                
                #row binding all the associations of i-th taxa to one data frame
                lmer_maaslin_ns <- rbind(lmer_maaslin_ns,
                                          NS_result %>% remove_rownames()
                )
                
        }
        
        lmer_maaslin_ns
        
        
        
        lmer_maaslin_spt <- data.frame()
        
        
        for(i in 1:155) {
                #Creating a data frame that includes CLR transformed data of i-th bug.
                
                #making differnt otu tables for each sample type
                otu_table_spt <- maaslin_spt_transformed_otu %>% column_to_rownames("feature")
                
                #tax table for different sample type
                taxa_data_spt <- otu_table_spt[i] %>% data.frame()
                
                #Making a merged dataframe having sample data and CLR transformed output
                lme_data_spt <- merge(taxa_data_spt, sample_data(physeq0_clr_spt), by = 0) %>% column_to_rownames("Row.names")
                
                #generating a character of formula.
                #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
                lme4_formula <- paste(names(taxa_data_spt), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
                
                #BAL stratified analysis
                Sputum_result <- lme4::lmer(formula = lme4_formula,
                                        data = lme_data_spt) %>% 
                        lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                        summary %>%
                        .$coefficients %>%
                        data.frame(check.names = F) %>% 
                        mutate(taxa = names(taxa_data_spt),
                               sample_type = "Sputum") %>% 
                        rownames_to_column("metadata")  %>%
                        subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
                
                
                #row binding all the associations of i-th taxa to one data frame
                lmer_maaslin_spt <- rbind(lmer_maaslin_spt,
                                         Sputum_result %>% remove_rownames()
                )
                
        }
        
        
        lmer_maaslin_bal <- lmer_maaslin_bal %>% 
                mutate(q_val = qvalue::qvalue(`Pr(>|t|)`)$qvalue)
        
        lmer_maaslin_ns <- lmer_maaslin_ns %>% 
                mutate(q_val = qvalue::qvalue(`Pr(>|t|)`)$qvalues)
        
        
        lmer_maaslin_spt <- lmer_maaslin_spt %>% 
                mutate(q_val = qvalue::qvalue(`Pr(>|t|)`)$qvalue)
        
#Still showing slightly higher q-value than maaslin.
        
        maaslin_ns %>% subset(., .$qval < 0.1) %>% .$feature %>% table
        lmer_maaslin_bal %>% subset(., .$q_val < 0.1) %>% .$taxa %>% table
        
        maaslin_ns %>% subset(., .$qval < 0.1) %>% .$feature %>% table
        lmer_maaslin_ns %>% subset(., .$q_val < 0.1) %>% .$taxa %>% table
        
        maaslin_sputum %>% subset(., .$qval < 0.1) %>% .$feature %>% table
        lmer_maaslin_spt %>% subset(., .$q_val < 0.1) %>% .$taxa %>% table
        
        
#They are quite similar!
        

# Amy's transform ---------------------------------------------------------

        
#Amy's suggestion about CLR normalization
        
        
        geometric_mean <- function(x, replace_zeros = "minimum") {
                
                # pseudocounts
                if (replace_zeros == "minimum") {
                        x[x == 0] <- min(x[x > 0]) # cheap pseudocount - take minimum
                } else if (replace_zeros == "half") {
                        x[x == 0] <- min(x[x > 0])/2 # cheap pseudocount - take half minimum
                } else if (is.numeric(replace_zeros)) {
                        if (replace_zeros <= 0) stop("You've told me to replace zeroes with zero, or a nonnegative value!")
                        x[x == 0] <- replace_zeros 
                }
                if (any(x == 0)) stop("There are zeros in x and you haven't told me what to do with them!")
                
                stopifnot(all(x >= 0))
                
                ## convert to relative abundance:
                pz <- x/sum(x)
                
                prod(pz)^(1/length(pz))
                
        }
        
        clr <- function(x, replace_zeros = "minimum") {
                
                stopifnot(all(x >= 0))
                
                # pseudocounts
                if (replace_zeros == "minimum") {
                        x[x == 0] <- min(x[x > 0]) # cheap pseudocount - take minimum
                } else if (replace_zeros == "half") {
                        x[x == 0] <- min(x[x > 0])/2 # cheap pseudocount - take half minimum
                } else if (is.numeric(replace_zeros)) {
                        if (replace_zeros <= 0) stop("You've told me to replace zeroes with zero, or a nonnegative value!")
                        x[x == 0] <- replace_zeros 
                }
                if (any(x == 0)) stop("There are zeros in x and you haven't told me what to do with them!")
                
                ## convert to relative abundance:
                pz <- x/sum(x)
                
                gz <- prod(pz)^(1/length(pz))
                
                log(pz/gz)
                
        }
        
        ### examples
        #clr(c(1))
        #clr(c(0.2, 0.1, 0.7, 0), replace_zeros=0.05) 
        #clr(c(0.2, 0.1, 0.7, 0), replace_zeros=0.01)
        #clr(c(0.2, 0.1, 0.7, 0))
        #clr(c(0.2, 0.1, 0.7, 0), replace_zeros=0.1)
        #clr(c(0.2, 0.1, 0.7, 0), replace_zeros=-0.1) # good
        
        clr(otu_table(phyloseq$phyloseq_rel_filtered) %>% data.frame %>% .$lyPMA_pos_1,
            replace_zeros = "minimum") 
        
        
        clr(otu_table(phyloseq$phyloseq_rel_filtered) %>% data.frame %>% .$lyPMA_pos_1,
            replace_zeros = "half") 
        
        clr(otu_table(phyloseq$phyloseq_rel_filtered) %>% data.frame %>% .$lyPMA_pos_1,
            replace_zeros = 0.000000000000001) 
        
        clr(otu_table(phyloseq$phyloseq_count_filtered) %>% data.frame %>% .$lyPMA_pos_1,
            replace_zeros = 1) 
        
        # For now it is resulting in Inf! Something need to be fixed.
        