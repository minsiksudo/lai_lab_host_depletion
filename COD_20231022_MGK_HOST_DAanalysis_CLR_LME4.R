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


dummy <- compositions::clr(otu_table(phyloseq$phyloseq_count_filtered)) %>% data.frame()
dummy2 <- otu_table(physeq0_clr_all) %>% data.frame()


# Using phyloseq$phyloseq_rel_filtered, do CLR transform 

physeq0_clr_all <- phyloseq$phyloseq_count_filtered %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        subset_samples(., sample_type %in% c("BAL", "Nasal", "Sputum")) %>% 
        microbiome::transform(., 'clr')

#for function,
physeq0_clr_all_f <- phyloseq$phyloseq_path_rpk_filtered %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        subset_samples(., sample_type %in% c("BAL", "Nasal", "Sputum")) %>% 
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

#for function


# CLR transformation was conducted after stratifying the data by sample type

physeq0_clr_spt_f <- subset_samples(phyloseq$phyloseq_path_rpk_filtered, sample_type == "Sputum") %>% 
        prune_taxa(taxa_sums(.) > 0, .) %>%
        microbiome::transform(., 'clr')

physeq0_clr_bal_f <- subset_samples(phyloseq$phyloseq_path_rpk_filtered, sample_type == "BAL") %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        microbiome::transform(., 'clr')

physeq0_clr_ns_f <- subset_samples(phyloseq$phyloseq_path_rpk_filtered, sample_type == "Nasal") %>%
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


#with all samples

lmer_taxa_all <- data.frame()

all <- physeq0_clr_all %>% taxa_sums() %>% length()

for(i in 1:all) {
        #Creating a data frame that includes CLR transformed data of i-th bug.
        #making differnt otu tables for each sample type
        otu_table <- physeq0_clr_all %>% otu_table
        
        #tax table for different sample type
        taxa_data <- otu_table[i] %>% t %>% data.frame()
        
        #Making a merged dataframe having sample data and CLR transformed output
        lme_data <- merge(taxa_data, sample_data(physeq0_clr_all), by = 0) %>% column_to_rownames("Row.names")
        
        #generating a character of formula.
        #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
        lme4_formula <- paste(names(taxa_data), "~", "sample_type + lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
        
        all_result <- lme4::lmer(formula = lme4_formula,
                                 data = lme_data) %>% 
                lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                summary %>%
                .$coefficients %>%
                data.frame(check.names = F) %>% 
                mutate(feature = names(taxa_data)) %>% 
                rownames_to_column("value")  %>%
                subset(., .$value %in% c("sample_typeBAL",
                                            "sample_typeNasal",
                                            "sample_typeSputum",
                                            "lypma", "benzonase", "host_zero", "molysis", "qiaamp")) %>%
                mutate(metadata = case_when(grepl("sample_type", `value`) ~ "sample_type",
                                            grepl("lypma", `value`) ~ "lypma",
                                            grepl("benzonase", `value`) ~ "benzonase",
                                            grepl("host_zero", `value`) ~ "host_zero",
                                            grepl("molysis", `value`) ~ "molysis",
                                            grepl("qiaamp", `value`) ~ "qiaamp",
                                            .default = "treatment"),
                       value = gsub("sample_type", "", value)
                )
                
        
        #row binding all the associations of i-th taxa to one data frame
        lmer_taxa_all <- rbind(lmer_taxa_all,
                               all_result) %>% 
                remove_rownames()
        
}


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



# Functional DA analysis --------------------------------------------------


# LME4 for CLR normalized data --------------------------------------------

#Running LME4 for each taxa


#with all samples

lmer_taxa_all_f <- data.frame()

all <- physeq0_clr_all_f %>% taxa_sums() %>% length()

for(i in 1:all) {
        #Creating a data frame that includes CLR transformed data of i-th bug.
        #making differnt otu tables for each sample type
        otu_table_f <- physeq0_clr_all_f %>% otu_table
        
        #tax table for different sample type
        taxa_data_f <- otu_table_f[i] %>% t %>% data.frame()
        
        #Making a merged dataframe having sample data and CLR transformed output
        lme_data_f <- merge(taxa_data_f, sample_data(physeq0_clr_all_f), by = 0) %>% column_to_rownames("Row.names")
        
        #generating a character of formula.
        #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
        lme4_formula_f <- paste(names(taxa_data_f), "~", "sample_type + lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
        
        #BAL stratified analysis
        all_result_f <- lme4::lmer(formula = lme4_formula_f,
                                 data = lme_data_f) %>% 
                lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                summary %>%
                .$coefficients %>%
                data.frame(check.names = F) %>% 
                mutate(feature = names(taxa_data_f)) %>% 
                rownames_to_column("value")  %>%
                subset(., .$value %in% c("sample_typeBAL",
                                         "sample_typeNasal",
                                         "sample_typeSputum",
                                         "lypma", "benzonase", "host_zero", "molysis", "qiaamp")) %>%
                mutate(metadata = case_when(grepl("sample_type", `value`) ~ "sample_type",
                                            grepl("lypma", `value`) ~ "lypma",
                                            grepl("benzonase", `value`) ~ "benzonase",
                                            grepl("host_zero", `value`) ~ "host_zero",
                                            grepl("molysis", `value`) ~ "molysis",
                                            grepl("qiaamp", `value`) ~ "qiaamp",
                                            .default = "treatment"),
                       value = gsub("sample_type", "", value)
                )
        
        
        #row binding all the associations of i-th taxa to one data frame
        lmer_taxa_all_f <- rbind(lmer_taxa_all_f,
                               all_result_f) %>% 
                remove_rownames()
        
}


#BAL


lmer_taxa_bal_f <- data.frame()

a <- physeq0_clr_bal_f %>% taxa_sums() %>% length()

for(i in 1:a) {
        #Creating a data frame that includes CLR transformed data of i-th bug.
        #making differnt otu tables for each sample type
        otu_table_bal_f <- physeq0_clr_bal_f %>% otu_table
        
        #tax table for different sample type
        taxa_data_bal_f <- otu_table_bal_f[i] %>% t %>% data.frame()
        
        #Making a merged dataframe having sample data and CLR transformed output
        lme_data_bal_f <- merge(taxa_data_bal_f, sample_data(physeq0_clr_bal_f), by = 0) %>% column_to_rownames("Row.names")
        
        #generating a character of formula.
        #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
        lme4_formula_f <- paste(names(taxa_data_bal_f), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
        
        #BAL stratified analysis
        BAL_result_f <- lme4::lmer(formula = lme4_formula_f,
                                 data = lme_data_bal_f) %>% 
                lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                summary %>%
                .$coefficients %>%
                data.frame(check.names = F) %>% 
                mutate(feature = names(taxa_data_bal_f),
                       sample_type = "BAL") %>% 
                rownames_to_column("metadata")  %>%
                subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        
        #row binding all the associations of i-th taxa to one data frame
        lmer_taxa_bal_f <- rbind(lmer_taxa_bal_f,
                               BAL_result_f) %>% 
                remove_rownames()
        
}

#Nasal

lmer_taxa_ns_f <- data.frame()

b <- physeq0_clr_ns_f %>% taxa_sums() %>% length()


for(i in 1:b) {
        #Creating a data frame that includes CLR transformed data of i-th bug.
        #making differnt otu tables for each sample type
        otu_table_ns_f <- physeq0_clr_ns_f %>% otu_table
        
        #tax table for different sample type
        taxa_data_ns_f <- otu_table_ns_f[i] %>% t %>% data.frame()
        
        #Making a merged dataframe having sample data and CLR transformed output
        lme_data_ns_f <- merge(taxa_data_ns_f, sample_data(physeq0_clr_ns_f), by = 0) %>% column_to_rownames("Row.names")
        
        #generating a character of formula.
        #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
        lme4_formula_f <- paste(names(taxa_data_ns_f), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
        
        #BAL stratified analysis
        Nasal_result_f <- lme4::lmer(formula = lme4_formula_f,
                                   data = lme_data_ns_f) %>% 
                lmerTest::as_lmerModLmerTest() %>% # p-value calculated by lmerTest
                summary %>%
                .$coefficients %>%
                data.frame(check.names = F) %>% 
                mutate(feature = names(taxa_data_ns_f),
                       sample_type = "Nasal") %>% 
                rownames_to_column("metadata")  %>%
                subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        
        #row binding all the associations of i-th taxa to one data frame
        lmer_taxa_ns_f <- rbind(lmer_taxa_ns_f,
                              Nasal_result_f) %>% 
                remove_rownames()
        
        
}

lmer_taxa_spt_f <- data.frame()

c <- physeq0_clr_spt_f %>% taxa_sums() %>% length()


for(i in 1:c) {
        #Creating a data frame that includes CLR transformed data of i-th bug.
        
        #making differnt otu tables for each sample type
        otu_table_spt_f <- physeq0_clr_spt_f %>% otu_table
        
        #tax table for different sample type
        taxa_data_spt_f <- otu_table_spt_f[i] %>% t %>% data.frame()
        
        #Making a merged dataframe having sample data and CLR transformed output
        lme_data_spt_f <- merge(taxa_data_spt_f, sample_data(physeq0_clr_spt_f), by = 0) %>% column_to_rownames("Row.names")
        
        #generating a character of formula.
        #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
        lme4_formula_f <- paste(names(taxa_data_spt_f), "~", "lypma + benzonase + host_zero + molysis + qiaamp + (1|subject_id)")
        
        #stratified analysis
        Sputum_result_f <- lme4::lmer(formula = lme4_formula_f,
                                    data = lme_data_spt_f) %>% 
                lmerTest::as_lmerModLmerTest() %>%
                summary %>%
                .$coefficients %>%
                data.frame(check.names = F) %>% 
                mutate(feature = names(taxa_data_spt_f),
                       sample_type = "Sputum") %>% 
                rownames_to_column("metadata")  %>%
                subset(., .$metadata %in% c("lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        
        #row binding all the associations of i-th taxa to one data frame
        lmer_taxa_spt_f <- rbind(lmer_taxa_spt_f,
                               Sputum_result_f) %>% 
                remove_rownames()
        
        
}

#Calculating q-value

lmer_taxa_all <- lmer_taxa_all %>%
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  

#lmer_taxa_all %>% subset(., .$qval < 0.1) %>% .$metadata %>% table

lmer_taxa_bal <- lmer_taxa_bal %>%
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  

lmer_taxa_ns <-  lmer_taxa_ns %>% 
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  

lmer_taxa_spt <-  lmer_taxa_spt %>% 
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  


#for function,


#Calculating q-value

lmer_taxa_all_f <- lmer_taxa_all_f %>%
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  

#lmer_taxa_all %>% subset(., .$qval < 0.1) %>% .$metadata %>% table

lmer_taxa_bal_f <- lmer_taxa_bal_f %>%
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  

lmer_taxa_ns_f <-  lmer_taxa_ns_f %>% 
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  

lmer_taxa_spt_f <-  lmer_taxa_spt_f %>% 
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH"))  








#Save LMER output


write.csv(lmer_taxa_all, "data/da_lmer_filt_maaslin_all.csv")
write.csv(lmer_taxa_bal, "data/da_lmer_filt_fit_data_bal.csv")
write.csv(lmer_taxa_spt, "data/da_lmer_filt_fit_data_spt.csv")
write.csv(lmer_taxa_ns, "data/da_lmer_filt_fit_data_ns.csv")


write.csv(lmer_taxa_all_f, "data/da_lmer_filt_f_maaslin_all.csv")
write.csv(lmer_taxa_bal_f, "data/da_lmer_filt_f_fit_data_bal.csv")
write.csv(lmer_taxa_spt_f, "data/da_lmer_filt_f_fit_data_spt.csv")
write.csv(lmer_taxa_ns_f, "data/da_lmer_filt_f_fit_data_ns.csv")




