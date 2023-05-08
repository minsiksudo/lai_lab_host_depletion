library(phyloseq)
library(Maaslin2)
library(tidyverse)


# Loading data ------------------------------------------------------------


phyloseq <- readRDS("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221129_MGK_host_tidy_tax.rds")

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
sample_data(phyloseq$phyloseq_path_rpkm) <- sample_data(alpha_diversity(phyloseq$phyloseq_path_rpkm))

# Filtering samples failed sequencing -------------------------------------

phyloseq$phyloseq_count <- subset_samples(phyloseq$phyloseq_count, S.obs != 0 & sample_type %in% c("BAL", "Nasal", "Sputum", "Mock"))
phyloseq$phyloseq_rel <- subset_samples(phyloseq$phyloseq_rel, S.obs != 0 & sample_type %in% c("BAL", "Nasal", "Sputum", "Mock"))


# Prevalence filtering  ---------------------------------------------------

taxa_qc <- data.frame("species" =  otu_table(subset_samples(phyloseq$phyloseq_count,
                                                            S.obs != 0 & sample_type %in% c("Mock", "BAL", "Nasal", "Sputum"))) %>%
                              t() %>% colnames(),
                      "prevalence" = ifelse(subset_samples(phyloseq$phyloseq_count, S.obs != 0 & sample_type %in% c("Mock", "BAL", "Nasal", "Sputum")) %>% otu_table() > 0, 1, 0) %>% t() %>% colSums(), #Prevalence of taxa
                      "mean_rel_abd" = subset_samples(phyloseq$phyloseq_count, S.obs != 0 & sample_type %in% c("Mock", "BAL", "Nasal", "Sputum")) %>%
                              otu_table() %>%
                              t() %>%
                              colMeans(na.rm = T) #mean relativ abundacne 
)

function_qc <- data.frame("function" =  otu_table(subset_samples(phyloseq$phyloseq_path_rpkm, S.obs != 0 & sample_type %in% c("Mock", "BAL", "Nasal", "Sputum"))) %>% t() %>% colnames(),
                          "prevalence" = ifelse(subset_samples(phyloseq$phyloseq_path_rpkm, S.obs != 0 & sample_type %in% c("Mock", "BAL", "Nasal", "Sputum")) %>% otu_table() > 0, 1, 0) %>% t() %>% colSums(), #Prevalence of taxa
                          "mean_rpkm" = subset_samples(phyloseq$phyloseq_path_rpkm, S.obs != 0 & sample_type %in% c("Mock", "BAL", "Nasal", "Sputum")) %>% otu_table() %>% t() %>% colMeans(na.rm = T) #mean relativ abundacne 
)


# Decontam ----------------------------------------------------------------





# Making a list of filters ------------------------------------------------


red_flag_taxa <- data.frame(species = taxa_qc$species,
                            red_flag_prev_abd = ifelse(taxa_qc$prevalence < otu_table(phyloseq_unfiltered$phyloseq_rel) %>%
                                                               t %>% rownames() %>%
                                                               length * 0.05 & taxa_qc$mean_rel_abd < quantile(taxa_qc$mean_rel_abd, 0.75), 1,0)) %>%
        mutate(red_flag_decontam_prev = species %in% (contaminant %>%
                                                              subset(., .$method == "prevalence" & .$sample_type != "all") %>%
                                                              .$contaminant %>% unique()))


red_flag_function <- data.frame(function. = function_qc$function., red_flag_prev_abd = ifelse(function_qc$prevalence < otu_table(phyloseq$phyloseq_path_rpkm) %>% t %>% rownames() %>% length * 0.05 & function_qc$mean_rpkm < quantile(function_qc$mean_rpkm, 0.75), 1, 0))



# Adding variables for MaAsLin --------------------------------------------

sample_data(phyloseq$phyloseq_rel)$log10.Final_reads <- log10(sample_data(phyloseq$phyloseq_rel)$Final_reads)
sample_data(phyloseq$phyloseq_rel)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_rel)$sample_type, sample_data(phyloseq$phyloseq_rel)$treatment, sep = ":")
sample_data(phyloseq$phyloseq_path_rpkm)$log10.Final_reads <- log10(sample_data(phyloseq$phyloseq_path_rpkm)$Final_reads)
sample_data(phyloseq$phyloseq_path_rpkm)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_path_rpkm)$sample_type, sample_data(phyloseq$phyloseq_path_rpkm)$treatment, sep = ":")

# Making new phyloseq object ----------------------------------------------

phyloseq$phyloseq_count_filtered <- prune_taxa(subset(red_flag_taxa,
                                                     red_flag_taxa$red_flag_prev_abd != 1 & !red_flag_taxa$red_flag_decontam_prev)$species,
                                              phyloseq$phyloseq_count)
phyloseq$phyloseq_path_rpkm_filtered <- prune_taxa(subset(red_flag_taxa, red_flag_function$red_flag_prev_abd != 1)$function., phyloseq$phyloseq_path_rpkm)
phyloseq$phyloseq_rel_filtered <- transform_sample_counts(phyloseq$phyloseq_count_filtered, function (x) {x/sum(x)})


#  MaAsLin --------------------------------------------------------

# Maaslin - # # y ~ log(final reads) + sample_type + treatment  -----------

#all samples
capture.output(maaslin_all2 = Maaslin2(input_data = otu_table(phyloseq$phyloseq_rel) %>% t %>% data.frame(), 
                input_metadata = phyloseq$phyloseq_rel %>% sample_data %>% data.frame(check.names = F), 
                output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                fixed_effects = c("sample_type", "log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                transform = "LOG", #default
                normalization = "TSS",
                random_effects = c("subject_id"), 
                reference = c("sample_type,BAL"), 
                plot_heatmap = F,
                plot_scatter = F))
maaslin_all <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#Mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel, sample_type == "Mock")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel, sample_type == "Mock")) %>% data.frame(), 
                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                       transform = "LOG", #default
                       normalization = "TSS",
                       plot_heatmap = F,
                       plot_scatter = F)
)
fit_data_pos <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


#NS
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel, sample_type == "Nasal")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel, sample_type == "Nasal")) %>% data.frame(), 
                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                       transform = "LOG", #default
                       normalization = "TSS",
                       random_effects = c("subject_id"), 
                       plot_heatmap = F,
                       plot_scatter = F)
)
fit_data_ns <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel, sample_type == "BAL")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel, sample_type == "BAL")) %>% data.frame(), 
                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                       transform = "LOG", #default
                       normalization = "TSS", # for RPKM, normalization is bit hard
                       random_effects = c("subject_id"), 
                       plot_heatmap = F,
                       plot_scatter = F)
)


fit_data_bal <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel, sample_type == "Sputum")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel, sample_type == "Sputum")) %>% data.frame(), 
                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                       transform = "LOG", #default
                       normalization = "TSS",
                       random_effects = c("subject_id"),
                       plot_heatmap = F,
                       plot_scatter = F)
)
fit_data_spt <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#NS + pos
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel, sample_type == "Nasal" | sample_type == "Mock")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel, sample_type == "Nasal" | sample_type == "Mock")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
fit_data_ns_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal + pos
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel, sample_type == "BAL" | sample_type == "Mock")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel, sample_type == "BAL" | sample_type == "Mock")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


fit_data_bal_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum + pos
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel, sample_type == "Sputum" | sample_type == "Mock")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel, sample_type == "Sputum" | sample_type == "Mock")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
fit_data_spt_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#Generating interaction term
sample_data(phyloseq$phyloseq_rel)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_rel)$sample_type, sample_data(phyloseq$phyloseq_rel)$treatment, sep = "*")


capture.output(maaslin_interaction = Maaslin2(input_data = otu_table(phyloseq$phyloseq_rel) %>% t %>% data.frame(), 
                                              input_metadata = phyloseq$phyloseq_rel %>% sample_data %>% data.frame(check.names = F), 
                                              output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                              fixed_effects = c("sample_type", "log10.Final_reads", "treatment", "sampletype_treatment"), 
                                              transform = "LOG", #default
                                              normalization = "TSS", 
                                              random_effects = c("subject_id"), 
                                              reference = c("sample_type,BAL", "treatment,Control", "sampletype_treatment,BAL*Control"), 
                                              plot_heatmap = F,
                                              plot_scatter = F))

maaslin_interaction <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


sample_data <- sample_data(phyloseq$phyloseq_path_rpkm) %>% data.frame(check.names = F) %>% subset(., !is.nan(.$simpson))
phyloseq$phyloseq_path_rpkm <- subset_samples(phyloseq$phyloseq_path_rpkm, S.obs != 0 & sample_type %in% c("BAL", "Nasal", "Sputum", "Mock", "Neg."))
sample_data(phyloseq$phyloseq_path_rpkm)$log10.Final_reads <- log10(sample_data(phyloseq$phyloseq_path_rpkm)$Final_reads)
sample_data(phyloseq$phyloseq_path_rpkm)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_path_rpkm)$sample_type, sample_data(phyloseq$phyloseq_path_rpkm)$treatment, sep = ":")




#all samples
capture.output(maaslin_all2 = Maaslin2(input_data = otu_table(phyloseq$phyloseq_path_rpkm) %>% t %>% data.frame(), 
                                       input_metadata = phyloseq$phyloseq_path_rpkm %>% sample_data %>% data.frame(check.names = F), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("sample_type", "log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       reference = c("sample_type,BAL"), 
                                       plot_heatmap = F,
                                       plot_scatter = F))
f_maaslin_all <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#Mock
# # y ~ log(final reads) + sample_type + treatment 


capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Mock")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Mock")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
f_fit_data_pos <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


#NS
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Nasal")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Nasal")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
f_fit_data_ns <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "BAL")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "BAL")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


f_fit_data_bal <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Sputum")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Sputum")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
f_fit_data_spt <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#NS + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Nasal" | sample_type == "Mock")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Nasal" | sample_type == "Mock")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
f_fit_data_ns_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "BAL" | sample_type == "Mock")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "BAL" | sample_type == "Mock")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


f_fit_data_bal_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Sputum" | sample_type == "Mock")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm, sample_type == "Sputum" | sample_type == "Mock")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
f_fit_data_spt_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#Generating interaction term
sample_data(phyloseq$phyloseq_path_rpkm)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_path_rpkm)$sample_type, sample_data(phyloseq$phyloseq_path_rpkm)$treatment, sep = "*")


capture.output(maaslin_interaction = Maaslin2(input_data = otu_table(phyloseq$phyloseq_path_rpkm) %>% t %>% data.frame(), 
                                              input_metadata = phyloseq$phyloseq_path_rpkm %>% sample_data %>% data.frame(check.names = F), 
                                              output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                              fixed_effects = c("sample_type", "log10.Final_reads", "treatment", "sampletype_treatment"), 
                                              transform = "LOG", #default
                                              normalization = "TSS", 
                                              random_effects = c("subject_id"), 
                                              reference = c("sample_type,BAL", "treatment,Control", "sampletype_treatment,BAL*Control"), 
                                              plot_heatmap = F,
                                              plot_scatter = F))

f_maaslin_interaction <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


write.csv(maaslin_all, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_all.csv")
write.csv(maaslin_interaction, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/maaslin_interaction.csv")
write.csv(fit_data_ns, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_ns.csv")
write.csv(fit_data_bal, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_bal.csv")
write.csv(fit_data_spt, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_spt.csv")
write.csv(fit_data_ns_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_ns_mock.csv")
write.csv(fit_data_bal_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_bal_mock.csv")
write.csv(fit_data_spt_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_spt_mock.csv")
write.csv(fit_data_pos, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_pos.csv")

write.csv(f_maaslin_all, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_maaslin_all.csv")
write.csv(f_maaslin_interaction, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_maaslin_interaction.csv")
write.csv(f_fit_data_bal, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_bal.csv")
write.csv(f_fit_data_spt, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_spt.csv")
write.csv(f_fit_data_ns, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_ns.csv")

write.csv(f_fit_data_bal_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_bal_mock.csv")
write.csv(f_fit_data_spt_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_spt_mock.csv")
write.csv(f_fit_data_ns_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_ns_mock.csv")

write.csv(f_fit_data_pos, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_pos.csv")




# Filtered MaAsLin --------------------------------------------------------



# Maaslin - # # y ~ log(final reads) + sample_type + treatment  -----------

#all samples
capture.output(maaslin_all2 = Maaslin2(input_data = otu_table(phyloseq$phyloseq_rel_filtered) %>% t %>% data.frame(), 
                                       input_metadata = phyloseq$phyloseq_rel_filtered %>% sample_data %>% data.frame(check.names = F), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("sample_type", "log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       reference = c("sample_type,BAL"), 
                                       plot_heatmap = F,
                                       plot_scatter = F))
filt_maaslin_all <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#Mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Mock")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Mock")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
filt_fit_data_pos <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


#NS
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Nasal")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Nasal")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
filt_fit_data_ns <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "BAL")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "BAL")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


filt_fit_data_bal <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Sputum")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Sputum")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
filt_fit_data_spt <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#NS + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Nasal" | sample_type == "Mock")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Nasal" | sample_type == "Mock")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
filt_fit_data_ns_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "BAL" | sample_type == "Mock")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "BAL" | sample_type == "Mock")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


filt_fit_data_bal_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Sputum" | sample_type == "Mock")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_rel_filtered, sample_type == "Sputum" | sample_type == "Mock")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
filt_fit_data_spt_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#Generating interaction term
sample_data(phyloseq$phyloseq_rel_filtered)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_rel_filtered)$sample_type, sample_data(phyloseq$phyloseq_rel_filtered)$treatment, sep = "*")


capture.output(maaslin_interaction = Maaslin2(input_data = otu_table(phyloseq$phyloseq_rel_filtered) %>% t %>% data.frame(), 
                                              input_metadata = phyloseq$phyloseq_rel_filtered %>% sample_data %>% data.frame(check.names = F), 
                                              output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                              fixed_effects = c("sample_type", "log10.Final_reads", "treatment", "sampletype_treatment"), 
                                              transform = "LOG", #default
                                              normalization = "TSS", 
                                              random_effects = c("subject_id"), 
                                              reference = c("sample_type,BAL", "treatment,Control", "sampletype_treatment,BAL*Control"), 
                                              plot_heatmap = F,
                                              plot_scatter = F))

filt_maaslin_interaction <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


#all samples
capture.output(maaslin_all2 = Maaslin2(input_data = otu_table(phyloseq$phyloseq_path_rpkm_filtered) %>% t %>% data.frame(), 
                                       input_metadata = phyloseq$phyloseq_path_rpkm_filtered %>% sample_data %>% data.frame(check.names = F), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("sample_type", "log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       reference = c("sample_type,BAL"), 
                                       plot_heatmap = F,
                                       plot_scatter = F))
filt_f_maaslin_all <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#Mock
# # y ~ log(final reads) + sample_type + treatment 


capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Mock")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Mock")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
filt_f_fit_data_pos <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#NS
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Nasal")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Nasal")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
filt_f_fit_data_ns <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "BAL")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "BAL")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


filt_f_fit_data_bal <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Sputum")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Sputum")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
filt_f_fit_data_spt <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#NS + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Nasal" | sample_type == "Mock")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Nasal" | sample_type == "Mock")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
filt_f_fit_data_ns_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "BAL" | sample_type == "Mock")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "BAL" | sample_type == "Mock")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


filt_f_fit_data_bal_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum + mock
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Sputum" | sample_type == "Mock")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq$phyloseq_path_rpkm_filtered, sample_type == "Sputum" | sample_type == "Mock")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
filt_f_fit_data_spt_mock <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#Generating interaction term
sample_data(phyloseq$phyloseq_path_rpkm_filtered)$sampletype_treatment <- paste(sample_data(phyloseq$phyloseq_path_rpkm_filtered)$sample_type, sample_data(phyloseq$phyloseq_path_rpkm_filtered)$treatment, sep = "*")


capture.output(maaslin_interaction = Maaslin2(input_data = otu_table(phyloseq$phyloseq_path_rpkm_filtered) %>% t %>% data.frame(), 
                                              input_metadata = phyloseq$phyloseq_path_rpkm_filtered %>% sample_data %>% data.frame(check.names = F), 
                                              output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                              fixed_effects = c("sample_type", "log10.Final_reads", "treatment", "sampletype_treatment"), 
                                              transform = "LOG", #default
                                              normalization = "TSS", 
                                              random_effects = c("subject_id"), 
                                              reference = c("sample_type,BAL", "treatment,Control", "sampletype_treatment,BAL*Control"), 
                                              plot_heatmap = F,
                                              plot_scatter = F))

filt_f_maaslin_interaction <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


write.csv(filt_maaslin_all, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_maaslin_all.csv")
write.csv(filt_maaslin_interaction, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_maaslin_interaction.csv")
write.csv(filt_fit_data_ns, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_fit_data_ns.csv")
write.csv(filt_fit_data_bal, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_fit_data_bal.csv")
write.csv(filt_fit_data_spt, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_fit_data_spt.csv")
write.csv(filt_fit_data_ns_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_fit_data_ns_mock.csv")
write.csv(filt_fit_data_bal_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_fit_data_bal_mock.csv")
write.csv(filt_fit_data_spt_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_fit_data_spt_mock.csv")
write.csv(filt_fit_data_pos, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_fit_data_pos.csv")

write.csv(filt_f_maaslin_all, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_maaslin_all.csv")
write.csv(filt_f_maaslin_interaction, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_maaslin_interaction.csv")
write.csv(filt_f_fit_data_bal, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_fit_data_bal.csv")
write.csv(filt_f_fit_data_spt, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_fit_data_spt.csv")
write.csv(filt_f_fit_data_ns, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_fit_data_ns.csv")

write.csv(filt_f_fit_data_bal_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_fit_data_bal_mock.csv")
write.csv(filt_f_fit_data_spt_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_fit_data_spt_mock.csv")
write.csv(filt_f_fit_data_ns_mock, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_fit_data_ns_mock.csv")

write.csv(filt_f_fit_data_pos, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/filt_f_fit_data_pos.csv")

