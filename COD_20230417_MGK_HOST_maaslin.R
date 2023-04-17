library(phyloseq)
library(Maaslin2)
library(tidyverse)

phyloseq <- readRDS("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221129_MGK_host_tidy_tax.rds")

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

phyloseq_rel_nz <- subset_samples(phyloseq$phyloseq_rel, S.obs != 0 & sample_type %in% c("BAL", "Nasal", "Sputum", "Mock", "Neg."))
sample_data(phyloseq_rel_nz)$log10.Final_reads <- log10(sample_data(phyloseq_rel_nz)$Final_reads)
sample_data(phyloseq_rel_nz)$sampletype_treatment <- paste(sample_data(phyloseq_rel_nz)$sample_type, sample_data(phyloseq_rel_nz)$treatment, sep = ":")



# Maaslin - # # y ~ log(final reads) + sample_type + treatment  -----------

#all samples
capture.output(maaslin_all2 = Maaslin2(input_data = otu_table(phyloseq_rel_nz) %>% t %>% data.frame(), 
                input_metadata = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), 
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

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Mock")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Mock")) %>% data.frame(), 
                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                       transform = "LOG", #default
                       normalization = "TSS",
                       plot_heatmap = F,
                       plot_scatter = F)
)
fit_data_pos <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


#Neg
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Neg.")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Neg.")) %>% data.frame(), 
                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                       transform = "LOG", #default
                       normalization = "TSS",
                       plot_heatmap = F,
                       plot_scatter = F)
)
fit_data_neg <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#NS
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Nasal")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Nasal")) %>% data.frame(), 
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

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "BAL")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "BAL")) %>% data.frame(), 
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

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Sputum")) %>% t %>% data.frame(),
                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Sputum")) %>% data.frame(), 
                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                       transform = "LOG", #default
                       normalization = "TSS",
                       random_effects = c("subject_id"),
                       plot_heatmap = F,
                       plot_scatter = F)
)
fit_data_spt <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#NS + neg
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Nasal" | sample_type == "Neg.")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Nasal" | sample_type == "Neg.")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
fit_data_ns_neg <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal + neg
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "BAL" | sample_type == "Neg.")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "BAL" | sample_type == "Neg.")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


fit_data_bal_neg <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum + neg
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Sputum" | sample_type == "Neg.")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Sputum" | sample_type == "Neg.")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "TSS",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
fit_data_spt_neg <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#Generating interaction term
sample_data(phyloseq_rel_nz)$sampletype_treatment <- paste(sample_data(phyloseq_rel_nz)$sample_type, sample_data(phyloseq_rel_nz)$treatment, sep = "*")


capture.output(maaslin_interaction = Maaslin2(input_data = otu_table(phyloseq_rel_nz) %>% t %>% data.frame(), 
                                              input_metadata = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), 
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
phyloseq_rel_nz <- subset_samples(phyloseq$phyloseq_path_rpkm, S.obs != 0 & sample_type %in% c("BAL", "Nasal", "Sputum", "Mock", "Neg."))
sample_data(phyloseq_rel_nz)$log10.Final_reads <- log10(sample_data(phyloseq_rel_nz)$Final_reads)
sample_data(phyloseq_rel_nz)$sampletype_treatment <- paste(sample_data(phyloseq_rel_nz)$sample_type, sample_data(phyloseq_rel_nz)$treatment, sep = ":")




#all samples
capture.output(maaslin_all2 = Maaslin2(input_data = otu_table(phyloseq_rel_nz) %>% t %>% data.frame(), 
                                       input_metadata = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("sample_type", "log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "NONE",
                                       random_effects = c("subject_id"), 
                                       reference = c("sample_type,BAL"), 
                                       plot_heatmap = F,
                                       plot_scatter = F))
f_maaslin_all <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#Mock
# # y ~ log(final reads) + sample_type + treatment 


capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Mock")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Mock")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "NONE",
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
f_fit_data_pos <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")


#Neg
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Neg.")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Neg.")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "TSS",
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
f_fit_data_neg <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#NS
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Nasal")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Nasal")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "NONE",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
f_fit_data_ns <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "BAL")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "BAL")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "NONE", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


f_fit_data_bal <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Sputum")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Sputum")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "NONE",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
f_fit_data_spt <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#NS + neg
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_ns2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Nasal" | sample_type == "Neg.")) %>% t %>% data.frame(),
                                       input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Nasal" | sample_type == "Neg.")) %>% data.frame(), 
                                       output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                       fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                       transform = "LOG", #default
                                       normalization = "NONE",
                                       random_effects = c("subject_id"), 
                                       plot_heatmap = F,
                                       plot_scatter = F)
)
f_fit_data_ns_neg <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#bal + neg
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_bal2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "BAL" | sample_type == "Neg.")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "BAL" | sample_type == "Neg.")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "NONE", # for RPKM, normalization is bit hard
                                        random_effects = c("subject_id"), 
                                        plot_heatmap = F,
                                        plot_scatter = F)
)


f_fit_data_bal_neg <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")

#sputum + neg
# # y ~ log(final reads) + sample_type + treatment 

capture.output(fit_data_spt2 = Maaslin2(input_data = otu_table(subset_samples(phyloseq_rel_nz, sample_type == "Sputum" | sample_type == "Neg.")) %>% t %>% data.frame(),
                                        input_metadata = sample_data(subset_samples(phyloseq_rel_nz, sample_type == "Sputum" | sample_type == "Neg.")) %>% data.frame(), 
                                        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
                                        transform = "LOG", #default
                                        normalization = "NONE",
                                        random_effects = c("subject_id"),
                                        plot_heatmap = F,
                                        plot_scatter = F)
)
f_fit_data_spt_neg <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output/all_results.tsv", sep = "\t")



#Generating interaction term
sample_data(phyloseq_rel_nz)$sampletype_treatment <- paste(sample_data(phyloseq_rel_nz)$sample_type, sample_data(phyloseq_rel_nz)$treatment, sep = "*")


capture.output(maaslin_interaction = Maaslin2(input_data = otu_table(phyloseq_rel_nz) %>% t %>% data.frame(), 
                                              input_metadata = phyloseq_rel_nz %>% sample_data %>% data.frame(check.names = F), 
                                              output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
                                              fixed_effects = c("sample_type", "log10.Final_reads", "treatment", "sampletype_treatment"), 
                                              transform = "LOG", #default
                                              normalization = "NONE", 
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
write.csv(fit_data_ns_neg, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_ns_neg.csv")
write.csv(fit_data_bal_neg, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_bal_neg.csv")
write.csv(fit_data_spt_neg, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_spt_neg.csv")
write.csv(fit_data_neg, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_neg.csv")
write.csv(fit_data_pos, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/fit_data_pos.csv")

write.csv(f_maaslin_all, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_maaslin_all.csv")
write.csv(f_maaslin_interaction, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_maaslin_interaction.csv")
write.csv(f_fit_data_bal, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_bal.csv")
write.csv(f_fit_data_spt, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_spt.csv")
write.csv(f_fit_data_ns, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_ns.csv")

write.csv(f_fit_data_bal_neg, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_bal_neg.csv")
write.csv(f_fit_data_spt_neg, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_spt_neg.csv")
write.csv(f_fit_data_ns_neg, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_ns_neg.csv")

write.csv(f_fit_data_neg, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_neg.csv")
write.csv(f_fit_data_pos, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/5_Scripts/MGK/Host_depletion_git/data/f_fit_data_pos.csv")