rm(list = ls())
path_working <- "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/"
path_library <- "/Library/Frameworks/R.framework/Resources/library"
str_libraries <- c("readxl", "ggplot2", "phyloseq", "vegan", "tidyverse", "microbiome")
seed <- 20221230
                               


# Initial setting for R ---------------------------------------------------
boot_r <- function() {
        .libPaths(path_library)
        
        require(pacman)
        knitr::opts_knit$set(root.dir = path_working)
        
        str_libraries |> unique() |> sort() -> str_libraries
        pacman::p_load(char = str_libraries)
        
        set.seed(seed)
}
boot_r()
 

# #Loading data -----------------------------------------------------------


        #Sample_data
                                
                                #Loading data
                                #Loading the data
                                sample_data1 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220324_MGK_Host_extraction.xlsx", sheet = 2)
                                sample_data1$depletion_date <- ifelse(sample_data1$host_zero == 1 | sample_data1$lypma == 1 | sample_data1$qiaamp == 1, 
                                                                      "20220317", 
                                                                      ifelse(sample_data1$molysis == 1, "20220318", 
                                                                             ifelse(sample_data1$benzonase == 1,  "20220321", NA)))
                                
                                #Cryopreserved data
                                sample_data2 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220404_MGK_Host_qPCR_cryopreserve.xlsx", sheet = 2)
                                sample_data2$depletion_date <- ifelse(sample_data2$host_zero == 1 | sample_data2$qiaamp == 1 | sample_data2$molysis == 1, 
                                                                      "20220330", 
                                                                      ifelse(sample_data2$lypma == 1 | sample_data2$benzonase == 1, "20220331", NA))
                                #BAL and sputum data
                                sample_data3 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220407_MGK_sputum_BAL_host_depletion.xlsx", sheet = 5)
                                sample_data3$extraction_date <- "20220412" #Modification of date
                                sample_data3$depletion_date <- ifelse(sample_data3$host_zero == 1 | sample_data3$qiaamp == 1 | sample_data3$molysis == 1, 
                                                                      "20220410", 
                                                                      ifelse(sample_data3$lypma == 1 | sample_data3$benzonase == 1, "20220411", NA))
                                
                                #Additional BAl and NS data
                                sample_data4 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220606_MGK_BAL_NB_host_depletion.xlsx", sheet = 2)
                                sample_data4$depletion_date  <- ifelse(sample_data4$control == 1, NA,  
                                                                       "20220609")
                                sample_data4$extraction_date <- "20220610"
                                #Laoding the qPCR data
                                #for sample_data1 and 2
                                #Loading data  - host
                                host1 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220324_MGK_Host_extraction.xlsx", sheet = 7)
                                names(host1)[3] <- "extraction_id"
                                host_simple1 <- host1[, c("extraction_id", "Quantity Mean")]
                                names(host_simple1)[2] <- "DNA_host"
                                host2 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220404_MGK_Host_qPCR_cryopreserve.xlsx", sheet = 6)
                                names(host2)[3] <- "extraction_id"
                                host_simple2 <- host2[, c("extraction_id", "Quantity Mean")]
                                names(host_simple2)[2] <- "DNA_host"
                                #Loading data  - Bac
                                bac1 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220324_MGK_Host_extraction.xlsx", sheet = 6)
                                names(bac1)[3] <- "extraction_id"
                                bac_simple1 <- bac1[, c("extraction_id", "Quantity Mean")]
                                names(bac_simple1)[2] <- "DNA_bac"
                                bac2 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220404_MGK_Host_qPCR_cryopreserve.xlsx", sheet = 5)
                                names(bac2)[3] <- "extraction_id"
                                bac_simple2 <- bac2[, c("extraction_id", "Quantity Mean")]
                                names(bac_simple2)[2] <- "DNA_bac"
                                
                                #Loading all sample data
                                sample_result1 <- merge(sample_data1, host_simple1, by = "extraction_id")
                                sample_result1 <- merge(sample_result1, bac_simple1, by = "extraction_id")
                                sample_result1 <- subset(sample_result1, !duplicated(sample_result1$extraction_id))
                                sample_result2 <- merge(sample_data2, host_simple2, by = "extraction_id")
                                sample_result2 <- merge(sample_result2, bac_simple2, by = "extraction_id")
                                sample_result2 <- subset(sample_result2, !duplicated(sample_result2$extraction_id))
                                
                                #for sampledata3     
                                
                                #Loading data  - host
                                host1 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220407_MGK_sputum_BAL_host_depletion.xlsx", sheet = 7)
                                host1$`Quantity Mean` <- ifelse(is.na(host1$`Quantity Mean`), host1$Quantity, host1$`Quantity Mean`)
                                names(host1)[3] <- "extraction_id"
                                host_simple1 <- host1[, c("extraction_id", "Quantity Mean")]
                                names(host_simple1)[2] <- "DNA_host"
                                #Loading data  - Bac
                                bac1 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220407_MGK_sputum_BAL_host_depletion.xlsx", sheet = 6)
                                names(bac1)[3] <- "extraction_id"
                                bac_simple1 <- bac1[, c("extraction_id", "Quantity Mean")]
                                names(bac_simple1)[2] <- "DNA_bac"
                                
                                #Loading all sample data
                                host_simple1 <- subset(host_simple1, !duplicated(host_simple1$extraction_id))
                                bac_simple1 <- subset(bac_simple1, !duplicated(bac_simple1$extraction_id))
                                
                                sample_result3 <- merge(sample_data3, host_simple1, by = "extraction_id")
                                sample_result3 <- merge(sample_result3, bac_simple1, by = "extraction_id")
                                
                                #for sampledata4     
                                
                                #Loading data  - host
                                host1 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220606_MGK_BAL_NB_host_depletion.xlsx", sheet = 6)
                                host1$`Quantity Mean` <- ifelse(is.na(host1$`Quantity Mean`), host1$Quantity, host1$`Quantity Mean`)
                                names(host1)[3] <- "extraction_id"
                                host_simple1 <- host1[, c("extraction_id", "Quantity Mean")]
                                names(host_simple1)[2] <- "DNA_host"
                                #Loading data  - Bac
                                bac1 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220606_MGK_BAL_NB_host_depletion.xlsx", sheet = 5)
                                names(bac1)[3] <- "extraction_id"
                                bac_simple1 <- bac1[, c("extraction_id", "Quantity Mean")]
                                names(bac_simple1)[2] <- "DNA_bac"
                                
                                #Loading all sample data
                                host_simple1 <- subset(host_simple1, !duplicated(host_simple1$extraction_id))
                                bac_simple1 <- subset(bac_simple1, !duplicated(bac_simple1$extraction_id))
                                
                                sample_result4 <- merge(sample_data4, host_simple1, by = "extraction_id")
                                sample_result4 <- merge(sample_result4, bac_simple1, by = "extraction_id")
                                
                                #Merging all the output                
                                
                                #merging data into one
                                sample_result1$pellet <- 0
                                sample_result2$pellet <- 0
                                sample_result3$pellet <- 0
                                
                                #Merging process
                                sample_result1$collection_date <- as.numeric(sample_result1$collection_date)
                                sample_result2$collection_date <- as.numeric(sample_result2$collection_date)
                                sample_result3$collection_date <- as.numeric(sample_result3$collection_date)
                                sample_result4$collection_date <- as.numeric(sample_result4$collection_date)
                                sample_result4
                                sample_result <- rbind(sample_result1, sample_result2, sample_result3, sample_result4)
                                
                                
                                #qPCR standards were 10, 0.1, 0.01.. etc. 
                                #Specification at Zymo was 20, 2, 0.2, ... per well 
                                sample_result$DNA_host <- sample_result$DNA_host/2
                                sample_result$DNA_bac <- sample_result$DNA_bac/2
                                
                                #qPCR was conducted for 10-fold diluted samples.
                                sample_result$DNA_host_nondil <- sample_result$DNA_host *10
                                sample_result$DNA_bac_nondil <- sample_result$DNA_bac *10
                                
                                sample_result$proportion <- sample_result$DNA_host/(sample_result$DNA_host + sample_result$DNA_bac)
                                
                                sample_result
                                
                                
                                #Methold labels
                                sample_result$treatment <- ifelse(sample_result$control==1, "control",
                                                                  ifelse(sample_result$lypma==1, "lyPMA",
                                                                         ifelse(sample_result$benzonase==1, "benzonase",
                                                                                ifelse(sample_result$qiaamp==1, "qiaamp",
                                                                                       ifelse(sample_result$host_zero==1, "host_zero",
                                                                                              ifelse(sample_result$molysis==1, "molysis", "mock"))))))
                                
                                
                                #Re-leveling factors before splitting data
                                sample_result$control <- factor(sample_result$control, levels = c(1, 0))
                                factor(sample_result$sample_type)
                                
                                
                                #subsetting control groups
                                sample_result_all <- sample_result
                                sample_result_all <- subset(sample_result_all, sample_result_all$sample_id != "N")
                                sample_result_all <- subset(sample_result_all, sample_result_all$sample_id != "P")
                                sample_result_all <- subset(sample_result_all, sample_result_all$sample_id != "Pos")
                                sample_result_all <- subset(sample_result_all, sample_result_all$sample_id != "Neg")
                                sample_result_all <- subset(sample_result_all, !is.na(sample_result_all$sample_id))
                                
                                sample_result <- subset(sample_result_all, sample_result_all$pellet == 0)
                                sample_result$sample_type <- factor(sample_result$sample_type, levels = c("contorl", "saliva", "nasal_blow", "nasal_swab", "Sputum", "BAL"))
                                
                                

# #Loading phyloseq object ------------------------------------------------

                
                
                                
                                path_rel <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20221102_PQ00331_Host_Depletion/humann3/PathAbundance.relab.metagenome.txt", sep = "\t")
                                path_rel$X..Pathway
                                #It seems this uses BioCyc pathway ids.....
                                #metacyc : refer to Jake
                                
                                #       Inquire Jake
                                path_rpk <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20221102_PQ00331_Host_Depletion/humann3/PathAbundance.relab.metagenome.txt", sep = "\t")
                                
                                kegg_count <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20221102_PQ00331_Host_Depletion/humann3/KeggPathwaysAbundance.metagenome.txt", sep = "\t", check.names = F)
                                kegg_count$`# Pathway`
                                #This seems to be using Kegg pathway ids
                                        #link at
                                        #kegg_url <- "https://www.genome.jp/kegg/pathway.html"

                                
                                gene_rpk <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20221102_PQ00331_Host_Depletion/humann3/GeneFamilies.KO.RPK.metagenome.txt", sep = "\t")
                                gene_rpk$X..Gene.Family
                                #Kegg ontholog
                                # --> irrelevant to the publication.. (maybe)
                
                
                #Adding Kegg pathway details 
                                
                        #Loaidng name of Keggs
                                kegg_path_url <- "https://www.genome.jp/kegg/pathway.html"
                                kegg_path_data <- readLines(con = kegg_path_url) %>% data.frame()
                                
                                kegg_path_data <- subset(kegg_path_data, grepl("<h4 id|<b id|<b>|/pathway/", kegg_path_data$.)) %>% data.frame()
                                
                                kegg_path_data <- gsub("\"small\"|</a></dd>|/a> <span class=\"new\">New!</span></dd>", "", kegg_path_data$.) %>% data.frame()
                                kegg_path_data$level1 <- ifelse(grepl("h4", kegg_path_data$.), kegg_path_data$., NA)
                                kegg_path_data$level2 <- ifelse(grepl("<b id|<b>", kegg_path_data$.), kegg_path_data$., NA)
                                #adding levels to the data
                                                for (row in c(2:length(kegg_path_data$level1))) {
                                                        kegg_path_data$level1[row] <- ifelse(is.na(kegg_path_data$level1[row]),
                                                                                             kegg_path_data$level1[row-1],
                                                                                             kegg_path_data$level1[row]
                                                        )
                                                }
                                                for (row in c(2:length(kegg_path_data$level2))) {
                                                        kegg_path_data$level2[row] <- ifelse(is.na(kegg_path_data$level2[row]),
                                                                                        kegg_path_data$level2[row-1],
                                                                                        kegg_path_data$level2[row]
                                                        )
                                                }
                                
                                
                                                kegg_path_data$level1 <- str_split_fixed(kegg_path_data$level1, " ", 3) %>% data.frame() %>% .[,3] %>% gsub("</h4>", "", .)
                                                kegg_path_data$level2 <- kegg_path_data$level2 %>% gsub("b id", "", .) %>% str_split_fixed(., " ", 2) %>% data.frame() %>% .[,2] %>% gsub("</b>", "", .)
                                                
                                                kegg_path_data$id <- str_split_fixed(kegg_path_data$., "\">", 2) %>% data.frame() %>% .[,1]
                                                kegg_path_data$pathway <- str_split_fixed(kegg_path_data$., "\">", 2) %>% data.frame() %>% .[,2]
                                                kegg_path_data$id <- substr(kegg_path_data$id, 9, 13) %>% paste("ko", ., sep = "")
                                                
                                                kegg_path_data$pathway <- str_split_fixed(kegg_path_data$pathway, "<", 2) %>% data.frame() %>% .[,1]
                                                kegg_path_data <- subset(kegg_path_data, !(kegg_path_data$id %>% substr(3,7) %>% as.numeric %>% is.na()))
                                                row.names(kegg_path_data) <- kegg_path_data$id
                                                kegg_path_data <- kegg_path_data[,c("level1", "level2", "pathway")]
                                                
                                #Kegg couint - removing sample names
                                row.names(kegg_count) <- kegg_count$`# Pathway`
                                kegg_count <- kegg_count[,-1]
                                kegg_path_data %>% row.names %>% length()
                                
                                
                                #`keggorthology` package have data(KeggOrthoDF)
                                #library(keggorthology)
                                #data(KeggOrtho)
                                # which have slightly fewer data than the webpage of kegg
                                
                                
                                
                        #QC of Keggs,
                                #unmapped
                                kegg_count %>% t() %>% data.frame()%>% .$UNMAPPED %>% log10 %>% hist
                                kegg_count %>% t() %>% data.frame()%>% .$UNINTEGRATED %>% log10 %>% hist
                                
                                
                        #Adding sample_id
                                sample_result$sample_id <- ifelse(sample_result$sample_type == "BAL",
                                                                  paste(sample_result$sample_id, sample_result$aliquot, sep = "_"),
                                                                  sample_result$sample_id)
                                
                                sample_result$sample_id <- ifelse(sample_result$cryo_preserve == 1,
                                                                  paste(sample_result$sample_id, "cryo", sep = "_"),
                                                                  sample_result$sample_id)
                                sample_result$sample_id <- ifelse(grepl( "NS_26|NS_37", sample_result$sample_id), 
                                                                  paste(sample_result$sample_id, sample_result$aliquot, sep = "_"),
                                                                  sample_result$sample_id)
                                sample_result$SampleID <- paste(sample_result$sample_id, sample_result$treatment, sep = "_")
                                
                                dummy <- kegg_count %>% t %>% .[,1:2] %>% data.frame()
                                dummy$SampleID <- row.names(dummy)
                                names(dummy) <- c("kegg_unmapped", "kegg_unintegrated", "SampleID")
                                dummy <- merge(dummy, subset(sample_result, sample_result$cryo_preserve == 0), by = "SampleID", all = T)
                                sample_data <- subset(dummy, dummy$SampleID %in% (t(kegg_count) %>% row.names()))
                                sample_data$treated <- rowSums(sample_data[, c("lypma", "benzonase", "molysis", "host_zero", "qiaamp")])
                                row.names(sample_data) <- sample_data$SampleID
                        
                        #constructing phylose object
                                #tax_table
                                        tax_table <- tax_table(kegg_path_data) 
                                        taxa_names(tax_table) <- row.names(kegg_path_data)
                                        kegg_path_data
                                #row.names(sample_data) <- sample_data$SampleID
                                phyloseq <- merge_phyloseq(otu_table(kegg_count, taxa_are_rows = T), tax_table, sample_data(sample_data))
                                colnames(tax_table(phyloseq)) <- c("level1", "level2", "level3")
                                
                                sample_data(phyloseq)$total_reads <- phyloseq %>% otu_table %>% colSums()
                                
                                
                                sample_data <- sample_data(phyloseq)
                                #sample_data$treatment <- factor(sample_data$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))
                                #sample_data$sample_type <- factor(sample_data$sample_type, levels = c("nasal_swab", "sputum", "bal", "pos_control", "neg_control"))
                                
                                #Q0. How many samples failed in sequencing
                                sample_data$total_reads %>% str
                                par(mfrow = c(3,1))
                                hist((log10(sample_data$total_reads + 1)), xlim = c(0,8), main = "log(Sum of otu table)", xlab = "") # one somaple showed 0 reads
                                kegg_count %>% t() %>% data.frame()%>% .$UNMAPPED %>% log10 %>% hist(xlim = c(0,8), main = "log(Unmapped)")
                                kegg_count %>% t() %>% data.frame()%>% .$UNINTEGRATED %>% log10 %>% hist(xlim = c(0,8), main = "log(Unintegrated)")
                                
                #No samples had 0 reads
                                
                                #how were the samples failed in library prep?
                                
                                #BAL
                                lib_fail_list <- c("NS_6_b_lyPMA", "NS_8_b_lyPMA", "NS_17_b_lyPMA", "NS_21_b_host_zero", "NS_22_b_host_zero", "NS_19_d_molysis", "NS_21_d_molysis", "NS_22_d_molysis", "NS_23_d_molysis", "NS_26_b_lyPMA", "BAL_073_b_host_zero", "BAL_073_c_molysis", "BAL_078_d_lyPMA")
                                sample_data(phyloseq)$lib_failed <- row.names(sample_data(phyloseq)) %in% lib_fail_list
                                
                                #I have no idea how Baylor decisded those sample as "library prep failed samples
                                
                                # Q0. How was the positive control and negative control? were there any contaminants? 
                                sample_data
                                phyloseq %>% sample_names()
                                phyloseq_control <- subset_samples(phyloseq, grepl("2022", SampleID))
                                phyloseq_sample <- subset_samples(phyloseq, !grepl("2022", SampleID))
                                
                                sample_data(phyloseq_control)$is.neg <- 1
                                sample_data(phyloseq_sample)$is.neg <- 0
                                
                                
                                phyloseq_control <- prune_taxa(taxa_sums(phyloseq_control) != 0, phyloseq_control) 
                                
                                phyloseq_control_rel <- transform_sample_counts(phyloseq_control, function(x){x / sum(x)})
                                
                                plot_bar(phyloseq_control_rel, fill="level1") + 
                                        facet_wrap (~ SampleID, scales= "free_x", nrow=1) +
                                        ylab("Relative abudnacne") +
                                        title("Species richenss") +
                                        theme_classic(base_size = 11, base_family = "serif")
                                plot_bar(phyloseq_control_rel, fill="level2") + 
                                        facet_wrap (~ SampleID, scales= "free_x", nrow=1) +
                                        ylab("Relative abudnacne") +
                                        title("Species richenss") +
                                        theme_classic(base_size = 11, base_family = "serif")
                                plot_bar(phyloseq_control_rel, fill="level3") + 
                                        facet_wrap (~ SampleID, scales= "free_x", nrow=1) +
                                        ylab("Relative abudnacne") +
                                        title("Species richenss") +
                                        theme_classic(base_size = 11, base_family = "serif")
                                
                                #issue with interpretation....it seems there is no change
                                
                #Decontam : not compatitive with functional analysis
                                # - omitted
                                
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
                                        alpha_diversity$SampleID <- row.names(alpha_diversity)
                                        sample_data <- merge(data.frame(sample_data), alpha_diversity, by = "SampleID")
                                        sample_data$sample_type <- as.character(sample_data$sample_type)
                                        row.names(sample_data) <- sample_data$SampleID
                                        sample_data
                                }
                                
                                sample_data(phyloseq) <- alpha_diversity(phyloseq) %>% sample_data
        
                                sample_data(phyloseq)$treatment <- tolower(sample_data(phyloseq)$treatment)
                                sample_data(phyloseq)$treatment <- factor(sample_data(phyloseq)$treatment, levels = c("control", "lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
                                sample_data(phyloseq)$sample_type <- factor(sample_data(phyloseq)$sample_type, levels = c("BAL", "Sputum", "nasal_swab"))

# Saving object -----------------------------------------------------------


                                write_rds(phyloseq, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221229_MGK_host_tidy_path.rds")

                                
                        
                                
