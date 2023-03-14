rm(list = ls())

#Required pacages

library(readxl)
library(phyloseq)
library(tidyverse)

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
                        
                        
                        
                        #control_data
                        sample_data5 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20230112_MGK_host_depletion_controls.xlsx", sheet = 2)
                        
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
                                                host3 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20230112_MGK_host_depletion_controls.xlsx", sheet = 8)
                                                names(host3)[2] <- "extraction_id"
                                                host_simple3 <- host3[, c("extraction_id", "Quantity Mean")]
                                                names(host_simple3)[2] <- "DNA_host"
                                                
                                                #Loading data  - Bac
                                                bac1 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220324_MGK_Host_extraction.xlsx", sheet = 6)
                                                names(bac1)[3] <- "extraction_id"
                                                bac_simple1 <- bac1[, c("extraction_id", "Quantity Mean")]
                                                names(bac_simple1)[2] <- "DNA_bac"
                                                bac2 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20220404_MGK_Host_qPCR_cryopreserve.xlsx", sheet = 5)
                                                names(bac2)[3] <- "extraction_id"
                                                bac_simple2 <- bac2[, c("extraction_id", "Quantity Mean")]
                                                names(bac_simple2)[2] <- "DNA_bac"
                                                bac3 <- read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/2_Protocols/host_depletion/SOP_20230112_MGK_host_depletion_controls.xlsx", sheet = 9)
                                                names(bac3)[2] <- "extraction_id"
                                                bac_simple3 <- bac3[, c("extraction_id", "Quantity Mean")]
                                                names(bac_simple3)[2] <- "DNA_bac"
                                                
                                                
                                                #Loading all sample data
                                                sample_result1 <- merge(sample_data1, host_simple1, by = "extraction_id")
                                                sample_result1 <- merge(sample_result1, bac_simple1, by = "extraction_id")
                                                sample_result1 <- subset(sample_result1, !duplicated(sample_result1$extraction_id))
                                                sample_result2 <- merge(sample_data2, host_simple2, by = "extraction_id")
                                                sample_result2 <- merge(sample_result2, bac_simple2, by = "extraction_id")
                                                sample_result2 <- subset(sample_result2, !duplicated(sample_result2$extraction_id))
                                                
                                                sample_result_control <- merge(sample_data5, host_simple3, by = "extraction_id")
                                                sample_result_control <- merge(sample_result_control, bac_simple3, by = "extraction_id")
                                                sample_result_control <- subset(sample_result_control, !duplicated(sample_result_control$extraction_id))
                                                
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
                                                sample_result_control$pellet <- 0
                                #Merging process
                                                sample_result1$collection_date <- as.numeric(sample_result1$collection_date)
                                                sample_result2$collection_date <- as.numeric(sample_result2$collection_date)
                                                sample_result3$collection_date <- as.numeric(sample_result3$collection_date)
                                                sample_result4$collection_date <- as.numeric(sample_result4$collection_date)
                                                sample_result_control$collection_date <- as.numeric(sample_result_control$collection_date)
                                                
                                sample_result <- rbind(sample_result1, sample_result2, sample_result3, sample_result4)
                                sample_result <- plyr::rbind.fill(sample_result, sample_result_control)
                                
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
                        sample_result$treatment <- ifelse(sample_result$treatment == "mock" | is.na(sample_result$treatment),
                                                    "control", sample_result$treatment)
                        sample_result$treatment
                        
                        #Re-leveling factors before splitting data
                                sample_result$control <- factor(sample_result$control, levels = c(1, 0))
                                factor(sample_result$sample_type)
                                
                                
                        #subsetting control groups
                                sample_result$sample_type <- ifelse(sample_result$extraction_id == "N",
                                                                    "neg_control_extraction",
                                                                    ifelse(sample_result$extraction_id == "P",
                                                                           "pos_control_extraction",
                                                                           sample_result$sample_type)
                                )
                                
                                controls <- subset(sample_result, grepl("control", sample_result$sample_type))
                                sample_result_all <- sample_result
                                sample_result_all <- subset(sample_result_all, sample_result_all$sample_id != "N")
                                sample_result_all <- subset(sample_result_all, sample_result_all$sample_id != "P")
                                sample_result_all <- subset(sample_result_all, sample_result_all$sample_id != "Pos")
                                sample_result_all <- subset(sample_result_all, sample_result_all$sample_id != "Neg")
                                sample_result_all <- subset(sample_result_all, !is.na(sample_result_all$sample_id))
                                
                                sample_result <- subset(sample_result_all, sample_result_all$pellet == 0)
                                sample_result$sample_type <- factor(sample_result$sample_type, levels = c("contorl", "saliva", "nasal_blow", "nasal_swab", "Sputum", "BAL"))
                                

                                
#Loading phyloseq object
                                
                                list.files("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20221102_PQ00331_HOST_Depletion/metaphlan3")
                                phyloseq <- microbiome::read_biom2phyloseq("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20221102_PQ00331_HOST_Depletion/metaphlan3/metaphlan3.allkingdoms.EstCount.biom")
                                readstat <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/4_Data/1_Raw/Baylor_Processed/20221102_PQ00331_HOST_Depletion/ReadStats.txt", sep = "\t")
                                readstat$SampleID
                                sample_names(phyloseq)
                                sample_result$cryo_preserve
                                
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
                                sample_result
                                
                                dummy <- merge(readstat, sample_result, by = "SampleID", all = T)
                                sample_data <- subset(dummy, !is.na(dummy$Host_mapped))
                                row.names(sample_data) <- sample_data$SampleID
                                sample_data
                                #Adding more data to the sampledata
                                sample_data$sequencing_host_prop <- sample_data$Host_mapped/sample_data$Raw_reads
                                
                                hist(sample_data$sequencing_host_prop)
                                sample_data$sample_type <- as.character(sample_data$sample_type)
                                sample_data$sample_type <- ifelse(grepl("Neg", sample_data$SampleID),
                                                                        "neg_control", sample_data$sample_type)
                                sample_data$sample_type <- ifelse(grepl("Pos", sample_data$SampleID),
                                                                  "pos_control", sample_data$sample_type)
                                otu_table <- otu_table(phyloseq) %>% data.frame(check.names = F)
                                total_reads <- otu_table %>% colSums() %>% data.frame()
                                total_reads
                                names(total_reads) <- "total_reads"
                                sample_data <- merge(sample_data, data.frame(total_reads), by = 0)
                                sample_data
                                sample_data$treated <- rowSums(sample_data[, c("lypma", "benzonase", "molysis", "host_zero", "qiaamp")])
                                
                                
                                
                                ggplot(sample_result_all, aes(x = sample_type, y = proportion, fill = treatment)) +
                                        geom_boxplot () 
                                        #xlab("Host DNA proportion by qPCR") +
                                        #ylab("(Host DNA / mapped bacterial DNA) by seqeuncing")
                                
                                #constructing phylose object
                                row.names(sample_data) <- sample_data$Row.names
                                phyloseq <- merge_phyloseq(phyloseq, sample_data(sample_data))
                                
                                

# Adding library prep results --------------------------------------------

                                lib_fail_list <- c("NS_6_b_lyPMA", "NS_8_b_lyPMA", "NS_17_b_lyPMA", "NS_21_b_host_zero", "NS_22_b_host_zero", "NS_19_d_molysis", "NS_21_d_molysis", "NS_22_d_molysis", "NS_23_d_molysis", "NS_26_b_lyPMA", "BAL_073_b_host_zero", "BAL_073_c_molysis", "BAL_078_d_lyPMA")
                                sample_data(phyloseq)$lib_failed <- row.names(sample_data(phyloseq)) %in% lib_fail_list
                                

# Factorize predicters ----------------------------------------------------
                                sample_data(phyloseq)$treatment <- tolower(sample_data(phyloseq)$treatment)
                                sample_data(phyloseq)$treatment <- factor(sample_data(phyloseq)$treatment, levels = c("control", "lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
                                sample_data(phyloseq)$sample_type <- factor(sample_data(phyloseq)$sample_type, levels = c("BAL", "nasal_swab", "Sputum"))

# Saving object -----------------------------------------------------------

                                         
                                write.csv(controls, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Tables/DAT_20230122_MGK_host_control_qPCR.csv")
                                write_rds(phyloseq, "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221129_MGK_host_tidy_tax.rds")
                                
