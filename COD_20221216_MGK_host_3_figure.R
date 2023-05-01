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



# Loading files -----------------------------------------------------------
#loading tidy phyloseq object
phyloseq <- read_rds("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221129_MGK_host_tidy_tax.rds")
phyloseq_path <- read_rds("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/4_Data/2_Tidy/Phyloseq/PHY_20221229_MGK_host_tidy_path.rds")
#control data
controls <- read.csv("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Tables/DAT_20230122_MGK_host_control_qPCR.csv")

#sample data loading
sample_data <- sample_data(phyloseq$phyloseq_count) %>% data.frame(check.names = F)
sample_data_path <- phyloseq$phyloseq_path_rpkm %>% sample_data %>% data.frame(check.names = F)

data_qPCR <- read_csv("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Tables/DAT_20230122_MGK_host_control_qPCR.csv")
#qPCR data of all the samples sent out sequencing
data_qPCR <- subset(data_qPCR, data_qPCR$baylor_other_id %in% c(sample_names(phyloseq$phyloseq_count)) | data_qPCR$baylor_other_id %in% c(read_excel("/Users/minsikkim/Dropbox (Partners HealthCare)/Project_Baylor/3_Documentation/Communications/2023-01-24_baylor_shipping_sicas2_nasal_host_depleted/CMMR_MetadataCapture_20230124_LaiP-PQ00430_SICAS2_NS.xlsx", skip = 27) %>% data.frame %>% .$`Optional..............secondary.ID`))

data_qPCR <- subset(data_qPCR, data_qPCR$sample_type %in% c("BAL", "nasal_swab", "Sputum", "neg_depletion", "pos_depletion"))
data_qPCR$sample_type <- factor(data_qPCR$sample_type, levels = c("BAL", "nasal_swab", "Sputum", "pos_depletion", "neg_depletion"),
                                labels = c("BAL", "Nasal swab", "Sputum", "P depletion", "N depletion"))
data_qPCR$treatment <- factor(data_qPCR$treatment, levels = c("control", "lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"),
                              labels = c("Control", "lyPMA", "Benzonase", "Host zero", "Molysis", "QIAamp"))



#making labe4ls for sample type
label <-  c("BAL","Nasal swab","Sputum")
names(label) <- c("BAL","Nasal swab","Sputum")


#Formattings
sample_data <- sample_data(phyloseq$phyloseq_count) %>% data.frame(check.names = F)


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




#phyloseq object
phyloseq_rel_nz = subset_samples(phyloseq$phyloseq_rel, S.obs != 0)
phyloseq_rel_nz = subset_samples(phyloseq_rel_nz, sample_type %in% c("BAL", "Nasal swab", "Sputum"))
phyloseq_path_rel = phyloseq$phyloseq_path_rpkm
phyloseq_path_rel_nz = subset_samples(phyloseq_path_rel, S.obs != 0)
phyloseq_path_rel_nz = subset_samples(phyloseq_path_rel_nz, sample_type %in% c("BAL", "Nasal swab", "Sputum"))


# FIgure 1 ----------------------------------------------------------------

#Figure 1. Overview of study design (like Birdy's Figure 1; make it simple and clear, or like Jake's Figure 1 for his CFB oxygen mSystems paper)

# FIgure 2 ----------------------------------------------------------------

#Figure 2. Host depletion effects (qPCR)


#2A: Change in total DNA (qPCR)


f2a <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = DNA_host_nondil + DNA_bac_nondil)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Total DNA by qPCR (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "A") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))

#2B: Change in human DNA (qPCR)


f2b <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = DNA_host_nondil)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Host DNA (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif")+ 
        labs(tag = "B") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))



#2C: Change in 16S DNA (qPCR)

f2c <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = DNA_bac_nondil)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridisscale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Bacterial DNA (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif")+ 
        labs(tag = "C") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))



#2D. Change in % host (qPCR)

f2d <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = host_proportion * 100)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Host DNA by qPCR (%)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "D") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))





f2a <- ggplot(data_qPCR, aes(x = sample_type, y = log10(DNA_host_nondil + DNA_bac_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("log<sub>10</sub>(qPCR total DNA)<br>(ng/μL)") +
        xlab("Sample type") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "A") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        theme(plot.tag = element_text(size = 15), axis.title.y = element_markdown()) +              # Plot title size
        guides(fill = guide_legend(nrow = 1, title = "Treatment"))


#2B: Change in human DNA (qPCR)
f2b <- ggplot(data_qPCR, aes(x = sample_type, y = log10(DNA_host_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")) +
        ylab("log<sub>10</sub>(qPCR host DNA)<br>(ng/μL)") +
        xlab("Sample type") +
        theme_classic (base_size = 12, base_family = "serif")+ 
        labs(tag = "B") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        theme(plot.tag = element_text(size = 15), axis.title.y = element_markdown()) +              # Plot title size
        guides(fill = guide_legend(nrow = 1, title = "Treatment"))
#2C: Change in 16S DNA (qPCR)
f2c <- ggplot(data_qPCR, aes(x = sample_type, y = log10(DNA_bac_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")) +
        ylab("log<sub>10</sub>(qPCR bacterial DNA)<br>(ng/μL)") +
        xlab("Sample type") +
        theme_classic (base_size = 12, base_family = "serif")+ 
        labs(tag = "C") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        theme(plot.tag = element_text(size = 15), axis.title.y = element_markdown()) +              # Plot title size
        guides(fill = guide_legend(nrow = 1, title = "Treatment"))

#2D. Change in % host (qPCR)
f2d <- ggplot(data_qPCR, aes(x = sample_type, y = log10(host_proportion * 100))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")) +
        ylab("log<sub>10</sub>(qPCR host DNA) (%)") +
        xlab("Sample type") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "D") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        theme(plot.tag = element_text(size = 15), axis.title.y = element_markdown()) +              # Plot title size
        guides(fill = guide_legend(nrow = 1, title = "Treatment"))

#output for markdown
ggarrange(f2a, f2b, f2c, f2d, common.legend = T , align = "hv")


png(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/Figure2.png",   # The directory you want to save the file in
    width = 183, # The width of the plot in inches
    height = 183, # The height of the plot in inches
    units = "mm",
    res = 600
) #fixing multiple page issue
ggarrange(f2a, f2b, f2c, f2d, common.legend = T , align = "hv")

dev.off()


# log - scale fig 2 -------------------------------------------------------


#2A: Change in total DNA (qPCR)


f2a_log <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = log10(DNA_host_nondil + DNA_bac_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("log10(Total DNA) (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "A") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))
f2a
#2B: Change in human DNA (qPCR)


f2b_log <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = log10(DNA_host_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("","Nasal swab","Sputum"))+
        ylab("log10(Host DNA) (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif")+ 
        labs(tag = "B") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))



#2C: Change in 16S DNA (qPCR)


f2c_log <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = log10(DNA_bac_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridisscale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("","Nasal swab","Sputum"))+
        ylab("log10(Bacterial DNA) (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif")+ 
        labs(tag = "C") +
        
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))



#2D. Change in % host (qPCR)

f2d <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = host_proportion * 100)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Host DNA by qPCR (%)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "D") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))



png(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/Figure2_log.png",   # The directory you want to save the file in
    width = 183, # The width of the plot in inches
    height = 183, # The height of the plot in inches
    units = "mm",
    res = 600
) #fixing multiple page issue
ggarrange(f2a_log, f2b_log, f2c_log, f2d, common.legend = T , align = "hv")

dev.off()


# Fig. SX - qPCR of controls ----------------------------------------------
controls$treatment <- factor(controls$treatment, levels = c("control","lyPMA", "benzonase", "host_zero", "molysis", "qiaamp"))

fSa <- ggplot(controls, aes(x = sample_type, y = DNA_host_nondil + DNA_bac_nondil)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         labels=c("HD_neg","HD_pos","EXT_neg", "EXT_pos"))+
        ylab("Total DNA by qPCR (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "A") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))

fSa_log <- ggplot(controls, aes(x = sample_type, y = log10(DNA_host_nondil + DNA_bac_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         labels=c("HD_neg","HD_pos","EXT_neg", "EXT_pos"))+
        ylab("log10(Total DNA) (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "A") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))
fSa
#2B: Change in human DNA (qPCR)


fSb <-ggplot(controls, aes(x = sample_type, y = DNA_host_nondil)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         labels=c("HD_neg","HD_pos","EXT_neg", "EXT_pos"))+
        ylab("Host DNA (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "B") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))

fSb_log <-ggplot(controls, aes(x = sample_type, y = log10(DNA_host_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         labels=c("HD_neg","HD_pos","EXT_neg", "EXT_pos"))+
        ylab("log10(Host DNA) (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "B") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))


#2C: Change in 16S DNA (qPCR)


fSc <- ggplot(controls, aes(x = sample_type, y = DNA_bac_nondil)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         labels=c("HD_neg","HD_pos","EXT_neg", "EXT_pos"))+
        ylab("Bacterial DNA (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "C") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))

fSc_log <- ggplot(controls, aes(x = sample_type, y = log10(DNA_bac_nondil))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         labels=c("HD_neg","HD_pos","EXT_neg", "EXT_pos"))+
        ylab("log10(Bacterial DNA) (ng/μL)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "C") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))

fSc

#2D. Change in % host (qPCR)

fSd <- ggplot(controls, aes(x = sample_type, y = host_proportion * 100)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         labels=c("HD_neg","HD_pos","EXT_neg", "EXT_pos"))+
        ylab("Host DNA (%)") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "D") +
        theme(plot.tag = element_text(size = 15), legend.position = "top") +              # Plot title size
        guides(fill = guide_legend(nrow = 1))


png(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/FigureSX.png",   # The directory you want to save the file in
    width = 183, # The width of the plot in inches
    height = 183, # The height of the plot in inches
    units = "mm",
    res = 600
) #fixing multiple page issue
ggarrange(fSa, fSb, fSc, fSd, common.legend = T , align = "hv")

dev.off()
png(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/FigureSX_log.png",   # The directory you want to save the file in
    width = 183, # The width of the plot in inches
    height = 183, # The height of the plot in inches
    units = "mm",
    res = 600
) #fixing multiple page issue
ggarrange(fSa_log, fSb_log, fSc_log, fSd, common.legend = T , align = "hv")

dev.off()



# FIgure 3 ----------------------------------------------------------------


#Figure 3. Host depletion effects (sequencing); may be better as a Table like Table 2 in Birdy's paper since the axes are problematic but the actual units are actually important
#	- Raw_reads
#	- Host_mapped
#	- % Host (we have used Host_mapped/Raw_reads in prior papers)
#	- Final_reads




#       - Raw_reads


f3a <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = Raw_reads)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Raw reads") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "A") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))

#	- Host_mapped


f3b <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = Host_mapped)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Host mapped reaeds") +
        theme_classic (base_size = 12, base_family = "serif")+ 
        labs(tag = "B") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))


#	- % Host (we have used Host_mapped/Raw_reads in prior papers)



f3d <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = Host_mapped / (Host_mapped + Metaphlan_mapped))) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Host mapped / total mapped (%)") +
        theme_classic (base_size = 12, base_family = "serif")+ 
        labs(tag = "D") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))


#	- Final_reads

f3c <- ggplot(subset(sample_data, sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(x = sample_type, y = Final_reads)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum"))+
        ylab("Final reads") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "C") +
        theme(plot.tag = element_text(size = 15)) +              # Plot title size
        guides(fill = guide_legend(nrow = 1))


png(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/Figure3.png",   # The directory you want to save the file in
    width = 183, # The width of the plot in inches
    height = 183, # The height of the plot in inches
    units = "mm",
    res = 600
    ) #fixing multiple page issue
ggarrange(f3a, f3b, f3c, f3d, common.legend = T, align = "hv")

dev.off()

# Figure 4 -----------------------------------------------------------------

# Figure 4. Alpha and beta diversity. I think we are going to have to highlight somehow the relationship between effective sequencing depth (ie Final_reads) and alpha diversity/beta diversity.
#Using the boxplots/PCoA you have shown previously alone will make the reader misinterpret the result without showing this relationship


facet_figure <- function(data, target) {
        sample_data <- sample_data(data) %>% data.frame(check.names = F)
        sample_data$treatment <- factor(sample_data$treatment, levels = c("control", "lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        sample_data$sample_type <- tolower(sample_data$sample_type)
        sample_data$`log10(Final_reads)` <- log10(sample_data$Final_reads)
        names(sample_data) <- ifelse(names(sample_data) == "S.obs", "species_richness", names(sample_data))
        names(sample_data) <- ifelse(names(sample_data) == "dbp", "berger_parker", names(sample_data))
        sample_data$sample_type <- factor(sample_data$sample_type, levels = c("bal", "Nasal swab", "sputum"))
        file1_long <- sample_data %>%
                filter(sample_type %in% c("bal", "Nasal swab", "sputum")) %>%
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


facet_figure(phyloseq, c("data_shannon", "data_invsimpson", "berger_parker")) +
#scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) +
        scale_x_discrete(name ="Sample type",
                         breaks=c("BAL","Nasal swab","Sputum"),
                         labels=c("BAL","Nasal swab","Sputum")) +
        ylab("") +
        theme_classic (base_size = 6, base_family = "serif") + 
        labs(tag = "D") +
        theme(plot.tag = element_text(size = 10)) 



label <-  c("BAL","Nasal swab","Sputum")
names(label) <- c("BAL","Nasal swab","Sputum")

f4a <-        ggplot(subset(sample_data(phyloseq), sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(y = S.obs)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Species richness") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "A") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~sample_type, labeller = labeller(sample_type = label)) + 
        guides(fill = guide_legend(nrow = 1))


f4b <-        ggplot(subset(sample_data(phyloseq), sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(y = data_shannon)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Shannon") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "B") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~sample_type, labeller = labeller(sample_type = label)) + 
        guides(fill = guide_legend(nrow = 1))

f4c <-        ggplot(subset(sample_data(phyloseq), sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(y = data_invsimpson)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Inverse simpson") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "C") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~sample_type, labeller = labeller(sample_type = label)) + 
        guides(fill = guide_legend(nrow = 1))

f4d <-        ggplot(subset(sample_data(phyloseq), sample_data$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(y = dbp)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Berger-Parker index") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "D") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~sample_type, labeller = labeller(sample_type = label)) + 
        guides(fill = guide_legend(nrow = 1))



f4ad <- ggarrange(f4a, f4b, f4c, f4d, common.legend = T, align = "hv") # alpha diversity plots


#Relative abundacnes for beta diveristy indices

pca_bray <- ordinate(phyloseq_rel_nz,  method = "PCoA", distance = "bray")

#ordination plot

f4e <- ordinate(phyloseq_rel_nz,  method = "PCoA", distance = "bray") %>%
        plot_ordination(phyloseq_rel_nz, ., col = "treatment", shape = "sample_type" ) + 
        #scale_color_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) +
        scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_shape(name = "Sample type", labels = c("BAL", "Nasal swab", "Sputum")) +
        geom_point(size = 3) +
        theme_classic (base_size = 12, base_family = "serif") +
        theme(plot.tag = element_text(size = 15), legend.spacing = unit(0, 'cm'), legend.key.height = unit(0.4, "cm")) + #legend.position = c(0.9, 0.4)
        labs(tag = "E")

f4e



#distances of betadiversity - boxplots
library(harrietr)

bray_dist_long <- distance(phyloseq_rel_nz, method="bray") %>% as.matrix() %>% melt_dist() #making long data of distance matrices
#Adding sample type and treatment name. 
#this can be also done by merging metadata into the `bray_dist_long`
names <- data.frame(str_split_fixed(bray_dist_long$iso1, "_", 3))
names2 <- data.frame(str_split_fixed(bray_dist_long$iso2, "_", 3))
bray_dist_long$sample_id_1 <- paste(names$X1, names$X2, sep = "_")
bray_dist_long$method_1 <- ifelse(grepl("control", bray_dist_long$iso1),"control", 
                                  ifelse(grepl("lyPMA", bray_dist_long$iso1),"lypma", 
                                         ifelse(grepl("benzonase", bray_dist_long$iso1),"benzonase", 
                                                ifelse(grepl("host", bray_dist_long$iso1),"host_zero", 
                                                       ifelse(grepl("qia", bray_dist_long$iso1),"qiaamp", 
                                                              ifelse(grepl("moly", bray_dist_long$iso1),"molysis", 
                                                                     NA))))))
#Adding data for iso 2 also should be done
bray_dist_long$sample_id_2 <- paste(names2$X1, names2$X2, sep = "_")
bray_dist_long$method_2 <-ifelse(grepl("control", bray_dist_long$iso2),"control", 
                                 ifelse(grepl("lyPMA", bray_dist_long$iso2),"lypma", 
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
bray_dist_long_within_sampleid_from_control$sample_type <- ifelse(grepl("NS", bray_dist_long_within_sampleid_from_control$iso1), "Nasal swab",
                                                                  ifelse(grepl("CFB", bray_dist_long_within_sampleid_from_control$iso1), "Sputum",
                                                                         ifelse(grepl("BAL", bray_dist_long_within_sampleid_from_control$iso1), "BAL", NA)))
bray_dist_long_within_sampleid_from_control
f4f <- ggplot(bray_dist_long_within_sampleid_from_control, aes(y = dist, fill = treatment)) +
        geom_boxplot() +
        #scale_fill_manual(values = c(viridis(6)[2:6])) +
        scale_fill_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Sample pair distances") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "F") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
        facet_wrap(~sample_type, labeller = labeller(sample_type = label))



png(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/Figure4.png",   # The directory you want to save the file in
    width = 183, # The width of the plot in inches
    height = 220, # The height of the plot in inches
    units = "mm",
    res = 600
) #fixing multiple page issue
ggarrange(f4ad, ggarrange(f4e, f4f, ncol = 1, align = "hv"), ncol = 1, align = "hv") +
        guides(fill = guide_legend(nrow = 1))
 # alpha diversity plots
#ggarrange(f4ad, ggarrange(f4e, f4f, ncol = 2),
#          ncol = 1) # alpha diversity plots


dev.off()


# Figure 5 ----------------------------------------------------------------



#Figure 5. Differential abundance. You have a few choices here.

#- Bubble plot (or heatmap) of relative abundance of top 20 bugs for each sample type stratified by host depletion method
#- See Figure 3D in https://www.frontiersin.org/articles/10.3389/fbinf.2021.774631/full
#- If bubble plot would have 5A be nasal, 5B be sputum, 5C be BAL
#- Something like Figure 6 in Birdy's paper where you plot bugs that were differentially abundant by treatment method
#- The problem with stacked bar charts is that you don't have enough colors for species level depiction and you lose too much granularity at the phylum or other higher level taxonomic ranks



#DA analysis - MaAslin

library(Maaslin2)
sample_data(phyloseq_rel_nz)$log10.Final_reads <- log10(sample_data(phyloseq_rel_nz)$Final_reads)
sample_data_sample <- sample_data(phyloseq_rel_nz) %>% data.frame()
otu_data_sample <- otu_table(phyloseq_rel_nz) %>% t %>% data.frame()
sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")


#Running MaAslin for all sample without decontam
#for taxa differentially abundant by host depletion method, look to see which ones overlap with potential contaminant taxa

# Maaslin - # # y ~ log(final reads) + sample_type + treatment  -----------



phyloseq_rel_nz
#all samples
sample_data_sample <- sample_data(phyloseq_rel_nz) %>% data.frame()
otu_data_sample <- otu_table(phyloseq_rel_nz) %>% t %>% data.frame()
sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")

fit_data = Maaslin2(
        input_data = otu_data_sample, 
        input_metadata = sample_data_sample, 
        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
        fixed_effects = c("sample_type", "log10.Final_reads", "lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
        transform = "LOG", #default
        normalization = "TSS", # default
        random_effects = c("original_sample"), 
        reference = c("sample_type,BAL"))

#Checking number of bugs for the main text of the manuscript 
fit_data$results %>% .$feature %>% subset(., !duplicated(.)) %>% length
fit_data$results %>% subset(., .$qval < 0.1) %>% .$feature %>% subset(., !duplicated(.)) %>% length
fit_data$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()

fit_data$results %>% .$coef %>% max
fit_data$results %>% subset(., .$metadata != "sample_type") %>% .$coef %>% max
fit_data$results %>% subset(., .$metadata != "sample_type") %>%  subset(., .$qval < 0.1) %>% .$coef %>% max
#NS
# # y ~ log(final reads) + sample_type + treatment 

phyloseq_rel_nz_ns <- subset_samples(phyloseq_rel_nz, sample_type == "Nasal swab")
sample_data_sample <- sample_data(phyloseq_rel_nz_ns) %>% data.frame()
otu_data_sample <- otu_table(phyloseq_rel_nz_ns) %>% t %>% data.frame()
sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")

fit_data_ns = Maaslin2(
        input_data = otu_data_sample, 
        input_metadata = sample_data_sample, 
        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
        transform = "LOG", #default
        normalization = "TSS", # default
        random_effects = c("original_sample"))

#Checking number of bugs for the main text of the manuscript 
fit_data_ns$results %>% .$feature %>% subset(., !duplicated(.)) %>% length
fit_data_ns$results %>% subset(., .$qval < 0.1) %>% .$feature %>% subset(., !duplicated(.)) %>% length
fit_data_ns$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()

fit_data_ns$results %>% .$coef %>% max
fit_data_ns$results %>% subset(., .$metadata != "sample_type") %>% .$coef %>% max
fit_data_ns$results %>% subset(., .$metadata != "sample_type") %>%  subset(., .$qval < 0.1) %>% .$coef %>% max

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
        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
        transform = "LOG", #default
        normalization = "TSS", # default
        random_effects = c("original_sample"))

fit_data_bal$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()


#sputum
# # y ~ log(final reads) + sample_type + treatment 
phyloseq_rel_nz_spt <- subset_samples(phyloseq_rel_nz, sample_type == "Sputum")
sample_data_sample <- sample_data(phyloseq_rel_nz_spt) %>% data.frame(check.names = F)
otu_data_sample <- otu_table(phyloseq_rel_nz_spt) %>% t %>% data.frame(check.names = F)
sample_data_sample$sampletype_treatment <- paste(sample_data_sample$sample_type, sample_data_sample$treatment, sep = ":")


fit_data_spt = Maaslin2(
        input_data = otu_data_sample, 
        input_metadata = sample_data_sample, 
        output = "/Users/minsikkim/Dropbox (Partners HealthCare)/@minsik/project_host_dna_depletion/Data/maaslin_output", 
        fixed_effects = c("log10.Final_reads","lypma", "benzonase", "host_zero", "molysis", "qiaamp"), 
        transform = "LOG", #default
        normalization = "TSS", # default
        random_effects = c("original_sample"))


#Checking number of bugs for the main text of the manuscript 

fit_data_spt$results %>% .$feature %>% subset(., !duplicated(.)) %>% length
fit_data_spt$results %>% subset(., .$qval < 0.1) %>% .$feature %>% subset(., !duplicated(.)) %>% length
fit_data_spt$results %>% subset(., .$qval < 0.1) %>% .$metadata %>% table()

fit_data_spt$results %>% .$coef %>% max
fit_data_spt$results %>% subset(., .$metadata != "sample_type") %>% .$coef %>% max
fit_data_spt$results %>% subset(., .$metadata != "sample_type") %>%  subset(., .$qval < 0.1) %>% .$coef %>% max




#Making significance table for figure

fit_data_spt_20 <- spread(fit_data_spt$results[order(fit_data_spt$results$qval),c("feature", "metadata", "qval")], metadata, qval)
fit_data_spt_20$min <- apply(fit_data_spt_20, 1, FUN = min)
fit_data_spt_20 <- fit_data_spt_20[order(fit_data_spt_20$min),] %>% .[1:20,]

fit_data_ns_20 <- spread(fit_data_ns$results[order(fit_data_ns$results$qval),c("feature", "metadata", "qval")], metadata, qval)
fit_data_ns_20$min <- apply(fit_data_ns_20, 1, FUN = min)
fit_data_ns_20 <- fit_data_ns_20[order(fit_data_ns_20$min),] %>% .[1:20,]
fit_data_ns_20$feature <-  ifelse(fit_data_ns_20$feature == "X.Collinsella._massiliensis",
                                  "[Collinsella]_massiliensis", fit_data_ns_20$feature)

fit_data_bal_20 <- spread(fit_data_bal$results[order(fit_data_bal$results$qval),c("feature", "metadata", "qval")], metadata, qval)
fit_data_bal_20$min <- apply(fit_data_bal_20, 1, FUN = min)
fit_data_bal_20 <- fit_data_bal_20[order(fit_data_bal_20$min),] %>% .[1:20,]
taxa_names(phyloseq_rel_nz_ns)[grepl("Collin", taxa_names(phyloseq_rel_nz_ns))]

phyloseq_rel_nz_spt_sig <- subset_taxa(phyloseq_rel_nz_spt, taxa_names(phyloseq_rel_nz_spt) %in% fit_data_spt_20$feature)
phyloseq_rel_nz_bal_sig <- subset_taxa(phyloseq_rel_nz_bal, taxa_names(phyloseq_rel_nz_bal) %in% fit_data_bal_20$feature)
phyloseq_rel_nz_ns_sig <- subset_taxa(phyloseq_rel_nz_ns, taxa_names(phyloseq_rel_nz_ns) %in% fit_data_ns_20$feature)

phyloseq_rel_nz_spt_sig_fig <- cbind(phyloseq_rel_nz_spt_sig %>% otu_table %>% t, phyloseq_rel_nz_spt_sig %>% sample_data) %>% group_by(treatment) %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% .[, 1:21] %>% column_to_rownames(., "treatment") %>% t ()
phyloseq_rel_nz_bal_sig_fig <- cbind(phyloseq_rel_nz_bal_sig %>% otu_table %>% t, phyloseq_rel_nz_bal_sig %>% sample_data) %>% group_by(treatment) %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% .[, 1:21] %>% column_to_rownames(., "treatment") %>% t ()
phyloseq_rel_nz_ns_sig_fig <- cbind(phyloseq_rel_nz_ns_sig %>% otu_table %>% t, phyloseq_rel_nz_ns_sig %>% sample_data) %>% group_by(treatment) %>% summarise_if(is.numeric, mean, na.rm = TRUE)%>% .[, 1:21] %>% column_to_rownames(., "treatment") %>% t ()


species_italic <- function(data) {
        names <- gsub("_", " ", rownames(data))
        names <- gsub("[]]|[[]", "", names)
        names <- gsub(" sp", " sp.", names)
        names <- gsub(" sp.", "* sp.", names)
        names <- gsub(" group", "* group.", names)
        names <- ifelse(grepl("[*]", names), 
                                                        paste("*", names, sep = ""), 
                                                        paste("*", names, "*", sep = ""))
        rownames(data) <- names
        data
}


phyloseq_rel_nz_spt_sig_fig <- species_italic(phyloseq_rel_nz_spt_sig_fig)
phyloseq_rel_nz_ns_sig_fig <- species_italic(phyloseq_rel_nz_ns_sig_fig)
phyloseq_rel_nz_bal_sig_fig <- species_italic(phyloseq_rel_nz_bal_sig_fig)
phyloseq_rel_nz_ns_sig_fig

#qval sig table


        sig_data_ns <- ifelse(fit_data_ns_20 %>% rownames_to_column(var = "-") %>% column_to_rownames(var = "feature") %>% species_italic %>% select(-c("-", "min", "log10.Final_reads")) < 0.1,
                              "*", NA)
        
        sig_data_bal <- ifelse(fit_data_bal_20 %>% rownames_to_column(var = "-") %>% column_to_rownames(var = "feature") %>% species_italic %>% select(-c("-", "min", "log10.Final_reads")) < 0.1,
                               "*", NA)
        
        sig_data_spt <- ifelse(fit_data_spt_20 %>% rownames_to_column(var = "-") %>% column_to_rownames(var = "feature") %>% species_italic %>% select(-c("-", "min", "log10.Final_reads")) < 0.1,
                               "*", NA)
        sig_data_spt <- gather(data.frame(sig_data_spt) %>% rownames_to_column(var = "feature"), treatment, significance, benzonase:qiaamp, factor_key=TRUE)
        sig_data_ns <- gather(data.frame(sig_data_ns) %>% rownames_to_column(var = "feature"), treatment, significance, benzonase:qiaamp, factor_key=TRUE)
        sig_data_bal <- gather(data.frame(sig_data_bal) %>% rownames_to_column(var = "feature"), treatment, significance, benzonase:qiaamp, factor_key=TRUE)

        sig_data_spt$treatment <- factor(sig_data_spt$treatment, levels =c("control","lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        sig_data_bal$treatment <- factor(sig_data_bal$treatment, levels =c("control","lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        sig_data_ns$treatment <- factor(sig_data_ns$treatment, levels =c("control","lypma", "benzonase", "host_zero", "molysis", "qiaamp"))
        
phyloseq_rel_nz_bal_sig_fig[phyloseq_rel_nz_bal_sig_fig == 0] <- NA
phyloseq_rel_nz_ns_sig_fig[phyloseq_rel_nz_ns_sig_fig == 0] <- NA
phyloseq_rel_nz_spt_sig_fig[phyloseq_rel_nz_spt_sig_fig == 0] <- NA
#Volcano plot


f5a <- ggplot(fit_data$results, aes(y = -log10(qval), x = coef, col = metadata)) +
        theme_classic(base_family = "serif") +
        labs(tag = "A") +
        geom_point(size = 2) +
        xlab("MaAslin coefficient") +
        ylab("-log<sub>10</sub>(*q*-value)") +
        geom_hline(yintercept = 1, col = "gray") +
        geom_vline(xintercept = 0, col = "gray") +
        geom_richtext(aes( 4, 8, label = "*q*-value = 0.1, fold-change = 0", vjust = -1, fontface = 1), col = "grey", size = 3, family = "serif") +
        theme(legend.position = "top", axis.title.y = ggtext::element_markdown()) +
        scale_color_manual(values = c("#4daf4a",  "#984ea3", "#f781bf", "#377eb8", "#ff7f00", "#ffff33", "#a65628"), labels = c("Benzonase", "Host zero", "log10 (Final reads)", "lyPMA", "Molysis",  "QIAaamp", "Sample type")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        guides(col = guide_legend(title = "Fixed effects", title.position = "top", nrow = 3))


#ffff33 qia

f5b <- ggballoonplot(phyloseq_rel_nz_ns_sig_fig, fill = "value") +
        theme_classic(base_family = "serif") +
        scale_fill_viridis() +
        xlab("Experimental group") +
        ylab("Species") +
        labs(tag = "B") +
        theme(panel.grid.major = element_line(colour = "grey"), legend.position = "top", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), axis.text.y = ggtext::element_markdown()) +
        scale_x_discrete(labels=c("control" = "Control", "lypma" = "lyPMA", "benzonase" = "Benzonase", "host_zero" = "Host-zero", "molysis" = "Molysis", "qiaamp" = "QIAamp")) +
        geom_text(data = data.frame(subset(sig_data_ns, !is.na(sig_data_ns$significance))), aes(y = feature, x = treatment, label = significance, col = significance), hjust = -0.5, vjust = 0.8, size = 5) +
        guides(col = guide_legend(nrow = 3, override.aes = aes(label = "*", size = 10, color = "red"), title="Significance", title.position = "top", order = 3), size = guide_legend(title = "Relative abundance", title.position = "top", order = 1, nrow = 2), fill = F) + 
        scale_color_manual(values = c("red", "black"), labels = c(expression(paste(italic("q"), "-value < 0.1", sep = "")), NA))
        
      


legend_1 <- get_legend(f5a)
legend_2 <- get_legend(f5b)

rm_legend <- function(p){p + theme(legend.position = "none")}


f5c <- ggballoonplot(phyloseq_rel_nz_spt_sig_fig, fill = "value") +
        theme_classic(base_family = "serif") +
        scale_fill_viridis() +
        xlab("Experimental group") +
        ylab("Species") +
        labs(tag = "C") +
        theme(panel.grid.major = element_line(colour = "grey"), legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), axis.text.y = ggtext::element_markdown())  +
        geom_text(data = data.frame(subset(sig_data_spt, !is.na(sig_data_spt$significance))), aes(y = feature, x = treatment, label = significance, col = significance), hjust = -0.5, vjust = 0.8, size = 5) +
        guides(fill = FALSE, shape = F, col = F) +
        scale_x_discrete(labels=c("control" = "Control", "lypma" = "lyPMA", "benzonase" = "Benzonase", "host_zero" = "Host-zero", "molysis" = "Molysis", "qiaamp" = "QIAamp"))


f5d <- ggballoonplot(phyloseq_rel_nz_bal_sig_fig, fill = "value") +
        theme_classic(base_family = "serif") +
        scale_fill_viridis() +
        xlab("Experimental group") +
        ylab("Species") +
        labs(tag = "D") +
        theme(panel.grid.major = element_line(colour = "grey"), legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), axis.text.y = ggtext::element_markdown()) +
        geom_text(data = data.frame(subset(sig_data_bal, !is.na(sig_data_bal$significance))), aes(y = feature, x = treatment, label = significance, col = significance), hjust = -0.5, vjust = 0.8, size = 5) +
        scale_x_discrete(labels=c("control" = "Control", "lypma" = "lyPMA", "benzonase" = "Benzonase", "host_zero" = "Host-zero", "molysis" = "Molysis", "qiaamp" = "QIAamp")) +
        guides(fill = FALSE, shape = F, col = F)





png(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/Figure5.png",   # The directory you want to save the file in
    width = 183, # The width of the plot in inches
    height = 220, # The height of the plot in inches
    units = "mm",
    res = 600
) #fixing multiple page issue

ggarrange(ggarrange(legend_1, legend_2, align = "hv"),
        ggarrange(rm_legend(f5a), rm_legend(f5b)),
        ggarrange(f5c, f5d), ncol = 1, heights = c(2, 5, 5))

dev.off()


# Figure 6 ----------------------------------------------------------------

f6a <-        ggplot(subset(sample_data(phyloseq_path), sample_data_path$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(y = S.obs)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Species richness") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "A") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~sample_type) + 
        guides(fill = guide_legend(nrow = 1))


f6b <-        ggplot(subset(sample_data(phyloseq_path), sample_data_path$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(y = data_shannon)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Shannon") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "B") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~sample_type) + 
        guides(fill = guide_legend(nrow = 1))

f6c <-        ggplot(subset(sample_data(phyloseq_path), sample_data_path$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(y = data_invsimpson)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Inverse simpson") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "C") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~sample_type) + 
        guides(fill = guide_legend(nrow = 1))

f6d <-        ggplot(subset(sample_data(phyloseq_path), sample_data_path$sample_type %in% c("Sputum", "Nasal swab", "BAL")), aes(y = dbp)) +
        geom_boxplot(aes(fill = treatment), lwd = 0.2) +
        #scale_fill_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + # color using viridis
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Berger-Parker index") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "D") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(~sample_type) + 
        guides(fill = guide_legend(nrow = 1))



f6ad <- ggarrange(f6a, f6b, f6c, f6d, common.legend = T, align = "hv") # alpha diversity plots


#Relative abundacnes for beta diveristy indices


pca_bray_path <- ordinate(phyloseq_path_rel_nz,  method = "PCoA", distance = "bray")

#ordination plot

f6e <- ordinate(phyloseq_path_rel_nz,  method = "PCoA", distance = "bray") %>%
        plot_ordination(phyloseq_path_rel_nz, ., col = "treatment", shape = "sample_type" ) + 
        #scale_color_viridis(discrete = 6, name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) +
        scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        scale_shape(name = "Sample type", labels = c("BAL", "Nasal swab", "Sputum")) +
        geom_point(size = 3) +
        theme_classic (base_size = 12, base_family = "serif") +
        theme(plot.tag = element_text(size = 15), legend.spacing = unit(0, 'cm'), legend.key.height = unit(0.4, "cm")) + #legend.position = c(0.9, 0.4)
        labs(tag = "E")


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
path_bray_dist_long_within_sampleid_from_control$sample_type <- ifelse(grepl("NS", path_bray_dist_long_within_sampleid_from_control$iso1), "Nasal swab",
                                                                  ifelse(grepl("CFB", path_bray_dist_long_within_sampleid_from_control$iso1), "Sputum",
                                                                         ifelse(grepl("BAL", path_bray_dist_long_within_sampleid_from_control$iso1), "BAL", NA)))


f6f <- ggplot(path_bray_dist_long_within_sampleid_from_control, aes(y = dist, fill = treatment)) +
        geom_boxplot() +
        #scale_fill_manual(values = c(viridis(6)[2:6])) +
        scale_fill_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), name = "Treatment", labels = c("Control","lyPMA", "Benzonase", "Host zero", "Molysis", "QIAaamp")) + #color using https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
        ylab("Sample pair distances") +
        theme_classic (base_size = 12, base_family = "serif") + 
        labs(tag = "F") +
        theme(plot.tag = element_text(size = 15),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
        facet_wrap(~sample_type, labeller = labeller(sample_type = label))



png(file = "/Users/minsikkim/Dropbox (Partners HealthCare)/Project_SICAS2_microbiome/7_Manuscripts/2022_MGK_Host_Depletion/Figures/Figure6.png",   # The directory you want to save the file in
    width = 183, # The width of the plot in inches
    height = 220, # The height of the plot in inches
    units = "mm",
    res = 600
) #fixing multiple page issue
ggarrange(f6ad, ggarrange(f6e, f6f, ncol = 1, align = "hv"), ncol = 1, align = "hv") +
        guides(fill = guide_legend(nrow = 1))
# alpha diversity plots
#ggarrange(f4ad, ggarrange(f4e, f4f, ncol = 2),
#          ncol = 1) # alpha diversity plots


dev.off()







