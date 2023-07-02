library(pacman)
pacman::p_load(tidyverse, plotly, grid, gridExtra, ggrepel, readxl)

setwd("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/2. Scripts/")
sheets <- readxl::excel_sheets("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. Data/01_raw_data/dsDNA_IMM_Sequencing.xlsx")
reads_data <- readxl::read_excel("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. Data/01_raw_data/dsDNA_IMM_Sequencing.xlsx", sheet = sheets[1])

#Batch 1 -----------------------------------------------------------------------

reads_data$Definition <- factor(reads_data$Definition)

ggplot(reads_data, aes(x = factor(Definition, levels = c("Homogenization", "Washing", "Swabbing","Control: Direct homogenization",
                                                         "Control: extracted bacteria","Control: Buffer",
                                                         "Control: DNA")), 
                       y = (prokaryote_total_bases/total_bases), color = Kit, shape=as.character(Zymo_community_uL))) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     labels = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     values = c("#E71D36", "#619CFF", "#00BA38")) +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Samples") + ylab("Proakyrotic bases/total bases") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 20, hjust = 1)) 
ggsave("r_figures/experiment_1/experiment_1_prokaryotic_proportion.png", units = "cm", height = 15, width = 22)

 ggplot(reads_data, aes(x = factor(Definition, levels = c("Homogenization", "Washing", "Swabbing","Control: Direct homogenization",
                                                        "Control: extracted bacteria","Control: Buffer",
                                                        "Control: DNA")),
                      y = prokaryote_total_bases,
                      shape = as.character(Zymo_community_uL))) +
  geom_point(alpha = 0.9, size = 3, aes(color = Kit)) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     labels = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     values = c("#E71D36", "#619CFF", "#00BA38")) +
  scale_y_log10() +
  geom_hline(yintercept = 3e8, linetype = "dashed", color = "black") +
  annotate("text", x = 2.5, y = 5e8, label = "100% coverage at 10x depth") +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Samples") +
  ylab("Prokaryotic bases") +
  theme_bw() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("r_figures/experiment_1/experiment_1_prokaryotic_bases.png", dpi=500, units ="cm", height = 15,width = 22)

ggplot(reads_data, aes(x = prokaryote_median_read_quality, y = prokaryote_total_bases, color = Kit, shape = as.character(Zymo_community_uL), label = Definition)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     labels = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     values = c("#E71D36", "#619CFF", "#00BA38")) +
  xlim(9, NA) +
  geom_text(hjust = 1.2, nudge_y = 0.01, size = 2) +
  scale_y_log10() +
  geom_hline(yintercept = 3e8, linetype = "dashed", color = "black") +
  annotate("text", x = 13, y = 5e8, size = 3, label = "100% coverage at 10x depth") +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Median Read quality [Q score]") +
  ylab("Prokaryotic bases") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("r_figures/experiment_1/experiment_1_prokaryotic_bases_vs_quality.png", dpi = 500, units = "cm", height = 15, width = 22)


ggplot(reads_data, aes(x = prokaryote_median_read_length / 1000, y = prokaryote_total_bases, color = Kit, shape = as.character(Zymo_community_uL), label = Definition)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     labels = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     values = c("#E71D36", "#619CFF", "#00BA38")) +
  xlim(NA, 3) +
  geom_text(hjust = -0.4, nudge_y = 0.01, size = 2) +
  scale_y_log10() +
  geom_hline(yintercept = 3e8, linetype = "dashed", color = "black") +
  annotate("text", x = 2, y = 5e8, size = 3, label = "100% coverage at 10x depth") +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Median read length [kbp]") +
  ylab("Prokaryotic bases") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("r_figures/experiment_1/experiment_1_prokaryotic_bases_vs_length.png", dpi = 500, units = "cm", height = 15, width = 22)


ggplot(reads_data, aes(x = prokaryote_median_read_length / 1000, y = prokaryote_median_read_quality, color = Kit, shape = as.character(Zymo_community_uL), label = Definition)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     labels = c("Cultured Cells", "Fecal Microbiome", "DNA Standard"),
                     values = c("#E71D36", "#619CFF", "#00BA38")) +
  xlim(NA, 3) +
  geom_text_repel(hjust = -0.4, nudge_y = 0.01, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Median prokaryotic read length [kbp]") +
  ylab("Median read quality of prokaryotic reads [Q score]") +theme_bw()
ggsave("r_figures/experiment_1/experiment_1_quality_vs_length.png", dpi = 500, units = "cm", height = 15, width = 24)



#------------------------------------------------------------------------------

sheets <- readxl::excel_sheets("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. Data/01_raw_data/dsDNA_IMM_Sequencing.xlsx")
reads_data <- readxl::read_excel("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. Data/01_raw_data/dsDNA_IMM_Sequencing.xlsx", sheet = sheets[2])

reads_data$Definition <- factor(reads_data$Definition)

ggplot(reads_data, aes(x = factor(Definition, levels = c("Homogenization", "Washing", "Swabbing", "Control: Direct homogenization", "Control: extracted zymo mock + in-house", "Control: extracted in-house strains", "Control: Buffer")), y = (prokaryote_total_bases / total_bases), color = Kit)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Fecal Microbiome", "Purefood Pathogen", "Qiagen DNA Microbiome"),
                     labels = c("Cultured Cells", "Fecal Microbiome", "Purefood Pathogen", "Qiagen DNA Microbiome"),
                     values = c("#E71D36", "#619CFF", "#C77CFF", "#FF9F1C")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Samples") +
  ylab("Prokaryotic bases/total bases") +theme_bw()

ggsave("r_figures/experiment_2/experiment_2_prokaryotic_proportion.png", units = "cm", height = 15, width = 22)


ggplot(reads_data, aes(x = factor(Definition, levels = c("Homogenization", "Washing", "Swabbing", "Control: Direct homogenization", "Control: extracted zymo mock + in-house", "Control: extracted in-house strains", "Control: Buffer")), y = prokaryote_total_bases, color = Kit)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Fecal Microbiome", "Purefood Pathogen", "Qiagen DNA Microbiome"),
                     labels = c("Cultured Cells", "Fecal Microbiome", "Purefood Pathogen", "Qiagen DNA Microbiome"),
                     values = c("#E71D36", "#619CFF", "#C77CFF", "#FF9F1C")) +
  scale_y_log10() +
  geom_hline(yintercept = 3.9e8, linetype = "dashed", color = "black") +
  annotate("text", x = 2, y = 5e8, label = "100% coverage at 10x depth") +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Samples") +
  ylab("Prokaryotic bases") + theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("r_figures/experiment_2/experiment_2_prokaryotic_bases.png", dpi = 500, units = "cm", height = 15, width = 22)

ggplot(reads_data, aes(x = prokaryote_median_read_length / 1000, y = prokaryote_median_read_quality, color = Kit, shape = as.character(Zymo_community_uL), label = Definition)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Fecal Microbiome", "Purefood Pathogen", "Qiagen DNA Microbiome"),
                     labels = c("Cultured Cells", "Fecal Microbiome", "Purefood Pathogen", "Qiagen DNA Microbiome"),
                     values = c("#E71D36", "#619CFF", "#C77CFF", "#FF9F1C")) +
  xlim(NA, 3) +
  geom_text_repel(hjust = -0.4, nudge_y = 0.01, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Median read length [kbp]") +
  ylab("Median read quality [Q score]") +theme_bw()
ggsave("r_figures/experiment_2/experiment_2_quality_vs_length.png", dpi = 500, units = "cm", height = 15, width = 24)

#------------------------------------------------------------------------------

sheets <- readxl::excel_sheets("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. Data/01_raw_data/dsDNA_IMM_Sequencing.xlsx")
reads_data <- readxl::read_excel("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. Data/01_raw_data/dsDNA_IMM_Sequencing.xlsx", sheet = sheets[3])

reads_data$Definition <- factor(reads_data$Definition)

ggplot(reads_data, aes(x = factor(Definition, levels = c("Host depleted sample","Control: Zymo mock community","Control: Klebsiella and Staphylococcus",               
                                                        "Control: Zymo mock + K. michiganensis and S. lugdunensis","Control: Direct extraction Zymo + K. michiganensis and S. lugdunensis",
                                                        "Control: Chicken","Control: Buffer", "Control: DNA", "Control: High molecular weight DNA")),
                       y = (prokaryote_total_bases / total_bases), color = Kit)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Purefood Pathogen", "DNA Standard"),
                     labels = c("Cultured Cells", "Purefood Pathogen", "DNA Standard"),
                     values = c("#E71D36", "#C77CFF", "#00BA38")) +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Samples") +
  ylab("Prokaryotic bases/total bases") +theme_bw() +theme(axis.text.x = element_text(angle = 50, hjust = 1)) + ylim(0,1)

ggsave("r_figures/experiment_3/experiment_3_prokaryotic_proportion.png", units = "cm", height = 25, width = 22)


ggplot(reads_data, aes(x = factor(Definition, levels = c("Host depleted sample","Control: Zymo mock community","Control: Klebsiella and Staphylococcus",               
                                                         "Control: Zymo mock + K. michiganensis and S. lugdunensis","Control: Direct extraction Zymo + K. michiganensis and S. lugdunensis",
                                                         "Control: Chicken","Control: Buffer", "Control: DNA", "Control: High molecular weight DNA")),
                       y = prokaryote_total_bases, color = Kit)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Purefood Pathogen", "DNA Standard"),
                     labels = c("Cultured Cells", "Purefood Pathogen", "DNA Standard"),
                     values = c("#E71D36", "#C77CFF", "#00BA38")) +
  scale_y_log10() +
  geom_hline(yintercept = 3.9e8, linetype = "dashed", color = "black") +
  annotate("text", x = 8, y = 5e8, label = "100% coverage at 10x depth") +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Samples") +
  ylab("Prokaryotic bases") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("r_figures/experiment_3/experiment_3_prokaryotic_bases.png", dpi = 500, units = "cm", height = 25, width = 22)

ggplot(reads_data, aes(x = prokaryote_median_read_length / 1000, y = prokaryote_median_read_quality, color = Kit, label = Definition)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells", "Purefood Pathogen", "DNA Standard"),
                     labels = c("Cultured Cells", "Purefood Pathogen", "DNA Standard"),
                     values = c("#E71D36", "#C77CFF", "#00BA38")) +
  xlim(NA, 3) +
  geom_text_repel(hjust = -0.4, nudge_y = 0.01, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(color = "Extraction Kit") +
  xlab("Median prokaryotic read length [kbp]") +
  ylab("Median read quality of prokaryotic reads [Q score]") +theme_bw()
ggsave("r_figures/experiment_3/experiment_3_quality_vs_length.png", dpi = 500, units = "cm", height = 15, width = 24)

#------------------------------------------------------------------------------

sheets <- readxl::excel_sheets("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. Data/01_raw_data/dsDNA_IMM_Sequencing.xlsx")
reads_data <- readxl::read_excel("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. Data/01_raw_data/dsDNA_IMM_Sequencing.xlsx", sheet = sheets[4])

reads_data$Definition <- factor(reads_data$Definition)

ggplot(reads_data, aes(x = factor(Definition, levels = c("Host depleted sample","Host depleted mock community",                        
                                                         "Direct extraction: mock community","Direct extraction: B. subtilis",         
                                                         "Direct extraction: P. aeruginosa","Direct extraction: E. coli",          
                                                         "Direct extraction: K. michiganensis","Direct extraction: S. enterica","Direct extraction: S. lugdunensis",
                                                         "Direct extraction: L. monocytogenes","Control: chicken only","Negative control: extraction",
                                                         "High molecular weight DNA standard","Negative control: Sequencing")),
                       y = (prokaryote_total_bases / total_bases), color = Kit)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction Kit",
                     breaks = c("Cultured Cells","DNA Standard","Sequencing buffer"),
                     labels = c("Cultured Cells","DNA Standard","Sequencing buffer"),
                     values = c("#E71D36", "#00BA38", "#808080")) +
  labs(color = "Extraction Kit") +
  xlab("Samples") +
  ylab("Prokaryotic bases/total bases") +theme_bw() +theme(axis.text.x = element_text(angle = 40, hjust = 1)) + ylim(0,1)

ggsave("r_figures/experiment_4/experiment_4_prokaryotic_proportion.png", units = "cm", height = 15, width = 22)


ggplot(reads_data, aes(x = factor(Definition, levels = c("Host depleted sample","Host depleted mock community",                        
                                                         "Direct extraction: mock community","Direct extraction: B. subtilis",         
                                                         "Direct extraction: P. aeruginosa","Direct extraction: E. coli",          
                                                         "Direct extraction: K. michiganensis","Direct extraction: S. enterica","Direct extraction: S. lugdunensis",
                                                         "Direct extraction: L. monocytogenes","Control: chicken only","Negative control: extraction",
                                                         "High molecular weight DNA standard","Negative control: Sequencing")),
                       y = prokaryote_total_bases, color = Kit)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction",
                     breaks = c("Cultured Cells","DNA Standard","Sequencing buffer"),
                     labels = c("Cultured Cells","DNA Standard","Sequencing buffer"),
                     values = c("#E71D36", "#00BA38", "#808080")) +
  scale_y_log10() +
  geom_hline(yintercept = 3.2e8, linetype = "dashed", color = "black") +
  annotate("text", x = 12, y = 5e8, label = "100% coverage at 10x depth") +
  labs(color = "Extraction Kit", shape = "Zymo mock community spiked [uL]") +
  xlab("Samples") +
  ylab("Prokaryotic bases") + theme_bw() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("r_figures/experiment_4/experiment_4_prokaryotic_bases.png", dpi = 500, units = "cm", height = 15, width = 22)

ggplot(reads_data, aes(x = prokaryote_median_read_length / 1000, y = prokaryote_median_read_quality, color = Kit, label = Definition)) +
  geom_point(alpha = 0.9, size = 3) +
  scale_color_manual(name = "Extraction",
                     breaks = c("Cultured Cells","DNA Standard","Sequencing buffer"),
                     labels = c("Cultured Cells","DNA Standard","Sequencing buffer"),
                     values = c("#E71D36", "#00BA38", "#808080")) +
  xlim(NA, 3) +
  geom_text_repel(hjust = -0.4, nudge_y = 0.01, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(color = "Extraction Kit") +
  xlab("Median prokaryotic read length [kbp]") +
  ylab("Median read quality of prokaryotic reads [Q score]") +theme_bw()
ggsave("r_figures/experiment_4/experiment_4_quality_vs_length.png", dpi = 500, units = "cm", height = 15, width = 24)

