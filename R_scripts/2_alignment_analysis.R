library(pacman)
pacman::p_load(tidyverse, plotly, grid, gridExtra, ggrepel, readxl)

setwd("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/2. Scripts/")

# Import and cleanup ----------------------------------------------------------

# import
file_path <- "/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. alignments/01_raw_alignments/Alignment_coverage_depth.tsv"
alignments <- read.table(file_path, sep = "\t", header = TRUE)

# Every sample possesses 26 rows in a very messy format
# we can add the name of the sample to every group of sequences iteratively 
group_size <- 26
num_groups <- ceiling(nrow(alignments) / group_size)
alignments$sample <- rep(1:num_groups, each = group_size, length.out = nrow(alignments))
sample_names <- unique(sub("_.*", "", alignments[seq(1, nrow(alignments), group_size), 1]))
alignments$sample <- sample_names[alignments$sample]

#Move the sample column to first place for easier visualization
alignments <- alignments %>%
  mutate(sample_col = sample) %>%
  select(-sample) %>% 
  relocate(sample_col, .before = 1) 

# Sample name is added to each column to avoid confusion
alignments <- alignments %>%
  group_by(sample_col) %>%
  mutate(group_prefix = str_extract(first(name), "^[^_]+")) %>%
  ungroup() %>%
  mutate(name = ifelse(!str_starts(name, group_prefix), paste(group_prefix, name, sep = "_"), name)) %>%
  select(-group_prefix)

# substitute alignmentsbase names with species names
alignments$name <- gsub("NZ_AP022547\\.1", "Klebsiella_michiganensis_chromosome", alignments$name)
alignments$name <- gsub("NZ_AP022548\\.1", "Klebsiella_michiganensis_plasmid", alignments$name)
alignments$name <- gsub("NZ_CP014022\\.1", "Staphylococcus_lugdunensis_chromosome", alignments$name)
alignments$name <- gsub("BS\\.pilon\\.polished\\.v3\\.ST170922", "Bacillus_subtilis_chromosome", alignments$name)

# Filter out all the headers inside the dataframe
alignments <- alignments %>%
  filter(!str_detect(startpos, "startpos"))

# Split the name string into different columns and rtemove the old one
alignments$name <- str_remove(alignments$name, "^[^_]+_")
alignments <- alignments %>%
  rowwise() %>%
  mutate(species = str_extract(name, "[^_]+_[^_]+"),
         genetic_element = str_extract(name, "[^_]+$"))
alignments$name <- NULL

for (i in 1:nrow(alignments)) {
  value <- alignments$sample_col[i]
  if (!grepl("^[0-9]", value)) {
    alignments$sample_col[i] <- paste0("1", value)
  }
}

# Reorder all rows based on sample name
alignments <- alignments[order(alignments$sample_col), ]

# Convert all numbers stored as strings to numbers
alignments$coverage <- as.numeric(alignments$coverage)
alignments$meandepth <- as.numeric(alignments$meandepth)

alignments <- alignments[!(grepl("^1", alignments$sample_col) & (alignments$species %in% c("Klebsiella_michiganensis", "Staphylococcus_lugdunensis"))), ]
alignments$genetic_element <- gsub("genome", "chromosome", alignments$genetic_element)

#write_csv(alignments, "/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. alignments/01_raw_alignments/Alignment_coverage_depth_processed.csv")

# -----------------------------------------------------------------------------

alignments <- read_csv("/Users/vbenvenga/Desktop/master_vanni/food_patogen_metagenomics/1. alignments/01_raw_alignments/Alignment_coverage_depth_processed.csv")

sample_col <- "7KBAK"
bacteria <- factor(rev(c("Bacillus_subtilis", "Listeria_monocytogenes", "Staphylococcus_aureus",
                     "Staphylococcus_lugdunensis", "Enterococcus_faecalis", "Lactobacillus_fermentum",
                     "Escherichia_coli", "Klebsiella_michiganensis", "Salmonella_enterica",
                     "Pseudomonas_aeruginosa")),
                   levels = rev(c("Bacillus_subtilis", "Listeria_monocytogenes", "Staphylococcus_aureus",
                              "Staphylococcus_lugdunensis", "Enterococcus_faecalis", "Lactobacillus_fermentum",
                              "Escherichia_coli", "Klebsiella_michiganensis", "Salmonella_enterica",
                              "Pseudomonas_aeruginosa")))

bacteria_colors <- c("#9dd58a", "#017bfc", "#39a30c", "#ff3c91", "#ffb337",
                     "#df9cff", "#9c2620", "#ffb4a6", "#8c2c6e", "#7d4237")

ggplot(alignments[alignments$sample_col == sample_col & grepl("chromosome", alignments$genetic_element),],
       aes(x = factor(species, levels = bacteria), y = coverage, width = meandepth/50, fill = species)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(meandepth, 2)), hjust = -0.2, colour = "black") +
  scale_fill_manual(name = "Species", breaks = bacteria, values = bacteria_colors) +
  theme_bw() +
  ylim(0, 105) +
  ylab("Coverage [%]") +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  geom_vline(xintercept = which(levels(bacteria) == "Lactobacillus_fermentum") - 0.5, linetype = "dashed", color = "black") +
  coord_flip() 

ggsave(paste0("r_figures/experiment_4/experiment_4_alignment_",sample_col,".png"), dpi = 500, units = "cm", height = 15, width = 22)


