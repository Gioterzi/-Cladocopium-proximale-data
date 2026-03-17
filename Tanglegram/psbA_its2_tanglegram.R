library(dendextend)
library(ggsci)
library(tidyverse)
library(dplyr)
library(forcats)
library(reshape2)
library(stringr)
library(tidyr)
library(tibble)
library(sangerseqR)
library(DECIPHER)
library(Biostrings)
library(phangorn)
library(ape)
library(ggplot2)
library(ggtree)
library(patchwork)
library(bioseq)
library(kmer)
library(GUniFrac)
library(seqinr)
library(vegan)
library(corrplot)
library(ggrepel)
library(ggmsa)
library(dendextend)
library(usedist)
library(readxl)
library(phytools)
library(RColorBrewer)
library(Polychrome)
library(svglite)



dna_to_DNAbin <- function (dna){
  DNAbin <- as_DNAbin(dna)
  names(DNAbin) <- names(dna)
  return(DNAbin)
}

read_fasta_df <- function (file = "") {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  tib <- tibble(name = gsub(">", "", fasta[ind]), sequence = seqs)
  return(tib)
}

setwd("/Users/bross/Desktop/AIMS/Analysis/ITS2:psbA_tanglegram")


datapsbA <- read.csv("/Users/bross/Desktop/AIMS/Analysis/NGS/All/feature_table_DIVs.csv")
data2 <- read.csv("/Users/bross/Desktop/AIMS/Analysis/ITS2_49:55_profile/SCF049_SCF055_ITS2.csv")
gal <- read.csv("/Users/bross/Desktop/AIMS/Analysis/ITS2_49:55_profile/galaxea_c1_samples.csv")
psa <- read_excel("post_RIC_data.xlsx")


#writing fasta file of first 5 ASVs

subset_data <- datapsbA[1:20, ]

fasta_file <- "/Users/bross/Desktop/AIMS/Analysis/NGS/All/first20_sequences.fasta"

# Write to FASTA format
writeLines(
  paste0(">", subset_data$ASV_code, "\n", subset_data$X),
  con = fasta_file
)



samples <- c("X","SCF049.1","SCF049.2","SCF049.3","SCF055.1","SCF055.4","SCF055.5","Gal.509", "ASV_code", "Psa.36", "Psa.44", "Psa.109", "Psa.113")

psbA_filt <- datapsbA %>%
  select(any_of(samples))

gal_filtered <- gal %>%
  filter(extraction_id == "509") %>%
  mutate(sample_name = paste0("Gal-",extraction_id ))


names(psa)[names(psa) == "A"] <- "sample_name"
psa <- psa[, !names(psa) %in% "symbiont"]

psa_filter <-  psa[grep("36|44|109|113", psa$sample_name), ]

seq_file <- list.files(pattern = "seqs.absolute.abund_and_meta.txt", recursive = FALSE, full.names = TRUE)
seqs <- read_tsv(seq_file) %>%
  filter(!(is.na(sample_name)))

seqs_long <- seqs %>%
  dplyr::select(sample_name, sample_type, host_genus,	host_species,	collection_latitude, collection_longitude, collection_date, collection_depth:last_col()) %>% # Select sample_names and the each column contain sequence count data
  pivot_longer(-sample_name:-collection_depth) %>% # make into long dataframe
  filter(value > 0) %>% # Remove zero values
  filter(!(is.na(sample_name)))

fasta_file <- list.files(pattern = ".fasta", recursive = FALSE, full.names = TRUE)  
fasta <- read_fasta_df(fasta_file)

seqs_long <- seqs_long %>% left_join(., fasta)

filtered_df <- seqs_long %>% 
  filter(sample_name %in% c("Ca1", "Ca2", "Ca3"))

common_cols <- intersect(names(data2), names(filtered_df))
df1_common <- data2[common_cols]
df2_common <- filtered_df[common_cols]
merged_df <- rbind(df1_common, df2_common)

merged_df <- merged_df %>% 
  filter(!sample_name %in% c("SCF055-C", "SCF049-C", "SCF055-1", "SCF055-4", "SCF055-5"))

merged_df$sample_name <- gsub("Ca1", "SCF055-1", merged_df$sample_name)
merged_df$sample_name <- gsub("Ca2", "SCF055-4", merged_df$sample_name)
merged_df$sample_name <- gsub("Ca3", "SCF055-5", merged_df$sample_name)


df1_subset <- gal_filtered[, names(gal_filtered) %in% names(merged_df)]
final_ITS2 <- rbind(df1_subset,merged_df )
subset2 <- final_ITS2[, names(final_ITS2) %in% names(psa_filter)]
final <- rbind(subset2,psa_filter )

final <- final %>% 
mutate(
  sample_name = str_replace_all(sample_name, "PS2_|PS4_", "Psa-")  # Replace '.' with '-'
)


#plot ITS2 profiles

final2 <- final %>%
  group_by(sample_name) %>%
  mutate(relative_abundance = value / sum(value)) %>%
  ungroup()


set.seed(5584)
P50 <- createPalette(100, c("#1F77B4", "#FFB74D", "#6A1B9A", "#FF4081", "#26A69A", 
                            "#FFEB3B", "#0288D1", "#E64A19", "#2E7D32", "#D32F2F"))
names(P50) <- NULL



final2$name <- factor(
  final2$name,
  levels = final2 %>%
    group_by(name) %>%
    summarise(value = sum(relative_abundance)) %>%
    arrange(desc(value)) %>%
    pull(name)
)

top_divs <- final2 %>%
  group_by(name) %>%
  summarise(value = sum(relative_abundance)) %>%
  arrange(desc(value)) %>%
  slice_head(n = 10) %>%
  pull(name)

ordered_divs <- final2 %>%
  group_by(name) %>%
  summarise(value = sum(relative_abundance)) %>%
  arrange(desc(value)) %>%
  pull(name)

all_DIV <- unique(final2$name)
color_mapping <- setNames(P50[1:length(ordered_divs)], ordered_divs)

sample_order <- c("Gal-509", "SCF049-3", "SCF049-2", "SCF049-1", 
                  "SCF055-4", "SCF055-1", "SCF055-5", 
                  "Psa-109", "Psa-113", "Psa-44", "Psa-36")

final2$sample_name <- factor(final2$sample_name, levels = sample_order)

# Create the plot with the top 10 DIVs in the legend
p1 <- ggplot(final2, aes(x = sample_name, y = relative_abundance, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "ITS2",
       x = "",
       y = "") +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank()   # Remove y-axis title
  ) +
  scale_fill_manual(values = c_pal, breaks = top_divs) +  # Use the top 10 DIVs
  guides(fill = guide_legend(title = "DIVs"))



set.seed(58584)
Pmega <- createPalette(2500,c("#FFB74D", "#6A1B9A", 
                              "#FFEB3B", "#E64A19"))
names(Pmega) <- NULL


feature_table_long <- read.csv("/Users/bross/Desktop/AIMS/Analysis/NGS/all/feature_table_long.csv")

feature_table_long <- feature_table_long %>%
  mutate(ASV_code = ifelse(ASV_code %in% c("ASV9", "ASV10"), "ASV1", ASV_code))

feature_table_long <- feature_table_long %>%
  mutate(ASV_code = ifelse(ASV_code %in% c("ASV11"), "ASV9", ASV_code)) 

feature_table_long <- feature_table_long %>%
  mutate(ASV_code = ifelse(ASV_code %in% c("ASV12"), "ASV10", ASV_code))

all_DIV2 <- unique(feature_table_long$ASV_code)
color_mapping2 <- setNames(Pmega[1:length(all_DIV2)], all_DIV2)

feature_table_long$ASV_code <- factor(
  feature_table_long$ASV_code,
  levels = feature_table_long %>%
    group_by(ASV_code) %>%
    summarise(total_abundance = sum(Relative_Abundance)) %>%
    arrange(desc(total_abundance)) %>%
    pull(ASV_code)
)

feature_table_long <- feature_table_long %>%
  filter(!Sample %in% c("G5-2", 
                     "G1-1",
                     "G1-2",
                     "G2-1",
                     "G2-2",
                     "G3-1",
                     "G3-2",
                     "G4-1",
                     "G4-2",
                     "G5-1")) 

feature_table_long$Sample <- factor(feature_table_long$Sample, levels = sample_order)

# Create the plot
p3 <- ggplot(feature_table_long, aes(x = Sample, y = Relative_Abundance, fill = ASV_code)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "psbA non-coding region",
       x = "",
       y = "") +
  guides(fill = guide_legend(title = "ASVs")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_fill_manual(
    values = color_mapping2,  # Use all colors in the plot
    breaks = unique(feature_table_long$ASV_code)[1:10]  # Limit the legend to the first 10 entries
  )

p <- p1/p3

ggsave("ITS2_psbA_profiles.pdf", plot = p)

ggsave("ITS2_psbA_profiles.svg", plot = p, device = "svg")












names(psbA_filt)[names(psbA_filt) == "X"] <- "sequence"

data1 <- psbA_filt %>%
  pivot_longer(
    cols = -c(sequence, ASV_code), # All columns except `sequence` and `ASV_code`
    names_to = "sample_name",        # New column to hold sample names
    values_to = "value"            # New column for the counts
  )%>%
  mutate(
    sample_name = str_replace_all(sample_name, "\\.", "-")  # Replace '.' with '-'
  )
names(data1)[names(data1) == "ASV_code"] <- "name"


#deframe getting fasta sequences

its2_fasta <- deframe(final %>% 
                        select(name, sequence)) %>%
  as_dna()

# Convert psbA to DNA sequence format
psba_fasta <- deframe(data1 %>% 
                        select(name, sequence)) %>%
  as_dna()


# 3. Calculate the k-distance matrix for both ITS2 and psbA sequences
kdist_its2 <- its2_fasta %>%
  dna_to_DNAbin() %>%  # Convert to binary format
  kdistance(k = 5, residues = "DNA", method = "edgar") %>%  # Calculate k-distance
  as.matrix()  # Convert to matrix format


kdist_psba <- psba_fasta %>%
  dna_to_DNAbin() %>%  # Convert to binary format
  kdistance(k = 5, residues = "DNA", method = "edgar") %>%  # Calculate k-distance
  as.matrix()  # Convert to matrix format


#Convert dataset in wide format, seqs_wide_its2 and seqs_wide_psba contain profiles for each sample, where rows are samples and columns are sequences.
seqs_wide_its2 <- final %>%
  select(sample_name, name, value) %>%  # Keep only relevant columns
  pivot_wider(
    names_from = name,          # Sequences become column names
    values_from = value,        # Values fill the columns
    values_fill = 0             # Missing values are set to 0
  ) %>%
  tibble::column_to_rownames(var = "sample_name")  # Use sample names as rownames


# same for psbA

seqs_wide_psba <- data1 %>%
  select(sample_name, name, value) %>%
  pivot_wider(
    names_from = name, 
    values_from = value, 
    values_fill = 0
  ) %>%
  tibble::column_to_rownames(var = "sample_name")


# Unifrac distances combining k-mer polygenetic distance matrix and sample profiles
#its2
k_tree_its2 <- kdist_its2 %>% phangorn::upgma()  # Build a tree from k-distance
k_unidist_its2 <- GUniFrac(seqs_wide_its2, k_tree_its2)  # Compute UniFrac distances
dist_its2 <- as.dist(k_unidist_its2$unifracs[, , "d_0.5"])  # Use alpha = 0.5

k_tree_psba <- kdist_psba %>% phangorn::upgma()  # Build a tree from k-distance
k_unidist_psba <- GUniFrac(seqs_wide_psba, k_tree_psba)  # Compute UniFrac distances
dist_psba <- as.dist(k_unidist_psba$unifracs[, , "d_0.5"])  # Use alpha = 0.5

hclust_samps_its2 <- upgma(as.matrix(dist_its2))  # Cluster samples
tree_its2 <- ggtree(hclust_samps_its2, size = 0.2) +
  geom_tiplab(angle = 270, size = 0.5) +
  theme(aspect.ratio = 0.3)

hclust_samps_psba <- upgma(as.matrix(dist_psba))  # Cluster samples
tree_psba <- ggtree(hclust_samps_psba, size = 0.2) +
  geom_tiplab(angle = 270, size = 0.5) +
  theme(aspect.ratio = 0.3)





# Convert ITS2 and psbA trees to `phylo` objects
phylo_its2 <- as.phylo(hclust_samps_its2)
phylo_psba <- as.phylo(hclust_samps_psba)


collapse_groups <- list(
  SCF049 = c("SCF049-1", "SCF049-2", "SCF049-3"),
  SCF055 = c("SCF055-1", "SCF055-4", "SCF055-5")
)

# Function to collapse groups in a phylo tree
collapse_samples <- function(tree, groups) {
  for (new_label in names(groups)) {
    samples_to_collapse <- groups[[new_label]]
    valid_samples <- intersect(tree$tip.label, samples_to_collapse)
    
    # Skip if less than 2 valid samples
    if (length(valid_samples) < 2) next
    
    # Find MRCA of the samples
    mrca_node <- getMRCA(tree, valid_samples)
    if (is.na(mrca_node)) {
      warning(paste("Skipping group", new_label, "- MRCA not found"))
      next
    }
    
    # Drop the old tips
    tree <- drop.tip(tree, valid_samples)
    
    # Add the new collapsed tip at the root
    tree <- bind.tip(tree, tip.label = new_label, where = length(tree$tip.label) + 1)
  }
  return(tree)
}


phylo_its2_collapsed <- collapse_samples(phylo_its2, collapse_groups)
phylo_psba_collapsed <- collapse_samples(phylo_psba, collapse_groups)


cophylo_trees <- cophylo(phylo_its2, phylo_psba)

cophylo_trees_collapsed <- cophylo(phylo_its2_collapsed, phylo_psba_collapsed)

pdf("ITS2_psbA_tanglegram.pdf")

plot.cophylo(
  cophylo_trees,
  link.type = "curved",
  link.lwd = 1,
  link.col = "black",
  fsize = 1,            # Font size for tip labels
  col = c("blue", "red"), # Color for ITS2 and psbA trees
  main = "Tanglegram: ITS2 vs psbA"
)

dev.off()

 pdf("ITS2_psbA_tanglegram_collapsed.pdf")

plot.cophylo(
  cophylo_trees_collapsed,
  link.type = "curved",
  link.lwd = 1,
  link.col = "gray",
  fsize = 0.7,            # Font size for tip labels
  col = c("blue", "red"), # Color for ITS2 and psbA trees
  main = "Tanglegram: ITS2 vs psbA"
)

dev.off()

