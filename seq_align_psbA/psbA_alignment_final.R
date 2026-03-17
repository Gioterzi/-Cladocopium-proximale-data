
library(tidyverse)
library(DECIPHER)
library(bioseq)
library(Biostrings)
library(ape)
library(sangerseqR)
library(ggtree)
library(reshape2)
library(readxl)
library(pheatmap)



read_fasta_df <- function (file = "") {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  tib <- tibble(label = gsub(">", "", fasta[ind]), sequence = seqs)
  return(tib)
}

read_fasta_df <- function (file = "") {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 
                                                       1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  tib <- tibble(label = gsub(">", "", fasta[ind]), sequence = seqs)
  return(tib)
}

write_fasta <- function(sequences, names, filename) {
  con <- file(filename, "w")
  for (i in seq_along(sequences)) {
    writeLines(paste0(">", names[i]), con)
    if (!is.character(sequences[[i]])) {
      sequences[[i]] <- as.character(sequences[[i]])
    }
    writeLines(sequences[[i]], con)
  }
  close(con)
}


write_fasta_df <- function (data, filename) 
{
  fastaLines = c()
  for (rowNum in 1:nrow(data)) {
    fastaLines = c(fastaLines, as.character(paste(">", 
                                                  data[rowNum, "label"], sep = "")))
    fastaLines = c(fastaLines, as.character(data[rowNum, 
                                                 "sequence"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

dna_to_DNAbin <- function (dna){
  DNAbin <- as_DNAbin(dna)
  names(DNAbin) <- names(dna)
  return(DNAbin)
}
dna_to_DNAStringset <- function(x) 
{
  bioseq:::check_dna(x)
  DNAstr <- DNAStringSet(paste(x))
  names(DNAstr) <- names(x)
  return(DNAstr)
}

DNAStringSet_to_dna <- function(x){
  x_dna <- as_dna(paste(x))
  names(x_dna) <- names(x)
  res <- tibble(label = names(x), sequence = x_dna)
  return(res)
}

# Convert DNAstringset to DNAbin
DNAStringSet_to_DNAbin <- function(DNAStringSet){
  DNAbin <- as.DNAbin(DNAStringSet)
  return(DNAbin)
}

replace_plus_with_minus <- function(seq) {
  return(chartr("+", "-", seq))
}

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2315-y
palette <- c("A" = "#46ff2d", 
             "G" = "#ffae01", 
             "C" = "#f24641", 
             "T" = "#4294fa", 
             "K" = "#8b4816",
             "M" = "#83831f",
             "R" = "#ffff81",
             "S" = "#ff9d80",
             "Y" = "#e381f2",
             "W" = "#80fff2",
             "V" = "#fde4b8",
             "B" = "#f9c1bf",
             "H" = "#c0d9f9",
             "D" = "#c7ffba",
             "U" = "#8989fb",
             "N" = "black", 
             "-" = "white",
             "+" = "deeppink")


pal_df <- data.frame(names = names(palette), col = palette)

read_ab1_df <- function(file) {
  ab1_data <- read.abif(file)
  if (inherits(ab1_data, "abif")) {
    seq <- as.character(sangerseqR::primarySeq(sangerseqR::sangerseq(ab1_data)))
    return(DNAStringSet(seq, names = basename(file)))
  } else {
    stop("The file does not seem to be a valid AB1 file.")
  }
}


setwd("/Users/bross/Desktop/AIMS/Analysis/Alignments/psbA_ASV_align/5asvs")

########### Align and pull together all sequences from different species based on marker
psba_all <- sort(list.files(path = "/Users/bross/Desktop/AIMS/Analysis/Alignments/old/first_markers_alignment/psbA", pattern = "psbA", full.names = TRUE))


fasta_files <- c("final_sequences_psbA.fasta", "all_markers.fasta")
combined_sequences <- DNAStringSet()
for (file in fasta_files) {
  combined_sequences <- c(combined_sequences, readDNAStringSet(file))
}

for (file in psba_all) {
  combined_sequences <- c(combined_sequences, readDNAStringSet(file))
}

keep_indices <- !grepl("SCF049-3|madreporum|patulum|sodalum|SCF055.01-psbA|goreaui", names(combined_sequences))
filtered_sequences <- combined_sequences[keep_indices]



marker_list <- c("psbA")

aligned_sequences_list <- list()
# Iterate over each marker
for (substring in marker_list) {
  # Initialize a DNAStringSet object to store sequences for the current marker
  marker_sequences <- DNAStringSet()
  
  # Iterate over combined sequences
  for (seq_index in 1:length(filtered_sequences)) {
    # Check if the sequence name contains the current marker
    if (grepl(substring, names(filtered_sequences)[seq_index])) {
      # If it contains the marker, add the sequence to marker_sequences
      marker_sequences <- c(marker_sequences, filtered_sequences[seq_index])
    }
  }
  
  # Perform multiple sequence alignment
  aligned_sequences <- marker_sequences %>%
    OrientNucleotides() %>%
    AlignSeqs() 
  
  # Store aligned sequences in the list
  aligned_sequences_list[[substring]] <- aligned_sequences
  
  # Write the aligned sequences to a file
  writeXStringSet(aligned_sequences, paste0("aligned_sequences_", substring, ".fasta"))
}

#######################
# View the alignments
aligned_df <- data.frame()
for(i in 1:length(aligned_sequences_list)){
  a_pair <- aligned_sequences_list[[i]] %>% writeXStringSet("temp_file.fasta")
  a_df <- read_fasta_df("temp_file.fasta")
  aligned_df <- rbind(aligned_df, a_df)
}

aligned_plotting <- aligned_df %>%
  mutate(sample_id = str_replace(label, "(psbA)", ""),
         label = str_extract(label, "(psbA)"))



name_map <- c("SCF049.01-" = "SCF049",
              "proliferum-055-" = "proliferum-055", 
              "ASV1-" =  "ASV1",
              "ASV2-" =  "ASV2", 
              "ASV3-" =  "ASV3", 
              "ASV4-" =  "ASV4", 
              "ASV5-" =  "ASV5", 
              "SCF049.01-newcons-" = "SCF049 - Sanger Internal",
              "SCF055.01-newcons-" = "proliferum - Sanger Internal",
             "proliferum_A02_71_C1_Galaxea_fasicularis-" = "proliferum_W4",
             "proliferum_A03_50_C1_Acropora_tenuis-" = "proliferum_W3",
             "proliferum_JPB08_77_C1_Siderastrea-" = "proliferum_W2",
             "proliferum_Zam03_70_C1_Astropora- " = "proliferum_W1", 
             "vulgare_A03_392_C1_Merulina-" = "vulgare_W6", 
             "vulgare_A03_311_C1_Turbinaria-" = "vulgare_W5",
             "vulgare_Pho09_239_C1_Acropora-" = "vulgare_W4", 
             "vulgare_A03_395_C1_Merulina-" = "vulgare_W3",
             "vulgare_Palau_offshore_C1_Pachyser_col3-" = "vulgare_W2", 
             "vulgare_NewCal16_51_C1_Pectinia-" = "vulgare_W1",
             "pacificum_VR19_LJ31-" = "pacificum_W4", 
             "pacificum_VR19_LJ05-" = "pacificum_W3",
             "pacificum_VR19_LJ13-" = "pacificum_W2", 
             "pacificum_VR19_LJ20-" = "pacificum_W1", 
             "latusorum_Haw09_014-" = "latusorum_W1",
             "latusorum_TW_117-" = "latusorum_W4",
             "latusorum_SML07_60-" = "latusorum_W2",
             "latusorum_Haw09_012-" = "latusorum_W3")

aligned_plotting$sample_id <- name_map[aligned_plotting$sample_id]

# Create profile key
key <- aligned_plotting %>%
  tibble::rownames_to_column(var = "id")


# Create long dataframe for ggplot
long_sequences <- str_split(aligned_plotting$sequence, "") %>%
  reshape2::melt() %>%
  group_by(L1) %>%
  mutate(x = row_number(),
         L1 = as.character(L1)) %>%
  left_join(., key, by = c("L1" = "id")) %>%
  ungroup()

long_sequences <- long_sequences %>%
  mutate(sample_id = factor(sample_id, levels = unique(key$sample_id)))


# Plot alignment
p1 <- ggplot(long_sequences, aes(y = sample_id, x = x)) +
  geom_tile(aes(fill = value), size = 1, name = "base") +
  facet_wrap(~ label, nrow = 3, scales = "free") +
  geom_vline(xintercept = seq(50, max(long_sequences$x), by = 50), color = "black", linetype = "dashed", size = 0.05) +
  scale_fill_manual(values = palette) +
  theme(aspect.ratio = 0.5,
        axis.title.y = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 10, 0)), # Add space below strip text
        strip.background = element_blank(), # Remove strip background
        strip.placement = "outside",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 8),
        axis.ticks.width = unit(0.1, "cm"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(min(long_sequences$x), max(long_sequences$x), by = 100)) +
  xlab("Position")
p1

ggsave("psbA_alignment_Cladocopium.pdf", plot = p1, device="pdf", width = 10, height = 4)


####### Concatenate based on species


fasta_files_prex <- c("aligned_sequences_psbA.fasta")

filtered_df <- fasta_files_prex %>%
  map(read_fasta_df) %>% # read in all the files
  purrr::reduce(rbind) # reduce with rbind into one dataframe

final_sequences <- filtered_df %>%
  mutate(sequence = case_when(str_detect(label, "psbA") ~ str_sub(sequence, start = 942, end = 1290)))

#psbA starting at 989 for gab and ambigous removal

aligned_plotting <- final_sequences %>%
  mutate(sample_id = str_replace(label, "(psbA)", ""),
         label = str_extract(label, "(psbA)"))


aligned_plotting$sample_id <- name_map[aligned_plotting$sample_id]
# Create profile key
key <- aligned_plotting %>%
  tibble::rownames_to_column(var = "id")


# Create long dataframe for ggplot
long_sequences <- str_split(aligned_plotting$sequence, "") %>%
  reshape2::melt() %>%
  group_by(L1) %>%
  mutate(x = row_number(),
         L1 = as.character(L1)) %>%
  left_join(., key, by = c("L1" = "id")) %>%
  ungroup()

long_sequences <- long_sequences %>%
  mutate(sample_id = factor(sample_id, levels = unique(key$sample_id)))

long_sequences$x <- long_sequences$x + 953

# long_sequences <- long_sequences %>%
#    filter(!str_detect(sample_id, "goreaui"))

# Plot alignment
p1 <- ggplot(long_sequences, aes(y = sample_id, x = x)) +
  geom_tile(aes(fill = value), size = 1, name = "base") +
  facet_wrap(~ label, nrow = 3, scales = "free") +
  geom_vline(xintercept = seq(50, max(long_sequences$x), by = 50), color = "black", linetype = "dashed", size = 0.05) +
  scale_fill_manual(values = palette) +
  theme(aspect.ratio = 0.5,
        axis.title.y = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 10, 0)), # Add space below strip text
        strip.background = element_blank(), # Remove strip background
        strip.placement = "outside",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 8),
        axis.ticks.width = unit(0.1, "cm"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(min(long_sequences$x), max(long_sequences$x), by = 100)) +
  xlab("Position") + coord_cartesian(xlim = c(953, max(long_sequences$x)))
p1

ggsave("psbA_alignment_Cladocopium_trimmed.svg", p1, width = 16, height = 20, dpi = 300)

  ggsave("psbA_alignment_Cladocopium_trimmed.pdf", plot = p1, device="pdf", width = 10, height = 4)


##########################################################################
#Including psbA, dropping Justin's SCF049 sequences, removed madreporum, patulum and sodalum
##########################################################################

species <- c("SCF049.01","proliferum-055", "ASV1", "ASV2", "ASV3", "ASV4", "ASV5","SCF049.01-newcons","SCF055.01-newcons",
             "proliferum_A02","proliferum_A03","proliferum_JPB08","proliferum_Zam03", 
             "vulgare_A03_392", "vulgare_A03_311", "vulgare_Pho09_239", "vulgare_A03_395", "vulgare_Palau", "vulgare_NewCal",
             "pacificum_VR19_LJ31", "pacificum_VR19_LJ05", "pacificum_VR19_LJ13", "pacificum_VR19_LJ20", 
             "latusorum_Haw09_014", "latusorum_TW_117", "latusorum_SML07_60", "latusorum_Haw09_012")
pattern <- paste(species, collapse = "|")

keep_indices <- grepl(pattern, final_sequences$label)
final_sequences <- final_sequences[keep_indices, ]


trimmed_fasta_filenames <- c("trimmed_psbA.fasta")

for (i in seq_along(trimmed_fasta_filenames)) {
  # Filter dataframe for sequences with corresponding labels
  trimmed_sequences <- final_sequences[grepl(gsub("trimmed_", "", gsub(".fasta", "", trimmed_fasta_filenames[i])), final_sequences$label), ]
  # Write sequences to trimmed FASTA file using your custom function
  write_fasta(trimmed_sequences$sequence, trimmed_sequences$label, trimmed_fasta_filenames[i])
}


compute_sequence_lengths <- function(fasta_file) {
  # Read the FASTA file
  fasta_sequences <- readDNAStringSet(fasta_file, format="fasta")
  
  # Compute lengths of sequences
  sequence_lengths <- width(fasta_sequences)
  
  # Return a named vector of sequence lengths
  return(sequence_lengths)
}

# Apply the function to each file and store the results in a list
sequence_lengths_list <- lapply(trimmed_fasta_filenames, compute_sequence_lengths)

# Display sequence lengths for each file
names(sequence_lengths_list) <- trimmed_fasta_filenames
sequence_lengths_list


# Write trimmed sequences to respective trimmed FASTA files


#fasta_files_prex2 <- c("trimmed_Cox1.fasta", "trimmed_LSU.fasta", "trimmed_cp23s.fasta", "trimmed_Cob.fasta") #, "aligned_sequences_psbA.fasta")
fasta_files_2 <- lapply(trimmed_fasta_filenames, readDNAStringSet)
# Define substrings to search for


concatenate_sequences <- function(species_name, fasta_list) {
  concatenated_seq <- DNAString("")
  for (fasta in fasta_list) {
    # Find the sequence for the given species in the current fasta file
    seq_index <- which(sapply(names(fasta), function(x) grepl(species_name, x)))
    if (length(seq_index) == 1) {
      concatenated_seq <- c(concatenated_seq, fasta[[seq_index]])
    } else {
      stop(paste("Sequence for species", species_name, "not found or duplicated in one of the FASTA files"))
    }
  }
  concatenated_seq
}

# Loop over each species and generate concatenated sequences
for (sp in species) {
  cat("Processing species:", sp, "\n")
  concat_seq <- concatenate_sequences(sp, fasta_files_2)
  
  # Create a DNAStringSet and set the name to the species name
  concat_seq_set <- DNAStringSet(concat_seq)
  names(concat_seq_set) <- sp
  
  # Write the concatenated sequence to a new FASTA file
  writeXStringSet(concat_seq_set, paste0(sp, "_concatenated_3.fasta"))
}

all_sequences <- list.files(pattern="concatenated_3.fasta", full.names = TRUE)
fasta_contents <- lapply(all_sequences, readLines)
all_sequences_content <- unlist(fasta_contents)
writeLines(all_sequences_content, "all_sequences_psbA.fasta")


clean <- readDNAStringSet("all_sequences_psbA.fasta")

clean <- lapply(clean, replace_plus_with_minus)
clean <- DNAStringSet(clean)
writeXStringSet(clean, "tree_psbA.fasta")


tree_sequence <- read.dna("tree_psbA.fasta", format = "fasta")
distance_matrix <- dist.dna(tree_sequence, model = "raw", as.matrix = TRUE) #calculate distance matrix
pairwise_differences <- dist.dna(tree_sequence, model = "raw")

nj_tree <- nj(pairwise_differences)
rooted_tree <- root(nj_tree, outgroup = "madreporum", resolve.root = TRUE)

pdf("tree_rooted_psbA.pdf", width = 30, heigh = 10)
plot(rooted_tree )
dev.off()

#add.scale.bar(length = mean(pairwise_differences), lwd = 1, col = "black", cex = 0.8)


#pairwise base substitution
num_sequences <- nrow(tree_sequence)
sequence_length <- ncol(tree_sequence)


num_diff_matrix <- matrix(0, nrow = num_sequences, ncol = num_sequences)
substitution_positions <- rep(FALSE, sequence_length)
labels <- rownames(tree_sequence)

for (i in 1:(num_sequences - 1)) {
  for (j in (i+1):num_sequences) {
    differences <- tree_sequence[i, ] != tree_sequence[j, ]
    num_diff <- sum(differences)
    num_diff_matrix[i, j] <- num_diff
    num_diff_matrix[j, i] <- num_diff  # Matrix is symmetric
    substitution_positions <- substitution_positions | differences
  }
}

rownames(num_diff_matrix) <- labels
colnames(num_diff_matrix) <- labels

nj_tree2 <- nj(num_diff_matrix)
plot(nj_tree2)
rooted_tree2 <- root(nj_tree2, outgroup = "madreporum", resolve.root = TRUE)

pdf("tree_rooted_psbA_extended.pdf", width = 30, heigh = 10)
plot(rooted_tree2)
dev.off()
heat <- pheatmap(num_diff_matrix, 
                 cluster_rows = FALSE, 
                 cluster_cols = FALSE, 
                 display_numbers = TRUE,  # Optional: to display the distance values
                 color = colorRampPalette(c("white", "blue", "red"))(100),  # Custom color scale
                 main = "Distance Heatmap",
                 number_color = "green",
                 fontsize_number = 10,  # Optional: font size for numbers
                 cellwidth = 50,  # Optional: cell width
                 cellheight = 50) 

pdf("distance_heatmap_psbA.pdf", width = 10, height = 10)
heat
dev.off()


############# only psbA tree

tree_sequence_psbA <- read.dna("trimmed_psbA.fasta", format = "fasta")
distance_matrix <- dist.dna(tree_sequence_psbA, model = "raw", as.matrix = TRUE) #calculate distance matrix
pairwise_differences <- dist.dna(tree_sequence_psbA, model = "raw")

nj_tree <- nj(pairwise_differences)
rooted_tree <- root(nj_tree, outgroup = "madreporum-psbA", resolve.root = TRUE)

pdf("tree_rooted_psbA_only.pdf", width = 30, heigh = 10)
plot(rooted_tree )
dev.off()





##########################
#plotting substitutions

substitution_indices <- which(substitution_positions)

# Create a matrix to store the bases at substitution positions
base_matrix <- tree_sequence[, substitution_indices]
colnames(base_matrix) <- substitution_indices
rownames(base_matrix) <- rownames(tree_sequence)

spe <- c("madreporum", "patulum", "proliferum", "SCF049-1","SCF049-2","SCF049-J", "SCF055-1","SCF055-4","SCF055-5", "sodalum", "vulgare")

# Convert the matrix to a data frame for easier plotting with ggplot2
base_df <- as.data.frame(t(base_matrix)) %>% 
  head(., 18)
rownames(base_df) <- substitution_indices
base_df$Position <- as.numeric(substitution_indices)


#export base_df in excel, manipulate it and reimport it
base_sub <- read_excel("/Users/bross/Desktop/AIMS/Analysis/Alignments/markers_alignment/seq_substitutions_position.xlsx", sheet = "Foglio9")

base_sub_long <- base_sub %>%
  pivot_longer(cols = -Position, names_to = "Species", values_to = "Base")

base_sub_long <- base_sub_long %>%
  mutate(Marker = case_when(
    Position >= 0 & Position <= 470 ~ "LSU",
    Position >= 470 & Position <= 987 ~ "cp23s",
    Position >= 987 & Position <= 1823 ~ "Cob",
    Position >= 1823 & Position <= 2724 ~ "Cox1",
    TRUE ~ "Unknown"  # In case any position falls outside the given ranges
  ))
base_sub_long$Base <- toupper(base_sub_long$Base)

marker_ranges <- tribble(
  ~Marker, ~start, ~end,
  "LSU", 0, 470,
  "cp23s", 470, 987,
  "Cob", 987, 1823,
  "Cox1", 1823, 2724)

p3 <- ggplot(base_sub_long, aes(x = Position, y = Species)) +
  geom_tile(aes(fill = Base), size = 4) +
  #geom_point(aes(color = base_sub_long$Base), size = 2) +
  scale_fill_manual(values = palette, name = "Base")+
  theme_minimal() +
  labs(title = "",
       x = "Position",
       y = "Species",
       color = "Base") +
  geom_segment(data = marker_ranges, aes(x = start, xend = end, y = Inf, yend = Inf, color = Marker), 
               size = 1) +
  geom_text(data = marker_ranges, aes(x = (start + end) / 2, y = Inf, label = Marker, color = Marker), 
            vjust = -0.5) +
  scale_color_manual(values = c("Cox1" = "red", "LSU" = "blue", "cp23s" = "green", 
                                "Cob" = "purple", "psbA" = "orange", "Unknown" = "gray")) +
  geom_hline(yintercept = c("SCF055-1","SCF055-4","SCF055-5", "SCF049-1","SCF049-2","SCF049-3", "proliferum"), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(0, max(base_sub_long$Position)), clip = 'off') +  # Set x-axis limits explicitly
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))  # Add margin to avoid clipping

ggsave("substitution_position.pdf", plot = p3, device="pdf", width = 10, height = 4)

##########################


       #psbA - now included in the all.inclusive code up here


##########################
psba_all <- sort(list.files(path = "/Users/bross/Desktop/AIMS/Analysis/markers_alignment/psbA", pattern = "psbA", full.names = TRUE))

psba_all_fasta <- psba_all %>%
  map(read_fasta_df) %>% # read in all the files
  purrr::reduce(rbind) # reduce with rbind into one dataframe

psba_strings <- psba_all_fasta %>%
  deframe() %>%
  as_dna() %>%
  dna_to_DNAStringset()

alligned_psba <- psba_strings %>%
  OrientNucleotides() %>%
  AlignSeqs()

writeXStringSet(alligned_psba, paste0("aligned_sequences_psbA.fasta"))

alligned_df <- read_fasta_df("aligned_sequences_psbA.fasta")

key <- alligned_df %>%
  tibble::rownames_to_column(var = "id")

long_sequences <- str_split(key$sequence, "") %>%
  reshape2::melt() %>%
  group_by(L1) %>%
  mutate(x = row_number(),
         L1 = as.character(L1)) %>%
  left_join(., key, by = c("L1" = "id")) %>%
  ungroup()

# Plot alignment
p2 <- ggplot(long_sequences, aes(y = label, x = x)) +
  geom_tile(aes(fill = value), size = 1, name = "base") +
  geom_vline(xintercept = seq(50, max(long_sequences$x), by = 50), color = "black", linetype = "dashed", size = 0.05) +
  scale_fill_manual(values = palette) +
  theme(aspect.ratio = 0.3,
        axis.title.y = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 10, 0)), # Add space below strip text
        strip.background = element_blank(), # Remove strip background
        strip.placement = "outside",
        axis.text.x = element_text(size = 0.5),
        axis.ticks.width = unit(0.1, "cm"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(min(long_sequences$x), max(long_sequences$x), by = 100)) +
  xlab("Position")
p2

trim <- alligned_df %>%
  mutate(sequence = case_when(str_detect(label, "psbA") ~ str_sub(sequence, start = 900, end = 1150),
                              TRUE ~ str_sub(sequence, start = 30, end = 600)))

key <- trim %>%
  tibble::rownames_to_column(var = "id")

long_sequences <- str_split(key$sequence, "") %>%
  reshape2::melt() %>%
  group_by(L1) %>%
  mutate(x = row_number(),
         L1 = as.character(L1)) %>%
  left_join(., key, by = c("L1" = "id")) %>%
  ungroup()

p2 <- ggplot(long_sequences, aes(y = label, x = x)) +
  geom_tile(aes(fill = value), size = 1, name = "base") +
  geom_vline(xintercept = seq(50, max(long_sequences$x), by = 50), color = "black", linetype = "dashed", size = 0.05) +
  scale_fill_manual(values = palette) +
  theme(aspect.ratio = 0.3,
        axis.title.y = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 10, 0)), # Add space below strip text
        strip.background = element_blank(), # Remove strip background
        strip.placement = "outside",
        axis.text.x = element_text(size = 0.5),
        axis.ticks.width = unit(0.1, "cm"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(min(long_sequences$x), max(long_sequences$x), by = 100)) +
  xlab("Position")
p2


ggsave("psbA_alignment_all_Cladocopium.pdf", plot = p2, device="pdf", width = 10, height = 4)

tree <- read.tree(text = "((madreporum-psbA:0.033327,patulum-psbA:0.019071):0.041766,((proliferum-psbA:0.000001,SCF049-psbA:0.022572):0.268198,vulgare-psbA:0.028671):0.032945,sodalum-psbA:0.030775);")
plot(tree)




