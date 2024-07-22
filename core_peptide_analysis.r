########################################
## Packages - don't package shame people
########################################

.libPaths()
.libPaths( c( "/data/san/data0/users/david/rstudio/packages" , .libPaths() ) )
newlib <- "/data/san/data0/users/david/rstudio/packages"
.libPaths()
library(janitor)
library(thacklr)
library(gggenomes)
library(dplyr)
library(remotes)
library(data.table)
library(stringr)
library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(taxonomizr)
library(tidyr)
library(ggsci)
library(Peptides)


########################################
## TREE of genomes / Taxonomy of genomes
########################################


tax_df = fread("/data/san/data1/users/david/mining_2023/class_IId_analysis/data/genomic_phylogenomic_tree/classify/classify/gtdbtk.bac120.summary.tsv") %>%
  clean_names()

assembly_tax_df = tax_df %>%
  as.data.frame() %>%
  separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(assembly = str_sub(user_genome, start = 1L, end = 15L)) 

tree_annot_df = tax_df %>% 
  column_to_rownames(var = "user_genome") %>%
  as.data.frame() %>%
  separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  select(order)

tree_annot_df$family = as.factor(tree_annot_df$family)

family_df <- tax_df %>%
  as.data.frame() %>%
  separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "merge") %>%
  select(user_genome, order, phylum, class)

order_colors <- c(pal_npg("nrc", alpha = 0.2)(9), pal_npg("nrc", alpha = 1)(9))
names(order_colors) <- unique(family_df$order)
family_df$color <- order_colors[family_df$order]

genome_tree = read.tree("/data/san/data1/users/david/mining_2023/class_IId_analysis/data/genomic_phylogenomic_tree/infer/intermediate_results/gtdbtk.unrooted.tree")
genome_tree = ape::root(genome_tree, outgroup = "GCF_018618955.1_ASM1861895v1")

genome_lead <- ggtree(genome_tree, size=0.4, layout = "circular") +
  theme_tree() +
  geom_hilight(node=704, # p__Actinomycetota
    fill="#F39B7FB2",
    alpha=0.5) +
  geom_hilight(node=588, # p__Bacillota
    fill="#91D1C2B2",
    alpha=0.5)

genome_lead$layers <- rev(genome_lead$layers)

circ2 <- genome_lead +
  geom_tippoint(pch=16, size=1, color="black") +
  geom_text(x=-0.8, y=-15, label="Bacillota", size = 6) +
  geom_text(x=2, y=20, label="Actinomycetota", size = 6) +
  theme(
    plot.margin = unit(c(1, 5, 1, 1), "cm")  # Top, Right, Bottom, Left margins
  )

ggsave(circ2, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/leaderless_genomes_tree.png",
	width = 20, height = 20, units = "cm", dpi=600)


########################################
## Top familes, genus and species with leaderless bacteriocins
########################################


tax_tables = tax_df %>% 
  as.data.frame() %>%
  separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  select(user_genome, phylum, family, genus, species)

top_10_family = tax_tables %>% 
  group_by(family) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(type = "Family") %>%
  dplyr::rename(id = family) %>%
  head(10)

top_10_genus = tax_tables %>%
  group_by(genus) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(type = "Genus") %>%
  dplyr::rename(id = genus) %>%
  head(10) 

top_10_species = tax_tables %>%
  group_by(species) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(type = "Species") %>%
  dplyr::rename(id = species) %>%
  head(10)

df_for_tax_bars = bind_rows(top_10_species, top_10_genus, top_10_family)
df_for_tax_bars$id <- gsub(".__", "", df_for_tax_bars$id)


tax_bar_plot = ggplot(df_for_tax_bars, aes(x = reorder(id, -count), y = count, fill = type)) +
  geom_bar(stat = "identity", color = "black", size = 1) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = NULL, y = "Count", title = element_blank()) +
  facet_wrap(~type, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "bold"),
    legend.position = "none")
ggsave(tax_bar_plot, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/tax_bar_plot.png",
  width = 20, height = 10, units = "cm", dpi=600)


# get what % of the total each phylum is
tax_tables %>%
  select(phylum) %>%
  group_by(phylum) %>%
  summarise(count = n()) %>%
  mutate(percent = count / sum(count) * 100) %>%
  arrange(desc(count)) %>%
  head(10)


########################################
## Gene clusters of arachnia and arcanobacterium (synthesized peptides)
########################################

bakta_dir = "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/mgnify/bakta"
bakta_files = list.files(bakta_dir, pattern = "*.tsv", recursive = TRUE, full.names = TRUE)
bakta_files = grep("hypo", bakta_files, value = TRUE, invert = TRUE)
bakta_df = do.call(rbind, lapply(bakta_files, fread, fill = TRUE, skip=5, sep="\t")) %>%
  janitor::clean_names() 

bakta_df$strand <- case_when(
  bakta_df$strand == "+" ~ "+",
  bakta_df$strand == "-" ~ "-",
  bakta_df$strand == "?" ~ NA,
  TRUE ~ bakta_df$strand
  )

bakta_df$type <- case_when(
  bakta_df$type == "cds" ~ "CDS",
  TRUE ~ bakta_df$type
  )

# focus on core
core = c("PKHINC_08475", "JMKLCA_05205")

gg_df = bakta_df %>%
  dplyr::rename(feat_id = locus_tag) %>%
  mutate(seq_id = number_sequence_id,
  	end = stop) %>%
  mutate(focus_point = case_when(feat_id %in% core ~ "focus",
                     TRUE ~ NA)) %>%
  mutate(color = case_when(
    feat_id == "JMKLCA_05195" ~ "#00A087B2", # membrane transpoer
    feat_id == "JMKLCA_05205" ~ "black",
    feat_id == "PKHINC_08475" ~ "black",
    feat_id == "JMKLCA_05210" ~ "#E64B35FF", # membrane transpoer
    feat_id == "JMKLCA_05215" ~ "#E64B35FF", # also membrane transport 
    feat_id == "JMKLCA_05220" ~ "#E64B35FF", # ATP ase transmembrane
    feat_id == "JMKLCA_05225" ~ "#8491B4FF", # bacteriocin immunity
    feat_id == "JMKLCA_05230" ~ "#8491B4FF",  # PH domain-containing protein
    feat_id == "JMKLCA_05235" ~ "#F39B7FB2", # nucleotide excision repair (NER)
    feat_id == "JMKLCA_05240" ~ "gray", # UDP-glucose 4-epimerase
    feat_id == "JMKLCA_05245" ~ "gray", #	type I glutamate--ammonia ligase 
    feat_id == "PKHINC_08480" ~ "gray",
    feat_id == "PKHINC_08485" ~ "#00A087B2", # response regulator
    feat_id == "PKHINC_08490" ~ "#E64B35FF", # abs transporter permease
    feat_id == "PKHINC_08495" ~ "#E64B35FF",  # ATP 
    feat_id == "PKHINC_08500" ~ "gray", 
    feat_id == "PKHINC_08505" ~ "gray", 
    feat_id == "PKHINC_08510" ~ "gray", # alpha/beta hydrolase-fold protein
    feat_id == "PKHINC_08515" ~ "gray", # VanZ protein
    feat_id == "PKHINC_08520" ~ "gray",
    TRUE ~ "grey90")
  ) %>%
  mutate(seq_id = case_when(number_sequence_id == "MGYG000000642_8" ~ "Arcanobacterium",
    number_sequence_id == "MGYG000298963_42-length-NA-cov-NA" ~ "Arachnia",
    TRUE ~ number_sequence_id))


plot1 = gggenomes(genes = gg_df) %>%
  focus(focus_point == "focus", 
    .expand = c(5e3, 9e3),
    .locus_id = str_glue("{seq_id}")) +
  geom_seq() +
  geom_gene(aes(fill = color),
        size = 3) +
  scale_fill_identity() +
  geom_seq_label(size = 1.5, vjust = 2, fontface = "bold.italic", family = "Arial") 

ggsave(plot1, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/bakta_gene_clusters.png",
  width = 8, height = 4, units = "cm", dpi=600)



########################################
## Core peptides properties
########################################


a58_core = readAAStringSet("/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/final_proteins.faa")
a58_core_df = as.data.frame(a58_core)

# initialise properties
a58_core_df$hydrophobicity <- NULL
a58_core_df$length <- NULL
a58_core_df$charge <- NULL
a58_core_df$iso <- NULL
a58_core_df$mw <- NULL
a58_core_df$aromatic <- NULL

# loop over each row in the dataframe and calculate properties for each sequence
for(i in 1:nrow(a58_core_df)) {
  a58_core_df$hydrophobicity[i] <- hydrophobicity(a58_core_df$x[i], scale = "KyteDoolittle")
  a58_core_df$length[i] <- nchar(a58_core_df$x[i])
  a58_core_df$charge[i] <- charge(a58_core_df$x[i], pKscale = "EMBOSS")
  a58_core_df$iso[i] <- pI(a58_core_df$x[i])
  a58_core_df$mw[i] <- mw(a58_core_df$x[i])
  # a58_core_df$aromatic[i] <- aromaticAA(a58_core_df$x[i])
}

fwrite(a58_core_df, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/a58_core_peptide_properties.csv",
  row.names = TRUE)

hydro_plot <- a58_core_df %>%
  arrange(desc(hydrophobicity)) %>%
  ggplot(aes(x = forcats::fct_inorder(x), y = hydrophobicity)) +
  geom_point(shape=21, color = "black", fill="white", alpha=0.1, size = 2) +
  geom_point(data = filter(a58_core_df, grepl("MGYG000298963", row.names(a58_core_df))),
     aes(fill = "#E64B35B2"), size = 2 , alpha = 1, shape=21, color = "black") +
  geom_point(data = filter(a58_core_df, grepl("MGYG000000642_01040", row.names(a58_core_df))),  
     aes(fill = "#4DBBD5B2"), size = 2, alpha = 1, shape=21, color = "black") +
  scale_color_identity() +  
  theme_classic() +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),  # Remove ticks on x-axis
    legend.position = "none") +  # Remove legend for colors
  xlab("") + 
  coord_cartesian(clip = "off")

length_plot <- a58_core_df %>%
  arrange(desc(length)) %>%
  ggplot(aes(x = forcats::fct_inorder(x), y = length)) +
geom_point(shape=21, color = "black", fill="white", alpha=0.1, size = 2) +
  geom_point(data = filter(a58_core_df, grepl("MGYG000298963", row.names(a58_core_df))),
             aes(fill = "#E64B35B2"), size = 2 , alpha = 1, shape=21, color = "black") +
  geom_point(data = filter(a58_core_df, grepl("MGYG000000642_01040", row.names(a58_core_df))),  
             aes(fill = "#4DBBD5B2"), size = 2, alpha = 1, shape=21, color = "black") +
  scale_color_identity() +  # Use actual color names
  theme_classic()+
  xlab("") +
    theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none") +  # Remove legend for colors
  ylab("Length") + 
  coord_cartesian(clip = "off")

charge_plot <- a58_core_df %>%
  arrange(desc(charge)) %>%
  ggplot(aes(x = forcats::fct_inorder(x), y = charge)) +
geom_point(shape=21, color = "black", fill="white", alpha=0.1, size = 2) +
  geom_point(data = filter(a58_core_df, grepl("MGYG000298963", row.names(a58_core_df))),
             aes(fill = "#E64B35B2"), size = 2 , alpha = 1, shape=21, color = "black") +
  geom_point(data = filter(a58_core_df, grepl("MGYG000000642_01040", row.names(a58_core_df))),  
             aes(fill = "#4DBBD5B2"), size = 2, alpha = 1, shape=21, color = "black") +
  scale_color_identity() +  # Use actual color names
  theme_classic()+
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none") +  # Remove legend for colors
  xlab("") +
  ylab("Charge") + 
  coord_cartesian(clip = "off")

iso_plot <- a58_core_df %>%
  arrange(desc(iso)) %>%
  ggplot(aes(x = forcats::fct_inorder(x), y = iso)) +
  geom_point(shape=21, color = "black", fill="white", alpha=0.1, size = 2) +
  geom_point(data = filter(a58_core_df, grepl("MGYG000298963", row.names(a58_core_df))),
             aes(fill = "#E64B35B2"), size = 2 , alpha = 1, shape=21, color = "black") +
  geom_point(data = filter(a58_core_df, grepl("MGYG000000642_01040", row.names(a58_core_df))),  
             aes(fill = "#4DBBD5B2"), size = 2, alpha = 1, shape=21, color = "black") +
  scale_color_identity() +  # Use actual color names
  theme_classic()+
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none") +  # Remove legend for colors
  xlab("") +
  ylab("pI") + 
  coord_cartesian(clip = "off")

mw_plot <- a58_core_df %>%
  arrange(desc(mw)) %>%
  ggplot(aes(x = forcats::fct_inorder(x), y = mw)) +
  geom_point(shape=21, color = "black", fill="white", alpha=0.1, size = 2) +
  geom_point(data = filter(a58_core_df, grepl("MGYG000298963", row.names(a58_core_df))),
             aes(fill = "#E64B35B2"), size = 2 , alpha = 1, shape=21, color = "black") +
  geom_point(data = filter(a58_core_df, grepl("MGYG000000642_01040", row.names(a58_core_df))),  
             aes(fill = "#4DBBD5B2"), size = 2, alpha = 1, shape=21, color = "black") +
  scale_color_identity() +  # Use actual color names
  theme_classic() +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none") +  # Remove legend for colors
  xlab("") +
  ylab("MW") + 
  coord_cartesian(clip = "off")

ggsave(hydro_plot, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/hydrophobicity.png", 
  width = 5, height = 5, units = "cm", dpi=600)
ggsave(length_plot, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/length.png", 
  width = 5, height = 5, units = "cm", dpi=600)
ggsave(charge_plot, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/charge.png",
  width = 5, height = 5, units = "cm", dpi=600)


library(gridExtra)

plot_f1_multi = arrangeGrob(
  hydro_plot,
  charge_plot, 
  mw_plot,
  iso_plot,
  length_plot,
  ncol = 3, nrow = 2)
ggsave(plot_f1_multi, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/figure_1_peptide_properties.png", 
  width = 14, height = 8, units = "cm", dpi = 600)





core_long <- a58_core_df %>%
  rownames_to_column(var = "name") %>%
  pivot_longer(cols = c(hydrophobicity, charge, iso, mw), names_to = "property", values_to = "value")


# Plot
pivot_plot <- core_long %>%
  ggplot(aes(x = length, y = value)) +
  geom_point(shape=21, color = "black", fill=NA, alpha=0.1, size = 2) +
  geom_point(data = filter(core_long, grepl("MGYG000298963", name)),
     aes(fill = "#E64B35B2"), size = 2 , alpha = 1, shape=21, color = "black") +
  geom_point(data = filter(core_long, grepl("MGYG000000642_01040", name)),  
     aes(fill = "#4DBBD5B2"), size = 2, alpha = 1, shape=21, color = "black") +
  scale_color_identity() +  # Use actual color names
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),  # Remove background ticks
    legend.position = "none",
    text = element_text(face = "bold", family = "Arial")) +  # Make text bold and Arial
  xlab("length") +
  ylab("") + 
  coord_cartesian(clip = "off") +
  facet_wrap(~property, ncol = 2, scales = "free")

ggsave(pivot_plot, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/figure_1_pivot_plot_peptide_properties.png", 
  width = 10, height = 10, units = "cm", dpi = 600)


########################################
## Split by Phylum
########################################

# make rownmame as column and remove everything after first white space
core_peptide_df = a58_core_df %>%
  rownames_to_column("id") %>%
  # remove everything after first whitespace 
  mutate(id = sub("\\s.*", "", id))

# Set the directory
bakta_dir_2 <- "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/bakta_out"

# List files, ensuring hypothetical files are excluded
bakta_files_2 <- list.files(bakta_dir_2, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
bakta_files_2 <- grep("hypo", bakta_files_2, value = TRUE, invert = TRUE)

# Read files and combine into a single dataframe, adding a 'file_name' column to track the file origin
bakta_df_2 <- rbindlist(lapply(bakta_files_2, function(x) {
  df <- fread(x, fill = TRUE, skip = 5, sep = "\t")
  # Remove ".tsv" from the basename to get the complete file name
  df[, assembly := str_remove(basename(x), "\\.tsv$")]  # This keeps the full basename including the version number
  return(df)
}), use.names = TRUE, fill = TRUE)

# Clean column names
bakta_df_2 <- clean_names(bakta_df_2)


bakta_df_2$strand <- case_when(
  bakta_df_2$strand == "+" ~ "+",
  bakta_df_2$strand == "-" ~ "-",
  bakta_df_2$strand == "?" ~ NA,
  TRUE ~ bakta_df_2$strand
  )

bakta_df_2$type <- case_when(
  bakta_df_2$type == "cds" ~ "CDS",
  TRUE ~ bakta_df_2$type
  )

peptide_phylum_properties = left_join(core_peptide_df, bakta_df_2, by = c("id" = "locus_tag")) %>%
  left_join(., assembly_tax_df, by = c("assembly" = "assembly")) %>%
  select(id, x, hydrophobicity, length, charge, iso, mw, assembly, phylum, class, order, family, genus, species) %>%
  distinct(x, length, hydrophobicity, species, .keep_all = TRUE)

peptide_phylum_properties <- peptide_phylum_properties %>%
  group_by(phylum) %>%
  filter(n() > 30)  # Adjust this number based on your data distribution

# Reshape the data to long format
peptide_long <- peptide_phylum_properties %>%
  pivot_longer(cols = c("hydrophobicity", "length", "charge", "iso"), names_to = "property", values_to = "value")


#########################################
# Plot the distribution of each property between
#########################################

# Define a function to perform Kruskal-Wallis test (as a common choice for non-parametric data)
perform_tests <- function(data) {
  test_result <- kruskal.test(value ~ phylum, data = data)
  tidy(test_result)  # Return a tidy version of the test results
}

# Apply tests for each property and collect results
test_results <- peptide_long %>%
  group_by(property) %>%
  summarise(test_details = list(perform_tests(cur_data())), .groups = 'drop')  # Use cur_data() to pass current group's data

# Extract p-values for annotation
p_values <- test_results %>%
  mutate(p_value = sapply(test_details, function(x) x$p.value)) %>%
  select(property, p_value)

# Merge p-values back into the long format data for plotting
peptide_long <- left_join(peptide_long, p_values, by = "property")

# Define colors for each phylum
phylum_colors <- c("#F39B7F", "#91D1C2") # Add more colors if needed
# Plotting with annotations
peptide_phylum_properties_phylum_plot <- ggplot(peptide_long, aes(x = phylum, y = value, fill = phylum)) +
  geom_boxplot() +
  scale_fill_manual(values = phylum_colors) +  # Assigning cus
  facet_wrap(~ property, scales = "free_y") +
  geom_text(data = subset(peptide_long, !duplicated(property)), 
    aes(label = ifelse(p_value < 0.05, paste("p =", format(p_value, digits = 3)), "")), 
      y = Inf, vjust = 1.5, size=2) +
  labs(title = "Comparison of Distinct Protein Properties Across Phyla",
       x = "Phylum", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw(base_size = 8)

ggsave(peptide_phylum_properties_phylum_plot, 
  file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/peptide_phylum_properties_phylum_plot.png",
  width = 14, height = 12, units = "cm", dpi = 600)


########################################
## write peptides in the phylum actinomycetotoa to a fasta file, x is the sequence and is is the fasta header
########################################
# Function to write filtered entries into FASTA format
write_filtered_fasta <- function(data, phylum_name, file_path) {
  # Filter data for the specified phylum
  filtered_data <- data %>%
    filter(phylum == phylum_name)
  
  # Write filtered entries into FASTA format
  if (nrow(filtered_data) > 0) {
    # Open the file connection
    file_conn <- file(file_path, "w")
    # Write sequence headers and sequences
    for (i in 1:nrow(filtered_data)) {
      writeLines(paste(">", filtered_data$id[i], sep = ""), file_conn)
      writeLines(filtered_data$x[i], file_conn)
    }
    
    # Manually add the additional sequences
    additional_sequences <- c(
      ">MGYG000298963_01692",
      "MALFSRLVSYAWNFGRKVVNWIWSHRNIILDWLRNGLAFDVIIVRVRRYLGI",
      ">MGYG000000642_01040",
      "MSAVAKIISIVGKYGSKAVHWCKNNVGTILNWINAGQTIEWIVNKIRRIVGV"
    )
    writeLines(additional_sequences, file_conn)
    
    # Close the file connection
    close(file_conn)
    cat("Filtered entries for phylum", phylum_name, "have been written to", file_path, "\n")
  } else {
    cat("No entries found for phylum", phylum_name, "\n")
  }
}

write_filtered_fasta(distinct(peptide_phylum_properties, x, .keep_all = TRUE), "p__Actinomycetota", "/home/david/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/actinomycetota/actinomycetota_core.faa")
write_filtered_fasta(distinct(peptide_phylum_properties, x, .keep_all = TRUE), "p__Bacillota", "/home/david/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/bacillota/bacillota_core.faa")

########################################
## Actinomycetota Core peptides alignment
########################################
writeFasta <- function(data, columnName, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    header = paste(">", data[rowNum, "id"], sep = "")
    sequence = as.character(data[rowNum, columnName])
    fastaLines = c(fastaLines, header, sequence)
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

actino_pep_df_nr = supplementary_table_1 %>%
  filter(phylum == "p__Actinomycetota")

writeFasta(actino_pep_df_nr, 
  "x", "/home/david/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/actinomycetota/actinomycetota_core.faa")


mean(actino_pep_df_nr$length)
range(actino_pep_df_nr$length)

actino_peps = readAAStringSet("/home/david/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/actinomycetota/actinomycetota_core.faa")
actino_msa = msa(actino_peps, type = "protein", method = "ClustalW")

msaPrettyPrint(actino_msa, file="/home/david/david/mining_2023/class_IId_analysis/figures/actino_msa.tex", output="tex",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90,
               verbose=FALSE, askForOverwrite=FALSE, paperWidth = 18, paperHeight = 20, alFile =  "/home/david/david/mining_2023/class_IId_analysis/figures/actino_msa.fasta")
texfile <- "/home/david/david/mining_2023/class_IId_analysis/figures/actino_msa.tex"
tinytex::pdflatex(texfile)

########################################
## Alignment of synthesized peptides
########################################

synt_peps = readAAStringSet("/home/david/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/synthesized_peptides.faa")
synt_msa = msa(synt_peps, type = "protein", method = "ClustalW")

msaPrettyPrint(synt_msa, file="/home/david/david/mining_2023/class_IId_analysis/figures/synt_msa.tex", output="tex",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90,
               verbose=FALSE, askForOverwrite=FALSE, paperWidth = 12, paperHeight = 4, alFile =  "/home/david/david/mining_2023/class_IId_analysis/figures/synt_msa.fasta")
texfile2 <- "/home/david/david/mining_2023/class_IId_analysis/figures/synt_msa.tex"
tinytex::pdflatex(texfile2)

########################################
## Phylogenetic tree of all core peptides and clusters 
########################################

peptide_clusters = fread("/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/clusters/leaderless_c50id.table") %>%
  arrange(cluster)
peptide_clusters = peptide_clusters %>% column_to_rownames("protein_acc") 

tree_lead <- read.tree("/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/reduced_tree/leaderless_cdhit90_aln.fasta.raxml.supportFBP")
tree_lead$tip.label %in% "AOHCFJ_11375"
tree_lead <- ape::root(tree_lead, outgroup = "Pediocin")
gg_lead <- ggtree(tree_lead, size=0.4, layout="circular") +
  geom_text(data = subset(gg_lead$data, label == "BLJHKG_36645"), 
            aes(label="A", x=x, y=y), 
            hjust=2, vjust=1, color="black", size=5, fontface="bold") +
  geom_text(data = subset(gg_lead$data, label == "LHJIIJ_00720"), 
            aes(label="B", x=x, y=y), 
            hjust=0.5, vjust=-1, color="black", size=5, fontface="bold") +
  geom_text(data = subset(gg_lead$data, label == "ALNPPI_03270"), 
            aes(label="C", x=x, y=y), 
            hjust=-1, vjust=-1, color="black", size=5, fontface="bold") +
  theme_tree2() 

tree_lead$tip.label
peptide_tax_annot_df = peptide_phylum_properties %>%
 column_to_rownames("id") %>%
 select(phylum) 

manual_df <- data.frame(
  row_names = c("MGYG000000642_01040", "EKOHGI_07665", "MGYG000320920_03099", "MGYG000049946_02918", "MGYG000320920_03506", "MGYG000298963_01692"),  # Manually named rows
  phylum = c("p__Actinomycetota", "p__Bacillota_A", "p__Bacillota_A", "p__Bacillota","p__Bacillota", "p__Actinomycetota")  # Manually input data for column 1
)

# Set the row names for the dataframe
rownames(manual_df) <- manual_df$row_names

# Remove the column used for row names
manual_df$row_names <- NULL

# Print the manually created dataframe
peptide_tax_annot_df = rbind(peptide_tax_annot_df, manual_df)

# A 33 is 39.5 close to NTN-A
# B 44 is 43.75 to Thucin A
# C 63 33.3 is Geobacillin 6
# singleton 79
gg_lead$data[gg_lead$data$node == 33, ]
gg_lead$data[gg_lead$data$node == 44, ]
gg_lead$data[gg_lead$data$node == 63, ]


# # Create the heatmap
p5 <- gheatmap(gg_lead, peptide_tax_annot_df,
  legend_title = "Phylum",
  color = "black",
  font.size = 2,
  colnames_angle = -90,
  colnames_offset_x = 0.3,
  colnames_offset_y = -0.1,
  offset=-0.1,
  width=0.1) +
  theme_tree() +
  scale_fill_manual(values = c("#F39B7F","#91D1C2","grey"), name = "Phylum") +
  theme(legend.text = element_text(face = "bold"),  # Make legend text bold
        legend.title = element_text(face = "bold")) # M

ggsave(p5, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figres/leaderless_core_cluster_plot_tree.png",
	width = 16, height = 16, units = "cm", dpi=600)


########################################
## COLOR EFI-EST by Taxonomy
########################################
ssn_phylum = left_join(core_peptide_df, bakta_df_2, by = c("id" = "locus_tag")) %>%
  left_join(., assembly_tax_df, by = c("assembly" = "assembly")) %>%
  select(id, x, hydrophobicity, length, charge, iso, mw, assembly, phylum, class, order, family, genus, species)

ssn_old = fread("/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/efi-est/node_table_dl.csv") %>%
 clean_names() %>%
 select(description, name, shared_name) %>%
 # mutate and remove all after first whte space
 mutate(id = sub("^([^ ]*) .*", "\\1", description)) %>% 
 left_join(., ssn_phylum, by="id") 

fwrite(ssn_old, "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/efi-est/node_table_annotated.csv")



#########################################
# Stats on the peptides
#########################################

peptide_df_stats = left_join(core_peptide_df, bakta_df_2, by = c("id" = "locus_tag")) %>%
  left_join(., assembly_tax_df, by = c("assembly" = "assembly")) %>%
  select(id, x, hydrophobicity, length, charge, iso, mw, assembly, phylum, class, order, family, genus, species) %>%
  mutate(species_component = sub("_.*", "", id))

# Count occurrences of species, genus, and phylum
count_table <- peptide_df_stats %>%
  group_by(species_component, species, genus, phylum) %>%
  summarise(count = n()) %>%
  ungroup()

# Categorize species based on the number of occurrences
count_table <- count_table %>%
  mutate(
    occurrence_category = case_when(
      count == 1 ~ "1",
      count == 2 ~ "2",
      count == 3 ~ "3",
      count == 4 ~ "4",
      count == 5 ~ "5",
      count > 5 ~ "More than 5"
    )
  ) 

stat_df_phylum = count_table %>%
  group_by(phylum, count) %>%
  summarise(n = n()) #%>%
  # filter(!is.na(phylum)) 

stat_df_phylum %>%
                 filter(count == 1) %>%
                 summarise(sum_n = sum(n)) %>%
                 pull(sum_n) %>% sum()


stat_df_genus = count_table %>%
  group_by(genus, count) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  filter(count > 1)

# Select specific columns for summary statistics
selected_columns <- c("charge", "length", "hydrophobicity", "mw", "iso")

# Get summary statistics for selected columns
summary_stats <- summary(peptide_df_stats[, selected_columns])

#########################################
# Assembly stats table, (assembly, contig, core_peptides, is multi?, genus, species, pH domain containing protein?)
#########################################

known_species_l_bac = c("s__Enterococus faecium",
  "s__Enterocoocus faecalis",
  "s__Staphylococcus capitis",
  "s__Staphylococcus coagulans",
  "s__Staphylococcus rattus",
  "s__Staphylococcus epidermidis",
  "s__Lactococcus lactis",
  "s__Lactobacillus garvieae",
  "s__Weissella hellenica",
  "s__Staphylococcus pseudintermedius",
  "s__Staphylococcus aureus",
  "s__Staphylococcus ratti",
  "s__Parageobacillus thermoglucosidasius",
  "s__Bacillus cereus",
  "s__Bacillus toyonensis",
  "s__Bacillus weidmannii")


supplementary_table_1 = left_join(core_peptide_df, bakta_df_2, by = c("id" = "locus_tag")) %>%
  left_join(., assembly_tax_df, by = c("assembly" = "assembly")) %>%
  select(assembly, id, x, hydrophobicity, length, charge, iso, mw, assembly, phylum, class, order, family, genus, species) %>%
  mutate(phylum = case_when(
    id == "MGYG000000642_01040" ~ "p__Actinomycetota",
    id == "MGYG000001470_00197" ~ "p__Bacillota",
    id == "MGYG000031742_00801" ~ "p__Bacillota", # epidermidis
    id == "MGYG000049946_02918" ~ "p__Bacillota",
    id == "MGYG000049946_02917" ~ "p__Bacillota",
    id == "MGYG000298963_01692" ~ "p__Actinomycetota",
    id == "MGYG000320920_03099" ~ "p__Bacillota_A",
    id == "MGYG000320920_03506" ~ "p__Bacillota",
    TRUE ~ phylum
  ))


fwrite(supplementary_table_1, "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/supplementary_table_1.tsv", sep = "\t")

novel_sp_not_known_to_encode = supplementary_table_1 %>% 
  filter(!species %in% known_species_l_bac) %>%
  select(assembly, genus,species) %>% 
  # make species 1 column and put all assemblies in another if they share the same species seperated by a comma
  group_by(genus, species) %>%
  summarise(assembly = paste(assembly, collapse = ", ")) %>%
  filter(!is.na(genus))
fwrite(novel_sp_not_known_to_encode, "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/novel_species_not_known_to_encode.tsv", sep = "\t")


########################################
## Plot operons from Actinomycetota that are interesting (ones with multiple core peptides)
########################################
core_ids = fread("/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/final_proteins.txt",
  header = FALSE)

core_focus = c("MNCMEN_27380", "DEKJBI_27380","EIDIJP_16820",
  "AHFNNK_09245","BDCJEE_20345","DEHBAO_05325",
  "LEFGKO_06925","LKOEGJ_10960", "ADHNEM_05635","EAHPIB_01905", "IJELBH_01450",
  "KEBEAJ_16320","LJDPIO_05525","HMOIEN_15240",
  "EKOHGI_07665","OPADLO_00475","EKFLIJ_03105","BCOKOO_01030", "HOBAFF_06210", "IKAFHA_08300")

contig_focus = c("NZ_CP036048.1","NZ_JACSNC010000005.1", 
    "NZ_JAPEMT010000002.1","NZ_VFQE01000001.1",
    "NZ_FIFQ01000008.1","NZ_MTAU01000004.1","NZ_AOCD01000018.1",
    "NZ_JAGMVT010000033.1",
    "NZ_FOPE01000006.1","NZ_PVRR01000016.1", "NZ_JADKSF010000001.1",
    "NZ_JAKVOD010000001.1", "NZ_CANLBI010000007.1")



int_df = bakta_df_2 %>%
  dplyr::rename(feat_id = locus_tag) %>%
  mutate(seq_id = number_sequence_id,
  	end = stop) %>%
  filter(seq_id %in% contig_focus) %>%
  mutate(focus_point = case_when(feat_id %in% core_focus ~ "focus",
                     TRUE ~ NA)) %>%
  mutate(color = case_when(
    feat_id == "JMKLCA_05195" ~ "#00A087B2", # membrane transpoer
    feat_id %in%  core_ids$V1 ~ "black",
    grepl("PF03703", db_xrefs) ~ "#8491B4FF",  # PH domain-containing protein
    feat_id %in% c("IKAFHA_08300","IKAFHA_08320") ~ "black", # membrane transpoer
    feat_id == "HOBAFF_06205" ~ "#E64B35FF", # membrane transpoer
    feat_id == "JMKLCA_05210" ~ "#E64B35FF", # membrane transpoer
    feat_id == "JMKLCA_05215" ~ "#E64B35FF", # also membrane transport 
    feat_id == "JMKLCA_05220" ~ "#E64B35FF", # ATP ase transmembrane
    feat_id == "JMKLCA_05225" ~ "#4DBBD5FF", # bacteriocin immunity
    product == "PH domain-containing protein" ~ "#8491B4FF",  # PH domain-containing protein
    feat_id == "HOBAFF_06190" ~ "#8491B4FF", # PH domain-containing protein
    product == "bPH-2 domain-containing protein" ~ "#8491B4FF", # bPH-2 domain-containing protein
    product == "BPH-2 domain-containing protein" ~ "#8491B4FF", # bPH-3 domain-containing protein
    product == "Putative membrane protein YdbT with pleckstrin-like domain" ~ "#8491B4FF", # membrane protein
    product == "Integral membrane protein" ~ "#8491B4FF", # membrane protein
    product == "LysR family transcriptional regulator" ~ "#00A087B2", # response regulator
    product == "Helix-turn-helix transcriptional regulator" ~ "#00A087B2", # response regulator
    product == "DNA-binding response regulator" ~ "#00A087B2", # response regulator
    product == "helix-turn-helix domain-containing protein" ~ "#00A087B2", # response regulator
    product == "Transcriptional regulator, LacI family" ~ "#00A087B2", # response regulator
    product == "LytTr DNA-binding domain-containing protein" ~ "#00A087B2", # response regulator
    product == "Positive transcriptional regulator, MutR family protein" ~ "#00A087B2", # response regulator
    product == "HTH lysR-type domain-containing protein" ~ "#00A087B2", # HTH lysR-type domain-containing protein
    product == "XRE family transcriptional regulator" ~ "#00A087B2", # HTH lysR-type domain-containing protein
    product == "ABC transporter ATP-binding protein" ~ "#E64B35FF", # abs transporter permease
    product == "ABC transporter permease subunit" ~ "#E64B35FF",  # ATP 
    product == "ABC transporter domain-containing protein" ~ "#E64B35FF", # ABC transporter permease
    product == "ABC-2 type transport system ATP-binding protein" ~ "#E64B35FF", # ABC-2 type transport system ATP-binding protein
    product == "Abc2" ~ "#E64B35FF", # ABC transporter permease
    product == "Ribose/galactoside ABC transporter permease" ~ "#E64B35FF", # ABC transporter permease
    product == "Ribose/galactoside ABC transporter ATP-binding protein" ~ "#E64B35FF", # ABC transporter permease
    product == "ABC transporter permease" ~ "#E64B35FF", # ABC transporter permease
    product == "Multidrug efflux pump subunit AcrA (membrane-fusion protein)" ~ "#E64B35FF", # multidrug efflux pump
    product == "ABC-type lipoprotein export system, ATPase component" ~ "#E64B35FF", # ABC-type lipoprotein export system, ATPase component
    product == "Macrolide ABC transporter permease" ~ "#E64B35FF", # Macrolide ABC transporter permease
    product == "Yip1 domain-containing protein" ~ "#343434", # Yip1 domain-containing protein
    product == "IS21 family transposase" ~ "orange", # IS21 family transposase
    TRUE ~ "grey90")
  ) %>%
  mutate(seq_id = case_when(
    seq_id == "NZ_CP036048.1" ~ "Bacillus mycoides BPN601",
    seq_id == "NZ_JACSNC010000005.1" ~ "Exiguobacterium indicum s191",
    seq_id == "NZ_JAPEMT010000002.1" ~ "Streptomyces sp. NBC_00237",
    seq_id == "NZ_VFQE01000001.1" ~ "Blastococcus colisei DSM 46837",
    seq_id == "NZ_FIFQ01000008.1" ~ "Streptococcus suis LSS13",
    seq_id == "NZ_MTAU01000004.1" ~ "Bacillus anthracis F34",
    seq_id == "NZ_AOCD01000018.1" ~ "Streptococcus ratti DSM 20564",
    seq_id == "NZ_JAGMVT010000033.1" ~ "Weissella diestrammenae DSM 27840",
    seq_id == "NZ_FOPE01000006.1" ~ "Lachnobacterium sp. C7",
    seq_id == "NZ_JADKSF010000001.1" ~ "Frondihabitans sp. VKM Ac-2883 n1",
    seq_id == "NZ_PVRR01000016.1" ~ "Bacillus wiedmannii NMSW16",
    seq_id ==  "NZ_JAKVOD010000001.1" ~ "Bifidobacterium tibiigranuli UW_MP_BIF2_1",
    seq_id == "NZ_CANLBI010000007.1" ~ "Brevilactibacter sinopodophylli 2B3-6",
    TRUE ~ seq_id)) 


plot_int = gggenomes(genes = int_df) %>%
  focus(focus_point == "focus", 
    .expand = c(5e3, 5e3),
    .locus_id = str_glue("{seq_id}")) +
  geom_seq() +
  geom_gene(aes(fill = color),
        size = 3) +
  scale_fill_identity() +
  geom_seq_label(size = 1.5, vjust = 2, fontface = "bold.italic", family = "Arial") 

ggsave(plot_int, file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/interesing_operons.png",
  width = 8, height = 10, units = "cm", dpi=600)

########################################
## align proteins from Clade A B and C and how structurs
########################################
output_dir = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/"

clade_rep <- c("BLJHKG_36645 aureocin A53 family class IId bacteriocin", 
  "LHJIIJ_00720 aureocin A53 family class IId bacteriocin", "ALNPPI_03270 aureocin A53 family class IId bacteriocin")
selected_sequences <- a58_core[ names(a58_core) %in% clade_rep]

clade_a 
clade_b
clade_c 

# align proteins
aligned_proteins = msa(selected_sequences, type = "protein", method = "ClustalW", 
  verbose= TRUE, order= "aligned" )


# Generate alignment output and convert to PDF
msaPrettyPrint(aligned_proteins, 
  file="/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/clade_protiens.tex", 
  output="tex", verbose=FALSE, askForOverwrite=FALSE,
  showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 100)

tinytex::pdflatex("/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/clade_protiens.tex")







########################################
## find most similar peptides to known
########################################
library(data.table)
library(tidyr)
library(dplyr)

# Read the data from the file without considering the first row as header
pid <- fread("/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/pid_to_known/known_and_new.txt", header = FALSE)

protein_names <- as.character(pid[[1]])
pid <- pid[, -1, with = FALSE]
colnames(pid) <- protein_names
pid_df <- as.data.frame(pid)
rownames(pid_df) <- protein_names

# make long
pid_long <- tibble::rownames_to_column(pid_df, "x")
pid_long <- pivot_longer(pid_long, -x, names_to = "y", values_to = "value")




# A BLJHKG_36645  
# B LHJIIJ_00720 
# C ALNPPI_03270
pid_long %>%
  filter(grepl("ALNPPI_03270", x)) %>%
  arrange(desc(value)) %>% View()

pid_long %>%
  filter(grepl("CCHAHL", x)) %>%
  arrange(desc(value)) %>% 
  filter(grepl("Thu", y)) %>% View()

# s suis ADHNEM_05635
pid_long %>%
  filter(grepl("ADHNEM_05635", x)) %>%
  arrange(desc(value)) %>% View() 


# s suis ADHNEM_05635
pid_long %>%
  filter(grepl("MGYG000000642_01040|MGYG000298963_01692", x)) %>%
  arrange(desc(value)) %>% View() 

  pid_long %>% filter(x == "MGYG000000642_01040" & y  == "MGYG000298963_01692")


IKAFHA_08300
pid_long %>%
  filter(grepl("IKAFHA_08320", x)) %>%
  arrange(desc(value)) %>%  View()
  filter(grepl("Thu", y)) %>% View()
########################################
## align interesting
########################################
bacillus8corefile = "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/interesting_alignments/8corebacillus.faa"
bacillus8core = readAAStringSet(bacillus8corefile)
msa_bacillus8core = msa(bacillus8core, type = "protein", method = "ClustalW")

msaPrettyPrint(msa_bacillus8core, file="/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/bacillus8core_msa.tex", output="tex",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90,
               verbose=FALSE, askForOverwrite=FALSE, paperWidth = 12, paperHeight = 4, alFile =  "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/bacillus8core_msa.fasta")

texfile3 <- "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/bacillus8core_msa.tex"
tinytex::pdflatex(texfile3)
########################################
corebacillusplotcorefile = "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/interesting_alignments/corebacillusplot.faa"
corebacillusplotcore = readAAStringSet(corebacillusplotcorefile)
msa_corebacillusplot = msa(corebacillusplotcore, type = "protein", method = "ClustalW")

msaPrettyPrint(msa_corebacillusplot, file="/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/corebacillusplot.tex", output="tex",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90,
               verbose=FALSE, askForOverwrite=FALSE, paperWidth = 12, paperHeight = 4, alFile =  "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/corebacillusplot_msa.fasta")

texfile4 <- "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/corebacillusplot.tex"
tinytex::pdflatex(texfile4)


########################################
# Define file paths
########################################

actinomycetota_core_df = peptide_df_stats %>% 
  filter(phylum == "p__Actinomycetota") %>% 
  select(id,x)
write_filtered_fasta()


input_file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/actinomycetota/actinomycetota_core_reduced.faa"
output_dir = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/"
proteins = readAAStringSet(input_file)

aligned_proteins = msa(proteins, type = "protein", method = "ClustalW", 
  verbose= TRUE, order= "aligned" )

# Generate output file path
output_base_name <- tools::file_path_sans_ext(basename(input_file))
output_tex = file.path(output_dir, paste0(output_base_name, "_alignment_output.tex"))

# Generate alignment output and convert to PDF
msaPrettyPrint(aligned_proteins, 
  file=output_tex, output="tex", verbose=FALSE, askForOverwrite=FALSE,
  showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90)

tinytex::pdflatex(output_tex)

########################################
# Define file paths
########################################
input_file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/actinomycetota/actinomycetota_core_reduced.faa"
output_dir = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/"
proteins = readAAStringSet(input_file)

aligned_proteins = msa(proteins, type = "protein", method = "ClustalW", 
  verbose= TRUE, order= "aligned" )

# Generate output file path
output_base_name <- tools::file_path_sans_ext(basename(input_file))
output_tex = file.path(output_dir, paste0(output_base_name, "_alignment_output.tex"))

# Generate alignment output and convert to PDF
msaPrettyPrint(aligned_proteins, 
  file=output_tex, output="tex", verbose=FALSE, askForOverwrite=FALSE,
  showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90)

tinytex::pdflatex(output_tex)


########################################
# Align proteins that were part of the search
########################################
input_file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/actinomycetota/actinomycetota_core_reduced.faa"
output_dir = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures/"
proteins = readAAStringSet(input_file)

aligned_proteins = msa(proteins, type = "protein", method = "ClustalW", 
  verbose= TRUE, order= "aligned" )

# Generate output file path
output_base_name <- tools::file_path_sans_ext(basename(input_file))
output_tex = file.path(output_dir, paste0(output_base_name, "_alignment_output.tex"))

# Generate alignment output and convert to PDF
msaPrettyPrint(aligned_proteins, 
  file=output_tex, output="tex", verbose=FALSE, askForOverwrite=FALSE,
  showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90)

tinytex::pdflatex(output_tex)


########################################
# Align proteins from arcanobacterium
########################################
input_file = "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis/actinomycetota/compre.faa"
output_dir = "/data/san/data1/users/david/mining_2023/class_IId_analysis/figures"
proteins = readAAStringSet(input_file)

aligned_proteins = msa(proteins, type = "protein", method = "ClustalW", 
  verbose= TRUE, order= "aligned" )

# Generate output file path
output_base_name <- tools::file_path_sans_ext(basename(input_file))
output_tex = file.path(output_dir, paste0(output_base_name, "_alignment_output.tex"))

# Generate alignment output and convert to PDF
msaPrettyPrint(aligned_proteins, 
  file=output_tex, output="tex", verbose=FALSE, askForOverwrite=FALSE,
  showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90)

tinytex::pdflatex(output_tex)
