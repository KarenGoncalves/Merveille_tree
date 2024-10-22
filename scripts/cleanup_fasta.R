#### Telecharger et ouvrir les librairies ####
devtools::source_gist("https://gist.github.com/KarenGoncalves/0db105bceff4ff69547ee25460dda978")

install_from_dif_sources(
  cran_packages = c("tidyverse", "tinytex", "patchwork"),
  bioconductor_packages = c("Biostrings", "msa", "treeio", "ggtree", "ape", "seqinr", "phangorn"),
  github_packages = "YuLab-SMU/ggmsa"
)

#### Ouvrir le fichier fasta  Amaryllidaceae####
# To open a fasta file with multiple sequences, use the function readDNAStringSet() or readAAStringSet()
fasta_input = "Inputs/BBE_Amaryllidaceae.fasta" %>%
  readAAStringSet(format = "fasta")

names(fasta_input)

#### Clean up names and get best sequences per gene ####
# Here, I am using regular expressions inside the "gsub" function. 
# It is a way to search and replace multiple texts at the same time.
# Whenever you see things inside [], it means "look for any of these characters": [A-Z] means "look for any upper-case letter
# \\w means any alpha numeric character (upper and lower case letters, _ and numbers)
# \\d means any number
# \\. and \\- means these exact symbols: . and -
# Text inside "()" (in gsub, in quotes) are patterns we want to refer to, you count the number of () to call it back
# \\1 means: first pattern; \\2: 2nd pattern; etc.+score=([\\d\\.\\-]+),[A-Z_\\d]+\\|.+"

prediction_scores <- 
  gsub("^([A-Za-z0-9_\\.\\-]+.p[0-9]+) .+score=([0-9\\.\\-]+)[ ,].+",
       "\\1 \\2",
       names(fasta_input)) %>%
  str_split(" ", simplify = T) %>% 
  data.frame() %>% 
  rename(X1="Protein_ID",
         X2="score") %>% 
  mutate(score = as.numeric(score),
         gene_ID = gsub("(_i\\d+)*.p\\d+" , "", Protein_ID))

# To get the length of the proteins, we can use the function nchar (count the number of letters in the text)
prediction_scores <- prediction_scores %>% 
  mutate(protein_seq = fasta_input %>% as.character(),
           protein_length = nchar(protein_seq))  %>% 
  group_by(gene_ID) %>% # now we group the proteins derived from the same gene
  # put the longest, with greatest prediction score and the protein ID at the top
  arrange(desc(protein_length), desc(score), Protein_ID) %>% 
  # Keep only the first protein for each gene
  slice_head()

#### Create a subset of the fasta file with only the selected proteins ####
# grep is like CTRL+F: searches for a text
selected_protein_IDs <- sapply(prediction_scores$Protein_ID, \(x) {
  grep(x, names(fasta_input))
})


selected_proteins <- 
  fasta_input[selected_protein_IDs]

# Cean up the sequence names
names(selected_proteins) <- 
  gsub("^([A-Za-z0-9_\\.\\-]+)(_i[0-9]+)*.p[0-9]+ .+score=[0-9\\.\\-]+,([A-Z_0-9]+)\\|.+", "\\3",
       names(selected_proteins)) %>%
  gsub("^([A-Za-z0-9_\\.]+)(_i[0-9]+)*.p[0-9]+ .+",
       "\\1\\2", .) %>% 
  gsub("_[A-Z]+$", "", .) %>%
  gsub("RETO", "Reticuline oxidase", .) %>%
  gsub("ASO", "O-acetylstemmadenine oxidase", .) %>%
  gsub("CASL1", "Cannabidiolic acid synthase-like 1", .) %>% 
  paste(prediction_scores$gene_ID)

writeXStringSet(selected_proteins, 
                filepath = "results/Amaryllidaceae_BBEs.fasta", 
                format = "fasta", width = 100000) # maximum width=200001


#### Ouvrir le fichier fasta  de Uniprot####
# To open a fasta file with multiple sequences, use the function readDNAStringSet() or readAAStringSet()
fasta_input = "Inputs/Characterized_seqs.fa" %>%
  readAAStringSet(format = "fasta")
names(fasta_input)

#### Clean up names and get best sequences per gene ####
names(fasta_input) <- 
  gsub("^sp\\|[A-Z0-9\\.]+\\|([A-Z0-9_]+) .+", "\\1", names(fasta_input))
#### Create a subset of the fasta file with only the selected proteins ####
# grep is like CTRL+F: searches for a text

writeXStringSet(fasta_input, 
                filepath = "results/Characterized_BBEs.fasta", format = "fasta", width = 10000)


#### Ouvrir les deux ensemble

fasta_amaryllidaceae = "results/Amaryllidaceae_BBEs.fasta" %>%
  readAAStringSet(format = "fasta")
fasta_input_characterized = "results/Characterized_BBEs.fasta" %>%
  readAAStringSet(format = "fasta")

writeXStringSet(c(fasta_amaryllidaceae, fasta_input_characterized), 
                filepath = "results/Combined_fasta.fa",
                format = "fasta", width=10000)
