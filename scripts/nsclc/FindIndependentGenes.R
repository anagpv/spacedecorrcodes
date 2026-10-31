#----------------------------------
#---- Identifying independent genes
#----------------------------------
library(tidyverse)

dir.data <- "/path/to/zenodo/download/"
# CosMx Lung data ------------------------------------------------------------
# Get genes in ST dataset
library(Giotto)
load(paste0(dir.data, "Lung_Giotto.RData"))
geneids <- unique(gem@feat_metadata$rna$feat_ID)

# Get ontology of genes
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

go_data <- getBM(
  attributes = c("hgnc_symbol", "go_id", "name_1006", "namespace_1003"),
  filters = "hgnc_symbol",
  values = geneids,
  mart = ensembl
)

# get entrez id
library(org.Hs.eg.db)
library(AnnotationDbi)
entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = geneids,
                                    keytype = "SYMBOL",
                                    columns = c("ENTREZID"))

# BioGRID -----------------------------------------------------------------
# Protein-protein interactions (PPI)- physical
# Genetic interactions
# Post-Translational Modifications (PTMs) - related to protein function or regulation
# chemical interactions
# CRISPR screens

biogrid <- data.table::fread(paste0(dir.data, "BIOGRID-ALL-4.4.245.tab3.zip"))
gene_pairs_biogrid <- biogrid %>%
  filter(`Official Symbol Interactor A` %in% geneids,
         `Official Symbol Interactor B` %in% geneids) %>%
  dplyr::select(`Official Symbol Interactor A`,`Official Symbol Interactor B`,
         `Experimental System Type`) %>%
  dplyr::rename(gene_to = `Official Symbol Interactor A`,
         gene_from = `Official Symbol Interactor B`) %>%
  mutate(biogrid = 1)

gene_pairs_biogrid_t <- gene_pairs_biogrid
colnames(gene_pairs_biogrid_t) <- c("gene_from","gene_to","Experimental System Type","biogrid")

gene_pairs_biogrid <- rbind(gene_pairs_biogrid, gene_pairs_biogrid_t[,c(2,1,3,4)])

# STRING ------------------------------------------------------------
# Get genes without known protein interactions
library(STRINGdb)
string_db <- STRINGdb$new(version = "12", species = 9606,
                          score_threshold = 0)  # Human species code

# Map to string ID
genes <- data.frame(gene = geneids)
mapped <- string_db$map( genes,"gene", removeUnmappedRows = FALSE )

# Get protein interactionsppi_multiscore <- string_db$get_multiscore(mapped$STRING_id)
ppi_network <- string_db$get_interactions(mapped$STRING_id)

# Change back to gene symbol
ppi_network <- merge(ppi_network, mapped[, c("STRING_id", "gene")],by.x = "from", by.y = "STRING_id", all.x = TRUE)
ppi_network <- merge(ppi_network, mapped[, c("STRING_id", "gene")],by.x = "to", by.y = "STRING_id", all.x = TRUE,
                     suffixes = c("_from", "_to"))
ppi_network <- ppi_network[, c("gene_from", "gene_to", "combined_score")]

# Find the ones that do not show interaction
gene_pairs <- expand.grid(mapped$gene, mapped$gene)
gene_pairs <- gene_pairs[gene_pairs$Var1 != gene_pairs$Var2, ]  # Remove self-interactions
colnames(gene_pairs) <- c("gene_from", "gene_to")

# Merge to check which pairs exist in the STRING interaction data
gene_pairs_stringdp <- merge(gene_pairs, ppi_network, by = c("gene_from", "gene_to"), all.x = TRUE)
gene_pairs_stringdp <- gene_pairs_stringdp %>%
  mutate(stringdp_score = ifelse(is.na(combined_score),0,combined_score)) %>%
  dplyr::select(!combined_score)

gene_pairs_stringdp_t <- gene_pairs_stringdp
colnames(gene_pairs_stringdp_t) <- c("gene_from","gene_to","stringdp_score")
gene_pairs_stringdp <- rbind(gene_pairs_stringdp, gene_pairs_stringdp_t[,c(2,1,3)])

# Filter out pairs that exist in STRING (i.e., keep only NA values in 'combined_score')
#non_interacting_pairs <- non_interacting_pairs[is.na(non_interacting_pairs$combined_score), c("gene_from", "gene_to")]


# TRRUST ------------------------------------------------------------------
trrust <- read.delim(paste0(dir.data,"trrust_rawdata.human.tsv"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(trrust) <- c("TF", "Target", "Mode", "PubMed_ID")

# Filter only gene pairs that are in your list
gene_pairs_trrust <- trrust %>%
  filter(TF %in% geneids, Target %in% geneids) %>%
  dplyr::select(TF, Target, Mode) %>%
  dplyr::rename(gene_from = TF,
                gene_to = Target) %>%
  mutate(trrust = 1)

# Add the reverse direction (symmetric interactions)
gene_pairs_trrust_t <- gene_pairs_trrust
colnames(gene_pairs_trrust_t) <- c("gene_to", "gene_from", "Mode", "trrust")

# Combine original and transposed
gene_pairs_trrust <- rbind(gene_pairs_trrust, gene_pairs_trrust_t[, c("gene_from", "gene_to", "Mode", "trrust")])


# DoRothEA ----------------------------------------------------------------
library(dorothea)
data(dorothea_hs)  # Human TF-target interactions

gene_pairs_dorothea <- dorothea_hs %>%
  filter(tf %in% geneids, target %in% geneids) %>%
  dplyr::select(tf, target, confidence) %>%
  dplyr::rename(gene_from = tf,
                gene_to = target) %>%
  mutate(dorothea = 1)

gene_pairs_dorothea_t <- gene_pairs_dorothea
colnames(gene_pairs_dorothea_t) <- c("gene_to", "gene_from", "confidence", "dorothea")

# Combine both directions
gene_pairs_dorothea <- rbind(gene_pairs_dorothea, gene_pairs_dorothea_t[, c("gene_from", "gene_to", "confidence", "dorothea")])


# CHEA --------------------------------------------------------------------
library(GSEABase)
library(dplyr)
process_chea_gmt <- function(gmt_path, geneids, source_label) {
  gmt <- getGmt(gmt_path)

  # Flatten into a TF–target data frame
  chea_df <- lapply(gmt, function(set) {
    data.frame(tf = setName(set), target = geneIds(set), stringsAsFactors = FALSE)
  }) %>% bind_rows()

  # Filter for TF–target pairs within your gene list
  gene_pairs <- chea_df %>%
    filter(tf %in% geneids, target %in% geneids) %>%
    dplyr::rename(gene_from = tf, gene_to = target) %>%
    mutate(!!source_label := 1)

  # Add reverse/symmetric interactions
  gene_pairs_t <- gene_pairs
  colnames(gene_pairs_t) <- c("gene_to", "gene_from", source_label)

  bind_rows(
    gene_pairs,
    gene_pairs_t[, c("gene_from", "gene_to", source_label)]
  )
}

gene_pairs_remap <- process_chea_gmt(paste0(dir.data,"ReMap_ChIP-seq.gmt"), geneids, "remap")
gene_pairs_encode <- process_chea_gmt(paste0(dir.data,"ENCODE_ChIP-seq.gmt"), geneids, "encode")


# Shared regulators -------------------------------------------------------
# Pairs of genes that are targets of the same TF (curated databases only).
# Unlike the blocks above, the TF does not need to be in the gene list; only
# the targets do. ReMap is not used here: its ChIP-based target sets are so
# broad that almost every pair of genes shares some bound TF.
regulons_curated <- bind_rows(
  trrust %>% dplyr::select(tf = TF, target = Target),
  dorothea_hs %>% filter(confidence %in% c("A", "B", "C")) %>% dplyr::select(tf, target)
) %>%
  filter(target %in% geneids) %>%
  distinct()


# every pair of targets that share a TF (both directions are included)
gene_pairs_sharedtf <- regulons_curated %>%
  inner_join(regulons_curated, by = "tf", relationship = "many-to-many") %>%
  filter(target.x != target.y) %>%
  dplyr::rename(gene_from = target.x, gene_to = target.y) %>%
  group_by(gene_from, gene_to) %>%
  summarise(n_shared_tf = n_distinct(tf),
            shared_tfs = paste(sort(unique(tf)), collapse = ";"),
            .groups = "drop") %>%
  mutate(shared_tf = 1)

# Pathways ----------------------------------------------------------------
# Check number of pathways in common betwen genes
library(msigdbr)

# Get genes and pathways
msigdb_data_string <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gene_symbol %in% geneids,
         !gs_collection %in% c("C8")) #remove cell type and immune pathways
gene_pathways_string <- msigdb_data_string %>% dplyr::select(gene_symbol, gs_name) %>% distinct()

# Create a list of pathways per gene
gene_pathway_list_string <- split(gene_pathways_string$gs_name, gene_pathways_string$gene_symbol)

# Create a binary gene-pathway matrix
gene_pathway_matrix_string <- gene_pathways_string %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = gs_name, values_from = value, values_fill = 0) %>%
  column_to_rownames("gene_symbol")

colSums( gene_pathway_matrix_string)[colSums( gene_pathway_matrix_string) > 1]

# Jaccard distance - 1-(# intersection)/(# union)
# distance = 1 no pathways in common, distance = 0 all pathways in common
distance_matrix_string <- as.matrix(dist(gene_pathway_matrix_string, method = "binary"))
colnames(distance_matrix_string) <- rownames(gene_pathway_matrix_string)
rownames(distance_matrix_string) <- rownames(gene_pathway_matrix_string)

gene_pairs_pathway <- reshape2::melt(distance_matrix_string)
colnames(gene_pairs_pathway) <- c("gene_from","gene_to","pathway_score")

gene_pairs_pathway_t <- gene_pairs_pathway
colnames(gene_pairs_pathway_t) <- c("gene_from","gene_to","pathway_score")
gene_pairs_pathway <- rbind(gene_pairs_pathway, gene_pairs_pathway_t[,c(2,1,3)])


# Combining  -------------------------------------------
dependenceGenes <- gene_pairs_stringdp %>%
  full_join(gene_pairs_pathway, by = c("gene_from", "gene_to")) %>%
  full_join(gene_pairs_biogrid, by = c("gene_from", "gene_to")) %>%
  full_join(gene_pairs_trrust, by = c("gene_from", "gene_to")) %>%
  full_join(gene_pairs_dorothea, by = c("gene_from", "gene_to")) %>%
  full_join(gene_pairs_remap, by = c("gene_from", "gene_to")) %>%
  full_join(gene_pairs_sharedtf, by = c("gene_from", "gene_to"))

dependenceGenes %>%
  mutate(gene_from = as.character(gene_from),gene_to = as.character(gene_to),
    pair = paste(pmin(gene_from, gene_to), pmax(gene_from, gene_to), sep = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  dplyr::select(-pair) %>%
  dplyr::rename(pathway_distance = pathway_score) %>%
  write.csv(file = paste0(dir.data,"KnownInteractionsDatabases_Lung.csv"))
