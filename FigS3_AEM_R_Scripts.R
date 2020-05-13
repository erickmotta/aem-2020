### DIVERSITY ANALYSES:

### Step 1: Load packages
library(ggplot2)
library(vegan)
library(phyloseq)
library(DESeq2)
library(grid)
library(gridExtra)

### Step 2: Import files needed for plotting diversity
# Use ASV table used for other analyses in QIIME2, converted into JSON format with the script biom convert.
asv_table_wre <- import_biom("tables/diversity_tables/WRE/feature-table-sorted_w_tax_json.biom")

# Use the tree.nwk tree from QIIME2
# Obtain the rooted-tree.qza (command:qiime phylogeny align-to-tree-mafft-fasttree) and use it as input to obtain the tree.nwk (command: qiime tools export)
tree_file_wre <- read_tree_greengenes("tables/diversity_tables/WRE/tree.nwk")

#Before importing the mapping file, make a new copy of the file and remove the # from the header
mapping_file_wre <- import_qiime_sample_data("tables/diversity_tables/WRE/sample_map.txt")

### Step 3: Merge data into one phyloseq-class object and check contents
data_wre <- merge_phyloseq(asv_table_wre,mapping_file_wre,tree_file_wre)

# Check taxonomy rank names of phyloseq object
rank_names(data_wre)

# Use colnames to assign the correct names to the columns
tax_names <- data.frame(tax_table(data_wre))
colnames(tax_names) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_table(data_wre) <- tax_table(as.matrix(tax_names))

# Check taxonomy rank names again
rank_names(data_wre)

#Check what genera are represented in dataset
get_taxa_unique(data_wre,"Genus")

#Check number of taxa and samples in dataset
ntaxa(data_wre)
nsamples(data_wre)

### Step 4: Scale reads to even depth
### matround function
# Better rounding function than R's base round
matround <- function(x){trunc(x+0.5)}

### scale_reads function
# Function to scale reads (written by my prev lab mate)
# Modified from code written by Michelle Berry, available at http://deneflab.github.io/MicrobeMiseq/ 
# Scales reads by 
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding 
# Default for n is the minimum sample size in your library
# Default for round is floor
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "round") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, 
                                          function(x) {(n * x/sum(x))}
  )
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- round(otu_table(physeq.scale))
  } else if (round == "matround"){
    otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

### Run the commands
data_wre_scale <- scale_reads(physeq = data_wre,
                              n = min(sample_sums(data_wre)),
                              round = "matround")

#See OTU table:
otu_table(data_wre_scale)

#Extract abundance matrix from the phyloseq object
OTU1_scale = as(otu_table(data_wre_scale), "matrix")
# transpose if necessary
if(taxa_are_rows(data_wre_scale)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

sum_wre_scale <- colSums(OTU1_scale)

### Step 5: Plotting beta diversity metrics
### Bray-Curtis
library(proto)
library(devtools)
library(plotly)

# Make ordination, specifying the type of ordination (PCoA), the metric (bray-curtis).
bray_curtis_ordination_wre <- ordinate(data_wre_scale,"PCoA",distance="bray")

# Plot ordination
pcoa_bray_curtis_wre = plot_ordination(data_wre_scale,bray_curtis_ordination_wre,color="Condition",title="Bray-Curtis")

PCoA_Bray_WRE <- pcoa_bray_curtis_wre + 
  scale_color_manual(values=c('#CDCDCD','#FFCDA1','#A8A8A8','#F7921E','#7A7A7A','#D35000')) +
  geom_point(aes(shape=Condition), size=3) + 
  scale_shape_manual(values = c(15,15,16,16,17,17)) +
  theme_bw() +
  ggtitle("PCoA plot of Bray-Curtis dissimilarity") +
  theme(plot.title = element_text(face = "plain", hjust = 0.5, size = 10)) +
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  xlab('PCo1 - 24.3%') +
  ylab('PCo2 - 14.4%')

### Weighted UniFrac
# Make ordination, specifying the type of ordination (PCoA), the metric (unifrac) and if it is weighted or not. 
weighted_unifrac_ordination_wre <- ordinate(data_wre_scale,"PCoA",distance="unifrac",weighted=TRUE)

# Plot ordination
pcoa_weighted_unifrac_wre = plot_ordination(data_wre_scale,weighted_unifrac_ordination_wre,color="Condition",title="Weighted UniFrac")

PCoA_weighted_WRE <- pcoa_weighted_unifrac_wre + 
  scale_color_manual(values=c('#CDCDCD','#FFCDA1','#A8A8A8','#F7921E','#7A7A7A','#D35000')) +
  geom_point(aes(shape=Condition), size=3) + 
  scale_shape_manual(values = c(15,15,16,16,17,17)) +
  theme_bw() +
  ggtitle("PCoA plot of weighted UniFrac dissimilarity") +
  theme(plot.title = element_text(face = "plain", hjust = 0.5, size = 10)) +
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  xlab('PCo1 - 39.2%') +
  ylab('PCo2 - 24.6%')

# Plotiing graphs together
plot_grid(PCoA_Bray_WRE, PCoA_weighted_WRE, labels = "AUTO", nrow = 2)
