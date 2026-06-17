library("phyloseq") # carico Phyloseq
library("dplyr")
library("tidyr") 
library("ggplot2")
library("cowplot")
library("rstatix") # PostHoc test (alpha-diversity)
library("vegan") # PERMANOVA test

############Analisi###############

# Path della working directory
setwd("C:/Users/alexa/OneDrive/Documenti/Documenti vari/2026/Lezioni/lm7_maggio/Biotecnologie Fitopatologiche - Materiale esercitazione R")

# Phyloseq----------------------------------------------------------------------------------------------------------------
#Carico l'ASV table e preparo i file
ASVtable <- read.csv(file = "16S_ASVTable.txt", sep = "\t")

colnames(ASVtable)
ASVtable_Counts <- ASVtable[,1:28]

rownames(ASVtable_Counts) <- ASVtable_Counts$ASVID
ASVtable_Counts$ASVID <- ordered(ASVtable_Counts$ASVID)

ASVtable_Counts$ASVID <- NULL

OTU <- otu_table(ASVtable_Counts, taxa_are_rows = TRUE)
OTU[is.na(OTU)] <- 0
otu_table(OTU) <- otu_table(OTU, taxa_are_rows = taxa_are_rows(OTU))

colnames(ASVtable)
ASVtable_Taxa <- ASVtable[,c(1,29:34)]
rownames(ASVtable_Taxa) <- ASVtable_Taxa$ASVID
ASVtable_Taxa$ASVID <- ordered(ASVtable_Taxa$ASVID)
ASVtable_Taxa$ASVID <- NULL
ASVtable_Taxa <- as.matrix(ASVtable_Taxa)

TAXA <- tax_table(ASVtable_Taxa)

MAP <- read.csv(file = "metadata_16s.txt", sep = "\t")
#MAP$Terminal <- NULL

rownames(MAP) <- MAP$ASVID
MAP$X.ASV.ID <- NULL
MAP <- sample_data(MAP)

data_ps <- merge_phyloseq(OTU, TAXA, MAP)

# Filtraggio: Voglio vedere i mitocondri e cloroplasti
data_ps_CLORO_MITO <- subset_taxa(data_ps, Order == "Chloroplast" | Family == "Mitochondria")
data_ps_CLORO_MITO

# Filtraggio: Considero tutto ciò che è stato identificato come batterio, scartando Mitochondria e Chloroplast
data_ps <- subset_taxa(data_ps, Order != "Chloroplast" & Family != "Mitochondria")
data_ps

# Ordino il nome delle variabili "Compartment"
data_ps@sam_data$Compartment <- factor(
  data_ps@sam_data$Compartment,
  levels = c("Rhizosphere", "Endorhizosphere", "Xylem")
)

##############1. Abundances-----------------------------------------------------------------------------------------
### 1.1 Absolute###
data_ps_lev2 <- tax_glom(data_ps, taxrank = 'Phylum')
data_ps_lev3 <- tax_glom(data_ps, taxrank = 'Class')
data_ps_lev4 <- tax_glom(data_ps, taxrank = 'Order')
data_ps_lev5 <- tax_glom(data_ps, taxrank = 'Family')
data_ps_lev6 <- tax_glom(data_ps, taxrank = 'Genus')

### Example: save this object as a tabular dataframe###
counts_lev2 <- as.data.frame(data_ps_lev2@otu_table)
taxa_lev2 <- as.data.frame(data_ps_lev2@tax_table)

counts_lev2$ASVID <- row.names(taxa_lev2)
taxa_lev2$ASVID <- row.names(taxa_lev2)

finalTable_lev2 <- merge(counts_lev2,taxa_lev2, by = "ASVID")
finalTable_lev2$ASVID <- NULL
finalTable_lev2 <- finalTable_lev2[,c(29,1:28)]
finalTable_lev2

finalTable_lev2$somma <- rowSums(finalTable_lev2[,2:7])
finalTable_lev2 <- finalTable_lev2 %>%
  arrange(desc(somma))

### 1.2. Relative phylum
data_ps_lev2_REL <- (transform_sample_counts(data_ps_lev2, function(x) 100 * x/sum(x)))

counts_lev2_REL <- as.data.frame(data_ps_lev2_REL@otu_table)
taxa_lev2_REL <- as.data.frame(data_ps_lev2_REL@tax_table)

counts_lev2_REL$ASVID <- row.names(taxa_lev2)
taxa_lev2_REL$ASVID <- row.names(taxa_lev2)

finalTable_lev2_REL <- merge(counts_lev2_REL, taxa_lev2_REL, by = "ASVID")
finalTable_lev2_REL$ASVID <- NULL

colnames(finalTable_lev2_REL)

finalTable_lev2_REL <- finalTable_lev2_REL[,c(29,1:28)]
finalTable_lev2_REL$somma <- rowSums(finalTable_lev2_REL[,2:27])
finalTable_lev2_REL <- finalTable_lev2_REL %>%
  arrange(desc(somma))

finalTable_lev2_REL
#Relative family
data_ps_lev5_REL <- (transform_sample_counts(data_ps_lev5, function(x) 100 * x/sum(x)))

counts_lev5_REL <- as.data.frame(data_ps_lev5_REL@otu_table)
taxa_lev5_REL <- as.data.frame(data_ps_lev5_REL@tax_table)

counts_lev5_REL$ASVID <- row.names(taxa_lev5_REL)
taxa_lev5_REL$ASVID <- row.names(taxa_lev5_REL)

finalTable_lev5_REL <- merge(counts_lev5_REL, taxa_lev5_REL, by = "ASVID")
finalTable_lev5_REL$ASVID <- NULL

colnames(finalTable_lev5_REL)

finalTable_lev5_REL <- finalTable_lev5_REL[,c(32,1:31)]
finalTable_lev5_REL$somma <- rowSums(finalTable_lev5_REL[,2:28])
finalTable_lev5_REL <- finalTable_lev5_REL %>%
  arrange(desc(somma))

write.table(finalTable_lev5_REL, file = "16Data_Family.txt", sep = "\t",
            row.names = FALSE)
#Relative genus
data_ps_lev6_REL <- (transform_sample_counts(data_ps_lev6, function(x) 100 * x/sum(x)))

counts_lev6_REL <- as.data.frame(data_ps_lev6_REL@otu_table)
taxa_lev6_REL <- as.data.frame(data_ps_lev6_REL@tax_table)

counts_lev6_REL$ASVID <- row.names(taxa_lev6_REL)
taxa_lev6_REL$ASVID <- row.names(taxa_lev6_REL)

finalTable_lev6_REL <- merge(counts_lev6_REL, taxa_lev6_REL, by = "ASVID")
finalTable_lev6_REL$ASVID <- NULL

colnames(finalTable_lev6_REL)

finalTable_lev6_REL <- finalTable_lev6_REL[,c(33,1:32)]
finalTable_lev6_REL$somma <- rowSums(finalTable_lev6_REL[,2:28])
finalTable_lev6_REL <- finalTable_lev6_REL %>%
  arrange(desc(somma))

write.table(finalTable_lev6_REL, file = "16Data_Genus.txt", sep = "\t",
            row.names = FALSE)


### 1.3. Plot - Mean relative abundances ###

# 1.3.1. PHYLUM
df_phylum_media <- data_ps_lev2_REL %>%
  psmelt() %>%                                         # Converte l'oggetto phyloseq in dataframe
  group_by(Treatment_Compartment, Phylum) %>%                           # Raggruppa per tesi allo studio e phylum
  summarize(Mean_Relative_Abundance = mean(Abundance), .groups = 'drop') %>%
  mutate(Phylum = ifelse(Mean_Relative_Abundance < 1, "Other taxa (< 1%)", Phylum)) #Setto il valore soglia
  #summarize(Mean_Relative_Abundance = sum(Mean_Relative_Abundance), .groups = 'drop') %>%
  #pivot_wider(names_from = Phylum, values_from = Mean_Relative_Abundance, values_fill = 0)


#Ordering the feature table according the treatment variables
df_phylum_media$Treatment_Compartment <- factor(
  df_phylum_media$Treatment_Compartment,
  levels = c("Control_Rhizosphere", "Leaf_Rhizosphere", "Roots_Rhizosphere",
             "Control_Endorhizosphere", "Leaf_Endorhizosphere", "Roots_Endorhizosphere",
             "Control_Xylem", "Leaf_Xylem", "Roots_Xylem")
)

#Plot  
grafico_16S_phylum <- ggplot(df_phylum_media, aes(x = Treatment_Compartment, y = Mean_Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(y = "Rel. abundance (%)", x = "Samples") +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

grafico_16S_phylum

# 1.3.2 FAMIGLIA
df_family_media <- data_ps_lev5_REL %>%
  psmelt() %>%                                         # Converte l'oggetto phyloseq in dataframe
  group_by(Treatment_Compartment, Family) %>%                           # Raggruppa per tesi e phylum
  summarize(Mean_Relative_Abundance = mean(Abundance), .groups = 'drop') %>%
  mutate(Family = ifelse(Mean_Relative_Abundance < 5, "Other taxa (< 5%)", Family))


df_family_media$Treatment_Compartment <- factor(
  df_family_media$Treatment_Compartment,
  levels = c("Control_Rhizosphere", "Leaf_Rhizosphere", "Roots_Rhizosphere",
             "Control_Endorhizosphere", "Leaf_Endorhizosphere", "Roots_Endorhizosphere",
             "Control_Xylem", "Leaf_Xylem", "Roots_Xylem")
)

#Plot
grafico_16S_family <- ggplot(df_family_media, aes(x = Treatment_Compartment, y = Mean_Relative_Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(y = "Rel. abundance (%)", x = "Samples") +
  scale_fill_viridis_d(option = "H") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

grafico_16S_family

grafici_16S_MERGED <- plot_grid(grafico_16S_phylum, grafico_16S_family,
                               rel_widths = c(1, 1), labels = "AUTO")

grafici_16S_MERGED

ggsave2(plot = grafici_16S_MERGED,
       filename = "grafico_16S_RELATIVE.jpeg",
       width = 14,
       height = 8.50,
       dpi = 300)

########################2. Alfa diversità-------------------------------------------------------------------------

#2.1. Risultato tabulare

data_alpha <- estimate_richness(data_ps)

write.table(data_alpha, file = "output_alpha.txt", sep = "\t")

#2.2. Grafico

#Opzione A
data_alpha_grafico_OBS <- plot_richness(data_ps, measures=c("Observed"),
                                    x="Compartment", color="Treatment") +
  facet_grid("Treatment")

data_alpha_grafico_OBS

data_alpha_grafico_SHAN <- plot_richness(data_ps, measures=c("Shannon"),
                                    x="Compartment", color="Treatment") +
  facet_grid("Treatment")

data_alpha_grafico_SHAN
#Opzione B
# Boxplot
metadata <- read.csv(file = "metadata_16s.txt", sep = "\t")
colnames(metadata)[1] <- "Samples"

data_alpha$Samples <- rownames(data_alpha)

data_alpha_metadata <- merge(data_alpha, metadata, by = "Samples")

data_alpha_metadata$Compartment <- factor(
  data_alpha_metadata$Compartment,
  levels = c("Rhizosphere", "Endorhizosphere", "Xylem")
)


data_alpha_ggplot2_OBSERVED <- ggplot(data = data_alpha_metadata,
                                      aes(x = Compartment, y = Observed, fill = Treatment)) +
  geom_boxplot()

data_alpha_ggplot2_OBSERVED

data_alpha_ggplot2_SHANNON <- ggplot(data = data_alpha_metadata,
                                     aes(x = Compartment, y = Shannon,
                                         fill = Treatment)) +
  geom_boxplot()

data_alpha_ggplot2_SHANNON

#UNIAMO I DUE GRAFICI
data_alpha_ggplot2_OBSERVED_SHANNON <- plot_grid(data_alpha_ggplot2_OBSERVED,
                                                 data_alpha_ggplot2_SHANNON,
                                                 labels = "AUTO",
                                                 label_size = 12)
data_alpha_ggplot2_OBSERVED_SHANNON

########################3. Beta########################
ord <- ordinate(data_ps, "PCoA", "bray", weighted = TRUE)

plot_beta <- plot_ordination(data_ps, ord, color = "Compartment", shape = "Treatment") + #, +
  #shape = "Treatment") +
  geom_point(size = 4, alpha = 0.5) +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(plot.subtitle = element_text(hjust=0.5)) +
  theme(legend.position = "right") +
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  scale_x_continuous(labels = ~sub("-", "\u2212", .x))

plot_beta <- plot_beta + geom_point(size=5, alpha=2) + 
  theme(text=element_text(size=14)) #+ scale_color_manual(values = c("#89CEEF", "#DC4F4D"))

plot_beta

ggsave(filename = "beta_div.jpeg",
       plot_beta,
       dpi = 300,
       width = 7.8,
       height = 4.09)

#############STATISTICA############

#1. Alpha Diversity
#Esempio basato sull'indice Observed
shapiro.test(data_alpha_metadata$Observed)
kruskal.test(Observed ~ Treatment_Compartment, data = data_alpha_metadata)

dunn <- dunn_test(data_alpha_metadata, Observed ~ Treatment_Compartment, 
                  p.adjust.method = "fdr", detailed = TRUE)

#2. Beta Diversity
dist_matrix <- phyloseq::distance(data_ps, method = "bray")

beta_stats <- adonis2(dist_matrix ~ Compartment*Treatment,
                             data= as.data.frame(as.matrix(sample_data(data_ps))),
                             permutations = 999,
                           by = "terms")

beta_stats$Variables <- rownames(beta_stats)
rownames(beta_stats) <- NULL
colnames(beta_stats)
beta_stats <- beta_stats[,c(6,1:5)] 

write.table(beta_stats,
            file = "beta_stats.txt",
            sep = "\t",
            row.names = FALSE)

#######4. Differential abundance analysis################

######Pairwise comparisons at the genus level
#Root inoculation vs control - Rhizosphere
data_ps_RRvsCK <- subset_samples(data_ps_lev6, Treatment_Compartment == "Control_Rhizosphere" | Treatment_Compartment == "Roots_Rhizosphere")

data_ps_RRvsCK <- prune_taxa(taxa_sums(data_ps_RRvsCK) > 0, data_ps_RRvsCK)
#RRvsCK
diagdds <- phyloseq_to_deseq2(data_ps_RRvsCK, ~ Treatment_Compartment)
diagdds$Treatment_Compartment <- relevel(diagdds$Treatment_Compartment, ref = "Control_Rhizosphere")

diagdds <- DESeq(diagdds, test="Wald", fitType="local")

res <- results(diagdds,
               cooksCutoff = FALSE)
alpha <- 0.1
sigtab <- res[which(res$padj < alpha), ]

sigtab <- cbind(as(sigtab, "data.frame"))#,
head(res)

sigtab$AsvId <- rownames(sigtab)
sigtab <- sigtab[order(sigtab$padj), ]

data_DA_otu <- as.data.frame(data_ps_RRvsCK@otu_table)
data_DA_otu$AsvId <- rownames(data_DA_otu)

data_DA_taxa <- as.data.frame(data_ps_RRvsCK@tax_table)
data_DA_taxa$AsvId <- rownames(data_DA_taxa)

data_DA_def <- merge(sigtab, data_DA_taxa, by = "AsvId")
data_DA_def <- merge(data_DA_def, data_DA_otu, by = "AsvId")

data_DA_def <- data_DA_def %>%
  arrange(padj, desc(log2FoldChange))

data_DA_def <- data_DA_def[,c(13,2:12,14:19)]
dir.create("output_differential/")
write.table(data_DA_def, file = "output_differential/diff_RRvsControl.tsv",
            sep = "\t",row.names = F)