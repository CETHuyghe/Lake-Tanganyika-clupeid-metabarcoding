#######################################################################################
# Clupeid Lake Tanganyika 16S Microbiome Analysis Huyghe & Aerts et al. (under review)
#######################################################################################

# Set working directory
setwd("Working_directory")

# Upload count matrix
counts <-read.csv(file="Count_matrix_16S_clupeid_17042023.csv",sep=",",check.names=FALSE,row.names = 1)

# Upload metadata of sardine host
metadata <-read.csv(file="COI_metadata_R_all_18042023.csv",sep=",",check.names=FALSE,row.names = 1)

# Upload OTU classifications
OTUs <-read.csv(file="OTU_ID_16S_Clupeid_17042023.csv",sep=",",check.names=FALSE,row.names = 1)

# Remove OTUs present in the blanks
# Remove OTUs of which more than 5% of the reads are present in blank
tot_reads_OTUs <- rowSums(counts) # Total nr of reads per OTU
rel_reads_OTUs <- as.data.frame((counts) / tot_reads_OTUs) # Calculate relative values
blank_otu <- rownames(rel_reads_OTUs[which(rel_reads_OTUs$BlankPCR > 0.05),]) # Make list with OTUs of which more than 5% present in the PCR blank 
counts <- counts[ !rownames(counts) %in% blank_otu, ] # Remove these
blank_otu2 <-rownames(rel_reads_OTUs[which(rel_reads_OTUs$BlankRedo > 0.05),]) # Make list with OTUs of which more than 5% present in the redone PCR blank
counts <- counts[ !rownames(counts) %in% blank_otu2, ] # Remove these

# Check with number of reads present in blank, OTUs with more than 100 reads in blank will be removed
rownames(counts[which(counts$BlankRedo > 100),]) #none
rownames(counts[which(counts$BlankPCR > 100),]) #none

# Remove Blank samples from dataset
counts <- subset(counts, select = -c(BlankPCR,BlankRedo))

#Remove samples with low counts, lower than 100 total read number
counts = counts[,colSums(counts) > 100]

#Remove OTUs with low counts, lower than 50 total read number
counts = counts[rowSums(counts) > 50,]

###################
# 1. Phylum level #
###################

## Make subsets for Phylum level ID
###################################

library(phyloseq)

# Order phylum count matrix
counts_phylum_1 <- as.matrix(counts)
counts_phylum_1 <- counts_phylum_1[ order((row.names(counts_phylum_1))), ]

# Order OTU identification dataset and keep only OTUs present in count matrix
otus_phylum <- OTUs[ order((row.names(OTUs))), ] # Order
OTUs_keep <-rownames(counts_phylum_1) # Make list with OTUs are present in count matrix
otus_phylum <- otus_phylum[ rownames(otus_phylum) %in% OTUs_keep, ] # Keep these OTUs in the count matrix

# Remove unassigned OTUs
unnas_OTUs <- rownames(otus_phylum[otus_phylum$Phylum == "",  ]) # Make list with OTUs that were not identified to Phylum level from OTU dataset
otus_phylum <- otus_phylum[ !rownames(otus_phylum) %in% unnas_OTUs, ] # Remove these from OTU dataset
counts_phylum_1 <- counts_phylum_1[ !rownames(counts_phylum_1) %in% unnas_OTUs, ] # Remove these from count matrix

# Make an phyloseq object otu table 
otus_phylum$Phylum <- as.factor(otus_phylum$Phylum)
counts_phylum_table <- otu_table(as.matrix(counts_phylum_1), taxa_are_rows=F) # Make otu table

# Order tables
counts_phylum_table <- counts_phylum_table[ order((row.names(counts_phylum_table))), ]
otus_phylum <- otus_phylum[ order((row.names(otus_phylum))), ]

# Merge reads of OTUs from the same Phylum 
counts_phylum_merged <- merge_samples(counts_phylum_table, t(otus_phylum$Phylum)) 
counts_phylum_merged2 <- as(otu_table(counts_phylum_merged), "matrix")

# Remove samples low counts, below 10 reads
counts_phylum_merged2 = counts_phylum_merged2[,colSums(counts_phylum_merged2) > 100] #none

# Remove samples low counts, below 10 reads
counts_phylum_merged2 = counts_phylum_merged2[rowSums(counts_phylum_merged2) > 100, ] #one

#Take relative values
counts_phylum_2 <- counts_phylum_merged2
tot_reads_samples <- colSums(counts_phylum_2) # Calculate total nr reads
counts_phylum_rel <- t(t(counts_phylum_2) / tot_reads_samples) #Calculate relative values

#Remove bacterial phyla with relative abundance <0,01
counts_phylum_rel_1 = counts_phylum_rel[rowSums(counts_phylum_rel) > 0.01,] #one

#Remove bacterial phyla removed from relative count matrix from absolute count matrix
phyla_keep <- rownames(counts_phylum_rel_1) # Make list with Phyla to keep
counts_phylum_2 <- counts_phylum_2[ rownames(counts_phylum_2) %in% phyla_keep, ] # Keep only these in count matrix

# Optimise metadata
###################

# Add "type"/category of host
metadata2 <- metadata
metadata2 <- metadata2 %>%
  unite("type", Species,Season)
metadata2 <- metadata2 %>%
  unite("type", type,Location)
metadata2 <- metadata2 %>%
  unite("type", type,Location_detailed)
metadata2 <- as.data.frame.matrix(metadata2)
metadata2$type <- as.factor(metadata2$type)

# Make metadataset with only samples occuring in counts matrix
samples_phylumtable <- colnames(counts_phylum_2) # Make list with samples of count matrix
metadata3 <- metadata2[ rownames(metadata2) %in% samples_phylumtable, ] # Keep these in metadata

# Add other variables to metadata
metadata3$Species <- with(metadata, Species[match(rownames(metadata3),rownames(metadata))])
metadata3$Season <- with(metadata, Season[match(rownames(metadata3),rownames(metadata))])
metadata3$Location <- with(metadata, Location[match(rownames(metadata3),rownames(metadata))])
metadata3$Location_detailed <- with(metadata, Location_detailed[match(rownames(metadata3),rownames(metadata))])

# Order the datasets
counts_phylum_3 <- counts_phylum_2[ ,order((colnames(counts_phylum_2))) ] # Count matrix
metadata3 <- metadata3[ order((row.names(metadata3))), ] # Metadata

# Add total number of reads per sample to metadata
reads_phyla <- as.data.frame(rowSums(t(counts_phylum_3))) # Calculate total nr reads count matrix per sample
colnames(reads_phyla) <- 'reads'
metadata3$Total_reads <- with(reads_phyla, reads[match(rownames(metadata3),rownames(reads_phyla))]) # Add these to metadata

## Stacked relative prey abundance plots
########################################

# Order samples by listing the order required
Type_ordering <- c("Limnothrissa_miodon_Dry_North_Uvira","Limnothrissa_miodon_Wet_North_Uvira", "Limnothrissa_miodon_Wet_North_Bujumbura","Limnothrissa_miodon_Wet_Central_Kalemie","Limnothrissa_miodon_Wet_South_Sumbu","Limnothrissa_miodon_Wet_South_Mpulungu","Stolothrissa_tanganicae_Dry_North_Uvira","Stolothrissa_tanganicae_Wet_North_Bujumbura","Stolothrissa_tanganicae_Wet_Central_Kalemie","Stolothrissa_tanganicae_Wet_South_Sumbu","Stolothrissa_tanganicae_Wet_South_Mpulungu") # List with order types for downstream analysis
metadata3$type <- factor(metadata3$type, levels=Type_ordering) # Add required order levels
metadata3 <- metadata3[order(metadata3$type),] # Order metadata by type as specified

# Melt the count matrix for stacked plot
library(reshape2)
counts_phylum_3_melt <- melt(counts_phylum_3, id.vars = "Sample", variable.name = "Phyla") # Melts samples for abundance bars

# Plot stacked barplot absolute values
library(ggplot2)
Order_cols = c("darkolivegreen","darkgoldenrod","lightsalmon2","chocolate4","darkcyan","lightskyblue1","dodgerblue3","firebrick4","darkolivegreen2","paleturquoise4","darkorange3","rosybrown2","darkseagreen1","darkslateblue","maroon","darkorchid4") # Define colours
ggplot(counts_phylum_3_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=Order_cols) + 
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(size = 8, angle = -80, vjust = 0.5, hjust=1))

# Order Phyla in relative count matrix from highest mean value to lowest
mean_values <- rowMeans(counts_phylum_rel_1) # Calculate mean relative values per phylum
mean_values_1 <- as.data.frame(mean_values) # As dataframe
counts_phylum_rel_1 <- as.data.frame(counts_phylum_rel_1) 
counts_phylum_rel_1$mean_values = with(mean_values_1, mean_values[match(rownames(counts_phylum_rel_1),rownames(mean_values_1))]) # Add mean values to dataset
counts_phylum_rel_1 <- counts_phylum_rel_1 %>% arrange(desc(mean_values)) # Order dataset by mean value
counts_phylum_rel_1 <- subset(counts_phylum_rel_1, select = -c(mean_values))
counts_phylum_rel_1 <- as.matrix(counts_phylum_rel_1)

# Plot stacked barplot relative values
counts_phylum_rel_melt <- melt(counts_phylum_rel_1, id.vars = "Sample", variable.name = "Phyla") # Melt relative count matrix
ggplot(counts_phylum_rel_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=Order_cols) + 
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(size = 8, angle = -80, vjust = 0.5, hjust=1))

# Add metadata to melted count matrix
counts_phylum_rel_melt$Species <-with(metadata, Species[match(counts_phylum_rel_melt$Var2,rownames(metadata))])
counts_phylum_rel_melt$Season <-with(metadata, Season[match(counts_phylum_rel_melt$Var2,rownames(metadata))])
counts_phylum_rel_melt$Location <-with(metadata, Location[match(counts_phylum_rel_melt$Var2,rownames(metadata))])
counts_phylum_rel_melt$Location_detailed <-with(metadata, Location_detailed[match(counts_phylum_rel_melt$Var2,rownames(metadata))])
counts_phylum_rel_melt$Total_reads <-with(metadata3, Total_reads[match(counts_phylum_rel_melt$Var2,rownames(metadata3))])
counts_phylum_rel_melt$type <-with(metadata3, type[match(counts_phylum_rel_melt$Var2,rownames(metadata3))])

library(tidyverse)
counts_phylum_rel_melt <- counts_phylum_rel_melt %>%
  unite("Loc_Season", Location_detailed:Season)

# Order by host species
counts_phylum_rel_melt$Species<-factor(counts_phylum_rel_melt$Species)
counts_phylum_rel_melt <- counts_phylum_rel_melt[order(counts_phylum_rel_melt$Species), ]
levels(counts_phylum_rel_melt$Species)
counts_phylum_rel_melt <- counts_phylum_rel_melt %>% arrange(Species)

# Plot relative stacked barplots per location and season
ggplot(counts_phylum_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~Loc_Season, scales=("free"),nrow=3) + 
  scale_fill_manual(values=Order_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  guides(fill=guide_legend(title="Microbiome Phyla")) +
  theme(axis.text.x = element_text(color = "black", size = 8, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Change host species names
SpeciesID=as.matrix(counts_phylum_rel_melt$Species)
names(SpeciesID)=counts_phylum_rel_melt$Var2
SpeciesID[SpeciesID == "Limnothrissa_miodon"] <- "L. miodon"
SpeciesID[SpeciesID == "Stolothrissa_tanganicae"] <- "S. tanganicae"

# Order sample order within grid
Sample_ordering <- rownames(metadata3)
counts_phylum_rel_melt$Var2 <- factor(counts_phylum_rel_melt$Var2, levels=Sample_ordering)
counts_phylum_rel_melt <- counts_phylum_rel_melt[order(counts_phylum_rel_melt$Var2),]

# Stacked barplot per type/category
ggplot(data=counts_phylum_rel_melt, aes(x=Var2, y=value, group=type)) +
  xlab("Sample") + ylab("Relative abundance") + #adds title axes
  geom_bar(aes(fill=Var1),stat = "identity") + #adds catagory TissueID by colouring bars
  theme_minimal() + 
  scale_fill_manual(values=Order_cols)+ #choose colour fill
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  ggtitle("Relative Prey abundance (phylum level)") 

# Nr reads per sample barplot
ggplot(data=metadata3, aes(x=row.names(metadata3), y=log10(Total_reads))) +
  geom_bar(aes(fill="Black"),stat = "identity") + 
  xlab("Sample") + ylab("Log10 number of reads") + #adds title axes
  theme_minimal() + 
  scale_fill_manual(values="black")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Log10 number of reads per sample")


# All relative abundance stacked barplots in different small plots
ggplot(counts_phylum_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_grid(Species~Loc_Season, scales="free_x", space="free") + 
  scale_fill_manual(values=Order_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  theme_bw() +
  guides(fill=guide_legend(title="Microbiome Phyla")) +
  theme(axis.text.x = element_text(color = "black", size = 7, face = "plain",angle=-80),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

## Bray-Curtis dissimilarity
############################

library(vegan)
library(ecodist)

# Calculate Bray-Curtis
# beta-diversity index, quantifies compositional dissimilarity between two samples or groups
dist.BrCu.all<-vegdist(t(counts_phylum_2),method="bray")
clust.BrCu.all<-hclust(dist.BrCu.all, method="ward.D") # agglomerative clustering using complete linkage

# Plot dendogram based on Bray-Curtis
plot(clust.BrCu.all)

#write csv for combined analysis
dist.BrCu.all_1 <- as.matrix(dist.BrCu.all)
write.csv(dist.BrCu.all_1, file = "16S_BrCu_table.csv")

# Make a PCoA based on Bray-curtis dissimilarities (BCD)
# Not done with count matrix since bc the data must be linear. BCD is commonly used in MB literature
endo.pco.bray<-wcmdscale(dist.BrCu.all)
plot(endo.pco.bray)

# Plot PCoA plot with ggplot
endo.pco.bray <- as.data.frame(endo.pco.bray)
ggplot(endo.pco.bray, aes(x=V1, y=V2, colour=metadata3$Location_detailed, shape=metadata3$Species)) +
  geom_point(size=5) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller

# Use BCD in Non Multidimensional scaling (NMDS) plot, works better for datasets with many zeroes
BC.nmds <- metaMDS(t(counts_phylum_rel_1), distance="bray", autotransform = F)
envfit <- envfit(BC.nmds, metadata3, permutations = 999, na.rm = TRUE) # Fits environmental vectors
sample.fit <- envfit(BC.nmds, t(counts_phylum_rel_1), permutations = 999) # Fits species vectors

# New dataset with point coordinates of NMDS results for ggplot plotting, with metadata
site.scrs <- as.data.frame(vegan::scores(BC.nmds, display = "sites"))
site.scrs$Species <- with(metadata3, Species[match(rownames(site.scrs),rownames(metadata3))]) 
site.scrs$Location_detailed <- with(metadata3, Location_detailed[match(rownames(site.scrs), rownames(metadata3))]) 
site.scrs$Season <- with(metadata3, Season[match(rownames(site.scrs), rownames(metadata3))])
site.scrs$SL <- with(metadata3, SL[match(rownames(site.scrs), rownames(metadata3))])
site.scrs$Total_reads <- with(metadata3, Total_reads[match(rownames(site.scrs), rownames(metadata3))])
site.scrs$Loc_season <- as.factor(paste(site.scrs$Location_detailed, site.scrs$Season, sep = ""))

# Make dataset with microbiome phylum vectors, to add to ggplot
spp.scrs <- as.data.frame(vegan::scores(sample.fit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) # Add phylum names
spp.scrs <- cbind(spp.scrs, pval = sample.fit$vectors$pvals) # Add pvalues
sig.spp.scrs <- subset(spp.scrs, pval<=0.01) # Shows phyla with significant values
head(spp.scrs)

# Make dataset to show extrinsic variables
env.scores<- as.data.frame(vegan::scores(envfit, display = "vectors")) # From envfit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) # Add name
env.scores <- cbind(env.scores, pval = envfit$vectors$pvals)# Add pvalue
sig.env.scrs <- subset(env.scores, pval<=0.01) # Significant
head(env.scores)

# ggplot of NMDS
Locations_shape <- c(24,23,25,6,14,2) # Define shapes
nmds.plot <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, fill = factor(Species),colour = factor(Species), shape = factor(Loc_season)), size = 4)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  scale_shape_manual(values=Locations_shape) +
  scale_fill_manual(values=c("coral4","lightblue3")) +
  scale_colour_manual(values=c("coral4","lightblue")) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Species", shape = "Tissue")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

nmds.plot + labs(title = "Basic ordination plot") #displays plot

# Add vectors significant phylum
nmds.plot+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 4, direction = "both", segment.size = 0.25,  max.overlaps=0)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Non-Metric Multidimensional Scaling")

# Check for nr reads
ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, fill = log10(Total_reads),colour=log10(Total_reads),  size = 6,shape = factor(Species)))+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

site.scrs$SL <- as.numeric(site.scrs$SL)
ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, fill = SL, colour=SL,  size = 6,shape = factor(Species)))+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

# Make heatmap based on BCD values
library(pheatmap)

BCD_matrix <- as.matrix(dist.BrCu.all)
pheatmap::pheatmap((1-BCD_matrix), symm = TRUE,clustering_method="ward.D",show_colnames=TRUE)
pheatmap::pheatmap((1-BCD_matrix), symm = TRUE,clustering_method="average",show_colnames=TRUE)

## Calculate diversity values
#############################

# Calculate Shannon diversity and add to metadata
Shannon <- as.data.frame(diversity(t(counts_phylum_2), index = "shannon")) # Calculate Shannon
Shannon$div <- diversity(t(counts_phylum_2), index = "shannon") # Make Shannon dataset
metadata3$Shannon <- with(Shannon, div[match(rownames(metadata3),rownames(Shannon))]) # Add to metadata

# Optimise metadata
metadata4 <- metadata3
metadata4 <- metadata4 %>%
  unite("LocSeas",Location,Season)
metadata3$LocSeas <- with(metadata4, LocSeas[match(rownames(metadata3),rownames(metadata4))])

ggplot(metadata3, aes(x=LocSeas, y=Shannon,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

#write csv
Shanno1 <- as.matrix(Shannon)
write.csv(Shanno1, file = "Shannon.csv")

# Calculate Chao1 diversity and add to metadata
library(fossil)
chao1div <-as.data.frame(t(estimateR(t(counts_phylum_2))))
chao1div$Chao1 <- as.numeric(chao1div[,2])
metadata4$Chao1 <- with(chao1div, Chao1[match(rownames(metadata3),rownames(chao1div))])

ggplot(metadata4, aes(x=LocSeas, y=Chao1,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

## Look for host factors influencing microbiome composition
###########################################################

# PERMANOVA based on BCD
# All samples
BCD_PERMA<-adonis2(dist.BrCu.all~Species+Location_detailed+Season, data=metadata3, permutations = 999, method="bray") 
adonis2(dist.BrCu.all~Species, data=metadata3, permutations = 10000, p.adjust.m = "bonferroni")
adonis2(dist.BrCu.all~Species+Location_detailed+Season, data=metadata3, permutations = 10000, p.adjust.m = "bonferroni")

# PERMANOVA SL and Total reads
adonis2(dist.BrCu.all~Total_reads, data=metadata3, permutations = 10000, p.adjust.m = "bonferroni")

no_SL_samples <- rownames(metadata3[which(is.na(metadata3$SL)),])
metadata_SL <- metadata3[ !rownames(metadata3) %in% no_SL_samples, ]
counts_SL <- counts_phylum_2[ ,!colnames(counts_phylum_2) %in% no_SL_samples ]
BrCu.SL<-vegdist(t(counts_SL),method="bray")
adonis2(BrCu.SL~SL, data=metadata_SL, permutations = 10000, p.adjust.m = "bonferroni")

# PERMANOVA Samples Uvira
Uvira_samples <- rownames(metadata3[which(metadata3$Location_detailed == "Uvira"),]) 
metadata_uvira <- metadata3[ rownames(metadata3) %in% Uvira_samples, ] 
counts_uvira <- counts_phylum_2[ ,colnames(counts_phylum_2) %in% Uvira_samples ]

BrCu.uvira<-vegdist(t(counts_uvira),method="bray")
adonis2(BrCu.uvira~Species*Season, data=metadata_uvira, permutations = 10000, p.adjust.m = "bonferroni")

# Uvira L.miodon
Uvira_LM_samples <- rownames(metadata_uvira[which(metadata_uvira$Species == "L_miodon"),]) 
metadata_uvira_LM <- metadata_uvira[ rownames(metadata_uvira) %in% Uvira_LM_samples, ] 
counts_uvira_LM <- counts_uvira[ ,colnames(counts_uvira) %in% Uvira_LM_samples ]

BrCu.uvira_LM<-vegdist(t(counts_uvira_LM),method="bray")
adonis2(BrCu.uvira_LM~Season, data=metadata_uvira_LM, permutations = 10000, p.adjust.m = "bonferroni")

# Uvira S.tanganicae not for microbiome
Uvira_ST_samples <- rownames(metadata_uvira[which(metadata_uvira$Species == "Stolothrissa_tanganicae"),]) 
metadata_uvira_ST <- metadata_uvira[ rownames(metadata_uvira) %in% Uvira_ST_samples, ] 
counts_uvira_ST <- counts_uvira[ ,colnames(counts_uvira) %in% Uvira_ST_samples ]

BrCu.uvira_ST<-vegdist(t(counts_uvira_ST),method="bray")
adonis2(BrCu.uvira_ST~Season, data=metadata_uvira_ST, permutations = 10000, p.adjust.m = "bonferroni")

# PERMANOVA without Uvira Wet
Uvira_Dry_samples <- rownames(metadata_uvira[which(metadata_uvira$Season == "Dry"),]) 
counts_wet <- counts_phylum_2[ ,! colnames(counts_phylum_2) %in% Uvira_Dry_samples ] 
metadata_wet <- metadata3[ ! rownames(metadata3) %in% Uvira_Dry_samples, ]

BrCu.wet<-vegdist(t(counts_wet),method="bray")
adonis2(BrCu.wet~Location_detailed+Species, data=metadata_wet, permutations = 999, method="bray") 
adonis2(BrCu.wet~Species+Location_detailed, data=metadata_wet, permutations = 999, method="bray")

adonis2(BrCu.wet~Species*Location_detailed, data=metadata_wet, permutations = 10000, p.adjust.m = "bonferroni")

# Wet L.miodon
Wet_LM_samples <- rownames(metadata_wet[which(metadata_wet$Species == "L_miodon"),]) 
metadata_wet_LM <- metadata_wet[ rownames(metadata_wet) %in% Wet_LM_samples, ] 
counts_wet_LM <- counts_wet[ ,colnames(counts_wet) %in% Wet_LM_samples ]

BrCu.wet_LM<-vegdist(t(counts_wet_LM),method="bray")
adonis2(BrCu.wet_LM~Location_detailed, data=metadata_wet_LM,permutations = 10000, p.adjust.m = "bonferroni")

# Wet S.tanganicae
Wet_ST_samples <- rownames(metadata_wet[which(metadata_wet$Species == "S_tanganicae"),]) 
metadata_wet_ST <- metadata_wet[ rownames(metadata_wet) %in% Wet_ST_samples, ] 
counts_wet_ST <- counts_wet[ ,colnames(counts_wet) %in% Wet_ST_samples ]

BrCu.wet_ST<-vegdist(t(counts_wet_ST),method="bray")
adonis2(BrCu.wet_ST~Location_detailed, data=metadata_wet_ST, permutations = 10000, p.adjust.m = "bonferroni")


## Additional general calculations
##################################

# Percentage abundance prey per sardine species
LM_samples <- rownames(metadata3[which(metadata3$Species == "L_miodon"),])
counts_LM <- counts_phylum_rel_1[ ,colnames(counts_phylum_rel_1) %in% LM_samples ]
LM_taxa <- as.matrix(rowSums(counts_LM))
LM_taxa_rel <- LM_taxa/colSums(LM_taxa)

ST_samples <- rownames(metadata3[which(metadata3$Species == "S_tanganicae"),])
counts_ST <- counts_phylum_rel_1[ ,colnames(counts_phylum_rel_1) %in% ST_samples ]
ST_taxa <- as.matrix(rowSums(counts_ST))
ST_taxa_rel <- ST_taxa/colSums(ST_taxa)

# Nr Reads
mean(metadata3$Total_reads)
median(metadata3$Total_reads)
min(metadata3$Total_reads)
sd(metadata3$Total_reads)
max(metadata3$Total_reads)

# Mean diversity
metadata_ST <- metadata3[ rownames(metadata3) %in% ST_samples , ]
mean(metadata_ST$Shannon)
mean(metadata_ST$Chao1)

metadata_LM <- metadata3[ rownames(metadata3) %in% LM_samples , ]
mean(metadata_LM$Shannon)
mean(metadata_LM$Chao1)

