################################################################################
# Clupeid Lake Tanganyika COI Prey Analysis Huyghe & Aerts et al. (under review)
################################################################################

# Set working directory
setwd("Working_directory") 

# Upload OTU classifications (not needed for MicroBiome MB)
OTUs <-read.csv(file="Sardine_selected_OTUs_07102022.csv",sep=",",check.names=FALSE,row.names = 1)

# Upload count matrix 
counts0 <-read.csv(file="Sardine_complete_count_matrix_07102022.csv",sep=",",check.names=FALSE,row.names = 1)

# Upload metadata of host individuals
metadata <-read.csv(file="COI_metadata_R_all_18042023.csv",sep=",",check.names=FALSE,row.names = 1)

# Keep only OTU's that were identified by BOLDdigger in count matrix, the ones that are present in "OTUs" dataset
otus_to_keep <- rownames(OTUs)
counts <- counts0[ rownames(counts0) %in% otus_to_keep, ]

# Remove OTUs present in the blanks
# Blanks are B_R, B_1 and B_2

#See whether any of the OTU's are present in the blank
blank_otu <- rownames(counts[which(counts$B_R > 0),]) #one
counts <- counts[ !rownames(counts) %in% blank_otu, ]
OTUs <- OTUs[ !rownames(OTUs) %in% blank_otu, ]
rownames(counts[which(counts$B_1 > 0),]) #none
rownames(counts[which(counts$B_2 > 0),]) #none

# Remove Blanks
counts <- subset(counts, select = -c(B_1,B_2,B_R))

#Remove samples with low counts, below 10 reads
counts = counts[,colSums(counts) > 10]

################################################################################
# 1. Order level #
################################################################################

# Make subsets with OTU's identified to order level
library(phyloseq)
counts_order_1 <- as.matrix(counts)
counts_order_1 <- counts_order_1[ order((row.names(counts_order_1))), ]
otus_order <- OTUs[ order((row.names(OTUs))), ]

# Remove OTUs which were not identified till order level
noorder_otu <- rownames(otus_order[which(otus_order$Similarity < 85),]) # 85 comes from the similarity threshold for order by Bolddigger
otus_order <- otus_order[ !rownames(otus_order) %in% noorder_otu, ] # Remove these from dataset
counts_order_1 <- counts_order_1[ !rownames(counts_order_1) %in% noorder_otu, ]
# Make subset of samples that have OTU's that were identified to order level
rmSamplesOrder <- colnames(counts_order_1[,colSums(counts_order_1) == 0]) 
counts_order_1 <- counts_order_1[ ,!colnames(counts_order_1) %in% rmSamplesOrder ]

## Merge abundances of same orders
##################################

# Make an phyloseq object otu table 
otus_order$Order <- as.factor(otus_order$Order)
counts_order_table <- otu_table(as.matrix(counts_order_1), taxa_are_rows=F)

# Order tables
counts_order_table <- counts_order_table[ order((row.names(counts_order_table))), ]
otus_order <- otus_order[ order((row.names(otus_order))), ]

# Merge reads of OTUs from the same order 
counts_order_merged <- merge_samples(counts_order_table, t(otus_order$Order))
counts_order_merged2 <- as(otu_table(counts_order_merged), "matrix")

# Remove samples low counts, below 10 reads
counts_order_merged2 = counts_order_merged2[,colSums(counts_order_merged2) > 10]

# Remove Orders with low count, below 5 reads
counts_order_merged2 = counts_order_merged2[rowSums(counts_order_merged2) > 5,]

#Take relative values per sample
tot_reads_samples <- colSums(counts_order_merged2)
order_samples <- colnames(counts_order_merged2)

counts_order_merged_rel <- t(t(counts_order_merged2) / tot_reads_samples)

## Optimise metadata for downstream analysis
############################################

library(dplyr)
library(tidyr)
metadata2 <- metadata

# Make category "type" (categories) where species, season, location and location_detailed or combined
metadata2 <- metadata2 %>%
  unite("type", Species,Season)
metadata2 <- metadata2 %>%
  unite("type", type,Location)
metadata2 <- metadata2 %>%
  unite("type", type,Location_detailed)
metadata2 <- as.data.frame.matrix(metadata2)
metadata2$type <- as.factor(metadata2$type)

# Order the metadata and count matrix
counts_order_merged2 <- counts_order_merged2[ order((row.names(counts_order_merged2))), ]
counts_order_merged2 <- counts_order_merged2[ ,order((colnames(counts_order_merged2))) ]
metadata2 <- metadata2[ order((row.names(metadata2))), ]
order_sample_table <- otu_table(t(counts_order_merged2), taxa_are_rows=F)

# Make metadata with only the samples that are present in our order count matrix
samples_ordertable <- colnames(counts_order_merged2)
metadata3 <- metadata2[ rownames(metadata2) %in% samples_ordertable, ]
metadata3$Species <- with(metadata, Species[match(rownames(metadata3),rownames(metadata))])
metadata3$Season <- with(metadata, Season[match(rownames(metadata3),rownames(metadata))])
metadata3$Location <- with(metadata, Location[match(rownames(metadata3),rownames(metadata))])
metadata3$Location_detailed <- with(metadata, Location_detailed[match(rownames(metadata3),rownames(metadata))])
metadata3 <- metadata3[ order((metadata3$type)), ]

# Add number of reads per sample to metadata
reads_order <- as.data.frame(rowSums(t(counts_order_merged2)))
colnames(reads_order) <- 'reads'
metadata3$Total_reads <- with(reads_order, reads[match(rownames(metadata3),rownames(reads_order))])

## Stacked relative prey abundance barplots
###########################################

# Order samples metadata for plot
Type_ordering <- c("Limnothrissa_miodon_Dry_North_Uvira","Limnothrissa_miodon_Wet_North_Uvira","Limnothrissa_miodon_Wet_North_Bujumbura","Limnothrissa_miodon_Wet_Center_Kalemie","Limnothrissa_miodon_Wet_South_Sumbu","Limnothrissa_miodon_Wet_South_Mpulungu","Stolothrissa_tanganicae_Dry_North_Uvira","Stolothrissa_tanganicae_Wet_North_Uvira","Stolothrissa_tanganicae_Wet_Center_Kalemie","Stolothrissa_tanganicae_Wet_South_Sumbu","Stolothrissa_tanganicae_Wet_South_Mpulungu")
metadata3$type <- factor(metadata3$type, levels=Type_ordering)
metadata3 <- metadata3[order(metadata3$type),]

# Order prey count matrix for plot
Prey_ordering <- c("Cichliformes","Perciformes","Cyprinodontiformes","Siluriformes","Decapoda","Calanoida","Cyclopoida","Diptera","Lepidoptera","Ephemeroptera","Limnomedusae")
counts_order_merged3 <- as.data.frame(counts_order_merged2)
counts_order_merged3$Prey <- factor(rownames(counts_order_merged3))
counts_order_merged3$Prey <- factor(counts_order_merged3$Prey, levels=Prey_ordering)
counts_order_merged3 <- counts_order_merged3[order(counts_order_merged3$Prey),]
counts_order_merged3 <- subset(counts_order_merged3, select = -Prey )

# Order prey relative count matrix for plot
counts_order_merged_rel2 <- as.data.frame(counts_order_merged_rel)
counts_order_merged_rel2$Prey <- factor(rownames(counts_order_merged_rel2), levels=Prey_ordering)
counts_order_merged_rel2 <- counts_order_merged_rel2[order(counts_order_merged_rel2$Prey),]
counts_order_merged_rel2 <- subset(counts_order_merged_rel2, select = -Prey )

# Melt the count matrix for stacked plot
library(reshape2)
counts_order_merged2_melt <- melt(counts_order_merged2, id.vars = "Sample", variable.name = "Phyla") # Melts samples for abundance bars

# Plot stacked barplot with absolute number of reads
library(ggplot2)
Order_cols = c("chocolate","cyan4","darkorange","dodgerblue4","firebrick","darkolivegreen3","darkolivegreen1","darkseagreen3","lightsteelblue1","lightskyblue4","paleturquoise3")

ggplot(counts_order_merged2_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=Order_cols) + 
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Melt the relative count matrix for stacked plot and order the prey as specified before
counts_order_merged_rel_melt <- melt(counts_order_merged_rel, id.vars = "Sample", variable.name = "Phyla")
counts_order_merged_rel_melt$Var1 <- factor(counts_order_merged_rel_melt$Var1, levels=Prey_ordering)
counts_order_merged_rel_melt <- counts_order_merged_rel_melt[order(counts_order_merged_rel_melt$Var1),]

# Plot stacked barplot with relative number of reads
Order_cols = c("cyan4","lightskyblue4","dodgerblue4","paleturquoise3","firebrick","chocolate","darkorange","darkolivegreen3","darkolivegreen4","darkolivegreen1","lightsteelblue1") # Other order because of different order taxa
ggplot(counts_order_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=Order_cols) + 
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Make grid for each host category
# Add metadata to merged melted
counts_order_merged_rel_melt$Species <-with(metadata, Species[match(counts_order_merged_rel_melt$Var2,rownames(metadata))])
counts_order_merged_rel_melt$Season <-with(metadata, Season[match(counts_order_merged_rel_melt$Var2,rownames(metadata))])
counts_order_merged_rel_melt$Location <-with(metadata, Location[match(counts_order_merged_rel_melt$Var2,rownames(metadata))])
counts_order_merged_rel_melt$Location_detailed <-with(metadata, Location_detailed[match(counts_order_merged_rel_melt$Var2,rownames(metadata))])
counts_order_merged_rel_melt$Total_reads <-with(metadata3, Total_reads[match(counts_order_merged_rel_melt$Var2,rownames(metadata3))])
counts_order_merged_rel_melt$type <-with(metadata3, type[match(counts_order_merged_rel_melt$Var2,rownames(metadata3))])
library(tidyverse)
counts_order_merged_rel_melt <- counts_order_merged_rel_melt %>%
  unite("Loc_Season", Location_detailed:Season)

# Order by species
counts_order_merged_rel_melt$Species<-factor(counts_order_merged_rel_melt$Species)
counts_order_merged_rel_melt <- counts_order_merged_rel_melt[order(counts_order_merged_rel_melt$Species), ]
levels(counts_order_merged_rel_melt$Species)
counts_order_merged_rel_melt <- counts_order_merged_rel_melt %>% arrange(Species)

# Relative value stacked barplot per detailed location and season
ggplot(counts_order_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~Loc_Season, scales=("free"),nrow=3) + 
  scale_fill_manual(values=Order_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  guides(fill=guide_legend(title="Prey order")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Change species names for plot
SpeciesID=as.matrix(counts_order_merged_rel_melt$Species)
names(SpeciesID)=counts_order_merged_rel_melt$Var2
SpeciesID[SpeciesID == "Limnothrissa_miodon"] <- "L. miodon"
SpeciesID[SpeciesID == "Stolothrissa_tanganicae"] <- "S. tanganicae"

# Order by host category
Sample_ordering <- rownames(metadata3)
counts_order_merged_rel_melt$Var2 <- factor(counts_order_merged_rel_melt$Var2, levels=Sample_ordering)
counts_order_merged_rel_melt <- counts_order_merged_rel_melt[order(counts_order_merged_rel_melt$Var2),]

# All types arranged in different small plots
ggplot(counts_order_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_grid(Species~Loc_Season, scales="free_x", space="free") + 
  scale_fill_manual(values=Order_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  theme_bw() +
  guides(fill=guide_legend(title="Prey order")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Stacked barplot arranged per host category (type)
ggplot(data=counts_order_merged_rel_melt, aes(x=Var2, y=value, group=type)) +
  xlab("Sample") + ylab("Relative abundance") + #adds title axes
  scale_fill_manual(values=Order_cols) +
  geom_bar(aes(fill=Var1),stat = "identity") + #adds catagory TissueID by colouring bars
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  theme_classic() +
  ggtitle("Relative Prey abundance (order level)") 

# Nr reads per sample barplot to add
ggplot(data=metadata3, aes(x=row.names(metadata3), y=log10(Total_reads))) +
  geom_bar(aes(fill="Black"),stat = "identity") + 
  xlab("Sample") + ylab("Log10 number of reads") + #adds title axes
  theme_classic() +
  scale_fill_manual(values="black")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  ggtitle("Log10 number of reads per sample") 


## Bray-Curtis dissimilarity
############################

library(vegan)
library(ecodist)

# Calculate Bray-Curtis
dist.BrCu.all<-vegdist(t(counts_order_merged2),method="bray")
clust.BrCu.all<-hclust(dist.BrCu.all, method="ward.D") #agglomerative clustering using complete linkage

# Plot dendogram
plot(clust.BrCu.all)

# Make a PCoA based on Bray-curtis dissimilarities (BCD)
# Not done with count matrix since because the data must be linear.
endo.pco.bray<-wcmdscale(dist.BrCu.all)
plot(endo.pco.bray)

# Plot with ggplot
endo.pco.bray <- as.data.frame(endo.pco.bray)
ggplot(endo.pco.bray, aes(x=V1, y=V2, colour=metadata3$Location_detailed, shape=metadata3$Species)) +
  geom_point(size=5) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller

# Use BCD in Non Multidimensional scaling (NMDS) plot, works better for datasets with many missing values
BC.nmds <- metaMDS(t(counts_order_merged_rel), distance="bray", autotransform = F)
envfit <- envfit(BC.nmds, metadata3, permutations = 999,na.rm = TRUE) # fits environmental vectors
sample.fit <- envfit(BC.nmds, t(counts_order_merged_rel), permutations = 999) # fits species vectors

# New dataset with point coordinates of NMDS results for ggplot plotting
site.scrs <- as.data.frame(vegan::scores(BC.nmds, display = "sites"))

# Add metadata
site.scrs$Species <-with(metadata3, Species[match(rownames(site.scrs),rownames(metadata3))])
site.scrs$Location_detailed <-with(metadata3, Location_detailed[match(rownames(site.scrs),rownames(metadata3))])
site.scrs$Season <-with(metadata3, Season[match(rownames(site.scrs),rownames(metadata3))])
site.scrs$SL <-with(metadata3, SL[match(rownames(site.scrs),rownames(metadata3))])
site.scrs$Total_reads <-with(metadata3, Total_reads[match(rownames(site.scrs),rownames(metadata3))])
site.scrs$Loc_season <- as.factor(paste(site.scrs$Location_detailed, site.scrs$Season, sep = ""))
head(site.scrs)

# Make dataset with prey vectors, to add to ggplot
spp.scrs <- as.data.frame(vegan::scores(sample.fit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) # add prey names
spp.scrs <- cbind(spp.scrs, pval = sample.fit$vectors$pvals) # add pvalues
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) # shows prey with significant values
head(spp.scrs)

# Make dataset to show extrinsic variables
env.scores<- as.data.frame(vegan::scores(envfit, display = "vectors")) # From envfit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) # Add name
env.scores <- cbind(env.scores, pval = envfit$vectors$pvals) # Add pvalues
sig.env.scrs <- subset(env.scores, pval<=0.05) # significant
head(env.scores)

# ggplot of NMDS
Locations_shape <- c(24,23,25,6,14,2) # specify shapes plot points
nmds.plot <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = factor(Species),colour = factor(Species), shape = factor(Loc_season)), size = 6)+ 
  coord_fixed()+
  scale_shape_manual(values=Locations_shape) +
  scale_fill_manual(values=c("coral4","lightblue3")) +
  scale_colour_manual(values=c("coral4","lightblue3")) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Species", shape = "Tissue")+ 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))

nmds.plot + labs(title = "Basic ordination plot") 

# Add vectors significant prey
nmds.plot+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + # add vectors
 labs(title = "Non-Metric Multidimensional Scaling")

# Make same plot with log10 total nr of reads to check for correlation
ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = log10(Total_reads),colour=log10(Total_reads),  size = 6))+
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

# Make same plot with standard length to check for correlation
ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = SL,colour=SL,  size = 6,shape = factor(Species)))+ 
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

# Make heatmap based on BCD values
library(pheatmap)

BCD_matrix <- as.matrix(dist.BrCu.all)
pheatmap::pheatmap((1-BCD_matrix), symm = TRUE,clustering_method="ward.D",show_colnames=TRUE)

## Calculate diversity values
#############################

# Calculate Shannon diversity and add to metadata
diversity(counts_order_merged2, index = "shannon") 
Shannon <- as.data.frame(diversity(t(counts_order_merged2), index = "shannon"))
Shannon$div <- diversity(t(counts_order_merged2), index = "shannon")
metadata$Shannon <- with(Shannon, div[match(rownames(metadata),rownames(Shannon))])
metadata3$Shannon <- with(Shannon, div[match(rownames(metadata3),rownames(Shannon))])

# Make new metadata
metadata4 <- metadata3
metadata4 <- metadata4 %>%
  unite("LocSeas",Location,Season)
metadata3$LocSeas <- with(metadata4, LocSeas[match(rownames(metadata3),rownames(metadata4))])

# Plot diversity barplots
ggplot(metadata3, aes(x=LocSeas, y=Shannon,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

# Calculate Chao1 diversity and add to metadata
chao1div <-as.data.frame(t(estimateR(t(counts_order_merged2))))
chao1div$Chao1 <- as.numeric(chao1div[,2])
metadata3$Chao1 <- with(chao1div, Chao1[match(rownames(metadata3),rownames(chao1div))])

# Plot diversity barplots
ggplot(metadata3, aes(x=LocSeas, y=Chao1,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

## Look for host factors influencing prey composition
#####################################################

# PERMANOVA based on BCD
adonis2(dist.BrCu.all~Species, data=metadata3, permutations = 999, p.adjust.m = "bonferroni")
adonis2(dist.BrCu.all~Location_detailed, data=metadata3, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA Total reads
adonis2(dist.BrCu.all~Total_reads, data=metadata3, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA Standard length
no_SL_samples <- rownames(metadata3[which(is.na(metadata3$SL)),])
metadata_SL <- metadata3[ !rownames(metadata3) %in% no_SL_samples, ]
counts_SL <- counts_order_merged2[ ,!colnames(counts_order_merged2) %in% no_SL_samples ]
BrCu.SL<-vegdist(t(counts_SL),method="bray")

adonis2(BrCu.SL~SL, data=metadata_SL, permutations = 999, p.adjust.m = "bonferroni")

# SL Limothrissa miodon
LM_SL_samples <- rownames(metadata_SL[which(metadata_SL$Species == "Limnothrissa_miodon"),])
metadata_SL_LM <- metadata_SL[ rownames(metadata_SL) %in% LM_SL_samples, ]
counts_SL_LM <- counts_SL[ ,colnames(counts_SL) %in% LM_SL_samples ]
BrCu.SL_LM<-vegdist(t(counts_SL_LM),method="bray")

adonis2(BrCu.SL_LM~SL, data=metadata_SL_LM, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.SL_LM~Total_reads, data=metadata_SL_LM, permutations = 999, p.adjust.m = "bonferroni")

# SL Stolothrissa tanganicae
ST_SL_samples <- rownames(metadata_SL[which(metadata_SL$Species == "Stolothrissa_tanganicae"),])
metadata_SL_ST <- metadata_SL[ rownames(metadata_SL) %in% ST_SL_samples, ]
counts_SL_ST <- counts_SL[ ,colnames(counts_SL) %in% ST_SL_samples ]
BrCu.SL_ST<-vegdist(t(counts_SL_ST),method="bray")

adonis2(BrCu.SL_ST~SL, data=metadata_SL_ST, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.SL_ST~Total_reads, data=metadata_SL_ST, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA Samples Uvira
Uvira_samples <- rownames(metadata3[which(metadata3$Location_detailed == "Uvira"),])
metadata_uvira <- metadata3[ rownames(metadata3) %in% Uvira_samples, ]
counts_uvira <- counts_order_merged2[ ,colnames(counts_order_merged2) %in% Uvira_samples ]
BrCu.uvira<-vegdist(t(counts_uvira),method="bray")

adonis2(BrCu.uvira~Species, data=metadata_uvira, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.uvira~Season, data=metadata_uvira, permutations = 999, p.adjust.m = "bonferroni")

# Uvira L.miodon
Uvira_LM_samples <- rownames(metadata_uvira[which(metadata_uvira$Species == "Limnothrissa_miodon"),])
metadata_uvira_LM <- metadata_uvira[ rownames(metadata_uvira) %in% Uvira_LM_samples, ]
counts_uvira_LM <- counts_uvira[ ,colnames(counts_uvira) %in% Uvira_LM_samples ]
BrCu.uvira_LM<-vegdist(t(counts_uvira_LM),method="bray")

adonis2(BrCu.uvira_LM~Season, data=metadata_uvira_LM, permutations = 999, p.adjust.m = "bonferroni")

# Uvira S.tanganicae
Uvira_ST_samples <- rownames(metadata_uvira[which(metadata_uvira$Species == "Stolothrissa_tanganicae"),])
metadata_uvira_ST <- metadata_uvira[ rownames(metadata_uvira) %in% Uvira_ST_samples, ]
counts_uvira_ST <- counts_uvira[ ,colnames(counts_uvira) %in% Uvira_ST_samples ]
BrCu.uvira_ST<-vegdist(t(counts_uvira_ST),method="bray")

adonis2(BrCu.uvira_ST~Season, data=metadata_uvira_ST, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA without Uvira Dry
Uvira_Dry_samples <- rownames(metadata_uvira[which(metadata_uvira$Season == "Dry"),])
counts_wet <- counts_order_merged2[ ,! colnames(counts_order_merged2) %in% Uvira_Dry_samples ]
metadata_wet <- metadata3[ ! rownames(metadata3) %in% Uvira_Dry_samples, ]
BrCu.wet<-vegdist(t(counts_wet),method="bray")

adonis2(BrCu.wet~Location_detailed, data=metadata_wet, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.wet~Location, data=metadata_wet, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.wet~Species, data=metadata_wet, permutations = 999, p.adjust.m = "bonferroni")

# Wet L.miodon
Wet_LM_samples <- rownames(metadata_wet[which(metadata_wet$Species == "Limnothrissa_miodon"),])
metadata_wet_LM <- metadata_wet[ rownames(metadata_wet) %in% Wet_LM_samples, ]
counts_wet_LM <- counts_wet[ ,colnames(counts_wet) %in% Wet_LM_samples ]
BrCu.wet_LM<-vegdist(t(counts_wet_LM),method="bray")

adonis2(BrCu.wet_LM~Location_detailed, data=metadata_wet_LM,permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.wet_LM~Location, data=metadata_wet_LM,permutations = 999, p.adjust.m = "bonferroni")

# Wet L.miodon remove low sample locations
Wet_LM_low_samples <- c("B1_L_8","S1_S_9")
metadata_wet_LM_2 <- metadata_wet_LM[ !rownames(metadata_wet_LM) %in% Wet_LM_low_samples, ]
counts_wet_LM_2 <- counts_wet_LM[ ,!colnames(counts_wet_LM) %in% Wet_LM_low_samples ]
BrCu.wet_LM_2<-vegdist(t(counts_wet_LM_2),method="bray")

adonis2(BrCu.wet_LM_2~Location_detailed, data=metadata_wet_LM_2,permutations = 999, p.adjust.m = "bonferroni")

# Wet S.tanganicae
Wet_ST_samples <- rownames(metadata_wet[which(metadata_wet$Species == "Stolothrissa_tanganicae"),])
metadata_wet_ST <- metadata_wet[ rownames(metadata_wet) %in% Wet_ST_samples, ]
counts_wet_ST <- counts_wet[ ,colnames(counts_wet) %in% Wet_ST_samples ]
BrCu.wet_ST<-vegdist(t(counts_wet_ST),method="bray")

adonis2(BrCu.wet_ST~Location_detailed, data=metadata_wet_ST, permutations = 999, p.adjust.m = "bonferroni")

## Additional general calculations
##################################

# Percentage abundance prey per sardine species
LM_samples <- rownames(metadata3[which(metadata3$Species == "Limnothrissa_miodon"),])
counts_LM <- counts_order_merged_rel[ ,colnames(counts_order_merged_rel) %in% LM_samples ]
LM_taxa <- as.matrix(rowSums(counts_LM))
LM_taxa_rel <- LM_taxa/colSums(LM_taxa)
LM_taxa_rel

ST_samples <- rownames(metadata3[which(metadata3$Species == "Stolothrissa_tanganicae"),])
counts_ST <- counts_order_merged_rel[ ,colnames(counts_order_merged_rel) %in% ST_samples ]
ST_taxa <- as.matrix(rowSums(counts_ST))
ST_taxa_rel <- ST_taxa/colSums(ST_taxa)
ST_taxa_rel

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

## MANTEL Test with 16S
#######################

BrCu_16S0 <-read.csv(file="16S_BrCu_table.csv",sep=",",check.names=FALSE)
rownames(BrCu_16S0) <-  BrCu_16S0[,1]
BrCu_16S <- subset(BrCu_16S0, select = -c(1) )

BrCu_16S_samples <- row.names(BrCu_16S)
BrCu_COI_samples <- row.names(BCD_matrix)

Samples_mantel <- BrCu_COI_samples[BrCu_COI_samples %in% BrCu_16S_samples]

BrCu_16S_mantel <- BrCu_16S[ rownames(BrCu_16S) %in% Samples_mantel, colnames(BrCu_16S) %in% Samples_mantel]
BrCu_COI_mantel <- BCD_matrix[ rownames(BCD_matrix) %in% Samples_mantel, colnames(BCD_matrix) %in% Samples_mantel]

BrCu_COI_mantel_dist <- as.dist(BrCu_COI_mantel)
BrCu_16S_mantel_dist <- as.dist(BrCu_16S_mantel)

library(ade4)
mantel.rtest(BrCu_COI_mantel_dist, BrCu_16S_mantel_dist, nrepet = 9999)
plot(BrCu_COI_mantel_dist,  BrCu_16S_mantel_dist, ylab="Distance matrix of the microbiome composition(16S rRNA)", xlab="Distance matrix of the prey item composition (COI)")
abline(lm(BrCu_16S_mantel_dist ~ BrCu_COI_mantel_dist))
reg <- (lm(BrCu_16S_mantel_dist ~ BrCu_COI_mantel_dist))
coeff=coefficients(reg)
eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))

# Mantel L.miodon
LM_samples
Samples_mantel_LM <- Samples_mantel[Samples_mantel %in% LM_samples]

BrCu_16S_mantel_LM <- BrCu_16S[ rownames(BrCu_16S) %in% Samples_mantel_LM, colnames(BrCu_16S) %in% Samples_mantel_LM]
BrCu_COI_mantel_LM <- BCD_matrix[ rownames(BCD_matrix) %in% Samples_mantel_LM, colnames(BCD_matrix) %in% Samples_mantel_LM]

BrCu_COI_mantel_LM_dist <- as.dist(BrCu_COI_mantel_LM)
BrCu_16S_mantel_LM_dist <- as.dist(BrCu_16S_mantel_LM)

mantel.rtest(BrCu_COI_mantel_LM_dist, BrCu_16S_mantel_LM_dist, nrepet = 9999)
plot(BrCu_COI_mantel_LM_dist,  BrCu_16S_mantel_LM_dist, ylab="Distance matrix of the microbiome composition(16S rRNA)", xlab="Distance matrix of the prey item composition (COI)")
abline(lm(BrCu_16S_mantel_LM_dist ~ BrCu_COI_mantel_LM_dist))
reg <- (lm(BrCu_16S_mantel_LM_dist ~ BrCu_COI_mantel_LM_dist))
coeff=coefficients(reg)
eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))

# Mantel S.tanganicae
ST_samples
Samples_mantel_ST <- Samples_mantel[Samples_mantel %in% ST_samples]

BrCu_16S_mantel_ST <- BrCu_16S[ rownames(BrCu_16S) %in% Samples_mantel_ST, colnames(BrCu_16S) %in% Samples_mantel_ST]
BrCu_COI_mantel_ST <- BCD_matrix[ rownames(BCD_matrix) %in% Samples_mantel_ST, colnames(BCD_matrix) %in% Samples_mantel_ST]

BrCu_COI_mantel_ST_dist <- as.dist(BrCu_COI_mantel_ST)
BrCu_16S_mantel_ST_dist <- as.dist(BrCu_16S_mantel_ST)

mantel.rtest(BrCu_COI_mantel_ST_dist, BrCu_16S_mantel_ST_dist, nrepet = 9999)
plot(BrCu_COI_mantel_ST_dist,  BrCu_16S_mantel_ST_dist, ylab="Distance matrix of the microbiome composition(16S rRNA)", xlab="Distance matrix of the prey item composition (COI)")
abline(lm(BrCu_16S_mantel_ST_dist ~ BrCu_COI_mantel_ST_dist))
reg <- (lm(BrCu_16S_mantel_ST_dist ~ BrCu_COI_mantel_ST_dist))
coeff=coefficients(reg)
eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))



################################################################################
# 2. Class #
################################################################################

# Make subsets for order level ID
counts_class_1 <- as.matrix(counts)
counts_class_1 <- counts_class_1[ order((row.names(counts_class_1))), ]
otus_class <- OTUs[ order((row.names(OTUs))), ]

# Remove OTU's which were not identified till family level
noclass_otu <- rownames(otus_class[which(otus_class$Class == ""),])
otus_class <- otus_class[ !rownames(otus_class) %in% noclass_otu, ]
counts_class_1 <- counts_class_1[ !rownames(counts_class_1) %in% noclass_otu, ]
rmSamplesclass <- colnames(counts_class_1[,colSums(counts_class_1) == 0])
counts_class_1 <- counts_class_1[ ,!colnames(counts_class_1) %in% rmSamplesclass ]

# Merge samples
otus_class$Class <- as.factor(otus_class$Class)
counts_class_table <- otu_table(as.matrix(counts_class_1), taxa_are_rows=F)
counts_class_merged <- merge_samples(counts_class_table, t(otus_class$Class))
counts_class_merged2 <- as(otu_table(counts_class_merged), "matrix")

# Remove samples and Classes with low counts
counts_class_merged2 = counts_class_merged2[,colSums(counts_class_merged2) > 10]
counts_class_merged2 = counts_class_merged2[rowSums(counts_class_merged2) > 10,]

# Remove these samples from metadata
samples_Class_keep <- colnames(counts_class_merged2)
metadata_class0 <- metadata
metadata_class0 <- metadata_class0[ rownames(metadata_class0) %in% samples_Class_keep, ]

#Take relative values
tot_reads_class_samples <- colSums(counts_class_merged2)
class_samples <- colnames(counts_class_merged2)
counts_class_merged_rel <- t(t(counts_class_merged2) / tot_reads_class_samples)

# Total counts per sample
reads_class <- as.data.frame(rowSums(t(counts_class_merged2)))
colnames(reads_class) <- 'reads'
metadata_class0$Total_reads <- with(reads_class, reads[match(rownames(metadata_class0),rownames(reads_class))])

#Take relative values per sample
tot_reads_class_samples <- colSums(counts_class_merged2)
class_samples <- colnames(counts_class_merged2)
counts_class_merged_rel <- t(t(counts_class_merged2) / tot_reads_class_samples)

## Optimise metadata for downstream analysis
############################################

library(dplyr)
library(tidyr)
metadata_class <- metadata_class0

# Make category "type" (categories) where species, season, location and location_detailed or combined
metadata_class <- metadata_class %>%
  unite("type", Species,Season)
metadata_class <- metadata_class %>%
  unite("type", type,Location)
metadata_class <- metadata_class %>%
  unite("type", type,Location_detailed)
metadata_class <- as.data.frame.matrix(metadata_class)
metadata_class$type <- as.factor(metadata_class$type)

# class the metadata and count matrix
counts_class_merged2 <- counts_class_merged2[ order((row.names(counts_class_merged2))), ]
counts_class_merged2 <- counts_class_merged2[ ,order((colnames(counts_class_merged2))) ]
metadata_class <- metadata_class[ order((row.names(metadata_class))), ]
class_sample_table <- otu_table(t(counts_class_merged2), taxa_are_rows=F)

# Make metadata with only the samples that are present in our class count matrix
samples_classtable <- colnames(counts_class_merged2)
metadata_class2 <- metadata_class[ rownames(metadata_class) %in% samples_classtable, ]
metadata_class2$Species <- with(metadata, Species[match(rownames(metadata_class2),rownames(metadata))])
metadata_class2$Season <- with(metadata, Season[match(rownames(metadata_class2),rownames(metadata))])
metadata_class2$Location <- with(metadata, Location[match(rownames(metadata_class2),rownames(metadata))])
metadata_class2$Location_detailed <- with(metadata, Location_detailed[match(rownames(metadata_class2),rownames(metadata))])
metadata_class2$Total_reads <- with(metadata_class0, Total_reads[match(rownames(metadata_class2),rownames(metadata_class0))])
metadata_class2 <- metadata_class2[ order((metadata_class2$type)), ]

## Stacked relative prey abundance plots
########################################

# Order samples metadata for plot
Type_classing <- c("Limnothrissa_miodon_Dry_North_Uvira","Limnothrissa_miodon_Wet_North_Uvira","Limnothrissa_miodon_Wet_North_Bujumbura","Limnothrissa_miodon_Wet_Center_Kalemie","Limnothrissa_miodon_Wet_South_Sumbu","Limnothrissa_miodon_Wet_South_Mpulungu","Stolothrissa_tanganicae_Dry_North_Uvira","Stolothrissa_tanganicae_Wet_North_Uvira","Stolothrissa_tanganicae_Wet_Center_Kalemie","Stolothrissa_tanganicae_Wet_South_Sumbu","Stolothrissa_tanganicae_Wet_South_Mpulungu")
metadata_class2$type <- factor(metadata_class2$type, levels=Type_classing)
metadata_class2 <- metadata_class2[order(metadata_class2$type),]

# Order prey count matrix for plot
Prey_classing <- c("Actinopterygii","Malacostraca","Copepoda","Insecta","Hydrozoa")
counts_class_merged3 <- as.data.frame(counts_class_merged2)
counts_class_merged3$Prey <- factor(rownames(counts_class_merged3))
counts_class_merged3$Prey <- factor(counts_class_merged3$Prey, levels=Prey_classing)
counts_class_merged3 <- counts_class_merged3[order(counts_class_merged3$Prey),]
counts_class_merged3 <- subset(counts_class_merged3, select = -Prey )

# Order prey relative count matrix for plot
counts_class_merged_rel2 <- as.data.frame(counts_class_merged_rel)
counts_class_merged_rel2$Prey <- factor(rownames(counts_class_merged_rel2), levels=Prey_classing)
counts_class_merged_rel2 <- counts_class_merged_rel2[order(counts_class_merged_rel2$Prey),]
counts_class_merged_rel2 <- subset(counts_class_merged_rel2, select = -Prey )

# Melt the count matrix for stacked plot
library(reshape2)
counts_class_merged2_melt <- melt(counts_class_merged2, id.vars = "Sample", variable.name = "Phyla") # Melts samples for abundance bars

# Plot stacked barplot with absolute number of reads
library(ggplot2)
class_cols = c("cyan4","chocolate","lightsteelblue1","darkolivegreen3","firebrick") 

ggplot(counts_class_merged2_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=class_cols) + 
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Melt the relative count matrix for stacked plot and class the prey as specified before
counts_class_merged_rel_melt <- melt(counts_class_merged_rel, id.vars = "Sample", variable.name = "Phyla")
counts_class_merged_rel_melt$Var1 <- factor(counts_class_merged_rel_melt$Var1, levels=Prey_classing)
counts_class_merged_rel_melt <- counts_class_merged_rel_melt[order(counts_class_merged_rel_melt$Var1),]

# Plot stacked barplot with relative number of reads
class_cols = c("cyan4","firebrick","chocolate","darkolivegreen3","lightsteelblue1") 
ggplot(counts_class_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=class_cols) + 
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Make grid for each host category
# Add metadata to merged melted
counts_class_merged_rel_melt$Species <-with(metadata_class0, Species[match(counts_class_merged_rel_melt$Var2,rownames(metadata_class0))])
counts_class_merged_rel_melt$Season <-with(metadata_class0, Season[match(counts_class_merged_rel_melt$Var2,rownames(metadata_class0))])
counts_class_merged_rel_melt$Location <-with(metadata_class0, Location[match(counts_class_merged_rel_melt$Var2,rownames(metadata_class0))])
counts_class_merged_rel_melt$Location_detailed <-with(metadata_class0, Location_detailed[match(counts_class_merged_rel_melt$Var2,rownames(metadata_class0))])
counts_class_merged_rel_melt$Total_reads <-with(metadata_class, Total_reads[match(counts_class_merged_rel_melt$Var2,rownames(metadata_class))])
counts_class_merged_rel_melt$type <-with(metadata_class, type[match(counts_class_merged_rel_melt$Var2,rownames(metadata_class))])
library(tidyverse)
counts_class_merged_rel_melt <- counts_class_merged_rel_melt %>%
  unite("Loc_Season", Location_detailed:Season)

# Order by species
counts_class_merged_rel_melt$Species<-factor(counts_class_merged_rel_melt$Species)
counts_class_merged_rel_melt <- counts_class_merged_rel_melt[order(counts_class_merged_rel_melt$Species), ]
levels(counts_class_merged_rel_melt$Species)
counts_class_merged_rel_melt <- counts_class_merged_rel_melt %>% arrange(Species)

# Relative value stacked barplot per detailed location and season
ggplot(counts_class_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~Loc_Season, scales=("free"),nrow=3) + 
  scale_fill_manual(values=class_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  guides(fill=guide_legend(title="Prey class")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Change species names for plot
SpeciesID_class=as.matrix(counts_class_merged_rel_melt$Species)
names(SpeciesID_class)=counts_class_merged_rel_melt$Var2
SpeciesID_class[SpeciesID_class == "Limnothrissa_miodon"] <- "L. miodon"
SpeciesID_class[SpeciesID_class == "Stolothrissa_tanganicae"] <- "S. tanganicae"

# Order by host category
Sample_classing <- factor(rownames(metadata_class2), levels = rownames(metadata_class2))
counts_class_merged_rel_melt$Var2 <- factor(counts_class_merged_rel_melt$Var2, levels=Sample_classing)
counts_class_merged_rel_melt <- counts_class_merged_rel_melt[order(counts_class_merged_rel_melt$Var2),]
counts_class_merged_rel_melt$type <- factor(counts_class_merged_rel_melt$type, levels=Type_classing)

# All types arranged in different small plots
ggplot(counts_class_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_grid(Species~Loc_Season, scales="free_x", space="free") + 
  scale_fill_manual(values=class_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  theme_bw() +
  guides(fill=guide_legend(title="Prey class")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Stacked barplot arranged per host category (type)
ggplot(data=counts_class_merged_rel_melt, aes(x=Var2, y=value, group=type)) +
  xlab("Sample") + ylab("Relative abundance") + #adds title axes
  scale_fill_manual(values=class_cols) +
  geom_bar(aes(fill=Var1),stat = "identity") + #adds catagory TissueID by colouring bars
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  theme_classic() +
  ggtitle("Relative Prey abundance (class level)") +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Nr reads per sample barplot to add
ggplot(data=metadata_class2, aes(x=row.names(metadata_class2), y=log10(Total_reads))) +
  geom_bar(aes(fill="Black"),stat = "identity") + 
  xlab("Sample") + ylab("Log10 number of reads") + #adds title axes
  theme_classic() +
  scale_fill_manual(values="black")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  ggtitle("Log10 number of reads per sample") 


## Bray-Curtis dissimilarity
############################

library(vegan)
library(ecodist)

# Calculate Bray-Curtis
dist.BrCu.all_class<-vegdist(t(counts_class_merged2),method="bray")
clust.BrCu.all_class<-hclust(dist.BrCu.all_class, method="ward.D") #agglomerative clustering using complete linkage

# Plot dendogram
plot(clust.BrCu.all_class)

# Make a PCoA based on Bray-curtis dissimilarities (BCD)
# Not done with count matrix since because the data must be linear.
endo.pco.bray_class<-wcmdscale(dist.BrCu.all_class)
plot(endo.pco.bray_class)

# Plot with ggplot
ggplot(endo.pco.bray_class, aes(x=V1, y=V2)) +
  geom_point(size=5) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller

# Use BCD in Non Multidimensional scaling (NMDS) plot, works better for datasets with many zeroes
BC.nmds_class <- metaMDS(t(counts_class_merged_rel), distance="bray", autotransform = F)
envfit_class <- envfit(BC.nmds_class, metadata_class2, permutations = 999,na.rm = TRUE) # fits environmental vectors
sample.fit_class <- envfit(BC.nmds_class, t(counts_class_merged_rel), permutations = 999) # fits species vectors

# New dataset with point coordinates of NMDS results for ggplot plotting
site.scrs_class <- as.data.frame(vegan::scores(BC.nmds_class, display = "sites"))

# Add metadata
site.scrs_class$Species <-with(metadata_class2, Species[match(rownames(site.scrs_class),rownames(metadata_class2))])
site.scrs_class$Location_detailed <-with(metadata_class2, Location_detailed[match(rownames(site.scrs_class),rownames(metadata_class2))])
site.scrs_class$Season <-with(metadata_class2, Season[match(rownames(site.scrs_class),rownames(metadata_class2))])
site.scrs_class$SL <-with(metadata_class2, SL[match(rownames(site.scrs_class),rownames(metadata_class2))])
site.scrs_class$Total_reads <-with(metadata_class2, Total_reads[match(rownames(site.scrs_class),rownames(metadata_class2))])
site.scrs_class$Loc_season <- as.factor(paste(site.scrs_class$Location_detailed, site.scrs_class$Season, sep = ""))
head(site.scrs_class)

# Make dataset with prey vectors, to add to ggplot
spp.scrs_class <- as.data.frame(vegan::scores(sample.fit_class, display = "vectors"))
spp.scrs_class <- cbind(spp.scrs_class, Species = rownames(spp.scrs_class)) # add prey names
spp.scrs_class <- cbind(spp.scrs_class, pval = sample.fit_class$vectors$pvals) # add pvalues
sig.spp.scrs_class_class <- subset(spp.scrs_class, pval<=0.05) # shows prey with significant values
head(spp.scrs_class)

# Make dataset to show extrinsic variables
env.scores_class<- as.data.frame(vegan::scores(envfit_class, display = "vectors")) #extracts relevant scores from envifit
env.scores_class <- cbind(env.scores_class, env.variables = rownames(env.scores_class)) #and then gives them their names
env.scores_class <- cbind(env.scores_class, pval = envfit_class$vectors$pvals) # add pvalues
sig.env.scrs_class <- subset(env.scores_class, pval<=0.05) # significant
head(env.scores_class)

# ggplot of NMDS
Locations_shape <- c(24,23,25,6,14,2) # specify shapes plot points
nmds.plot_class <- ggplot(site.scrs_class, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = factor(Species),colour = factor(Species), shape = factor(Loc_season)), size = 6)+ 
  coord_fixed()+
  scale_shape_manual(values=Locations_shape) +
  scale_fill_manual(values=c("coral4","lightblue3")) +
  scale_colour_manual(values=c("coral4","lightblue3")) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Species", shape = "Tissue")+ 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))

nmds.plot_class + labs(title = "Basic ordination plot") 

# Add vectors significant prey
nmds.plot_class+
  geom_segment(data = sig.spp.scrs_class_class, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + # add vectors
  labs(title = "Non-Metric Multidimensional Scaling")

# Make same plot with log10 total nr of reads to check for correlation
ggplot(site.scrs_class, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = log10(Total_reads),colour=log10(Total_reads),  size = 6, shape=Species))+
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

# Make same plot with standard length to check for correlation
ggplot(site.scrs_class, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = SL,colour=SL,  size = 6, shape=Species))+ 
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot


# Make heatmap based on BCD values
library(pheatmap)

BCD_matrix_class <- as.matrix(dist.BrCu.all_class)
pheatmap::pheatmap((1-BCD_matrix_class), symm = TRUE,clustering_method="ward.D",show_colnames=TRUE)

## Calculate diversity values
#############################

# Calculate Shannon_class diversity and add to metadata
diversity(counts_class_merged2, index = "shannon") 
Shannon_class <- as.data.frame(diversity(t(counts_class_merged2), index = "shannon"))
Shannon_class$div <- diversity(t(counts_class_merged2), index = "shannon")
metadata$Shannon_class <- with(Shannon_class, div[match(rownames(metadata),rownames(Shannon_class))])
metadata_class2$Shannon_class <- with(Shannon_class, div[match(rownames(metadata_class2),rownames(Shannon_class))])

# Make new metadata
metadata_class3 <- metadata_class2
metadata_class3 <- metadata_class3 %>%
  unite("LocSeas",Location,Season)
metadata_class2$LocSeas <- with(metadata_class3, LocSeas[match(rownames(metadata_class2),rownames(metadata_class3))])

# Plot diversity barplots
ggplot(metadata_class2, aes(x=LocSeas, y=Shannon_class,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

# Calculate Chao1 diversity and add to metadata
chao1div_class <-as.data.frame(t(estimateR(t(counts_class_merged2))))
chao1div_class$Chao1 <- as.numeric(chao1div_class[,2])
metadata_class2$Chao1 <- with(chao1div_class, Chao1[match(rownames(metadata_class2),rownames(chao1div_class))])

# Plot diversity barplots
ggplot(metadata_class2, aes(x=LocSeas, y=Chao1,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

## Look for host factors influencing prey composition
#####################################################

# PERMANOVA based on BCD
BCD_PERMA_class<-adonis2(dist.BrCu.all_class~Species+Location_detailed+Season, data=metadata_class2, permutations = 999, method="bray")
BCD_PERMA_class
adonis2(dist.BrCu.all_class~Species, data=metadata_class2, permutations = 999, p.adjust.m = "bonferroni")
adonis2(dist.BrCu.all_class~Species+Location_detailed, data=metadata_class2, permutations = 999, p.adjust.m = "bonferroni")
adonis2(dist.BrCu.all_class~Location_detailed, data=metadata_class2, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA SL and Total reads
adonis2(dist.BrCu.all_class~Total_reads, data=metadata_class2, permutations = 999, p.adjust.m = "bonferroni")

no_SL_samples_class <- rownames(metadata_class2[which(is.na(metadata_class2$SL)),])
metadata_SL_class <- metadata_class2[ !rownames(metadata_class2) %in% no_SL_samples_class, ]
counts_SL_class <- counts_class_merged2[ ,!colnames(counts_class_merged2) %in% no_SL_samples_class ]
BrCu.SL_class<-vegdist(t(counts_SL_class),method="bray")
adonis2(BrCu.SL_class~SL, data=metadata_SL_class, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA Samples Uvira
Uvira_samples_class <- rownames(metadata_class2[which(metadata_class2$Location_detailed == "Uvira"),])
metadata_uvira_class <- metadata_class2[ rownames(metadata_class2) %in% Uvira_samples_class, ]
counts_uvira_class <- counts_class_merged2[ ,colnames(counts_class_merged2) %in% Uvira_samples_class ]

BrCu.uvira_class<-vegdist(t(counts_uvira_class),method="bray")
adonis2(BrCu.uvira_class~Species+Season, data=metadata_uvira_class, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.uvira_class~Species, data=metadata_uvira_class, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.uvira_class~Season, data=metadata_uvira_class, permutations = 999, p.adjust.m = "bonferroni")

# Uvira L.miodon
Uvira_LM_samples_class <- rownames(metadata_uvira_class[which(metadata_uvira_class$Species == "Limnothrissa_miodon"),])
metadata_uvira_LM_class <- metadata_uvira_class[ rownames(metadata_uvira_class) %in% Uvira_LM_samples_class, ]
counts_uvira_LM_class <- counts_uvira_class[ ,colnames(counts_uvira_class) %in% Uvira_LM_samples_class ]

BrCu.uvira_LM_class<-vegdist(t(counts_uvira_LM_class),method="bray")
adonis2(BrCu.uvira_LM_class~Season, data=metadata_uvira_LM_class, permutations = 999, p.adjust.m = "bonferroni")

# Uvira S.tanganicae
Uvira_ST_samples_class <- rownames(metadata_uvira_class[which(metadata_uvira_class$Species == "Stolothrissa_tanganicae"),])
metadata_uvira_ST_class <- metadata_uvira_class[ rownames(metadata_uvira_class) %in% Uvira_ST_samples_class, ]
counts_uvira_ST_class <- counts_uvira_class[ ,colnames(counts_uvira_class) %in% Uvira_ST_samples_class ]

BrCu.uvira_ST_class<-vegdist(t(counts_uvira_ST_class),method="bray")
adonis2(BrCu.uvira_ST_class~Season, data=metadata_uvira_ST_class, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA without Uvira Dry
Uvira_Dry_samples_class <- rownames(metadata_uvira_class[which(metadata_uvira_class$Season == "Dry"),])
counts_wet_class <- counts_class_merged2[ ,! colnames(counts_class_merged2) %in% Uvira_Dry_samples_class ]
metadata_wet_class <- metadata_class2[ ! rownames(metadata_class2) %in% Uvira_Dry_samples_class, ]

BrCu.wet_class<-vegdist(t(counts_wet_class),method="bray")
adonis2(BrCu.wet_class~Location_detailed, data=metadata_wet_class, permutations = 999, method="bonferroni")
adonis2(BrCu.wet_class~Species+Location_detailed, data=metadata_wet_class, permutations = 999, method="bonferroni")

adonis2(BrCu.wet_class~Species+Location_detailed, data=metadata_wet_class, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.wet_class~Species, data=metadata_wet_class, permutations = 999, p.adjust.m = "bonferroni")

# Wet L.miodon
Wet_LM_samples_class <- rownames(metadata_wet_class[which(metadata_wet_class$Species == "Limnothrissa_miodon"),])
metadata_wet_LM_class <- metadata_wet_class[ rownames(metadata_wet_class) %in% Wet_LM_samples_class, ]
counts_wet_LM_class <- counts_wet_class[ ,colnames(counts_wet_class) %in% Wet_LM_samples_class ]

BrCu.wet_LM_class<-vegdist(t(counts_wet_LM_class),method="bray")
adonis2(BrCu.wet_LM_class~Location_detailed, data=metadata_wet_LM_class,permutations = 999, p.adjust.m = "bonferroni")

# Wet S.tanganicae
Wet_ST_samples_class <- rownames(metadata_wet_class[which(metadata_wet_class$Species == "Stolothrissa_tanganicae"),])
metadata_wet_ST_class <- metadata_wet_class[ rownames(metadata_wet_class) %in% Wet_ST_samples_class, ]
counts_wet_ST_class <- counts_wet_class[ ,colnames(counts_wet_class) %in% Wet_ST_samples_class ]

BrCu.wet_ST_class<-vegdist(t(counts_wet_ST_class),method="bray")
adonis2(BrCu.wet_ST_class~Location_detailed, data=metadata_wet_ST_class, permutations = 999, p.adjust.m = "bonferroni")

## Additional general calculations
##################################

# Percentage abundance prey per sardine species
LM_samples_class <- rownames(metadata_class2[which(metadata_class2$Species == "Limnothrissa_miodon"),])
counts_LM_class <- counts_class_merged_rel[ ,colnames(counts_class_merged_rel) %in% LM_samples_class ]
LM_taxa_class <- as.matrix(rowSums(counts_LM_class))
LM_taxa_rel_class <- LM_taxa_class/colSums(LM_taxa_class)

ST_samples_class <- rownames(metadata_class2[which(metadata_class2$Species == "Stolothrissa_tanganicae"),])
counts_ST_class <- counts_class_merged_rel[ ,colnames(counts_class_merged_rel) %in% ST_samples_class ]
ST_taxa_class <- as.matrix(rowSums(counts_ST_class))
ST_taxa_rel_class <- ST_taxa_class/colSums(ST_taxa_class)

# Nr Reads

mean(metadata_class2$Total_reads)
median(metadata_class2$Total_reads)
min(metadata_class2$Total_reads)
sd(metadata_class2$Total_reads)
max(metadata_class2$Total_reads)

# Mean diversity
metadata_ST_class <- metadata_class2[ rownames(metadata_class2) %in% ST_samples_class , ]
mean(metadata_ST_class$Shannon)
mean(metadata_ST_class$Chao1)

metadata_LM_class <- metadata_class2[ rownames(metadata_class2) %in% LM_samples_class , ]
mean(metadata_LM_class$Shannon)
mean(metadata_LM_class$Chao1)



################################################################################
# Family
################################################################################

# Make subsets for order level ID
counts_family_1 <- as.matrix(counts)
counts_family_1 <- counts_family_1[ order((row.names(counts_family_1))), ]
otus_family <- OTUs[ order((row.names(OTUs))), ]

# Remove OTU's which were not identified till family level
nofamily_otu <- rownames(otus_family[which(otus_family$Similarity < 90 ),])
otus_family <- otus_family[ !rownames(otus_family) %in% nofamily_otu, ]
counts_family_1 <- counts_family_1[ !rownames(counts_family_1) %in% nofamily_otu, ]
rmSamplesFamily <- colnames(counts_family_1[,colSums(counts_family_1) == 0])
counts_family_1 <- counts_family_1[ ,!colnames(counts_family_1) %in% rmSamplesFamily ]

# Merge samples
otus_family$Family <- as.factor(otus_family$Family)
counts_family_table <- otu_table(as.matrix(counts_family_1), taxa_are_rows=F)
counts_family_merged <- merge_samples(counts_family_table, t(otus_family$Family))
counts_family_merged2 <- as(otu_table(counts_family_merged), "matrix")

# Remove samples low counts
counts_family_merged2 = counts_family_merged2[,colSums(counts_family_merged2) > 10]
counts_family_merged2 = counts_family_merged2[rowSums(counts_family_merged2) > 10,]

# Remove these samples from metadata
samples_family_keep <- colnames(counts_family_merged2)
metadata_family0 <- metadata
metadata_family0 <- metadata_family0[ rownames(metadata_family0) %in% samples_family_keep, ]

#Take relative values
tot_reads_family_samples <- colSums(counts_family_merged2)
family_samples <- colnames(counts_family_merged2)
counts_family_merged_rel <- t(t(counts_family_merged2) / tot_reads_family_samples)

# Total counts per sample
reads_family <- as.data.frame(rowSums(t(counts_family_merged2)))
colnames(reads_family) <- 'reads'
metadata_family0$Total_reads <- with(reads_family, reads[match(rownames(metadata_family0),rownames(reads_family))])

#Take relative values per sample
tot_reads_family_samples <- colSums(counts_family_merged2)
family_samples <- colnames(counts_family_merged2)
counts_family_merged_rel <- t(t(counts_family_merged2) / tot_reads_family_samples)

## Optimise metadata for downstream analysis
############################################

library(dplyr)
library(tidyr)
metadata_family <- metadata_family0

# Make category "type" (categories) where species, season, location and location_detailed or combined
metadata_family <- metadata_family %>%
  unite("type", Species,Season)
metadata_family <- metadata_family %>%
  unite("type", type,Location)
metadata_family <- metadata_family %>%
  unite("type", type,Location_detailed)
metadata_family <- as.data.frame.matrix(metadata_family)
metadata_family$type <- as.factor(metadata_family$type)

# family the metadata and count matrix
counts_family_merged2 <- counts_family_merged2[ order((row.names(counts_family_merged2))), ]
counts_family_merged2 <- counts_family_merged2[ ,order((colnames(counts_family_merged2))) ]
metadata_family <- metadata_family[ order((row.names(metadata_family))), ]
family_sample_table <- otu_table(t(counts_family_merged2), taxa_are_rows=F)

# Make metadata with only the samples that are present in our family count matrix
samples_familytable <- colnames(counts_family_merged2)
metadata_family2 <- metadata_family[ rownames(metadata_family) %in% samples_familytable, ]
metadata_family2$Species <- with(metadata, Species[match(rownames(metadata_family2),rownames(metadata))])
metadata_family2$Season <- with(metadata, Season[match(rownames(metadata_family2),rownames(metadata))])
metadata_family2$Location <- with(metadata, Location[match(rownames(metadata_family2),rownames(metadata))])
metadata_family2$Location_detailed <- with(metadata, Location_detailed[match(rownames(metadata_family2),rownames(metadata))])
metadata_family2$Total_reads <- with(metadata_family0, Total_reads[match(rownames(metadata_family2),rownames(metadata_family0))])
metadata_family2 <- metadata_family2[ order((metadata_family2$type)), ]

## Stacked relative prey abundance plots
########################################

# Order samples metadata for plot
Type_familying <- c("Limnothrissa_miodon_Wet_North_Uvira","Limnothrissa_miodon_Wet_North_Bujumbura","Limnothrissa_miodon_Wet_Center_Kalemie","Limnothrissa_miodon_Wet_South_Sumbu","Limnothrissa_miodon_Wet_South_Mpulungu","Stolothrissa_tanganicae_Dry_North_Uvira","Stolothrissa_tanganicae_Wet_North_Uvira","Stolothrissa_tanganicae_Wet_Center_Kalemie","Stolothrissa_tanganicae_Wet_South_Sumbu","Stolothrissa_tanganicae_Wet_South_Mpulungu")
metadata_family2$type <- factor(metadata_family2$type, levels=Type_familying)
metadata_family2 <- metadata_family2[order(metadata_family2$type),]

# Order prey count matrix for plot
Prey_familying <- c("Cichlidae","Latidae","Procatopodidae","Claroteidae","Diaptomidae","Chironomidae","Chloropidae","Dolichopodidae","Sphaeroceridae","Syrphidae","Olindiidae")
counts_family_merged3 <- as.data.frame(counts_family_merged2)
counts_family_merged3$Prey <- factor(rownames(counts_family_merged3))
counts_family_merged3$Prey <- factor(counts_family_merged3$Prey, levels=Prey_familying)
counts_family_merged3 <- counts_family_merged3[order(counts_family_merged3$Prey),]
counts_family_merged3 <- subset(counts_family_merged3, select = -Prey )

# Order prey relative count matrix for plot
counts_family_merged_rel2 <- as.data.frame(counts_family_merged_rel)
counts_family_merged_rel2$Prey <- factor(rownames(counts_family_merged_rel2), levels=Prey_familying)
counts_family_merged_rel2 <- counts_family_merged_rel2[order(counts_family_merged_rel2$Prey),]
counts_family_merged_rel2 <- subset(counts_family_merged_rel2, select = -Prey )

# Melt the count matrix for stacked plot
library(reshape2)
counts_family_merged2_melt <- melt(counts_family_merged2, id.vars = "Sample", variable.name = "Phyla") # Melts samples for abundance bars

# Plot stacked barplot with absolute number of reads
library(ggplot2)
family_cols = c("cyan4","lightskyblue4","dodgerblue4","paleturquoise3","chocolate","darkolivegreen3","darkolivegreen4","darkolivegreen1","darkolivegreen2","darkolivegreen","lightsteelblue1") 

ggplot(counts_family_merged2_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Melt the relative count matrix for stacked plot and family the prey as specified before
counts_family_merged_rel_melt <- melt(counts_family_merged_rel, id.vars = "Sample", variable.name = "Phyla")
counts_family_merged_rel_melt$Var1 <- factor(counts_family_merged_rel_melt$Var1, levels=Prey_familying)
counts_family_merged_rel_melt <- counts_family_merged_rel_melt[order(counts_family_merged_rel_melt$Var1),]

# Plot stacked barplot with relative number of reads
family_cols = c("cyan4","lightskyblue4","dodgerblue4","paleturquoise3","chocolate","darkolivegreen3","darkolivegreen4","darkolivegreen1","darkolivegreen2","darkolivegreen","lightsteelblue1") 

ggplot(counts_family_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=family_cols) + 
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Make grid for each host category
# Add metadata to merged melted
counts_family_merged_rel_melt$Species <-with(metadata_family0, Species[match(counts_family_merged_rel_melt$Var2,rownames(metadata_family0))])
counts_family_merged_rel_melt$Season <-with(metadata_family0, Season[match(counts_family_merged_rel_melt$Var2,rownames(metadata_family0))])
counts_family_merged_rel_melt$Location <-with(metadata_family0, Location[match(counts_family_merged_rel_melt$Var2,rownames(metadata_family0))])
counts_family_merged_rel_melt$Location_detailed <-with(metadata_family0, Location_detailed[match(counts_family_merged_rel_melt$Var2,rownames(metadata_family0))])
counts_family_merged_rel_melt$Total_reads <-with(metadata_family, Total_reads[match(counts_family_merged_rel_melt$Var2,rownames(metadata_family))])
counts_family_merged_rel_melt$type <-with(metadata_family, type[match(counts_family_merged_rel_melt$Var2,rownames(metadata_family))])
library(tidyverse)
counts_family_merged_rel_melt <- counts_family_merged_rel_melt %>%
  unite("Loc_Season", Location_detailed:Season)

# Order by species
counts_family_merged_rel_melt$Species<-factor(counts_family_merged_rel_melt$Species)
counts_family_merged_rel_melt <- counts_family_merged_rel_melt[order(counts_family_merged_rel_melt$Species), ]
levels(counts_family_merged_rel_melt$Species)
counts_family_merged_rel_melt <- counts_family_merged_rel_melt %>% arrange(Species)

# Relative value stacked barplot per detailed location and season
ggplot(counts_family_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~Loc_Season, scales=("free"),nrow=3) + 
  scale_fill_manual(values=family_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  guides(fill=guide_legend(title="Prey family")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Change species names for plot
SpeciesID_family=as.matrix(counts_family_merged_rel_melt$Species)
names(SpeciesID_family)=counts_family_merged_rel_melt$Var2
SpeciesID_family[SpeciesID_family == "Limnothrissa_miodon"] <- "L. miodon"
SpeciesID_family[SpeciesID_family == "Stolothrissa_tanganicae"] <- "S. tanganicae"

# Order by host category
Sample_familying <- factor(rownames(metadata_family2), levels = rownames(metadata_family2))
counts_family_merged_rel_melt$Var2 <- factor(counts_family_merged_rel_melt$Var2, levels=Sample_familying)
counts_family_merged_rel_melt <- counts_family_merged_rel_melt[order(counts_family_merged_rel_melt$Var2),]
counts_family_merged_rel_melt$type <- factor(counts_family_merged_rel_melt$type, levels=Type_familying)

# All types arranged in different small plots
ggplot(counts_family_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_grid(Species~Loc_Season, scales="free_x", space="free") + 
  scale_fill_manual(values=family_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  theme_bw() +
  guides(fill=guide_legend(title="Prey family")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Stacked barplot arranged per host category (type)
ggplot(data=counts_family_merged_rel_melt, aes(x=Var2, y=value, group=type)) +
  xlab("Sample") + ylab("Relative abundance") + #adds title axes
  scale_fill_manual(values=family_cols) +
  geom_bar(aes(fill=Var1),stat = "identity") + #adds catagory TissueID by colouring bars
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  theme_classic() +
  ggtitle("Relative Prey abundance (family level)") +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Nr reads per sample barplot to add
ggplot(data=metadata_family2, aes(x=row.names(metadata_family2), y=log10(Total_reads))) +
  geom_bar(aes(fill="Black"),stat = "identity") + 
  xlab("Sample") + ylab("Log10 number of reads") + #adds title axes
  theme_classic() +
  scale_fill_manual(values="black")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  ggtitle("Log10 number of reads per sample") 

## Bray-Curtis dissimilarity
############################

library(vegan)
library(ecodist)

# Calculate Bray-Curtis
dist.BrCu.all_family<-vegdist(t(counts_family_merged2),method="bray")
clust.BrCu.all_family<-hclust(dist.BrCu.all_family, method="ward.D") #agglomerative clustering using complete linkage

# Plot dendogram
plot(clust.BrCu.all_family)

# Make a PCoA based on Bray-curtis dissimilarities (BCD)
# Not done with count matrix since because the data must be linear.
endo.pco.bray_family<-wcmdscale(dist.BrCu.all_family)
plot(endo.pco.bray_family)

# Plot with ggplot
endo.pco.bray_family_df <- as.data.frame(endo.pco.bray_family)
ggplot(endo.pco.bray_family_df, aes(x=V1, y=V2)) +
  geom_point(size=5) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller

# Use BCD in Non Multidimensional scaling (NMDS) plot, works better for datasets with many zeroes
BC.nmds_family <- metaMDS(t(counts_family_merged_rel), distance="bray", autotransform = F)
envfit_family <- envfit(BC.nmds_family, metadata_family2, permutations = 999,na.rm = TRUE) # fits environmental vectors
sample.fit_family <- envfit(BC.nmds_family, t(counts_family_merged_rel), permutations = 999) # fits species vectors

# New dataset with point coordinates of NMDS results for ggplot plotting
site.scrs_family <- as.data.frame(vegan::scores(BC.nmds_family, display = "sites"))

# Add metadata
site.scrs_family$Species <-with(metadata_family2, Species[match(rownames(site.scrs_family),rownames(metadata_family2))])
site.scrs_family$Location_detailed <-with(metadata_family2, Location_detailed[match(rownames(site.scrs_family),rownames(metadata_family2))])
site.scrs_family$Season <-with(metadata_family2, Season[match(rownames(site.scrs_family),rownames(metadata_family2))])
site.scrs_family$SL <-with(metadata_family2, SL[match(rownames(site.scrs_family),rownames(metadata_family2))])
site.scrs_family$Total_reads <-with(metadata_family2, Total_reads[match(rownames(site.scrs_family),rownames(metadata_family2))])
site.scrs_family$Loc_season <- as.factor(paste(site.scrs_family$Location_detailed, site.scrs_family$Season, sep = ""))
head(site.scrs_family)

# Make dataset with prey vectors, to add to ggplot
spp.scrs_family <- as.data.frame(vegan::scores(sample.fit_family, display = "vectors"))
spp.scrs_family <- cbind(spp.scrs_family, Species = rownames(spp.scrs_family)) # add prey names
spp.scrs_family <- cbind(spp.scrs_family, pval = sample.fit_family$vectors$pvals) # add pvalues
sig.spp.scrs_family_family <- subset(spp.scrs_family, pval<=0.05) # shows prey with significant values
head(spp.scrs_family)

# Make dataset to show extrinsic variables
env.scores_family<- as.data.frame(vegan::scores(envfit_family, display = "vectors")) #extracts relevant scores from envifit
env.scores_family <- cbind(env.scores_family, env.variables = rownames(env.scores_family)) #and then gives them their names
env.scores_family <- cbind(env.scores_family, pval = envfit_family$vectors$pvals) # add pvalues
sig.env.scrs_family <- subset(env.scores_family, pval<=0.05) # significant
head(env.scores_family)

# ggplot of NMDS
Locations_shape <- c(24,23,25,6,14,2) # specify shapes plot points
nmds.plot_family <- ggplot(site.scrs_family, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = factor(Species),colour = factor(Species), shape = factor(Loc_season)), size = 6)+ 
  coord_fixed()+
  scale_shape_manual(values=Locations_shape) +
  scale_fill_manual(values=c("coral4","lightblue3")) +
  scale_colour_manual(values=c("coral4","lightblue3")) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Species", shape = "Tissue")+ 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))

nmds.plot_family + labs(title = "Basic ordination plot") 

# Add vectors significant prey
nmds.plot_family+
  geom_segment(data = sig.spp.scrs_family_family, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + # add vectors
  labs(title = "Non-Metric Multidimensional Scaling")

# Make same plot with log10 total nr of reads to check for correlation
ggplot(site.scrs_family, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = log10(Total_reads),colour=log10(Total_reads),  size = 6, shape=Species))+
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

# Make same plot with standard length to check for correlation
ggplot(site.scrs_family, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = SL,colour=SL,  size = 6, shape=Species))+ 
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot


# Make heatmap based on BCD values
library(pheatmap)

BCD_matrix_family <- as.matrix(dist.BrCu.all_family)
pheatmap::pheatmap((1-BCD_matrix_family), symm = TRUE,clustering_method="ward.D",show_colnames=TRUE)

## Calculate diversity values
#############################

# Calculate Shannon_family diversity and add to metadata
diversity(counts_family_merged2, index = "shannon") 
Shannon_family <- as.data.frame(diversity(t(counts_family_merged2), index = "shannon"))
Shannon_family$div <- diversity(t(counts_family_merged2), index = "shannon")
metadata$Shannon <- with(Shannon_family, div[match(rownames(metadata),rownames(Shannon_family))])
metadata_family2$Shannon <- with(Shannon_family, div[match(rownames(metadata_family2),rownames(Shannon_family))])

# Make new metadata
metadata_family3 <- metadata_family2
metadata_family3 <- metadata_family3 %>%
  unite("LocSeas",Location,Season)
metadata_family2$LocSeas <- with(metadata_family3, LocSeas[match(rownames(metadata_family2),rownames(metadata_family3))])

# Plot diversity barplots
ggplot(metadata_family2, aes(x=LocSeas, y=Shannon,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

# Calculate Chao1 diversity and add to metadata
chao1div_family <-as.data.frame(t(estimateR(t(counts_family_merged2))))
chao1div_family$Chao1 <- as.numeric(chao1div_family[,2])
metadata_family2$Chao1 <- with(chao1div_family, Chao1[match(rownames(metadata_family2),rownames(chao1div_family))])

# Plot diversity barplots
ggplot(metadata_family2, aes(x=LocSeas, y=Chao1,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

## Look for host factors influencing prey composition
#####################################################

# PERMANOVA based on BCD
BCD_PERMA_family<-adonis2(dist.BrCu.all_family~Species+Location_detailed+Season, data=metadata_family2, permutations = 999, method="bray")
BCD_PERMA_family
adonis2(dist.BrCu.all_family~Species, data=metadata_family2, permutations = 999, p.adjust.m = "bonferroni")
adonis2(dist.BrCu.all_family~Species+Location_detailed, data=metadata_family2, permutations = 999, p.adjust.m = "bonferroni")
adonis2(dist.BrCu.all_family~Location_detailed, data=metadata_family2, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA SL and Total reads
adonis2(dist.BrCu.all_family~Total_reads, data=metadata_family2, permutations = 999, p.adjust.m = "bonferroni")

no_SL_samples_family <- rownames(metadata_family2[which(is.na(metadata_family2$SL)),])
metadata_SL_family <- metadata_family2[ !rownames(metadata_family2) %in% no_SL_samples_family, ]
counts_SL_family <- counts_family_merged2[ ,!colnames(counts_family_merged2) %in% no_SL_samples_family ]
BrCu.SL_family<-vegdist(t(counts_SL_family),method="bray")
adonis2(BrCu.SL_family~SL, data=metadata_SL_family, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA Samples Uvira
Uvira_samples_family <- rownames(metadata_family2[which(metadata_family2$Location_detailed == "Uvira"),])
metadata_uvira_family <- metadata_family2[ rownames(metadata_family2) %in% Uvira_samples_family, ]
counts_uvira_family <- counts_family_merged2[ ,colnames(counts_family_merged2) %in% Uvira_samples_family ]

BrCu.uvira_family<-vegdist(t(counts_uvira_family),method="bray")
adonis2(BrCu.uvira_family~Species+Season, data=metadata_uvira_family, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.uvira_family~Species, data=metadata_uvira_family, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.uvira_family~Season, data=metadata_uvira_family, permutations = 999, p.adjust.m = "bonferroni")

# Uvira S.tanganicae
Uvira_ST_samples_family <- rownames(metadata_uvira_family[which(metadata_uvira_family$Species == "Stolothrissa_tanganicae"),])
metadata_uvira_ST_family <- metadata_uvira_family[ rownames(metadata_uvira_family) %in% Uvira_ST_samples_family, ]
counts_uvira_ST_family <- counts_uvira_family[ ,colnames(counts_uvira_family) %in% Uvira_ST_samples_family ]

BrCu.uvira_ST_family<-vegdist(t(counts_uvira_ST_family),method="bray")
adonis2(BrCu.uvira_ST_family~Season, data=metadata_uvira_ST_family, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA without Uvira Dry
Uvira_Dry_samples_family <- rownames(metadata_uvira_family[which(metadata_uvira_family$Season == "Dry"),])
counts_wet_family <- counts_family_merged2[ ,! colnames(counts_family_merged2) %in% Uvira_Dry_samples_family ]
metadata_wet_family <- metadata_family2[ ! rownames(metadata_family2) %in% Uvira_Dry_samples_family, ]

BrCu.wet_family<-vegdist(t(counts_wet_family),method="bray")
adonis2(BrCu.wet_family~Location_detailed, data=metadata_wet_family, permutations = 999, method="bonferroni")
adonis2(BrCu.wet_family~Species+Location_detailed, data=metadata_wet_family, permutations = 999, method="bonferroni")

adonis2(BrCu.wet_family~Species+Location_detailed, data=metadata_wet_family, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.wet_family~Species, data=metadata_wet_family, permutations = 999, p.adjust.m = "bonferroni")

# Wet L.miodon
Wet_LM_samples_family <- rownames(metadata_wet_family[which(metadata_wet_family$Species == "Limnothrissa_miodon"),])
metadata_wet_LM_family <- metadata_wet_family[ rownames(metadata_wet_family) %in% Wet_LM_samples_family, ]
counts_wet_LM_family <- counts_wet_family[ ,colnames(counts_wet_family) %in% Wet_LM_samples_family ]

BrCu.wet_LM_family<-vegdist(t(counts_wet_LM_family),method="bray")
adonis2(BrCu.wet_LM_family~Location_detailed, data=metadata_wet_LM_family,permutations = 999, p.adjust.m = "bonferroni")

# Wet S.tanganicae
Wet_ST_samples_family <- rownames(metadata_wet_family[which(metadata_wet_family$Species == "Stolothrissa_tanganicae"),])
metadata_wet_ST_family <- metadata_wet_family[ rownames(metadata_wet_family) %in% Wet_ST_samples_family, ]
counts_wet_ST_family <- counts_wet_family[ ,colnames(counts_wet_family) %in% Wet_ST_samples_family ]

BrCu.wet_ST_family<-vegdist(t(counts_wet_ST_family),method="bray")
adonis2(BrCu.wet_ST_family~Location_detailed, data=metadata_wet_ST_family, permutations = 999, p.adjust.m = "bonferroni")

# Percentage abundance prey per sardine species
LM_samples_family <- rownames(metadata_family2[which(metadata_family2$Species == "Limnothrissa_miodon"),])
counts_LM_family <- counts_family_merged_rel[ ,colnames(counts_family_merged_rel) %in% LM_samples_family ]
LM_taxa_family <- as.matrix(rowSums(counts_LM_family))
LM_taxa_rel_family <- LM_taxa_family/colSums(LM_taxa_family)

ST_samples_family <- rownames(metadata_family2[which(metadata_family2$Species == "Stolothrissa_tanganicae"),])
counts_ST_family <- counts_family_merged_rel[ ,colnames(counts_family_merged_rel) %in% ST_samples_family ]
ST_taxa_family <- as.matrix(rowSums(counts_ST_family))
ST_taxa_rel_family <- ST_taxa_family/colSums(ST_taxa_family)

## Additional general calculations
##################################

# Percentage abundance prey per sardine species
LM_samples_family <- rownames(metadata_family2[which(metadata_family2$Species == "Limnothrissa_miodon"),])
counts_LM_family <- counts_family_merged_rel[ ,colnames(counts_family_merged_rel) %in% LM_samples_family ]
LM_taxa_family <- as.matrix(rowSums(counts_LM_family))
LM_taxa_rel_family <- LM_taxa_family/colSums(LM_taxa_family)
LM_taxa_rel_family

ST_samples_family <- rownames(metadata_family2[which(metadata_family2$Species == "Stolothrissa_tanganicae"),])
counts_ST_family <- counts_family_merged_rel[ ,colnames(counts_family_merged_rel) %in% ST_samples_family ]
ST_taxa_family <- as.matrix(rowSums(counts_ST_family))
ST_taxa_rel_family <- ST_taxa_family/colSums(ST_taxa_family)
ST_taxa_rel_family

# Nr Reads
mean(metadata_family2$Total_reads)
median(metadata_family2$Total_reads)
min(metadata_family2$Total_reads)
sd(metadata_family2$Total_reads)
max(metadata_family2$Total_reads)

# Mean diversity
metadata_ST_family <- metadata_family2[ rownames(metadata_family2) %in% ST_samples_family , ]
mean(metadata_ST_family$Shannon)
mean(metadata_ST_family$Chao1)

metadata_LM_family <- metadata_family2[ rownames(metadata_family2) %in% LM_samples_family , ]
mean(metadata_LM_family$Shannon)
mean(metadata_LM_family$Chao1)


################################################################################
# Genus
################################################################################

# Make subsets for Species level ID
counts_genus <- as.matrix(counts)
counts_genus <- counts_genus[ order((row.names(counts_genus))), ]
otus_genus <- OTUs[ order((row.names(OTUs))), ]

# Remove OTU's which were not identified till species level
nogenus_otu <- rownames(otus_genus[which(otus_genus$Similarity < 95 ),])
otus_genus <- otus_genus[ !rownames(otus_genus) %in% nogenus_otu, ]
counts_genus <- counts_genus[ !rownames(counts_genus) %in% nogenus_otu, ]
rmSamplesGenus <- colnames(counts_genus[,colSums(counts_genus) == 0])
counts_genus <- counts_genus[ ,!colnames(counts_genus) %in% rmSamplesGenus ]

# Merge samples
otus_genus$Genus <- as.factor(otus_genus$Genus)
counts_genus_table <- otu_table(as.matrix(counts_genus), taxa_are_rows=F)
counts_genus_merged <- merge_samples(counts_genus_table, t(otus_genus$Genus))
counts_genus_merged2 <- as(otu_table(counts_genus_merged), "matrix")

# Remove samples low counts
counts_genus_merged2 = counts_genus_merged2[,colSums(counts_genus_merged2) > 10]
counts_genus_merged2 = counts_genus_merged2[rowSums(counts_genus_merged2) > 10,]

# Remove these samples from metadata
samples_genus_keep <- colnames(counts_genus_merged2)
metadata_genus0 <- metadata
metadata_genus0 <- metadata_genus0[ rownames(metadata_genus0) %in% samples_genus_keep, ]

#Take relative values
tot_reads_genus_samples <- colSums(counts_genus_merged2)
genus_samples <- colnames(counts_genus_merged2)
counts_genus_merged_rel <- t(t(counts_genus_merged2) / tot_reads_genus_samples)

# Total counts per sample
reads_genus <- as.data.frame(rowSums(t(counts_genus_merged2)))
colnames(reads_genus) <- 'reads'
metadata_genus0$Total_reads <- with(reads_genus, reads[match(rownames(metadata_genus0),rownames(reads_genus))])

#Take relative values per sample
tot_reads_genus_samples <- colSums(counts_genus_merged2)
genus_samples <- colnames(counts_genus_merged2)
counts_genus_merged_rel <- t(t(counts_genus_merged2) / tot_reads_genus_samples)

## Optimise metadata for downstream analysis
############################################

library(dplyr)
library(tidyr)
metadata_genus <- metadata_genus0

# Make category "type" (categories) where species, season, location and location_detailed or combined
metadata_genus <- metadata_genus %>%
  unite("type", Species,Season)
metadata_genus <- metadata_genus %>%
  unite("type", type,Location)
metadata_genus <- metadata_genus %>%
  unite("type", type,Location_detailed)
metadata_genus <- as.data.frame.matrix(metadata_genus)
metadata_genus$type <- as.factor(metadata_genus$type)

# genus the metadata and count matrix
counts_genus_merged2 <- counts_genus_merged2[ order((row.names(counts_genus_merged2))), ]
counts_genus_merged2 <- counts_genus_merged2[ ,order((colnames(counts_genus_merged2))) ]
metadata_genus <- metadata_genus[ order((row.names(metadata_genus))), ]
genus_sample_table <- otu_table(t(counts_genus_merged2), taxa_are_rows=F)

# Make metadata with only the samples that are present in our genus count matrix
samples_genustable <- colnames(counts_genus_merged2)
metadata_genus2 <- metadata_genus[ rownames(metadata_genus) %in% samples_genustable, ]
metadata_genus2$Species <- with(metadata, Species[match(rownames(metadata_genus2),rownames(metadata))])
metadata_genus2$Season <- with(metadata, Season[match(rownames(metadata_genus2),rownames(metadata))])
metadata_genus2$Location <- with(metadata, Location[match(rownames(metadata_genus2),rownames(metadata))])
metadata_genus2$Location_detailed <- with(metadata, Location_detailed[match(rownames(metadata_genus2),rownames(metadata))])
metadata_genus2$Total_reads <- with(metadata_genus0, Total_reads[match(rownames(metadata_genus2),rownames(metadata_genus0))])
metadata_genus2 <- metadata_genus2[ order((metadata_genus2$type)), ]

## Stacked relative prey abundance plots
########################################

# Order samples metadata for plot
Type_genusing <- c("Limnothrissa_miodon_Wet_North_Uvira","Limnothrissa_miodon_Wet_North_Bujumbura","Limnothrissa_miodon_Wet_Center_Kalemie","Limnothrissa_miodon_Wet_South_Sumbu","Limnothrissa_miodon_Wet_South_Mpulungu","Stolothrissa_tanganicae_Dry_North_Uvira","Stolothrissa_tanganicae_Wet_North_Uvira","Stolothrissa_tanganicae_Wet_Center_Kalemie","Stolothrissa_tanganicae_Wet_South_Sumbu","Stolothrissa_tanganicae_Wet_South_Mpulungu")
metadata_genus2$type <- factor(metadata_genus2$type, levels=Type_genusing)
metadata_genus2 <- metadata_genus2[order(metadata_genus2$type),]

# Order prey count matrix for plot
Prey_genusing <- c("Altolamprologus","Bathybates","Benthochromis","Callochromis","Chalinochromis","Cyprichromis","Eretmodus","Haplochromis","Julidochromis","Lamprologus","Lepidiolamprologus","Limnotilapia","Neolamprologus","Ophthalmotilapia","Petrochromis","Telmatochromis","Xenotilapia","Lates","Lamprichthys","Auchenoglanis","Chrysichthys","Tropodiaptomus","Limnocnida")

counts_genus_merged_rel2 <- as.data.frame(counts_genus_merged_rel)
counts_genus_merged_rel2$Prey <- factor(rownames(counts_genus_merged_rel2))
counts_genus_merged_rel2$Prey <- factor(rownames(counts_genus_merged_rel2), levels=Prey_genusing)
counts_genus_merged_rel2 <- counts_genus_merged_rel2[order(counts_genus_merged_rel2$Prey),]
counts_genus_merged_rel2 <- subset(counts_genus_merged_rel2, select = -Prey )

# Melt the count matrix for stacked plot
library(reshape2)
counts_genus_merged2_melt <- melt(counts_genus_merged2, id.vars = "Sample", variable.name = "Phyla") # Melts samples for abundance bars

# Plot stacked barplot with absolute number of reads
library(ggplot2)
genus_cols = c("maroon","slateblue1","mediumpurple4","palevioletred","thistle1","mediumorchid3","maroon2","darkorange4","tan1","gold","peachpuff4","coral2","peachpuff","firebrick3","yellow2","khaki1","goldenrod","olivedrab4","olivedrab1","forestgreen","darkolivegreen3","chocolate3","lightsteelblue1")

ggplot(counts_genus_merged2_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Melt the relative count matrix for stacked plot and genus the prey as specified before
counts_genus_merged_rel_melt <- melt(counts_genus_merged_rel, id.vars = "Sample", variable.name = "Phyla")
counts_genus_merged_rel_melt$Var1 <- factor(counts_genus_merged_rel_melt$Var1, levels=Prey_genusing)
counts_genus_merged_rel_melt <- counts_genus_merged_rel_melt[order(counts_genus_merged_rel_melt$Var1),]

# Plot stacked barplot with relative number of reads
genus_cols = c("maroon","slateblue1","mediumpurple4","palevioletred","thistle1","mediumorchid3","maroon2","darkorange4","tan1","gold","peachpuff4","coral2","peachpuff","firebrick3","yellow2","khaki1","goldenrod","olivedrab4","olivedrab1","forestgreen","darkolivegreen3","chocolate3","lightsteelblue1") 

ggplot(counts_genus_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=genus_cols) + 
  labs(y="Rarefied abundance of microbiome phylum (nr of reads)") + 
  guides(fill=guide_legend(title="Microbiome phylum")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Make grid for each host category
# Add metadata to merged melted
counts_genus_merged_rel_melt$Species <-with(metadata_genus0, Species[match(counts_genus_merged_rel_melt$Var2,rownames(metadata_genus0))])
counts_genus_merged_rel_melt$Season <-with(metadata_genus0, Season[match(counts_genus_merged_rel_melt$Var2,rownames(metadata_genus0))])
counts_genus_merged_rel_melt$Location <-with(metadata_genus0, Location[match(counts_genus_merged_rel_melt$Var2,rownames(metadata_genus0))])
counts_genus_merged_rel_melt$Location_detailed <-with(metadata_genus0, Location_detailed[match(counts_genus_merged_rel_melt$Var2,rownames(metadata_genus0))])
counts_genus_merged_rel_melt$Total_reads <-with(metadata_genus, Total_reads[match(counts_genus_merged_rel_melt$Var2,rownames(metadata_genus))])
counts_genus_merged_rel_melt$type <-with(metadata_genus, type[match(counts_genus_merged_rel_melt$Var2,rownames(metadata_genus))])
library(tidyverse)
counts_genus_merged_rel_melt <- counts_genus_merged_rel_melt %>%
  unite("Loc_Season", Location_detailed:Season)

# Order by species
counts_genus_merged_rel_melt$Species<-factor(counts_genus_merged_rel_melt$Species)
counts_genus_merged_rel_melt <- counts_genus_merged_rel_melt[order(counts_genus_merged_rel_melt$Species), ]
levels(counts_genus_merged_rel_melt$Species)
counts_genus_merged_rel_melt <- counts_genus_merged_rel_melt %>% arrange(Species)

# Relative value stacked barplot per detailed location and season
ggplot(counts_genus_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~Loc_Season, scales=("free"),nrow=3) + 
  scale_fill_manual(values=genus_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  guides(fill=guide_legend(title="Prey genus")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Change species names for plot
SpeciesID_genus=as.matrix(counts_genus_merged_rel_melt$Species)
names(SpeciesID_genus)=counts_genus_merged_rel_melt$Var2
SpeciesID_genus[SpeciesID_genus == "Limnothrissa_miodon"] <- "L. miodon"
SpeciesID_genus[SpeciesID_genus == "Stolothrissa_tanganicae"] <- "S. tanganicae"

# Order by host category
Sample_genusing <- factor(rownames(metadata_genus2), levels = rownames(metadata_genus2))
counts_genus_merged_rel_melt$Var2 <- factor(counts_genus_merged_rel_melt$Var2, levels=Sample_genusing)
counts_genus_merged_rel_melt <- counts_genus_merged_rel_melt[order(counts_genus_merged_rel_melt$Var2),]
counts_genus_merged_rel_melt$type <- factor(counts_genus_merged_rel_melt$type, levels=Type_genusing)

# All types arranged in different small plots
ggplot(counts_genus_merged_rel_melt, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity") +
  facet_grid(Species~Loc_Season, scales="free_x", space="free") + 
  scale_fill_manual(values=genus_cols) + 
  labs(y="Relative nr of reads", x= "Specimen") + 
  theme_bw() +
  guides(fill=guide_legend(title="Prey genus")) +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Stacked barplot arranged per host category (type)
ggplot(data=counts_genus_merged_rel_melt, aes(x=Var2, y=value, group=type)) +
  xlab("Sample") + ylab("Relative abundance") + #adds title axes
  scale_fill_manual(values=genus_cols) +
  geom_bar(aes(fill=Var1),stat = "identity") + #adds catagory TissueID by colouring bars
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  theme_classic() +
  ggtitle("Relative Prey abundance (genus level)") +
  theme(axis.text.x = element_text(color = "black", size = 11, face = "plain",angle=-90),
        axis.text.y = element_text(color = "black", size = 11),  
        axis.title.x = element_text(color = "black", size = 15, hjust = .5, vjust = 0),
        axis.title.y = element_text(color = "black", size = 15, vjust = 2),legend.title = element_text(size=15), legend.text = element_text(size=12))

# Nr reads per sample barplot to add
ggplot(data=metadata_genus2, aes(x=row.names(metadata_genus2), y=log10(Total_reads))) +
  geom_bar(aes(fill="Black"),stat = "identity") + 
  xlab("Sample") + ylab("Log10 number of reads") + #adds title axes
  theme_classic() +
  scale_fill_manual(values="black")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~type, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  ggtitle("Log10 number of reads per sample") 

## Bray-Curtis dissimilarity
############################

library(vegan)
library(ecodist)

# Calculate Bray-Curtis
dist.BrCu.all_genus<-vegdist(t(counts_genus_merged2),method="bray")
clust.BrCu.all_genus<-hclust(dist.BrCu.all_genus, method="ward.D") #agglomerative clustering using complete linkage

# Plot dendogram
plot(clust.BrCu.all_genus)

# Make a PCoA based on Bray-curtis dissimilarities (BCD)
# Not done with count matrix since because the data must be linear.
endo.pco.bray_genus<-wcmdscale(dist.BrCu.all_genus)
plot(endo.pco.bray_genus)

# Plot with ggplot
endo.pco.bray_genus_df <- as.data.frame(endo.pco.bray_genus)
ggplot(endo.pco.bray_genus_df, aes(x=V1, y=V2)) +
  geom_point(size=5) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller

# Use BCD in Non Multidimensional scaling (NMDS) plot, works better for datasets with many zeroes
BC.nmds_genus <- metaMDS(t(counts_genus_merged_rel), distance="bray", autotransform = F)
envfit_genus <- envfit(BC.nmds_genus, metadata_genus2, permutations = 999,na.rm = TRUE) # fits environmental vectors
sample.fit_genus <- envfit(BC.nmds_genus, t(counts_genus_merged_rel), permutations = 999) # fits species vectors

# New dataset with point coordinates of NMDS results for ggplot plotting
site.scrs_genus <- as.data.frame(vegan::scores(BC.nmds_genus, display = "sites"))

# Add metadata
site.scrs_genus$Species <-with(metadata_genus2, Species[match(rownames(site.scrs_genus),rownames(metadata_genus2))])
site.scrs_genus$Location_detailed <-with(metadata_genus2, Location_detailed[match(rownames(site.scrs_genus),rownames(metadata_genus2))])
site.scrs_genus$Season <-with(metadata_genus2, Season[match(rownames(site.scrs_genus),rownames(metadata_genus2))])
site.scrs_genus$SL <-with(metadata_genus2, SL[match(rownames(site.scrs_genus),rownames(metadata_genus2))])
site.scrs_genus$Total_reads <-with(metadata_genus2, Total_reads[match(rownames(site.scrs_genus),rownames(metadata_genus2))])
site.scrs_genus$Loc_season <- as.factor(paste(site.scrs_genus$Location_detailed, site.scrs_genus$Season, sep = ""))
head(site.scrs_genus)

# Make dataset with prey vectors, to add to ggplot
spp.scrs_genus <- as.data.frame(vegan::scores(sample.fit_genus, display = "vectors"))
spp.scrs_genus <- cbind(spp.scrs_genus, Species = rownames(spp.scrs_genus)) # add prey names
spp.scrs_genus <- cbind(spp.scrs_genus, pval = sample.fit_genus$vectors$pvals) # add pvalues
sig.spp.scrs_genus_genus <- subset(spp.scrs_genus, pval<=0.05) # shows prey with significant values
head(spp.scrs_genus)

# Make dataset to show extrinsic variables
env.scores_genus<- as.data.frame(vegan::scores(envfit_genus, display = "vectors")) #extracts relevant scores from envifit
env.scores_genus <- cbind(env.scores_genus, env.variables = rownames(env.scores_genus)) #and then gives them their names
env.scores_genus <- cbind(env.scores_genus, pval = envfit_genus$vectors$pvals) # add pvalues
sig.env.scrs_genus <- subset(env.scores_genus, pval<=0.05) # significant
head(env.scores_genus)

# ggplot of NMDS
Locations_shape <- c(24,23,25,6,14,2) # specify shapes plot points
nmds.plot_genus <- ggplot(site.scrs_genus, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = factor(Species),colour = factor(Species), shape = factor(Loc_season)), size = 6)+ 
  coord_fixed()+
  scale_shape_manual(values=Locations_shape) +
  scale_fill_manual(values=c("coral4","lightblue3")) +
  scale_colour_manual(values=c("coral4","lightblue3")) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Species", shape = "Tissue")+ 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))

nmds.plot_genus + labs(title = "Basic ordination plot") 

# Add vectors significant prey
nmds.plot_genus+
  geom_segment(data = sig.spp.scrs_genus_genus, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + # add vectors
  labs(title = "Non-Metric Multidimensional Scaling")

# Make same plot with log10 total nr of reads to check for correlation
ggplot(site.scrs_genus, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = log10(Total_reads),colour=log10(Total_reads),  size = 6, shape=Species))+
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

# Make same plot with standard length to check for correlation
ggplot(site.scrs_genus, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, fill = SL,colour=SL,  size = 6, shape=Species))+ 
  coord_fixed()+
  scale_color_gradientn(colours = hcl.colors(10)) +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot


# Make heatmap based on BCD values
library(pheatmap)

BCD_matrix_genus <- as.matrix(dist.BrCu.all_genus)
pheatmap::pheatmap((1-BCD_matrix_genus), symm = TRUE,clustering_method="ward.D",show_colnames=TRUE)

## Calculate diversity values
#############################

# Calculate Shannon_genus diversity and add to metadata
diversity(counts_genus_merged2, index = "shannon") 
Shannon_genus <- as.data.frame(diversity(t(counts_genus_merged2), index = "shannon"))
Shannon_genus$div <- diversity(t(counts_genus_merged2), index = "shannon")
metadata$Shannon <- with(Shannon_genus, div[match(rownames(metadata),rownames(Shannon_genus))])
metadata_genus2$Shannon <- with(Shannon_genus, div[match(rownames(metadata_genus2),rownames(Shannon_genus))])

# Make new metadata
metadata_genus3 <- metadata_genus2
metadata_genus3 <- metadata_genus3 %>%
  unite("LocSeas",Location,Season)
metadata_genus2$LocSeas <- with(metadata_genus3, LocSeas[match(rownames(metadata_genus2),rownames(metadata_genus3))])

# Plot diversity barplots
ggplot(metadata_genus2, aes(x=LocSeas, y=Shannon,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

# Calculate Chao1 diversity and add to metadata
chao1div_genus <-as.data.frame(t(estimateR(t(counts_genus_merged2))))
chao1div_genus$Chao1 <- as.numeric(chao1div_genus[,2])
metadata_genus2$Chao1 <- with(chao1div_genus, Chao1[match(rownames(metadata_genus2),rownames(chao1div_genus))])

# Plot diversity barplots
ggplot(metadata_genus2, aes(x=LocSeas, y=Chao1,fill=Species)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(0.75),dotsize=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("coral4","lightblue3"))

## Look for host factors influencing prey composition
#####################################################

# PERMANOVA based on BCD
BCD_PERMA_genus<-adonis2(dist.BrCu.all_genus~Species+Location_detailed+Season, data=metadata_genus2, permutations = 999, method="bray")
BCD_PERMA_genus
adonis2(dist.BrCu.all_genus~Species, data=metadata_genus2, permutations = 999, p.adjust.m = "bonferroni")
adonis2(dist.BrCu.all_genus~Species+Location_detailed, data=metadata_genus2, permutations = 999, p.adjust.m = "bonferroni")
adonis2(dist.BrCu.all_genus~Location_detailed, data=metadata_genus2, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA SL and Total reads
adonis2(dist.BrCu.all_genus~Total_reads, data=metadata_genus2, permutations = 999, p.adjust.m = "bonferroni")

no_SL_samples_genus <- rownames(metadata_genus2[which(is.na(metadata_genus2$SL)),])
metadata_SL_genus <- metadata_genus2[ !rownames(metadata_genus2) %in% no_SL_samples_genus, ]
counts_SL_genus <- counts_genus_merged2[ ,!colnames(counts_genus_merged2) %in% no_SL_samples_genus ]
BrCu.SL_genus<-vegdist(t(counts_SL_genus),method="bray")
adonis2(BrCu.SL_genus~SL, data=metadata_SL_genus, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA Samples Uvira
Uvira_samples_genus <- rownames(metadata_genus2[which(metadata_genus2$Location_detailed == "Uvira"),])
metadata_uvira_genus <- metadata_genus2[ rownames(metadata_genus2) %in% Uvira_samples_genus, ]
counts_uvira_genus <- counts_genus_merged2[ ,colnames(counts_genus_merged2) %in% Uvira_samples_genus ]

BrCu.uvira_genus<-vegdist(t(counts_uvira_genus),method="bray")
adonis2(BrCu.uvira_genus~Species+Season, data=metadata_uvira_genus, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.uvira_genus~Species, data=metadata_uvira_genus, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.uvira_genus~Season, data=metadata_uvira_genus, permutations = 999, p.adjust.m = "bonferroni")

# Uvira S.tanganicae
Uvira_ST_samples_genus <- rownames(metadata_uvira_genus[which(metadata_uvira_genus$Species == "Stolothrissa_tanganicae"),])
metadata_uvira_ST_genus <- metadata_uvira_genus[ rownames(metadata_uvira_genus) %in% Uvira_ST_samples_genus, ]
counts_uvira_ST_genus <- counts_uvira_genus[ ,colnames(counts_uvira_genus) %in% Uvira_ST_samples_genus ]

BrCu.uvira_ST_genus<-vegdist(t(counts_uvira_ST_genus),method="bray")
adonis2(BrCu.uvira_ST_genus~Season, data=metadata_uvira_ST_genus, permutations = 999, p.adjust.m = "bonferroni")

# PERMANOVA without Uvira Dry
Uvira_Dry_samples_genus <- rownames(metadata_uvira_genus[which(metadata_uvira_genus$Season == "Dry"),])
counts_wet_genus <- counts_genus_merged2[ ,! colnames(counts_genus_merged2) %in% Uvira_Dry_samples_genus ]
metadata_wet_genus <- metadata_genus2[ ! rownames(metadata_genus2) %in% Uvira_Dry_samples_genus, ]

BrCu.wet_genus<-vegdist(t(counts_wet_genus),method="bray")
adonis2(BrCu.wet_genus~Location_detailed, data=metadata_wet_genus, permutations = 999, method="bonferroni")
adonis2(BrCu.wet_genus~Species+Location_detailed, data=metadata_wet_genus, permutations = 999, method="bonferroni")

adonis2(BrCu.wet_genus~Species+Location_detailed, data=metadata_wet_genus, permutations = 999, p.adjust.m = "bonferroni")
adonis2(BrCu.wet_genus~Species, data=metadata_wet_genus, permutations = 999, p.adjust.m = "bonferroni")

# Wet L.miodon
Wet_LM_samples_genus <- rownames(metadata_wet_genus[which(metadata_wet_genus$Species == "Limnothrissa_miodon"),])
metadata_wet_LM_genus <- metadata_wet_genus[ rownames(metadata_wet_genus) %in% Wet_LM_samples_genus, ]
counts_wet_LM_genus <- counts_wet_genus[ ,colnames(counts_wet_genus) %in% Wet_LM_samples_genus ]

BrCu.wet_LM_genus<-vegdist(t(counts_wet_LM_genus),method="bray")
adonis2(BrCu.wet_LM_genus~Location_detailed, data=metadata_wet_LM_genus,permutations = 999, p.adjust.m = "bonferroni")

# Wet S.tanganicae
Wet_ST_samples_genus <- rownames(metadata_wet_genus[which(metadata_wet_genus$Species == "Stolothrissa_tanganicae"),])
metadata_wet_ST_genus <- metadata_wet_genus[ rownames(metadata_wet_genus) %in% Wet_ST_samples_genus, ]
counts_wet_ST_genus <- counts_wet_genus[ ,colnames(counts_wet_genus) %in% Wet_ST_samples_genus ]

BrCu.wet_ST_genus<-vegdist(t(counts_wet_ST_genus),method="bray")
adonis2(BrCu.wet_ST_genus~Location_detailed, data=metadata_wet_ST_genus, permutations = 999, p.adjust.m = "bonferroni")

# Percentage abundance prey per sardine species
LM_samples_genus <- rownames(metadata_genus2[which(metadata_genus2$Species == "Limnothrissa_miodon"),])
counts_LM_genus <- counts_genus_merged_rel[ ,colnames(counts_genus_merged_rel) %in% LM_samples_genus ]
LM_taxa_genus <- as.matrix(rowSums(counts_LM_genus))
LM_taxa_rel_genus <- LM_taxa_genus/colSums(LM_taxa_genus)

ST_samples_genus <- rownames(metadata_genus2[which(metadata_genus2$Species == "Stolothrissa_tanganicae"),])
counts_ST_genus <- counts_genus_merged_rel[ ,colnames(counts_genus_merged_rel) %in% ST_samples_genus ]
ST_taxa_genus <- as.matrix(rowSums(counts_ST_genus))
ST_taxa_rel_genus <- ST_taxa_genus/colSums(ST_taxa_genus)

## Additional general calculations
##################################

# Percentage abundance prey per sardine species
LM_samples_genus <- rownames(metadata_genus2[which(metadata_genus2$Species == "Limnothrissa_miodon"),])
counts_LM_genus <- counts_genus_merged_rel[ ,colnames(counts_genus_merged_rel) %in% LM_samples_genus ]
LM_taxa_genus <- as.matrix(rowSums(counts_LM_genus))
LM_taxa_rel_genus <- LM_taxa_genus/colSums(LM_taxa_genus)
LM_taxa_rel_genus

ST_samples_genus <- rownames(metadata_genus2[which(metadata_genus2$Species == "Stolothrissa_tanganicae"),])
counts_ST_genus <- counts_genus_merged_rel[ ,colnames(counts_genus_merged_rel) %in% ST_samples_genus ]
ST_taxa_genus <- as.matrix(rowSums(counts_ST_genus))
ST_taxa_rel_genus <- ST_taxa_genus/colSums(ST_taxa_genus)
ST_taxa_rel_genus

# Nr Reads
mean(metadata_genus2$Total_reads)
median(metadata_genus2$Total_reads)
min(metadata_genus2$Total_reads)
sd(metadata_genus2$Total_reads)
max(metadata_genus2$Total_reads)

# Mean diversity
metadata_ST_genus <- metadata_genus2[ rownames(metadata_genus2) %in% ST_samples_genus , ]
mean(metadata_ST_genus$Shannon)
mean(metadata_ST_genus$Chao1)

metadata_LM_genus <- metadata_genus2[ rownames(metadata_genus2) %in% LM_samples_genus , ]
mean(metadata_LM_genus$Shannon)
mean(metadata_LM_genus$Chao1)

#Take relative values
tot_reads_genus_samples <- colSums(counts_genus_merged2)
genus_samples <- colnames(counts_genus_merged2)

counts_genus_merged_rel <- t(t(counts_genus_merged2) / tot_reads_genus_samples)

# Total counts per sample
reads_genus <- as.data.frame(rowSums(t(counts_genus_merged2)))
colnames(reads_genus) <- 'reads'
metadata3$Total_reads_d_genus <- with(reads_genus, reads[match(rownames(metadata3),rownames(reads_genus))])

# Do same as with order level
counts_order_merged_rel <- counts_genus_merged_rel
counts_order_merged2 <- counts_genus_merged2




# Species
###################

# Make subsets for Species level ID
counts_species <- as.matrix(counts)
counts_species <- counts_species[ order((row.names(counts_species))), ]
otus_species <- OTUs[ order((row.names(OTUs))), ]

# Remove OTU's which were not identified till species level
nospecies_otu <- rownames(otus_species[which(otus_species$Similarity < 98 ),])
otus_species <- otus_species[ !rownames(otus_species) %in% nospecies_otu, ]
counts_species <- counts_species[ !rownames(counts_species) %in% nospecies_otu, ]
rmSamplesSpecies <- colnames(counts_species[,colSums(counts_species) == 0])
counts_species <- counts_species[ ,!colnames(counts_species) %in% rmSamplesSpecies ]

# Merge samples
otus_species$Species <- as.factor(otus_species$Species)
counts_spp_table <- otu_table(as.matrix(counts_species), taxa_are_rows=F)
counts_spp_merged <- merge_samples(counts_spp_table, t(otus_species$Species))
counts_spp_merged2 <- as(otu_table(counts_spp_merged), "matrix")

# Remove samples low counts
counts_spp_merged2 = counts_spp_merged2[,colSums(counts_spp_merged2) > 10]
counts_spp_merged2 = counts_spp_merged2[rowSums(counts_spp_merged2) > 10,]

#Take relative values
tot_reads_spp_samples <- colSums(counts_spp_merged2)
spp_samples <- colnames(counts_spp_merged2)

counts_spp_merged_rel <- t(t(counts_spp_merged2) / tot_reads_spp_samples)

# Total counts per sample
reads_spp <- as.data.frame(rowSums(t(counts_spp_merged2)))
colnames(reads_spp) <- 'reads'

# Check for L. tanganjicae and T. simplex
counts_spp_prey <- as.data.frame(t(counts_spp_merged2))
Trosim = counts_spp_prey[counts_spp_prey$simplex > 0,] # In 20 samples
Limtan = counts_spp_prey[counts_spp_prey$tanganjicae > 0,] # In 16 samples to species level
 
# Nr reads all 
##############
metadata3$Total_reads_a_class <- with(metadata_class3, Total_reads[match(rownames(metadata3),rownames(metadata_class3))])
metadata3$Total_reads_b_order <- metadata3$Total_reads
metadata3$Total_reads_c_family <- with(metadata_family3, Total_reads[match(rownames(metadata3),rownames(metadata_family3))])
metadata3$Total_reads_d_genus <- with(metadata_genus3, Total_reads[match(rownames(metadata3),rownames(metadata_genus3))])
metadata3$Total_reads_e_spp <- with(reads_spp, reads[match(rownames(metadata3),rownames(reads_spp))])

reads_counts <- as.data.frame(rowSums(t(counts0))) #counts has to be the original one
colnames(reads_counts) <- 'reads'
metadata3$Total_reads_all_otus <- with(reads_counts, reads[match(rownames(metadata3),rownames(reads_counts))])


Reads <- subset(metadata3, select = c(Total_reads_all_otus,Total_reads_a_class,Total_reads_b_order,Total_reads_c_family,Total_reads_d_genus,Total_reads_e_spp) )

Reads_table <- otu_table(as.matrix(Reads), taxa_are_rows=F)
Reads_melt <- melt(Reads_table, id.vars = "Sample", variable.name = "Rank") # Melts samples for abundance bars

reads_col <- c("gray","firebrick","chocolate","darkgoldenrod1","darkgreen","turquoise")

ggplot(data=Reads_melt, aes(x=Var2, y=log10(value), group=Var1)) +
  xlab("Sample") + ylab("Log10 number of reads") + #adds title axes
  scale_fill_manual(values=reads_col) + 
  geom_bar(aes(fill=Var2),stat = "identity") + #adds catagory TissueID by colouring bars
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3),axis.ticks.y = element_text(angle = 90, hjust = 1,size=7)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~Var1, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  theme_classic() +
  ggtitle("Number of total reads per rank")

ggplot(data=Reads_melt, aes(x=Var2, y=value, group=Var1)) +
  xlab("Sample") + ylab("Number of reads") + #adds title axes
  scale_fill_manual(values=reads_col) + 
  geom_bar(aes(fill=Var2),stat = "identity") + #adds catagory TissueID by colouring bars
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3),axis.ticks.y = element_text(angle = 90, hjust = 1,size=7)) + #signifies how you want the text on the x-axis, so it doesn't overlap
  facet_grid(.~Var1, scales = "free", switch = "x", space = "free_x") + #groups TissueTubeID by SpeciesID
  theme(strip.placement = "outside") +
  theme_classic() +
  ggtitle("Number of total reads per rank")








