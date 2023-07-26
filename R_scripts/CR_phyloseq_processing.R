library(phyloseq)
library(stringr)
library(dplyr)
library(tidyverse)


###############load in data, prepare for phyloseq
data <- read.csv("CR_taxmatrix.csv", row.names=1 , header=TRUE)
write.csv(data, "CR_counts_for_phyloseq.csv", row.names = TRUE)


#first tax table
tax1 <- read.csv("rep_bins_autometa_qa.csv")
colnames(tax1)[1] <- "bin"
tax1[c("classification","UID")] <- str_split_fixed(tax1$Markerlineage , ",", 2)
tax1$UID <- gsub("[(]","",as.character(tax1$UID))
tax1$UID <- gsub("[)]","",as.character(tax1$UID))
tax1 <- tax1[,c(1,2,30,31)]

#second tax table
tax2 <- read.table("out.BAT.bin2classification.official_names.txt", sep="\t", header=TRUE)
tax2$bin <- gsub(".fasta","",as.character(tax2$bin))
tax2 <- tax2[,c(1,6:12)]

#merge to use tax as key for later
tax <- merge(tax1, tax2, by = "bin")
tax <- tax[,c(1,5:11)]
rownames(tax) <- tax$bin
tax <- tax[,-1]

tax <- tax %>%
  rename_with(str_to_title)

tax$Superkingdom <- gsub("[0-9]+", "", tax$Superkingdom)
tax$Superkingdom <- gsub(":", "", tax$Superkingdom)
tax$Superkingdom <- gsub("[.]", "", tax$Superkingdom)
tax$Superkingdom <- gsub("[[:space:]]", "", tax$Superkingdom)

tax$Phylum <- gsub("[0-9]+", "", tax$Phylum)
tax$Phylum <- gsub(":", "", tax$Phylum)
tax$Phylum <- gsub("[.]", "", tax$Phylum)
tax$Phylum <- gsub("[[:space:]]", "", tax$Phylum)

tax$Class <- gsub("[0-9]+", "", tax$Class)
tax$Class <- gsub(":", "", tax$Class)
tax$Class <- gsub("[.]", "", tax$Class)
tax$Class <- gsub("[[:space:]]", "", tax$Class)

tax$Order <- gsub("[0-9]+", "", tax$Order)
tax$Order <- gsub(":", "", tax$Order)
tax$Order <- gsub("[.]", "", tax$Order)
tax$Order <- gsub("[[:space:]]", "", tax$Order)

tax$Family <- gsub("[0-9]+", "", tax$Family)
tax$Family <- gsub(":", "", tax$Family)
tax$Family <- gsub("[.]", "", tax$Family)
tax$Family <- gsub("[[:space:]]", "", tax$Family)

tax$Genus <- gsub("[0-9]+", "", tax$Genus)
tax$Genus <- gsub(":", "", tax$Genus)
tax$Genus <- gsub("[.]", "", tax$Genus)
tax$Genus <- gsub("[[:space:]]", "", tax$Genus)

tax$Species <- gsub("[0-9]+", "", tax$Species)
tax$Species <- gsub("[.]", "", tax$Species)

#save and edit the rest in excel. removing whtiespace from species can be confusing for later, so do it in excel.
write.csv(tax, "CR_taxa_for_phyloseq.csv", row.names = TRUE)

#create metadata table

metadata <- as.data.frame(colnames(data))
metadata <- as.data.frame(metadata[-1,])
colnames(metadata)[1] <- "Sample"
metadata[c("Treatment","Day","Replicate")] <- str_split_fixed(metadata$Sample , "_", 3)
rownames(metadata) <-metadata$Sample
metadata <- metadata[,-1]
write.csv(metadata, "CR_metadata_for_phyloseq.csv", row.names = TRUE)

