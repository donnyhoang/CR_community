library(phyloseq)
library(vegan)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(ggpubr)
library(viridis)


#read in files
tax <- read.csv("CR_taxa_for_phyloseq_nosupportedit.csv", header=TRUE, row.names=1)
counts <- read.csv("CR_counts_for_phyloseq.csv", header=TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("CR_metadata_for_phyloseq.csv", header=TRUE, row.names= 1)
tax$Bin <- row.names(tax)

#create phyloseq
count.ps <- otu_table(as.matrix(counts),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps <- tax_table(as.matrix(tax)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps <- sample_data(metadata) #define our metadata as sample information
ps <- phyloseq(count.ps, tax.ps, meta.ps) #put them all together into a phyloseq object


#transform seq counts into %
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))


#write out table from ps so I can more easily wrangle


new_counts <- as.data.frame(t(otu_table(ps)))
#new_counts$variable <- row.names(new_counts)
#row.names(new_counts) <- NULL
#metadata$variable <- row.names(metadata)
#row.names(metadata) <- NULL
#new_data <- merge(new_counts, metadata, by = "variable")
#new_data <- new_data[,c(1, 56:58,2:55)]

#ordination
attach(metadata)
#ordinate in 3 dimensions bc it cannot find solution with 2
ordination.model <- metaMDS(new_counts, distance='bray', k=3)

#extract coordinates so I can plot with ggplot2
nmds_coords <- ordination.model$points

data <- cbind(metadata, nmds_coords)
row.names(data) <- NULL
data$Treatment<-factor(data$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
palette <-c( "#56B4E9", "#999999","#E69F00","#0072B2", "#F0E442","#CC79A7")

pord <- ggplot() +
  geom_point(data=data,
             aes(x=MDS1, y=MDS2, color=Treatment, bg=Treatment),
             size=4) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  coord_fixed(ratio=1) +
  ggtitle("Bray-Curtis Dissimilarity of Community Expression") +
  scale_color_manual(values=palette) +
  scale_fill_manual(values=palette) 
pord

##############stats
dist <- vegdist(new_counts) #create distance matrix
attach(metadata) #attach metadata


#anosim
anoTreatment <- anosim(dist, Treatment) #similarity analysis
summary(anoTreatment) # see summary
plot(anoTreatment)
anoDay <- anosim(dist, Day)
anoRep <-anosim(dist, Replicate)


summary(anoTreatment)
summary(anoDay)
summary(anoRep)

#adonis
adTreatment <- adonis2(dist ~ Treatment, data=metadata, permutations=99, method="bray")
adDay <- adonis2(dist ~ Day, data=metadata,  permutations=99, method="bray")
adReplicate <- adonis2(dist ~ Replicate, data=metadata,  permutations=99, method="bray")

adTreatment
adDay
adReplicate

############# diversity #########

data <- cbind(new_counts, metadata)
data <- data[,c(55:57, 1:54)]
data$Sample_Name <- row.names(data)
row.names(data) <- NULL
data <- data[,c(58,1:3,4:57)]


###ALPHA DIVERSITY###

#observed
count <- ddply(data,~Sample_Name,function(x){
  data.frame(Count=sum(x[,5:58]>0))
})
#count_list<-list(metadata,richness)



menhinick <- function(x) {
  sum(x>0)/sqrt(sum(x))
}


richness <- ddply(data,~Sample_Name,function(x){
  data.frame(Richness=menhinick(x[,5:58]>0))
})
#abundance


shan <- ddply(data,~Sample_Name,function(x){
  data.frame(Shannon=diversity(x[,5:58],index="shannon"))
})
shan
#shan_list<-list(metadata,shan)



invsimpson <- ddply(data,~Sample_Name,function(x){
  data.frame(InverseSimpson=diversity(x[,5:58],index="invsimpson"))
})
invsimpson
#in_list<-list(metadata,invsimpson)

#evenness
pielou <-ddply(data,~Sample_Name,function(x){
  data.frame(Pielou=exp(diversity(x[,5:58], index="shannon"))/sum(x[,5:58]>0))
})
pielou
#pie_list<-list(metadata,pielou)

metadata <- data[,c(1:4)]
metadata <- merge(metadata, richness, by="Sample_Name")
metadata <- merge(metadata, shan, by="Sample_Name")
metadata <- merge(metadata, invsimpson, by="Sample_Name")
metadata <- merge(metadata, pielou, by="Sample_Name")
metadata <- merge(metadata, count, by="Sample_Name")


metadata$Treatment<-factor(metadata$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))

pshan <-ggboxplot(metadata, x="Day",y="Shannon", fill="Treatment", lwd=1) +
  theme_classic(base_size=22)+
  #theme(aspect.ratio=1) +
  scale_fill_manual(values=palette) + 
  labs(y="Shannon Index", x="Day") +
  #stat_compare_means(method="t.test", aes(group=Treatment), label="p.signif") +
  ggtitle("Shannon Diversity")
pshan

psimp <-ggboxplot(metadata, x="Day",y="InverseSimpson", fill="Treatment", lwd=1) +
  theme_classic(base_size=22)+
  #theme(aspect.ratio=1) +
  scale_fill_manual(values=palette) + 
  labs(y="Inverse Simpson Index", x="Day") + 
  #stat_compare_means(method="t.test", aes(group=Treatment), label="p.signif") +
  ggtitle("Inverse Simpson Diversity")
psimp

ggarrange(pshan, psimp, ncol=1, common.legend = TRUE, legend = "bottom")

######### abundance plots

#relative expression
#remake new_counts without relative abundance so I can calc % myself
data_m <- melt(data, id=c("Sample_Name", "Treatment", "Replicate", "Day"))
colnames(data_m)[5] <- "Bin"
row.names(tax) <- NULL
data_m <- merge(data_m, tax, by = "Bin")

tax2 <- read.csv("CR_bin_tax_cluster.csv", header = TRUE)
tax2 <- tax2[,c(1,2)]
colnames(tax2)[1] <- "Bin"
data_m <- merge(data_m, tax2, by = "Bin")


counts2 <- fread("total_reads_mt.txt", sep="\t")
colnames(counts2) <- c("Sample_Name", "counts")
counts2$Sample_Name <- gsub("_mt.R1.fastq.cleaned.forward","", as.character(counts2$Sample_Name))

data_m <- merge(data_m, counts2, by = "Sample_Name")

data_m$percent <- ((data_m$value)/(data_m$counts))*100

#data_m$bin_variable <- paste(data_m$Bin, data_m$Sample_Name, sep = ".")
#data_m <- data_m[,c(16,17)]
#colnames(data_m)[1] <- "percent_mt"
#write.csv(data_m, "CR_perentmtabundance.csv", row.names = FALSE)

data_m$Treatment<-factor(data_m$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))

data_m$cluster_total <- as.factor(data_m$cluster_total)

data_m$Order <- factor(data_m$Order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                      "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                      "Sphingobacteriales","Xanthomonadales","no support"))

palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")




pbar <- ggplot(data=data_m, aes(x=Sample_Name, y = percent, color = Order, fill = Order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_blank()) +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg) +
  labs(y = "Percent Metatranscriptomtic Reads Mapped", x = "Samples Arranged by Time") +
  ggtitle("Relative Expression") +
  facet_grid (~ Treatment, scales = "free")
pbar

ptile <- ggplot(data=data_m, aes(x=Sample_Name, y = Bin)) +
  geom_tile(aes(fill=percent)) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_blank()) +
  scale_fill_viridis(option = "inferno") +
  labs(y = "Assigned order", x = "Samples Arranged by Time", fill = "Percent Reads Mapped") +
  ggtitle("Relative Expression of Bacterial Orders") +
  facet_grid (~ Treatment, scales = "free")

ptile

#relative abundance