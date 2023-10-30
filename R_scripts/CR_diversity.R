library(phyloseq)
library(vegan)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(ggpubr)
library(stringr)
library(viridis)
library(cowplot)
library(ggpmisc)

############ EXPRESSION ################## 

###ALPHA DIVERSITY###


menhinick <- function(x) {
  sum(x>0)/sqrt(sum(x))
}

############## CAZY EXPRESSION ###############
#read in files

counts <- read.csv("CR_dbCAN_funcmatrix_cazy.csv", header=TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("CR_metadata_for_phyloseq.csv", header=TRUE, row.names= 1)
tax <- as.data.frame(row.names(counts))
colnames(tax)[1] <- "cazy"
rownames(tax) <- tax$cazy


#create phyloseq
count.ps <- otu_table(as.matrix(counts),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps <- tax_table(as.matrix(tax)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps <- sample_data(metadata) #define our metadata as sample information
ps <- phyloseq(count.ps, tax.ps, meta.ps) #put them all together into a phyloseq object


#transform seq counts into %
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))


#write out table from ps so I can more easily wrangle

new_counts <- as.data.frame(t(otu_table(ps)))

data <- cbind(metadata, new_counts)
data$Sample_Name <- row.names(data)
row.names(data) <- NULL
data <- data[,c(519,1:3,4:518)]


menhinick_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Richness_func_MT=menhinick(x[,5:518]>0))
})


shan_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Shannon_func_MT=diversity(x[,5:518],index="shannon"))
})




invsimpson_func <- ddply(data,~Sample_Name,function(x){
  data.frame(InverseSimpson_func_MT=diversity(x[,5:518],index="invsimpson"))
})


#evenness
pielou_func <-ddply(data,~Sample_Name,function(x){
  data.frame(Pielou_func_MT=exp(diversity(x[,5:518], index="shannon"))/sum(x[,5:518]>0))
})





#### merge #####

metadata <- data[,c(1:4)]
metadata <- merge(metadata, menhinick_func, by="Sample_Name")
metadata <- merge(metadata, shan_func, by="Sample_Name")
metadata <- merge(metadata, invsimpson_func, by="Sample_Name")
metadata <- merge(metadata, pielou_func, by="Sample_Name")

#write.csv(metadata, "CR_dbcan_orgvfunc_diversity_mt.csv", row.names = FALSE)

metadata$Treatment<-factor(metadata$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
palette <-c( "#56B4E9", "#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")



#################### MT Reads UNFILTERED (ribosomal reads were kept) ##################

counts <- read.csv("CR_matrix_tax_mt_unfiltered.csv", header=TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("CR_metadata_for_phyloseq.csv", header=TRUE, row.names= 1)
tax <- as.data.frame(row.names(counts))
colnames(tax)[1] <- "bin"
rownames(tax) <- tax$bin


#create phyloseq
count.ps <- otu_table(as.matrix(counts),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps <- tax_table(as.matrix(tax)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps <- sample_data(metadata) #define our metadata as sample information
ps <- phyloseq(count.ps, tax.ps, meta.ps) #put them all together into a phyloseq object


#transform seq counts into %
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))


#write out table from ps so I can more easily wrangle

new_counts <- as.data.frame(t(otu_table(ps)))

data <- cbind(metadata, new_counts)
data$Sample_Name <- row.names(data)
row.names(data) <- NULL
data <- data[,c(58,1:3,4:57)]


menhinick_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Richness_MT=menhinick(x[,5:58]>0))
})



shan_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Shannon_MT=diversity(x[,5:58],index="shannon"))
})




invsimpson_func <- ddply(data,~Sample_Name,function(x){
  data.frame(InverseSimpson_MT=diversity(x[,5:58],index="invsimpson"))
})


#evenness
pielou_func <-ddply(data,~Sample_Name,function(x){
  data.frame(Pielou_MT=exp(diversity(x[,5:58], index="shannon"))/sum(x[,5:58]>0))
})

### read in cazy and merge
cazy <- read.csv("CR_dbcan_orgvfunc_diversity_mt.csv", header = TRUE)

metadata <- merge(cazy, menhinick_func, by="Sample_Name")
metadata <- merge(metadata, shan_func, by="Sample_Name")
metadata <- merge(metadata, invsimpson_func, by="Sample_Name")
metadata <- merge(metadata, pielou_func, by="Sample_Name")

metadata_mt <- metadata


############ CAZYME ABUNDANCE ################## 


############## CAZY
#read in files

counts <- read.csv("CR_dbCAN_funcmatrix_cazy_mg.csv", header=TRUE, row.names = 1, check.names = FALSE)
metadata <- as.data.frame(colnames(counts))
colnames(metadata)[1] <- "variable"
metadata[c("Treatment", "Day", "Replicate")] <- str_split_fixed(metadata $variable, "_", 3)
row.names(metadata) <- metadata$variable
tax <- as.data.frame(row.names(counts))
colnames(tax)[1] <- "cazy"
rownames(tax) <- tax$cazy



#create phyloseq
count.ps <- otu_table(as.matrix(counts),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps <- tax_table(as.matrix(tax)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps <- sample_data(metadata) #define our metadata as sample information
ps <- phyloseq(count.ps, tax.ps, meta.ps) #put them all together into a phyloseq object


#transform seq counts into %
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))


#write out table from ps so I can more easily wrangle

new_counts <- as.data.frame(t(otu_table(ps)))

data <- cbind(metadata, new_counts)
data$Sample_Name <- row.names(data)
row.names(data) <- NULL
data <- data[,c(522,1:4,5:521)]


menhinick_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Richness_func_MG=menhinick(x[,6:522]>0))
})



shan_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Shannon_func_MG=diversity(x[,6:522],index="shannon"))
})




invsimpson_func <- ddply(data,~Sample_Name,function(x){
  data.frame(InverseSimpson_func_MG=diversity(x[,6:522],index="invsimpson"))
})


#evenness
pielou_func <-ddply(data,~Sample_Name,function(x){
  data.frame(Pielou_func_MG=exp(diversity(x[,6:522], index="shannon"))/sum(x[,6:522]>0))
})





#### merge #####

metadata <- data[,c(1,3:5)]
metadata <- merge(metadata, menhinick_func, by="Sample_Name")
metadata <- merge(metadata, shan_func, by="Sample_Name")
metadata <- merge(metadata, invsimpson_func, by="Sample_Name")
metadata <- merge(metadata, pielou_func, by="Sample_Name")

#write.csv(metadata, "CR_dbcan_orgvfunc_diversity_mg.csv", row.names = FALSE)


#################### MG Reads ##################

counts <- t(read.csv("CR_matrix_tax_mg.csv", header=TRUE, row.names = 1, check.names = FALSE))
counts <- as.data.frame(counts)
metadata <- as.data.frame(colnames(counts))
colnames(metadata)[1] <- "variable"
metadata[c("Treatment", "Day", "Replicate")] <- str_split_fixed(metadata $variable, "_", 3)
row.names(metadata) <- metadata$variable
tax <- as.data.frame(row.names(counts))
colnames(tax)[1] <- "bin"
rownames(tax) <- tax$bin


#create phyloseq
count.ps <- otu_table(as.matrix(counts),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps <- tax_table(as.matrix(tax)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps <- sample_data(metadata) #define our metadata as sample information
ps <- phyloseq(count.ps, tax.ps, meta.ps) #put them all together into a phyloseq object


#transform seq counts into %
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))


#write out table from ps so I can more easily wrangle

new_counts <- as.data.frame(t(otu_table(ps)))

data <- cbind(metadata, new_counts)
data$Sample_Name <- row.names(data)
row.names(data) <- NULL
data <- data[,c(59,2:4,5:58)]


menhinick_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Richness_MG=menhinick(x[,5:58]>0))
})



shan_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Shannon_MG=diversity(x[,5:58],index="shannon"))
})




invsimpson_func <- ddply(data,~Sample_Name,function(x){
  data.frame(InverseSimpson_MG=diversity(x[,5:58],index="invsimpson"))
})


#evenness
pielou_func <-ddply(data,~Sample_Name,function(x){
  data.frame(Pielou_MG=exp(diversity(x[,5:58], index="shannon"))/sum(x[,5:58]>0))
})

### read in cazy and merge

cazy <- read.csv("CR_dbcan_orgvfunc_diversity_mg.csv", header = TRUE)

metadata <- merge(cazy, menhinick_func, by="Sample_Name")
metadata <- merge(metadata, shan_func, by="Sample_Name")
metadata <- merge(metadata, invsimpson_func, by="Sample_Name")
metadata <- merge(metadata, pielou_func, by="Sample_Name")

metadata_mg <- metadata

metadata_mg <- metadata_mg[,-c(2:4)]

metadata <- merge(metadata_mt, metadata_mg, by = "Sample_Name")

#write.csv(metadata, "CR_dbcan_diversity_calculations.csv", row.names = FALSE)

#Plot
metadata$Treatment<-factor(metadata$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
metadata$Day <- as.character(metadata$Day)
palette <-c( "#56B4E9", "#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")
palday <- c("#7DCEA0","#52BE80","#27AE60","#229954","#1E8449","#196F3D","#145A32")
#palday <- c("#bdbdbd","#969696","#737373","#525252","#252525","#111111", "#000000")


my.formula <- y ~ x

shan_plot_mt <- ggplot(metadata, aes(x= Shannon_MT, y = Shannon_func_MT)) +
  geom_smooth(method=lm, se = FALSE, col = 'black', linetype = "dashed", formula = my.formula) +
  #stat_poly_line(linetype = "dotted", se = F, color = "black") +
  stat_poly_eq(use_label(c("eq","R2")),
               vjust = 0.5) +
  geom_point(aes(color = Day)) +
  theme_classic(base_size=15) +
  ylim(2,5.3) +
  labs(x= "Shannon Diversity, MAGs", y = "Shannon Diversity, CAZymes") +
  scale_color_manual(values = palday) +
  facet_grid(~Treatment) +
  ggtitle("CAZY Expression Shannon Diversity vs. MAG Expression Shannon Diversity")
shan_plot_mt



shan_plot_mg <- ggplot(metadata, aes(x= Shannon_MG, y = Shannon_func_MG)) +
  geom_smooth(method=lm, se = FALSE, col = 'black', linetype = "dashed", formula = my.formula) +
  #stat_poly_line(linetype = "dotted", se = F, color = "black") +
  stat_poly_eq(use_label(c("eq","R2")),
               vjust = 0.5) +
  geom_point(aes(color = Day)) +
  theme_classic(base_size=15) +
  ylim(2,5.3) +
  labs(x= "Shannon Diversity, MAGs", y = "Shannon Diversity, CAZymes") +
  scale_color_manual(values = palday) +
  facet_grid(~Treatment) +
  ggtitle("CAZY Abundance Shannon Diversity vs. MAG Abundance Shannon Diversity")
shan_plot_mg




pshannon <- ggarrange(shan_plot_mg, shan_plot_mt, ncol = 1,
                      common.legend=TRUE, legend = "bottom",
                      labels = "AUTO", font.label = list(size=25))
pshannon

ggsave("pshannon.tiff", device = "tiff", dpi = 700)
ggsave("pshannon.png", device = "png", dpi = 700)
ggsave("pshannon.pdf", device = "pdf", dpi = 700)

##calc R2, intercept, slope
#https://rpubs.com/ftrecca/Intercept_and_slope

ph6 <- subset(metadata, Treatment == "pH6")
ph8 <- subset(metadata, Treatment == "pH8")
ph10 <- subset(metadata, Treatment == "pH10")
mes <- subset(metadata, Treatment == "MES")
g05 <- subset(metadata, Treatment == "G05")
g10 <- subset(metadata, Treatment == "G10")

sum_ph6_MT <- summary(lm(Shannon_func_MT ~ Shannon_MT, data = ph6))
sum_ph8_MT <- summary(lm(Shannon_func_MT ~ Shannon_MT, data = ph8))
sum_ph10_MT <- summary(lm(Shannon_func_MT ~ Shannon_MT, data = ph10))
sum_mes_MT <- summary(lm(Shannon_func_MT ~ Shannon_MT, data = mes))
sum_g05_MT <- summary(lm(Shannon_func_MT ~ Shannon_MT, data = g05))
sum_g10_MT <- summary(lm(Shannon_func_MT ~ Shannon_MT, data = g10))

sum_ph6_MT
sum_ph8_MT
sum_ph10_MT
sum_mes_MT
sum_g05_MT
sum_g10_MT


sum_ph6_MG <- summary(lm(Shannon_func_MG ~ Shannon_MG, data = ph6))
sum_ph8_MG <- summary(lm(Shannon_func_MG ~ Shannon_MG, data = ph8))
sum_ph10_MG <- summary(lm(Shannon_func_MG ~ Shannon_MG, data = ph10))
sum_mes_MG <- summary(lm(Shannon_func_MG ~ Shannon_MG, data = mes))
sum_g05_MG <- summary(lm(Shannon_func_MG ~ Shannon_MG, data = g05))
sum_g10_MG <- summary(lm(Shannon_func_MG ~ Shannon_MG, data = g10))

sum_ph6_MG
sum_ph8_MG
sum_ph10_MG
sum_mes_MG
sum_g05_MG
sum_g10_MG



######## hillR functional diversity ####

library(hillR)



data <- read.csv("CR_dbcan_total_table.csv", header = TRUE)

comm <- dcast(data, formula = variable ~ bin, fun.aggregate = sum, value.var = "percent")
row.names(comm) <- comm$variable
comm <- comm[,-c(1)]


traits <- dcast(data, formula = bin ~ CAZy, fun.aggregate = sum, value.var = "percent")
row.names(traits) <- traits$bin
traits <- traits[,-c(1)]


fd0_mt <- as.data.frame(t(hill_func(comm, traits, q = 0)))
fd0_mt$Hill <- "Hill0"
fd1_mt <- as.data.frame(t(hill_func(comm, traits, q = 1)))
fd1_mt$Hill <- "Hill1"
fd2_mt <- as.data.frame(t(hill_func(comm, traits, q = 2)))
fd2_mt$Hill <- "Hill2"


fd0_mt$variable <- row.names(fd0_mt)
row.names(fd0_mt) <- NULL
fd1_mt$variable <- row.names(fd1_mt)
row.names(fd1_mt) <- NULL
fd2_mt$variable <- row.names(fd2_mt)
row.names(fd2_mt) <- NULL

fd_mt <- rbind(fd0_mt, fd1_mt, fd2_mt)
fd_mt$Type <- "MT"
fd_mt[c("Treatment", "Day", "Replicate")] <- str_split_fixed(fd_mt $variable, "_", 3)

data2 <- read.csv("CR_dbcan_total_table_mg.csv", header = TRUE)


comm2 <- dcast(data2, formula = variable ~ bin, fun.aggregate = sum, value.var = "percent")
row.names(comm2) <- comm2$variable
comm2 <- comm2[,-c(1)]


traits2 <- dcast(data2, formula = bin ~ CAZy, fun.aggregate = sum, value.var = "percent")
row.names(traits2) <- traits2$bin
traits2 <- traits2[,-c(1)]


fd0_mg <- as.data.frame(t(hill_func(comm2, traits2, q = 0)))
fd0_mg$Hill <- "Hill0"
fd1_mg <- as.data.frame(t(hill_func(comm2, traits2, q = 1)))
fd1_mg$Hill <- "Hill1"
fd2_mg <- as.data.frame(t(hill_func(comm2, traits2, q = 2)))
fd2_mg$Hill <- "Hill2"

fd0_mg$variable <- row.names(fd0_mg)
row.names(fd0_mg) <- NULL
fd1_mg$variable <- row.names(fd1_mg)
row.names(fd1_mg) <- NULL
fd2_mg$variable <- row.names(fd2_mg)
row.names(fd2_mg) <- NULL


fd_mg <- rbind(fd0_mg, fd1_mg, fd2_mg)
fd_mg$Type <- "MG"
fd_mg[c("Treatment", "Day", "Replicate")] <- str_split_fixed(fd_mg $variable, "_", 3)


fd <- rbind(fd_mt, fd_mg)

#functional diversity = total functional diversity
p_fd <- ggplot(fd, aes(x = Treatment, y = FD_q, color = Type)) +
  geom_boxplot() +
  facet_grid(~Hill) +
  stat_compare_means(method = "t.test", label = "p.signif")
p_fd




hill1_mt <- as.data.frame(hill_taxa(comm, q = 1))
hill1_mt$variable <- row.names(hill1_mt)
rownames(hill1_mt) <- NULL
#hill1_mt[c("Treatment", "Day", "Replicate")] <- str_split_fixed(hill1_mt $variable, "_", 3)
hill1_mt$Type <- "MT"
colnames(hill1_mt)[1] <- "Hill1"

hill1_mg <- as.data.frame(hill_taxa(comm2, q = 1))
hill1_mg$variable <- row.names(hill1_mg)
rownames(hill1_mg) <- NULL
#hill1_mg[c("Treatment", "Day", "Replicate")] <- str_split_fixed(hill1_mg $variable, "_", 3)
hill1_mg$Type <- "MG"
colnames(hill1_mg)[1] <- "Hill1"

hill_tax <- rbind(hill1_mt, hill1_mg)
hill_tax$variable_Type <- paste(hill_tax$variable, hill_tax$Type, sep = "_")
hill_tax <- hill_tax[,c(1,4)]
hill1_fd <- subset(fd, Hill == "Hill1")
hill1_fd$variable_Type <- paste(hill1_fd$variable, hill1_fd$Type, sep = "_")

hill_total <- merge(hill_tax, hill1_fd, by = "variable_Type")


p_hill <- ggplot(hill_total, aes(x = Hill1, y = FD_q, color = Day)) +
  geom_point() +
  facet_grid(Type~Treatment) 
p_hill


############ SEED ################## 

############ organism based #########
#read in files
tax <- read.csv("CR_taxa_for_phyloseq_nosupportedit.csv", header=TRUE, row.names=1)
counts <- read.csv("CR_matrix_tax_mt.csv", header=TRUE, row.names = 1, check.names = FALSE)
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



############# diversity #########

data <- cbind(new_counts, metadata)
data <- data[,c(55:57, 1:54)]
data$Sample_Name <- row.names(data)
row.names(data) <- NULL
data <- data[,c(58,1:3,4:57)]


###ALPHA DIVERSITY###


menhinick_org <- ddply(data,~Sample_Name,function(x){
  data.frame(Richness_org=menhinick(x[,5:58]>0))
})



shan_org <- ddply(data,~Sample_Name,function(x){
  data.frame(Shannon_org=diversity(x[,5:58],index="shannon"))
})




invsimpson_org <- ddply(data,~Sample_Name,function(x){
  data.frame(InverseSimpson_org=diversity(x[,5:58],index="invsimpson"))
})


#evenness
pielou_org <-ddply(data,~Sample_Name,function(x){
  data.frame(Pielou_org=exp(diversity(x[,5:58], index="shannon"))/sum(x[,5:58]>0))
})





############## function diversity
#read in files

counts <- read.csv("CR_seed_matrix.csv", header=TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("CR_metadata_for_phyloseq.csv", header=TRUE, row.names= 1)
tax <- as.data.frame(row.names(counts))
colnames(tax)[1] <- "seed"
rownames(tax) <- tax$seed


#create phyloseq
count.ps <- otu_table(as.matrix(counts),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps <- tax_table(as.matrix(tax)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps <- sample_data(metadata) #define our metadata as sample information
ps <- phyloseq(count.ps, tax.ps, meta.ps) #put them all together into a phyloseq object


#transform seq counts into %
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))


#write out table from ps so I can more easily wrangle

new_counts <- as.data.frame(t(otu_table(ps)))


attach(metadata)
ordination.model <- metaMDS(new_counts, distance='bray', k=3)


#extract coordinates so I can plot with ggplot2
nmds_coords <- ordination.model$points

data1 <- cbind(metadata, nmds_coords)

#write.csv(data1, "CR_seed_func_NMDS_coord.csv", row.names = FALSE)


data <- cbind(metadata, new_counts)
data$Sample_Name <- row.names(data)
row.names(data) <- NULL
data <- data[,c(600292,1:3,4:600291)]


menhinick_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Richness_func=menhinick(x[,5:600291]>0))
})



shan_func <- ddply(data,~Sample_Name,function(x){
  data.frame(Shannon_func=diversity(x[,5:600291],index="shannon"))
})




invsimpson_func <- ddply(data,~Sample_Name,function(x){
  data.frame(InverseSimpson_func=diversity(x[,5:600291],index="invsimpson"))
})


#evenness
pielou_func <-ddply(data,~Sample_Name,function(x){
  data.frame(Pielou_func=exp(diversity(x[,5:600291], index="shannon"))/sum(x[,5:600291]>0))
})





#### merge #####

metadata <- data[,c(1:4)]
metadata <- merge(metadata, menhinick_org, by="Sample_Name")
metadata <- merge(metadata, shan_org, by="Sample_Name")
metadata <- merge(metadata, invsimpson_org, by="Sample_Name")
metadata <- merge(metadata, pielou_org, by="Sample_Name")
metadata <- merge(metadata, menhinick_func, by="Sample_Name")
metadata <- merge(metadata, shan_func, by="Sample_Name")
metadata <- merge(metadata, invsimpson_func, by="Sample_Name")
metadata <- merge(metadata, pielou_func, by="Sample_Name")

#write.csv(metadata, "CR_seed_func_diversity.csv", row.names = FALSE)

metadata$Treatment<-factor(metadata$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
palette <-c( "#56B4E9", "#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")

shan_plot <- ggplot(metadata, aes(x= Shannon_org, y = Shannon_func, color = Treatment)) +
  geom_point(size=3) +
  geom_smooth(method=lm, se = FALSE, col = 'black') +
  theme_classic(base_size=15) +
  labs(x= "Shannon Diversity, Bins", y = "Shannon Diversity, SEED Subsystems Database") +
  scale_color_manual(values = palette) +
  facet_grid(~Treatment)
shan_plot


isimpplot <- ggplot(metadata, aes(x= InverseSimpson_org, y = InverseSimpson_func, color = Treatment)) +
  geom_point(size=3) +
  geom_smooth(method=lm, se = FALSE, col = 'black') +
  theme_classic(base_size=15) +
  scale_color_manual(values = palette) +
  labs(x= "Inverse Simpson Diversity, Bins", y = "Inverse Simpson Diversity, SEED Subsystems Database") +
  facet_grid(~Treatment)
isimpplot

shan_plot_all <- ggplot(metadata, aes(x= Shannon_org, y = Shannon_func, color = Treatment)) +
  geom_point(size=3) +
  geom_smooth(method=lm, se = FALSE, col = 'black') +
  theme_classic(base_size=15) +
  labs(x= "Shannon Diversity, Bins", y = "Shannon Diversity, SEED Subsystems Database") +
  scale_color_manual(values = palette)
shan_plot_all


isimpplot_all <- ggplot(metadata, aes(x= InverseSimpson_org, y = InverseSimpson_func, color = Treatment)) +
  geom_point(size=3) +
  geom_smooth(method=lm, se = FALSE, col = 'black') +
  theme_classic(base_size=15) +
  scale_color_manual(values = palette) +
  labs(x= "Inverse Simpson Diversity, Bins", y = "Inverse Simpson Diversity, SEED Subsystems Database")
isimpplot_all



ggarrange(shan_plot, isimpplot, ncol = 1, common.legend=TRUE)
ggarrange(shan_plot_all, isimpplot_all, ncol = 1, common.legend=TRUE)

