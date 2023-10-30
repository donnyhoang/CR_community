#April 4, 2023

#load packages
library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

#load files into dataframe
files <- list.files(pattern="checkm.autometa.tsv")
temp <- lapply(files, fread, sep="\t")
data <- rbindlist(temp)
data <- as.data.frame(data)

#remove space from colnames
colnames(data) = gsub(" ","", colnames(data))


#read in rep bins
#list<-fread("rep_bins_autometa.list", sep="\t", header = FALSE)
#list2 <- list$V1
#rep_table <- filter(data, BinId %in% list2)
#write.csv(rep_table, "rep_bins_qa.csv", row.names = FALSE)
#length(unique(rep_table$Markerlineage))

#write table of highqual bins
#data <- subset(data, Completeness >= 80)
#data <- subset(data, Contamination < 40)
#highqual <- as.data.frame(data$BinId)
#write.csv(highqual, "bins_to_keep.csv", row.names = FALSE)


#trim down columns and rearrange columns so
#it looks more like standard checkm qa table, write out

#newdf <-data[,c(1:8,24:29)]
#newdf <- newdf[,c(1:5,9:14,6:8)]
#write.table(newdf, file="cr_autometa_bins_hq_checkm_table.tsv", row.names=FALSE, sep="\t", quote=FALSE)


#write table of unique classificaitons
#unique <- as.data.frame(unique(data$Markerlineage))
#unique[c('1','2')] <- str_split_fixed(unique$`unique(data$Markerlineage)`, "__", 2)
#unique[c('lineage','3')] <- str_split_fixed(unique$`2`, " ", 2)
#unique <- as.data.frame(unique$lineage)
#unique <- unique(unique)
#write.csv(unique, "unique_checkm_uid.csv", row.names = FALSE)

#make a column where value is same for everything
data$all <- "all"

#string split name for metadata
data[c('Treatment', 'Day','unusedstring')] <- str_split_fixed(data$BinId, "_", 3)

#evaluate completeness+contamination for high quality bins
#completeness, adjust cutoffs as you want
data <- data %>%
  mutate(comp = case_when(
    Completeness >= 80 ~ "high",
    Completeness <80 ~ "low"
  ))
#contamination, adjust cutoffs as you want
data <- data %>%
  mutate(cont = case_when(
    Contamination > 40 ~ "high",
    Contamination <=40 ~ "low"
  ))

#assign high quality based on completeness and contamination
data <- data %>%
  mutate(quality = case_when(
    comp=="high" & cont=="low" ~ "high",
    comp=="high" & cont=="high" ~ "low",
    comp=="low" & cont=="high" ~ "low",
    comp=="low" & cont=="low" ~ "low"
  ))

#get genomesize in Mb
data$GenomesizeMb <- (data$`Genomesize(bp)`)/1000000

#get number of bins and highquality bins
nrow(data)
highqual <- subset(data, data$quality == "high")
nrow(highqual)


##play around with high qual bins

highqual[c('Lineage', 'UID')] <- str_split_fixed(highqual$Markerlineage, " ", 2)
length(unique(highqual$UID))
length(unique(highqual$Markerlineage))
length(str_count(rep$is_rep, "yes"))

#write.csv(highqual, "cr_autometa_bins_to_keep.csv", row.names=FALSE)


#add in galah representative bin info
rep <- read.table("rep_bins_autometa.tsv", header=FALSE, sep ="\t")
colnames(rep)[1] = "representative"
colnames(rep)[2] = "BinId"
rep$representative <- gsub("hq_bins/","", as.character(rep$representative))
rep$representative <- gsub(".fasta","", as.character(rep$representative))
rep$BinId <- gsub("hq_bins/","", as.character(rep$BinId))
rep$BinId <- gsub(".fasta","", as.character(rep$BinId))
rep$is_rep <- ifelse(rep$representative==rep$BinId, "yes", "no")
sum(str_count(rep$is_rep, "yes")) #check that I have correct number of rep bins


data <- merge(data, rep, by = "BinId")

#get taxa of high quality bins
#subset each treatment
ph6 <- subset(data, data$quality == "high" & data$Treatment=="pH6")
mes <- subset(data, data$quality == "high" & data$Treatment=="MES")
g05 <- subset(data, data$quality == "high" & data$Treatment=="G05")
g10 <- subset(data, data$quality == "high" & data$Treatment=="G10")
ph8 <- subset(data, data$quality == "high" & data$Treatment=="pH8")
ph10 <- subset(data, data$quality == "high" & data$Treatment=="pH10")

#make into individual dataframe so I can get specific % of each tax classification
taxph6 <- as.data.frame(table(ph6$Markerlineage))
taxph6$percent <- (taxph6$Freq/nrow(ph6))*100
sum(taxph6$percent)
taxph6$Treatment <- "pH6"

taxmes <- as.data.frame(table(mes$Markerlineage))
taxmes$percent <- (taxmes$Freq/nrow(mes))*100
taxmes$Treatment <- "MES"

taxg05 <- as.data.frame(table(g05$Markerlineage))
taxg05$percent <- (taxg05$Freq/nrow(g05))*100
taxg05$Treatment <- "G05"

taxg10 <- as.data.frame(table(g10$Markerlineage))
taxg10$percent <- (taxg10$Freq/nrow(g10))*100
taxg10$Treatment <- "G10"

taxph8 <- as.data.frame(table(ph8$Markerlineage))
taxph8$percent <- (taxph8$Freq/nrow(ph8))*100
taxph8$Treatment <- "pH8"

taxph10 <- as.data.frame(table(ph10$Markerlineage))
taxph10$percent <- (taxph10$Freq/nrow(ph10))*100
taxph10$Treatment <- "pH10"

#merge for total tax table
tax_total <- rbind(taxph6, taxmes,taxg05,taxg10,taxph8,taxph10)
#tax_total <- tax_total[ , -2]

#get most abundant groups
tax_total$percent <- as.numeric(tax_total$percent)
tax_top <- tax_total %>%
  group_by(Var1) %>%
  summarise(sum=sum(percent))
tax_top <- tax_top %>%
  arrange(desc(sum))

#reorder based on tax_top
tax_total$Var1 <- factor(tax_total$Var1, levels = tax_top$Var1)

# get top 10. change number range for different "top x" list
top_10 <- tax_top[1:15 , ]
top_10 <- top_10$Var1


#filter out high abundant microbes
lowtax <- filter(tax_total, !Var1 %in% top_10)
lowtax <- lowtax %>%
  group_by(Treatment) %>%
  summarise(sum=sum(Freq),
            sump=sum(percent))
lowtax$Var1 <- "Other"
colnames(lowtax)[2] <- "Freq"
colnames(lowtax)[3] <- "percent"

#fitler out "low", readd add lowtax back to table, as "other"
tax_total <- filter(tax_total, Var1 %in% top_10)

tax_total <- rbind(tax_total, lowtax)

tax_total$Treatment <- factor(tax_total$Treatment, levels=c("pH6", "MES", "G05", "G10", "pH8","pH10"))

######graphs

palette <- c("black","red")

p_rep <- ggplot(data, aes(x=Contamination , y=Completeness, color=is_rep, size = is_rep, shape= is_rep)) +
  geom_jitter(alpha=1/2) +
  theme_classic(base_size=15)  +
  scale_color_manual(values=palette)  +
  scale_size_manual(values=c(1,3)) +
  scale_shape_manual(values=c(16,17)) +
  labs(y="Percent Completeness", x="Percent Contamination") +
  theme(aspect.ratio = 1) +
  guides(color=FALSE, size=FALSE, shape=guide_legend("Representative Bin"))
p_rep

completeness <- ggplot(data, aes(x=all , y=Completeness, color=is_rep, size = is_rep)) +
  geom_jitter(alpha=0.7) +
  theme_classic(base_size=15)  +
  scale_color_manual(values=palette)  +
  scale_size_manual(values=c(1,2)) +
  labs(y="Percent Completeness") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  guides(size=FALSE, color = guide_legend("Representative Bin"))

completeness


contamination <- ggplot(data, aes(x=all , y=Contamination, color=is_rep, size = is_rep)) +
  geom_jitter(alpha=0.7) +
  theme_classic(base_size=15) +
  scale_color_manual(values=palette)  +
  scale_size_manual(values=c(1,2)) +
  labs(y="Percent Contamination")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  guides(size=FALSE, color = guide_legend("Representative Bin"))
contamination



bin_qual <- ggarrange(completeness, contamination, common.legend = TRUE, legend = "right")
bin_qual

ggsave("bin_qual.png", heigh = 6, width = 8 , units = "in")


############graph of rep bin genome size
tax1 <- read.csv("rep_bins_autometa_qa.csv")
colnames(tax1)[1] <- "bin"
tax1[c("classification","UID")] <- str_split_fixed(tax1$Markerlineage , ",", 2)
tax1$UID <- gsub("[(]","",as.character(tax1$UID))
tax1$UID <- gsub("[)]","",as.character(tax1$UID))

#second tax table
tax2 <- read.table("out.BAT.bin2classification.official_names.txt", sep="\t", header=TRUE)
tax2$bin <- gsub(".fasta","",as.character(tax2$bin))

tax <- merge(tax1, tax2, by = "bin")
tax <- tax[,c(1,6:9,13:21,36:42)]

tax$GenomesizeMb <- (tax$`Genomesize.bp.`)/1000000

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

palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462",
            "#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
tax$Order<-factor(tax$Order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                    "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                    "Sphingobacteriales","Xanthomonadales","no support"))

write.csv(tax, "MAG_info.csv", row.names= FALSE)

p <-ggplot() + 
  geom_point(data=tax, aes(y=Genomesizemb,x=Order), size = 3) +
  geom_jitter() + 
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  labs(fill='Order', y='Genome Size (Mb)') +
  scale_fill_manual(values=palorg)
  #scale_color_manual(values=palorg) 
  #guides(color = "none", fill = guide_legend(ncol=10)) 

p

