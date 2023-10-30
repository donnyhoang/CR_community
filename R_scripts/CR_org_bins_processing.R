##mostly copied from
#https://github.com/transcript/samsa2/tree/master/R_scripts
## May 18, 2023

library(ggplot2)
library(scales)
library(reshape2)
library(knitr)
library(stringr)
library(RColorBrewer)
library(dplyr)

control_files <- list.files(pattern = "*_function.tsv", full.names = T, recursive = FALSE)
control_names = ""
for (name in control_files) {
  control_names <- c(control_names, unlist(strsplit(name, split='.', fixed=TRUE))[2])}
control_names <- control_names[-1]
control_names_trimmed = ""
for (name in control_names) {
  control_names_trimmed <- c(control_names_trimmed, unlist(strsplit(name, split='/', fixed=TRUE))[2])}
control_names_trimmed <- control_names_trimmed[-1]


# READ IN FILES
# loading the control table
y <- 0
for (x in control_files) {
  y <- y + 1
  if (y == 1) {
    control_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(control_table) = c("DELETE", x, "V3")
    control_table <- control_table[,c(2,3)]      }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    control_table <- merge(control_table, temp_table, by = "V3", all = T)  }
}
control_table[is.na(control_table)] <- 0
rownames(control_table) = control_table$V3
control_table_trimmed <- control_table[,-1]

#table loaded, melt
data_table <- melt(cbind(control_table_trimmed,
                         Gene = rownames(control_table_trimmed)))

#cleanup strings
data_table$variable <- gsub("./", "", as.character(data_table$variable))
data_table$variable <- gsub("_mt.all_annot_function.tsv", "", as.character(data_table$variable))


#string split Gene column
data_table_org <- data_table
data_table_org[c("bin","fasta")] <- str_split_fixed(data_table_org$Gene, "[.]", 2)
data_table_org$bin <- gsub(" ", "", as.character(data_table_org$bin))
data_table_org <- data_table_org[,c(2:4)]

data_table_org <- data_table_org %>%
  group_by(variable, bin) %>%
  summarise(value = sum(value))

#add raw counts mt
counts <- read.table("total_reads_mt.txt", header=FALSE, sep ="\t")
counts$V1 <- gsub("_mt.R1.fastq.cleaned.forward","",as.character(counts$V1))

colnames(counts)[1] = "variable"
colnames(counts)[2] = "counts"

#add raw counts mg

counts <- read.csv("CR_read_counts_merged_mg.csv", header = FALSE)
counts$V2 <- gsub("_mg","",as.character(counts$V2))
colnames(counts)[2] = "variable"
colnames(counts)[1] = "counts"

#merge with counts
data_table_total <- merge (data_table_org, counts, by = "variable")
total_table <- data_table_total %>%
  mutate(percent=(value/counts)*100)



#read bin tables (checkm, CATBAT) that have taxonomy info
tax <- read.csv("CR_bin_tax_cluster.csv")

#merge tables
total_table <- merge(total_table, tax, by ="bin")
length(unique(total_table$bin))
total_table[c("Treatment","Day","Replicate")] <- str_split_fixed(total_table$variable, "_", 3)


#write.csv(total_table, "CR_org_total_table.csv", row.names = FALSE)


#total_table_trim <- total_table[,c(1:5,12,16:18)]

#write.csv(total_table, "CR_org_total_table_cleaned.csv", row.names = FALSE)

#total_table <- read.csv("CR_org_total_table.csv", header=TRUE)
#tax_matrix <- dcast(data=total_table, formula = bin~variable, fun.aggregate = sum, value.var="value")
#write.csv(tax_matrix, "CR_taxmatrix.csv", row.names = FALSE)


total_table$superkingdom <- gsub("[0-9]+", "", total_table$superkingdom)
total_table$superkingdom <- gsub(": .", "", total_table$superkingdom)
total_table$phylum <- gsub("[0-9]+", "", total_table$phylum)
total_table$phylum <- gsub(": .", "", total_table$phylum)
total_table$class <- gsub("[0-9]+", "", total_table$class)
total_table$class <- gsub(": .", "", total_table$class)
total_table$order <- gsub("[0-9]+", "", total_table$order)
total_table$order <- gsub(": .", "", total_table$order)
total_table$family <- gsub("[0-9]+", "", total_table$family)
total_table$family <- gsub(": .", "", total_table$family)
total_table$genus <- gsub("[0-9]+", "", total_table$genus)
total_table$genus <- gsub(": .", "", total_table$genus)
total_table$species <- gsub("[0-9]+", "", total_table$species)
total_table$species <- gsub(": .", "", total_table$species)



summary <- total_table %>%
  group_by(bin, percent, Treatment, Day, Replicate) %>%
  summarise(percent=sum(percent),
            count=n())

write.csv(summary, "CR_org_total_table_summary", row.names=FALSE)


total_table$order<-factor(total_table$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                    "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                    "Sphingobacteriales","Xanthomonadales","no support"))
total_table$Treatment<-factor(total_table$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
total_table$Day <- as.character(total_table$Day)

palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462",
            "#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
palorg <- c("#FCCDE5")


pbox <- ggplot(total_table, aes(x=Day, y=sum, fill=genus)) +
  geom_boxplot() +
  scale_fill_manual(values=palorg) +
  #scale_color_manual(values=pal1) +
  theme_classic(base_size=20) +
  #ylim(0,19) +
  ylab("Percent of reads mapped to Timonella bin") + 
  xlab("Time (Day sampled)") +
  theme(axis.text.x=element_blank()) +
  theme(aspect.ratio=1)+
  facet_wrap(  ~ Treatment , scales = "free_x")
pbox

detach("package:stringr")
#data <- read.csv("CR_org_total_table.csv", header=TRUE)

#data <- subset(total_table_trim, Treatment == "G05")
data$Treatment<-factor(data$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
data$order <- gsub("[0-9]+", "", data$order)
data$order <- gsub(":", "", data$order)
data$order <- gsub("[.]", "", data$order)
data$order <- gsub("ales ", "ales", data$order)

length(unique(data$order))



#pal for organisms
palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462",
            "#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")


p <-ggplot(data, aes(y=percent,x=variable,fill=order, color=order)) + 
  geom_bar(stat="identity") + 
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "bottom") +
  labs(fill='Order', y='Percent relative expression') +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  guides(color = "none", fill = guide_legend(ncol=10)) +
  facet_grid( ~ Treatment, scales = "free")

p


data_trim <- data %>%
  group_by(Treatment,Day, Replicate) %>%
  summarise(sum=sum(percent))



palette <-c( "#56B4E9", "#999999","#E69F00","#0072B2", "#F0E442","#CC79A7")
data_trim$Day <- as.character(data_trim$Day)

pbox <- ggplot(data_trim, aes(x=Day, y=sum, fill=Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values=palette) +
  #scale_color_manual(values=pal1) +
  theme_classic(base_size=18) +
  #ylim(0,19) +
  ylab("Percent reads mapped") + 
  xlab("Time (Day sampled)") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(aspect.ratio=1)+
  facet_wrap(  ~ Treatment , scales = "free_x")
pbox


