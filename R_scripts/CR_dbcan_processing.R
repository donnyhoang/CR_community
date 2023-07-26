##mostly copied from
#https://github.com/transcript/samsa2/tree/master/R_scripts
## May 10, 2023

library(ggplot2)
library(reshape2)
library(knitr)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(data.table)

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
data_table$variable <- gsub("_mt.dbcan_annot_function.tsv", "", as.character(data_table$variable))
data_table$Gene <- gsub(",", "", as.character(data_table$Gene))

#string split Gene column
data_table_dbcan <- data_table
data_table_dbcan[c("bin","unused")] <- str_split_fixed(data_table_dbcan$Gene, "[.]", 2)
data_table_dbcan <- data_table_dbcan[,c(1:4)]
colnames(data_table_dbcan)[1] = "cazy"
data_table_dbcan$cazy <- gsub(" ", "", as.character(data_table_dbcan$cazy))
data_table_dbcan$bin <- gsub(" ", "", as.character(data_table_dbcan$bin))

#add raw counts
counts <- read.table("total_reads_mt.txt", header=FALSE, sep ="\t")
counts$V1 <- gsub("_mt.R1.fastq.cleaned.forward","",as.character(counts$V1))

colnames(counts)[1] = "variable"
colnames(counts)[2] = "counts"

#merge with counts
data_table_total <- merge (data_table_dbcan, counts, by = "variable")
data_table_total <- data_table_total %>%
  mutate(percent=(value/counts)*100)

#read hmmscan *domtblout files
files <- list.files(pattern=".hmmscanout.tsv")
temp <- lapply(files, fread, sep="\t")
hmm <- rbindlist(temp,)
hmm <- as.data.frame(hmm)
#oops this is ugly
hmm <- hmm[-(1:3), , drop=FALSE]
colnames(hmm)[1] = "CAZy"
colnames(hmm)[2] = "cazy"
colnames(hmm)[3] = "start"
colnames(hmm)[4] = "stop"
hmm$CAZy <- gsub(".hmm","", as.character(hmm$CAZy))
hmm$range <- apply(hmm[,c(3,4)], 1, paste , collapse ="-")
hmm$cazy <- apply(hmm[,c(2,5)], 1, paste , collapse =":")



#first tax table
tax <- read.csv("CR_bin_tax_cluster.csv")



#merge tables
total_table <- merge(hmm, data_table_total, by ="cazy")

total_table <- merge(total_table, tax, by ="bin")


length(unique(total_table$Markerlineage))
length(unique(total_table$UID))
length(unique(total_table$phylum))
length(unique(total_table$class))
length(unique(total_table$order))
length(unique(total_table$CAZy))
length(unique(total_table$cazy))


#sort into cazyme functions

cellulases <- c("GH5_","GH5_1", "GH5_2", "GH5_4", "GH5_5", "GH6", "GH5_7", "GH5_8", "GH5_11", "GH5_12", "GH5_13", "GH5_15", "GH5_16", "GH5_18", "GH5_19", "GH5_2","GH5_22", "GH5_23", "GH5_24", "GH5_25", "GH5_26", "GH5_27", "GH5_28", "GH5_29", "GH5_30", "GH5_31", "GH5_36", "GH5_37", "GH5_38", "GH5_39", "GH5_4", "GH5_40", "GH5_41", "GH5_43", "GH5_44", "GH5_45", "GH5_46", "GH5_47", "GH5_48", "GH5_49", "GH5_5", "GH5_50", "GH5_51", "GH5_53", "GH5_7", "GH5_8", "GH5_9",
                "GH6","GH7","GH8","GH9","GH12","GH16","GH26","GH44","GH45","GH48", "GH51","GH61","GH74","GH131", "AA9", "AA10")
xylanases <- c("GH10","GH11",
               "GH30","GH30_","GH30_2","GH30_3","GH30_5","GH30_7",
               "GH67","GH115","GH129")
ligninases <- c("AA1","AA1_","AA1_1","AA1_2","AA1_3","AA2")
other_oxidases <- c("AA4","AA5","AA5_1","AA5_2","AA6","AA7","AA8","AA12")
GMC_oxidoreductases <- c("AA3","AA3_","AA3_1","AA3_2","AA3_3")
amylases <- c("GH13_1","GH13_2","GH13_3","GH13_4","GH13_5","GH13_6","GH13_7","GH13_8","GH13_9","GH13_10","GH13_11","GH13_12","GH13_13","GH13_14","GH13_15","GH13_16","GH13_17","GH13_18","GH13_19","GH13_20","GH13_21","GH13_22","GH13_23","GH13_24","GH13_25","GH13_26","GH13_27","GH13_28","GH13_29","GH13_30","GH13_31","GH13_32","GH13_33","GH13_34","GH13_35","GH13_36","GH13_37","GH13_38","GH13_39","GH13_40",
              "GH14")
betaglucanases <- c("GH1","GH3")

cazy_all <- c(cellulases, xylanases, ligninases, other_oxidases, GMC_oxidoreductases, amylases, betaglucanases)




cazy_cellu <- subset(total_table, CAZy %in% cellulases)
cazy_cellu$func <-"cellulase"
cazy_xylan <- subset(total_table, CAZy %in% xylanases)
cazy_xylan$func <- "xylanase"
cazy_lignin <- subset(total_table, CAZy %in% ligninases)
cazy_lignin$func <- "ligninase"
cazy_other_oxidases <- subset(total_table, CAZy %in% other_oxidases)
cazy_other_oxidases$func <- "other oxidases"
cazy_GMC_oxidoreductases <- subset(total_table, CAZy %in% GMC_oxidoreductases)
cazy_GMC_oxidoreductases$func <- "GMC oxidoreductases"
cazy_amylase <- subset(total_table, CAZy %in% amylases)
cazy_amylase$func <- "amylase"
cazy_beta <- subset(total_table, CAZy %in% betaglucanases)
cazy_beta$func <- "betaglucanase"
cazy_other <- subset(total_table, !(CAZy %in% cazy_all))
cazy_other$func <- "other"
cazy_total <- rbind(cazy_cellu, cazy_xylan,cazy_lignin,cazy_other_oxidases,cazy_GMC_oxidoreductases,cazy_amylase,cazy_beta, cazy_other)
cazy_total[c("Treatment","Day","Replicate")] <- str_split_fixed(cazy_total$variable, "_", 3)


#write.csv(cazy_total, "CR_dbcan_total_table.csv", row.names = FALSE)





data_trim <- cazy_total %>%
  group_by(Treatment,Day, Replicate) %>%
  summarise(sum=sum(percent),
            count=n())
data_trim$Day <- as.character(data_trim$Day)



palette <-c( "#56B4E9", "#999999","#E69F00","#0072B2", "#F0E442","#CC79A7")
pbox <- ggplot(data_trim, aes(x=Day, y=sum, fill=Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values=palette) +
  #scale_color_manual(values=pal1) +
  theme_classic(base_size=20) +
  #ylim(0,19) +
  ylab("Percent reads mapped to dbCAN") + 
  xlab("Time (Day sampled)") +
  #theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(aspect.ratio=1)+
  facet_wrap(  ~ Treatment , scales = "free_x")
pbox


