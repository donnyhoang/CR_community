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

#cleanup strings, change as needed when doing mg vs mt
data_table$variable <- gsub("./", "", as.character(data_table$variable))
data_table$variable <- gsub("_mg.dbcan_annot_function.tsv", "", as.character(data_table$variable))
data_table$Gene <- gsub(",", "", as.character(data_table$Gene))

#string split Gene column
data_table_dbcan <- data_table
data_table_dbcan[c("bin","unused")] <- str_split_fixed(data_table_dbcan$Gene, "[.]", 2)
data_table_dbcan <- data_table_dbcan[,c(1:4)]
colnames(data_table_dbcan)[1] = "cazy"
data_table_dbcan$cazy <- gsub(" ", "", as.character(data_table_dbcan$cazy))
data_table_dbcan$bin <- gsub(" ", "", as.character(data_table_dbcan$bin))
data_table_total <- data_table_dbcan



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


length(unique(total_table$bin))
length(unique(total_table$order))



#sort into cazyme functions

#read in table
cazy_func <- read.csv("cazyme_list_referenece.csv", header = TRUE)
cazy_func <- cazy_func[,c(1:3)]

#make everything not in cazy_func into "other"
other_cazy <- as.data.frame(unique(total_table$CAZy))
colnames(other_cazy)[1] <- "CAZy"
cazy_list <- cazy_func$CAZy
other_cazy <- filter(other_cazy, !CAZy %in% cazy_list)
other_cazy$Activity <- "Other"
other_cazy$Substrate <- "Other"

colnames(other_cazy) <- colnames(cazy_func)

cazy_func <- rbind(cazy_func, other_cazy)


total_table <- merge(total_table, cazy_func, by = "CAZy")


#write out metadata from variable
total_table[c("Treatment", "Day", "Replicate")] <- str_split_fixed(total_table$variable, "_", 3)
total_table <- total_table[,-c(4:6)]
write.csv(total_table, "CR_dbcan_total_table_14m.csv", row.names = FALSE)



############### Get MAG data ################
hmm[c("bin", "unused")] <- str_split_fixed(hmm$cazy, "[.]", 2)
hmm <- hmm[,c(1,2,6)]
hmm <- hmm[!grepl("#", hmm$CAZy),]

#make everything not in cazy_func into "other"

other_cazy <- as.data.frame(unique(hmm$CAZy))
colnames(other_cazy)[1] <- "CAZy"

cazy_func <- read.csv("cazyme_list_referenece.csv", header = TRUE)
cazy_func <- cazy_func[,c(1:3)]

cazy_list <- cazy_func$CAZy
other_cazy <- filter(other_cazy, !CAZy %in% cazy_list)
other_cazy$Activity <- "Other"
other_cazy$Substrate <- "Other"

colnames(other_cazy) <- colnames(cazy_func)

cazy_func <- rbind(cazy_func, other_cazy)

#merge

merge1 <- merge(hmm, tax, by = "bin")
merge2 <-merge(merge1, cazy_func, by = "CAZy")

write.csv(merge2, "CR_dbcan_MAG_annotation_total.csv", row.names = FALSE)
