##mostly copied from
#https://github.com/transcript/samsa2/tree/master/R_scripts
## June 2, 2023


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
data_table$variable <- gsub("_mg.all_annot_function.tsv", "", as.character(data_table$variable))

#string split Gene column
data_table_org <- data_table
data_table_org[c("bin","fasta")] <- str_split_fixed(data_table_org$Gene, "[.]", 2)
data_table_org$bin <- gsub(" ", "", as.character(data_table_org$bin))
data_table_org <- data_table_org[,c(2:4)]



#add raw counts
counts <- read.csv("CR_read_counts_merged_mg.csv", header=FALSE)
counts$V2 <- gsub("_mg","",as.character(counts$V2))

colnames(counts)[1] = "counts"
colnames(counts)[2] = "variable"

#merge with counts
data_table_total <- merge (data_table_org, counts, by = "variable")
total_table <- data_table_total %>%
  mutate(percent=(value/counts)*100)


##
tax_matrix <- dcast(total_table, variable ~ bin, value.var = "value", fun.aggregate = sum, fill = 0)
tax_matrix_rel <- dcast(total_table, variable ~ bin, value.var = "percent", fun.aggregate = sum, fill = 0)
#write.csv(tax_matrix, "CR_org_bins_mg_matrix.csv", row.names = FALSE)
#write.csv(tax_matrix_rel, "CR_org_bins_mg_matrix_relative.csv", row.names = FALSE)

#read bin metadata
tax <- read.csv("CR_bin_tax_cluster.csv", header=TRUE)
tax$cluster_total <- as.factor(tax$clustercluster_total)
tax <- tax %>%
  mutate(across('clustercluster_total', str_replace, '0', 'unclustered'))

#merge tables
total_table <- merge(total_table, tax, by ="bin")

total_table[c("Treatment","Day","Replicate")] <- str_split_fixed(total_table$variable, "_", 3)
total_table$Treatment<-factor(total_table$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))

length(unique(total_table$bin))



summary <- total_table %>%
  group_by(bin, percent, Treatment, Day, Replicate) %>%
  summarise(percent=sum(percent),
            count=n())

#write.csv(total_table, "CR_org_total_table_cleaned_mg.csv", row.names = FALSE)
#write.csv(summary, "CR_org_total_table_summary_mg.csv", row.names=FALSE)

#######visualize
