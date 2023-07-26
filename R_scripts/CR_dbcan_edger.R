library(edgeR)
library(stringr)
library(data.table)
library(dplyr)
#https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
data <- read.csv("CR_dbcan_complete_table.csv", row.names=1)






completeCondition <- colnames(data)
completeCondition <- as.data.frame(completeCondition)
completeCondition[c('condition', 'time', 'replicate')] <- str_split_fixed(completeCondition$completeCondition, "_", 3)
completeCondition <- completeCondition[ -c(1,4) ]




completeCondition <- completeCondition %>%
  mutate(time2 = case_when(
    time == "1" ~ "Day1to3",
    time == "2" ~ "Day1to3",
    time == "3" ~ "Day1to3",
    time == "4" ~ "Day4to7",
    time == "5" ~ "Day4to7",
    time == "6" ~ "Day4to7",
    time == "7" ~ "Day4to7"
  ))




completeCondition$condition <- as.factor(completeCondition$condition)
completeCondition$time2 <- as.factor(completeCondition$time2)

treatment <- completeCondition$condition
time <- completeCondition$time2

groups <- completeCondition$condition
groups <- factor(paste(treatment, time, sep="_"))

d <- DGEList(counts=data, group=factor(groups))
d

dim(d)
d.full <- d
#head(d$counts)
#head(cpm(d))

#find mean library size -> ~35million
summary(d$samples$lib.size)

# 0.3 * 35 is about 10 reads per sample
keep <- rowSums(cpm(d) > 0.3) >=5

d <-d[keep,]
dim(d)

d$samples$lib.size <- colSums(d$counts)
d$samples

d <- calcNormFactors(d)
d



d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)

#look at all groups for comparisons
grouptable <- as.data.frame(unique(d1$samples$group))



########################## Compare all treatments and timepoints to pH 6 ###################################################

#first group listed is the baseline
etg05 <- exactTest(d1, pair=c(35,7))
etg10 <- exactTest(d1, pair=c(35,15))
etmes <- exactTest(d1, pair=c(35,21))
etph10 <- exactTest(d1, pair=c(35,28)) 
etph8 <- exactTest(d1, pair=c(35,42)) 

topTags(etg05, n=5)
g05 <- as.data.frame(etg05)
g05$cazy <- row.names(g05)
g05$cazy <- gsub(" ", "", g05$cazy)
g05$Comparison <- "G05"
#write.csv(g05, "CR_g05_edger.csv", row.names = FALSE)

topTags(etg10, n=5)
g10 <- as.data.frame(etg10)
g10$cazy <- row.names(g10)
g10$cazy <- gsub(" ", "", g10$cazy)
g10$Comparison <- "G10"
#write.csv(g10, "CR_g10_edger.csv", row.names = FALSE)

topTags(etmes, n=5)
mes <- as.data.frame(etmes)
mes$cazy <- row.names(mes)
mes$cazy <- gsub(" ", "", mes$cazy)
mes$Comparison <- "MES"
#write.csv(mes, "CR_mes_edger.csv", row.names = FALSE)


topTags(etph8, n=5)
ph8 <- as.data.frame(etph8)
ph8$cazy <- row.names(ph8)
ph8$cazy <- gsub(" ", "", ph8$cazy)
ph8$Comparison <- "pH8"
#write.csv(ph8, "CR_ph8_edger.csv", row.names = FALSE)


topTags(etph10, n=5)
ph10 <- as.data.frame(etph10)
ph10$cazy <- row.names(ph10)
ph10$cazy <- gsub(" ", "", ph10$cazy)
ph10$Comparison <- "pH10"
#write.csv(ph10, "CR_ph10_edger.csv", row.names = FALSE)

total_table_7 <- rbind(g05, g10, mes, ph8, ph10)
total_table_7$Comparison_Day <- "7"
rownames(total_table_7) <- NULL


#total_table <- rbind(total_table_1,total_table_2,total_table_3,total_table_4,total_table_5,total_table_6,total_table_7)
############################################################################

############### Compare days1-3 to days 4-7 within treatment ##################
etph6 <- exactTest(d1, pair=c(9,10))
etg05 <- exactTest(d1, pair=c(1,2))
etg10 <- exactTest(d1, pair=c(3,4))
etmes <- exactTest(d1, pair=c(5,6))
etph10 <- exactTest(d1, pair=c(7,8)) 
etph8 <- exactTest(d1, pair=c(11,12))

topTags(etph6, n=5)
ph6 <- as.data.frame(etph6)
ph6$cazy <- row.names(ph6)
ph6$cazy <- gsub(" ", "", ph6$cazy)
ph6$Comparison <- "pH6"


topTags(etg05, n=5)
g05 <- as.data.frame(etg05)
g05$cazy <- row.names(g05)
g05$cazy <- gsub(" ", "", g05$cazy)
g05$Comparison <- "G05"

topTags(etg10, n=5)
g10 <- as.data.frame(etg10)
g10$cazy <- row.names(g10)
g10$cazy <- gsub(" ", "", g10$cazy)
g10$Comparison <- "G10"

topTags(etmes, n=5)
mes <- as.data.frame(etmes)
mes$cazy <- row.names(mes)
mes$cazy <- gsub(" ", "", mes$cazy)
mes$Comparison <- "MES"


topTags(etph8, n=5)
ph8 <- as.data.frame(etph8)
ph8$cazy <- row.names(ph8)
ph8$cazy <- gsub(" ", "", ph8$cazy)
ph8$Comparison <- "pH8"


topTags(etph10, n=5)
ph10 <- as.data.frame(etph10)
ph10$cazy <- row.names(ph10)
ph10$cazy <- gsub(" ", "", ph10$cazy)
ph10$Comparison <- "pH10"

total_table <- rbind(ph6, g05, g10, mes, ph8, ph10)

rownames(total_table) <- NULL


############################################




#read in cazy info
files <- list.files(pattern=".hmmscanout.tsv")
temp <- lapply(files, fread, sep="\t")
hmm <- rbindlist(temp,)
hmm <- as.data.frame(hmm)
#oops this is ugly
hmm <- hmm[-(1:3), , drop=FALSE]
#hmm$name <- gsub(".faa.domtblout","",as.character(hmm$name))
colnames(hmm)[1] = "CAZy"
colnames(hmm)[2] = "cazy"
colnames(hmm)[3] = "start"
colnames(hmm)[4] = "stop"
hmm$CAZy <- gsub(".hmm","", as.character(hmm$CAZy))
hmm$range <- apply(hmm[,c(3,4)], 1, paste , collapse ="-")
hmm$cazy <- apply(hmm[,c(2,5)], 1, paste , collapse =":")
hmm$CAZy <- gsub(" ","", as.character(hmm$CAZy))
hmm$cazy <- gsub(" ","", as.character(hmm$cazy))
hmm[c("bin","discard")] <- str_split_fixed(hmm$cazy,"[.]",  2)
hmm$bin <- gsub(" ","", as.character(hmm$bin))
hmm <- hmm[,c(1,2,6)]

#read in bin info
tax <- read.csv("CR_bin_tax_cluster.csv")


#combine tables
total_table <- merge (total_table, hmm, by = "cazy")
total_table <- merge(total_table, tax, by = "bin")


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
#cazy_total$Comparison <- factor(cazy_total$Comparison, levels = c("pH6","MES","G05","G10","pH8","pH10"))

#write.csv(cazy_total, "CR_dbcan_edger_earlyvslate.csv", row.names = FALSE)
#write.csv(cazy_total, "CR_dbcan_edger.csv", row.names = FALSE)
