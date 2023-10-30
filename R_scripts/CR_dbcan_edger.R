library(edgeR)
library(stringr)
library(data.table)
library(dplyr)
#https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html


data <- read.csv("CR_dbcan_matrixforedger.csv", row.names=1)






completeCondition <- colnames(data)
completeCondition <- as.data.frame(completeCondition)
completeCondition[c('condition', 'time', 'replicate')] <- str_split_fixed(completeCondition$completeCondition, "_", 3)
completeCondition <- completeCondition[ -c(1,4) ]


###################### Ignore, kind of arbitrary grouping #####

#completeCondition <- completeCondition %>%
#  mutate(time2 = case_when(
#    time == "1" ~ "Day1to3",
#    time == "2" ~ "Day1to3",
#    time == "3" ~ "Day1to3",
#    time == "4" ~ "Day4to7",
#    time == "5" ~ "Day4to7",
#    time == "6" ~ "Day4to7",
#    time == "7" ~ "Day4to7"
#  ))

completeCondition$condition <- as.factor(completeCondition$condition)
completeCondition$time2 <- as.factor(completeCondition$time2)

treatment <- completeCondition$condition
time <- completeCondition$time2


##########################################
treatment <- completeCondition$condition
time <- completeCondition$time


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

####################### COMPARE TIMEPOINTS OF EACH TREATMENT TO THEIR RESPECTIVE DAY 1 ###########################

#first group listed is the baseline
et_2_G05 <- exactTest(d1, pair=c(1,2))
et_3_G05 <- exactTest(d1, pair=c(2,3))
et_4_G05 <- exactTest(d1, pair=c(3,4))
et_5_G05 <- exactTest(d1, pair=c(4,5))
et_6_G05 <- exactTest(d1, pair=c(5,6))
et_7_G05 <- exactTest(d1, pair=c(6,7))

#check top changes, doesn't mean much without being keyed though
topTags(et_2_G05, n=5)
topTags(et_3_G05, n=5)
topTags(et_4_G05, n=5)
topTags(et_5_G05, n=5)
topTags(et_6_G05, n=5)
topTags(et_7_G05, n=5)

#write out to df
G05_2 <- as.data.frame(et_2_G05)
G05_3 <- as.data.frame(et_3_G05)
G05_4 <- as.data.frame(et_4_G05)
G05_5 <- as.data.frame(et_5_G05)
G05_6 <- as.data.frame(et_6_G05)
G05_7 <- as.data.frame(et_7_G05)

#write cazy into column
G05_2$cazy <- rownames(G05_2)
G05_3$cazy <- rownames(G05_3)
G05_4$cazy <- rownames(G05_4)
G05_5$cazy <- rownames(G05_5)
G05_6$cazy <- rownames(G05_6)
G05_7$cazy <- rownames(G05_7)

#at time comparison
G05_2$Time <- "1_to_2"
G05_3$Time <- "2_to_3"
G05_4$Time <- "3_to_4"
G05_5$Time <- "4_to_5"
G05_6$Time <- "5_to_6"
G05_7$Time <- "6_to_7"

#combine, tidy up
G05_table <- rbind(G05_2, G05_3,G05_4,G05_5,G05_6,G05_7)
G05_table$cazy <- gsub(" ", "", G05_table$cazy)
G05_table$Treatment <- "G05"
rownames(G05_table) <- NULL

### G10
#first group listed is the baseline
et_2_G10 <- exactTest(d1, pair=c(8,9))
et_3_G10 <- exactTest(d1, pair=c(9,10))
et_4_G10 <- exactTest(d1, pair=c(10,11))
et_5_G10 <- exactTest(d1, pair=c(11,12))
et_6_G10 <- exactTest(d1, pair=c(12,13))
et_7_G10 <- exactTest(d1, pair=c(13,14))

topTags(et_2_G10, n=5)
topTags(et_3_G10, n=5)
topTags(et_4_G10, n=5)
topTags(et_5_G10, n=5)
topTags(et_6_G10, n=5)
topTags(et_7_G10, n=5)

G10_2 <- as.data.frame(et_2_G10)
G10_3 <- as.data.frame(et_3_G10)
G10_4 <- as.data.frame(et_4_G10)
G10_5 <- as.data.frame(et_5_G10)
G10_6 <- as.data.frame(et_6_G10)
G10_7 <- as.data.frame(et_7_G10)

G10_2$cazy <- rownames(G10_2)
G10_3$cazy <- rownames(G10_3)
G10_4$cazy <- rownames(G10_4)
G10_5$cazy <- rownames(G10_5)
G10_6$cazy <- rownames(G10_6)
G10_7$cazy <- rownames(G10_7)

G10_2$Time <- "1_to_2"
G10_3$Time <- "2_to_3"
G10_4$Time <- "3_to_4"
G10_5$Time <- "4_to_5"
G10_6$Time <- "5_to_6"
G10_7$Time <- "6_to_7"

G10_table <- rbind(G10_2, G10_3,G10_4,G10_5,G10_6,G10_7)
G10_table$cazy <- gsub(" ", "", G10_table$cazy)
G10_table$Treatment <- "G10"
rownames(G10_table) <- NULL

### MES
#first group listed is the baseline
et_2_MES <- exactTest(d1, pair=c(15,16))
et_3_MES <- exactTest(d1, pair=c(16,17))
et_4_MES <- exactTest(d1, pair=c(17,18))
et_5_MES <- exactTest(d1, pair=c(18,19))
et_6_MES <- exactTest(d1, pair=c(19,20))
et_7_MES <- exactTest(d1, pair=c(20,21))

topTags(et_2_MES, n=5)
topTags(et_3_MES, n=5)
topTags(et_4_MES, n=5)
topTags(et_5_MES, n=5)
topTags(et_6_MES, n=5)
topTags(et_7_MES, n=5)

MES_2 <- as.data.frame(et_2_MES)
MES_3 <- as.data.frame(et_3_MES)
MES_4 <- as.data.frame(et_4_MES)
MES_5 <- as.data.frame(et_5_MES)
MES_6 <- as.data.frame(et_6_MES)
MES_7 <- as.data.frame(et_7_MES)

MES_2$cazy <- rownames(MES_2)
MES_3$cazy <- rownames(MES_3)
MES_4$cazy <- rownames(MES_4)
MES_5$cazy <- rownames(MES_5)
MES_6$cazy <- rownames(MES_6)
MES_7$cazy <- rownames(MES_7)

MES_2$Time <- "1_to_2"
MES_3$Time <- "2_to_3"
MES_4$Time <- "3_to_4"
MES_5$Time <- "4_to_5"
MES_6$Time <- "5_to_6"
MES_7$Time <- "6_to_7"

MES_table <- rbind(MES_2, MES_3,MES_4,MES_5,MES_6,MES_7)
MES_table$cazy <- gsub(" ", "", MES_table$cazy)
MES_table$Treatment <- "MES"
rownames(MES_table) <- NULL

### pH10
#first group listed is the baseline
et_2_pH10 <- exactTest(d1, pair=c(22,23))
et_3_pH10 <- exactTest(d1, pair=c(23,24))
et_4_pH10 <- exactTest(d1, pair=c(24,25))
et_5_pH10 <- exactTest(d1, pair=c(25,26))
et_6_pH10 <- exactTest(d1, pair=c(26,27))
et_7_pH10 <- exactTest(d1, pair=c(27,28))

topTags(et_2_pH10, n=5)
topTags(et_3_pH10, n=5)
topTags(et_4_pH10, n=5)
topTags(et_5_pH10, n=5)
topTags(et_6_pH10, n=5)
topTags(et_7_pH10, n=5)

pH10_2 <- as.data.frame(et_2_pH10)
pH10_3 <- as.data.frame(et_3_pH10)
pH10_4 <- as.data.frame(et_4_pH10)
pH10_5 <- as.data.frame(et_5_pH10)
pH10_6 <- as.data.frame(et_6_pH10)
pH10_7 <- as.data.frame(et_7_pH10)

pH10_2$cazy <- rownames(pH10_2)
pH10_3$cazy <- rownames(pH10_3)
pH10_4$cazy <- rownames(pH10_4)
pH10_5$cazy <- rownames(pH10_5)
pH10_6$cazy <- rownames(pH10_6)
pH10_7$cazy <- rownames(pH10_7)

pH10_2$Time <- "1_to_2"
pH10_3$Time <- "2_to_3"
pH10_4$Time <- "3_to_4"
pH10_5$Time <- "4_to_5"
pH10_6$Time <- "5_to_6"
pH10_7$Time <- "6_to_7"

pH10_table <- rbind(pH10_2, pH10_3,pH10_4,pH10_5,pH10_6,pH10_7)
pH10_table$cazy <- gsub(" ", "", pH10_table$cazy)
pH10_table$Treatment <- "pH10"
rownames(pH10_table) <- NULL

### pH6
#first group listed is the baseline
et_2_pH6 <- exactTest(d1, pair=c(29,30))
et_3_pH6 <- exactTest(d1, pair=c(30,31))
et_4_pH6 <- exactTest(d1, pair=c(31,32))
et_5_pH6 <- exactTest(d1, pair=c(32,33))
et_6_pH6 <- exactTest(d1, pair=c(33,34))
et_7_pH6 <- exactTest(d1, pair=c(34,35))

topTags(et_2_pH6, n=5)
topTags(et_3_pH6, n=5)
topTags(et_4_pH6, n=5)
topTags(et_5_pH6, n=5)
topTags(et_6_pH6, n=5)
topTags(et_7_pH6, n=5)

pH6_2 <- as.data.frame(et_2_pH6)
pH6_3 <- as.data.frame(et_3_pH6)
pH6_4 <- as.data.frame(et_4_pH6)
pH6_5 <- as.data.frame(et_5_pH6)
pH6_6 <- as.data.frame(et_6_pH6)
pH6_7 <- as.data.frame(et_7_pH6)

pH6_2$cazy <- rownames(pH6_2)
pH6_3$cazy <- rownames(pH6_3)
pH6_4$cazy <- rownames(pH6_4)
pH6_5$cazy <- rownames(pH6_5)
pH6_6$cazy <- rownames(pH6_6)
pH6_7$cazy <- rownames(pH6_7)

pH6_2$Time <- "1_to_2"
pH6_3$Time <- "2_to_3"
pH6_4$Time <- "3_to_4"
pH6_5$Time <- "4_to_5"
pH6_6$Time <- "5_to_6"
pH6_7$Time <- "6_to_7"

pH6_table <- rbind(pH6_2, pH6_3,pH6_4,pH6_5,pH6_6,pH6_7)
pH6_table$cazy <- gsub(" ", "", pH6_table$cazy)
pH6_table$Treatment <- "pH6"
rownames(pH6_table) <- NULL

### pH8
#first group listed is the baseline
et_2_pH8 <- exactTest(d1, pair=c(36,37))
et_3_pH8 <- exactTest(d1, pair=c(37,38))
et_4_pH8 <- exactTest(d1, pair=c(38,39))
et_5_pH8 <- exactTest(d1, pair=c(39,40))
et_6_pH8 <- exactTest(d1, pair=c(40,41))
et_7_pH8 <- exactTest(d1, pair=c(41,42))

topTags(et_2_pH8, n=5)
topTags(et_3_pH8, n=5)
topTags(et_4_pH8, n=5)
topTags(et_5_pH8, n=5)
topTags(et_6_pH8, n=5)
topTags(et_7_pH8, n=5)

pH8_2 <- as.data.frame(et_2_pH8)
pH8_3 <- as.data.frame(et_3_pH8)
pH8_4 <- as.data.frame(et_4_pH8)
pH8_5 <- as.data.frame(et_5_pH8)
pH8_6 <- as.data.frame(et_6_pH8)
pH8_7 <- as.data.frame(et_7_pH8)

pH8_2$cazy <- rownames(pH8_2)
pH8_3$cazy <- rownames(pH8_3)
pH8_4$cazy <- rownames(pH8_4)
pH8_5$cazy <- rownames(pH8_5)
pH8_6$cazy <- rownames(pH8_6)
pH8_7$cazy <- rownames(pH8_7)

pH8_2$Time <- "1_to_2"
pH8_3$Time <- "2_to_3"
pH8_4$Time <- "3_to_4"
pH8_5$Time <- "4_to_5"
pH8_6$Time <- "5_to_6"
pH8_7$Time <- "6_to_7"

pH8_table <- rbind(pH8_2, pH8_3,pH8_4,pH8_5,pH8_6,pH8_7)
pH8_table$cazy <- gsub(" ", "", pH8_table$cazy)
pH8_table$Treatment <- "pH8"
rownames(pH8_table) <- NULL


### combine all into one table

total_table <- rbind(G05_table, G10_table, MES_table, pH10_table, pH6_table, pH8_table)
total_table <- filter(total_table, PValue<0.05)
total_table[c("bin", "unused")] <- str_split_fixed(total_table$cazy, "[.]", 2)
total_table <- total_table[,-c(8)]

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


#read in tax
tax <- read.csv("CR_bin_tax_cluster.csv")



#merge tables
total_table <- merge(hmm, total_table, by ="cazy")
total_table <- merge(total_table, tax, by ="bin")


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

total_table <- total_table[,-c(4:6)]
#write.csv(total_table, "CR_dbcan_edger_time_treatments.csv", row.names = FALSE)

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


cellulases <- c("GH5_","GH5_1", "GH5_2", "GH4_to_5", "GH5_5", "GH6", "GH5_7", "GH5_8", "GH5_11", "GH5_12", "GH5_13", "GH5_15", "GH5_16", "GH5_18", "GH5_19", "GH5_2","GH5_22", "GH5_23", "GH5_24", "GH5_25", "GH5_26", "GH5_27", "GH5_28", "GH5_29", "GH5_30", "GH5_31", "GH5_36", "GH5_37", "GH5_38", "GH5_39", "GH4_to_5", "GH4_to_50", "GH4_to_51", "GH4_to_53", "GH4_to_54", "GH4_to_55", "GH4_to_56", "GH4_to_57", "GH4_to_58", "GH4_to_59", "GH5_5", "GH5_50", "GH5_51", "GH5_53", "GH5_7", "GH5_8", "GH5_9",
                "GH6","GH7","GH8","GH9","GH12","GH16","GH26","GH44","GH45","GH48", "GH51","GH61","GH74","GH131", "AA9", "AA10")
xylanases <- c("GH10","GH11",
               "GH30","GH30_","GH30_2","GH30_3","GH30_5","GH30_7",
               "GH67","GH115","GH129")
ligninases <- c("AA1","AA1_","AA1_1","AA1_2","AA1_3","AA2")
other_oxidases <- c("AA4","AA5","AA5_1","AA5_2","AA6","AA7","AA8","AA12")
GMC_oxidoreductases <- c("AA3","AA3_","AA3_1","AA2_to_3","AA3_3")
amylases <- c("GH13_1","GH12_to_3","GH13_3","GH13_4","GH13_5","GH13_6","GH13_7","GH13_8","GH13_9","GH13_10","GH13_11","GH13_12","GH13_13","GH13_14","GH13_15","GH13_16","GH13_17","GH13_18","GH13_19","GH12_to_30","GH12_to_31","GH12_to_32","GH12_to_33","GH12_to_34","GH12_to_35","GH12_to_36","GH12_to_37","GH12_to_38","GH12_to_39","GH13_30","GH13_31","GH13_32","GH13_33","GH13_34","GH13_35","GH13_36","GH13_37","GH13_38","GH13_39","GH13_40",
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
