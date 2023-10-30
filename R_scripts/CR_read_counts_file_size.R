library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(plyr)
library(ggpubr)

#Read in files, rename columns for easy merging
mg_mt_mgm <- fread("all_reads_counts.tsv", sep="\t")
colnames(mg_mt_mgm) <- c("counts", "Sample")
mtm <- fread("total_reads_mt.txt", sep="\t")
colnames(mtm)<- c("Sample", "counts")

#combine files, make dataframe
data <- rbind(mg_mt_mgm,mtm)
data <- as.data.frame(data)


#string-split into useful data
data[c("Treatment","Day","Replicate","Sample_type")] <- str_split_fixed(data$Sample, "_", 4)

#replace file extensions into useful data
new_str = c("mg.fastq.1.gz.counts"="Metagenome",
            "mt.R1.fastq.cleaned.forward"="Metatranscriptome",
            "mg.pear.fastq.gz.counts"="MG_merged",
            "mt.R1.fastq.merged.ribodepleted.fastq.gz.counts"="MT_merged")
data$Sample_type <- str_replace_all(data$Sample_type, new_str)

#divide counts by 1 million for easier readability
data$counts2 <- data$counts/1000000

data$Treatment<-factor(data$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
data$Replicate<-factor(data$Replicate, levels = c("F","E","D","C","B","A"))


#summarise

data_summary <- data %>%
  group_by(Sample_type) %>%
  summarise(mean=mean(counts2),
            median=median(counts2),
           sd=sd(counts2),
           min=min(counts2),
           max=max(counts2))


hist <- ggplot(data, aes(x=counts2)) +
  geom_histogram(binwidth=1) +
  facet_grid(~Sample_type, scales = "free")
hist

mgm_list <- subset(data, Sample_type =="MG_merged" & counts2 > 14)
mtm_list <- subset(data, Sample_type =="MT_merged" & counts2 > 30)


#write.csv(mgm_list, "mg_14m_list.csv", row.names = FALSE)
#write.csv(mtm_list, "mt_30m_list.csv", row.names = FALSE)
#add missing data



rng = range(data$counts2)

data_mg <- subset(data, Sample_type == "Metagenome")

data_mg_summary <- data_mg %>%
  summarise(mean = mean(counts2),
            sd=sd(counts2))

#add missing data so I can fill NA
Treatment <- c("pH6","pH6", "G05","G05","G05","G05","G05","G10","pH8", "pH10","pH10","pH10","pH10","pH10","pH10")
Day <- c("6","7","2","4","7","7","7","7","7","3","5","5","5","5","5")
Replicate <- c("C","E","D","C","D","E","F","D","E","E","A","B","C","D","E")
temp <- data.frame(Treatment, Day, Replicate)

data_mg <- rbind.fill(data_mg, temp)
data_mg$Treatment<-factor(data_mg$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
data_mg$Replicate<-factor(data_mg$Replicate, levels = c("F","E","D","C","B","A"))



p1 <- ggplot(data_mg, aes(y=Replicate, x=Day, fill = counts2)) +
  geom_tile() +
  scale_fill_gradient2(na.value = "black", low="#13378c", mid="#43b7c4", high="#ffffcc", midpoint = mean(rng), breaks = seq(0,100, 25),limits=c(floor(rng[1]),ceiling(rng[2]))) +
  theme_classic(base_size = 15) +
  coord_fixed() +
  labs(fill="Read pairs per million") +
  ggtitle("Number of Metagenomic Read Pairs") +
  facet_grid (~Treatment)
p1

data_mt <- subset(data, Sample_type == "Metatranscriptome")
data_mt_summary <- data_mt %>%
  summarise(mean = mean(counts2),
            sd=sd(counts2))


#add missing data so I can fill NA
Treatment <- c("pH6","pH6","pH6",
               "MES","MES","MES","MES","MES","MES","MES","MES",
               "G05","G05","G05","G05",
               "G10","G10",
               "pH8","pH8","pH8","pH8","pH8","pH8",
               "pH10","pH10","pH10","pH10","pH10","pH10","pH10","pH10","pH10","pH10","pH10","pH10")
Day <- c("4","4","7",
         "1","1","1","1","3","4","4","6",
         "4","7","7","7",
         "7","7",
         "1","3","3","4","4","7",
         "1","3","3","3","3","4","4","4","4","4","6","7")
Replicate <- c("A","D","E",
               "A","B","D","F","E","D","E","F",
               "A","D","E","F",
               "A","D",
               "A","B","C","A","D","E",
               "B","A","B","D","E","A","B","C","D","E","D","A")
temp2 <- data.frame(Treatment, Day, Replicate)
data_mt <- rbind.fill(data_mt, temp2)
data_mt$Treatment<-factor(data_mt$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
data_mt$Replicate<-factor(data_mt$Replicate, levels = c("F","E","D","C","B","A"))

p2 <- ggplot(data_mt, aes(y=Replicate, x=Day, fill = counts2)) +
  geom_tile() +
  scale_fill_gradient2(na.value = "black", low="#13378c", mid="#43b7c4", high="#ffffcc", midpoint = mean(rng), breaks = seq(0,100, 25),limits=c(floor(rng[1]),ceiling(rng[2]))) +
  theme_classic(base_size = 15) +
  coord_fixed() +
  labs(fill="Read pairs per million") +
  ggtitle("Number of Metatranscriptomic Read Pairs") +
  facet_grid (~Treatment)
p2


data_mgm <- subset(data, Sample_type == "MG_merged")
data_mgm <- rbind.fill(data_mgm, temp)
data_mgm$Treatment<-factor(data_mgm$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
data_mgm$Replicate<-factor(data_mgm$Replicate, levels = c("F","E","D","C","B","A"))


p3 <- ggplot(data_mgm, aes(y=Replicate, x=Day, fill = counts2)) +
  geom_tile() +
  scale_fill_gradient2(na.value = "black", low="#250e35", mid="#b13657", high="#ffffcc", midpoint = mean(rng), breaks = seq(0,100, 25),limits=c(floor(rng[1]),ceiling(rng[2]))) +
  theme_classic(base_size = 15) +
  coord_fixed() +
  labs(fill="Reads per million") +
  ggtitle("Number of Merged Metagenomic reads") +
  facet_grid (~Treatment)
p3

data_mtm <- subset(data, Sample_type == "MT_merged")
data_mtm_summary <- data_mtm %>%
  summarise(mean = mean(counts2),
            sd=sd(counts2))


data_mtm <- rbind.fill(data_mtm, temp2)
data_mtm$Treatment<-factor(data_mtm$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
data_mtm$Replicate<-factor(data_mtm$Replicate, levels = c("F","E","D","C","B","A"))
p4 <- ggplot(data_mtm, aes(y=Replicate, x=Day, fill = counts2)) +
  geom_tile() +
  scale_fill_gradient2(na.value = "black",low="#250e35", mid="#b13657", high="#ffffcc", midpoint = mean(rng), breaks = seq(0,100, 25),limits=c(floor(rng[1]),ceiling(rng[2]))) +
  theme_classic(base_size = 15) +
  coord_fixed() +
  labs(fill="Reads per million") +
  ggtitle("Number of Merged Metatranscriptomic reads") +
  facet_grid (~Treatment)
p4




ggarrange(p1, p2, p3, p4,ncol=1, labels = "AUTO", font.label = list(size = 20))

ggsave("read_counts.tiff", device = "tiff", dpi = 700)
ggsave("read_counts.png", device = "png", dpi = 700)
ggsave("read_counts.pdf", device = "pdf", dpi = 700)



######### file size table ###


data$Sample <-str_replace(data$Sample, ".fastq.1.gz.counts", "")
data$Sample <-str_replace(data$Sample, ".R1.fastq.cleaned.forward", "")

data_keep <- subset(data, Sample_type %in% c("Metagenome", "Metatranscriptome"))

data_mg <- subset(data, Sample_type == "Metagenome")
data_mt <- subset(data, Sample_type == "Metatranscriptome")



d_size <- read.csv("CR_filesizes.csv", header = TRUE)

d_size_summ <- d_size %>%
  group_by(Dataset) %>%
  summarise(mean = mean(Size_GB),
            sd=sd(Size_GB),
            sum=sum(Size_GB))
