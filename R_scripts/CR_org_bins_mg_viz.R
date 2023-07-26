library(ggplot2)
library(scales)
library(reshape2)
library(knitr)
library(RColorBrewer)
library(dplyr)
library(stringr)


#read in summary file
otu <- read.csv("CR_org_bins_mg_matrix_relative.csv", header=TRUE)
cluster_dat <- read.csv("CR_bin_tax_cluster.csv", header=TRUE)

org_table <- melt(otu, id = c("variable"))
colnames(org_table)[2] <- "bin"

#by org

data <- merge(org_table, cluster_dat, by = "bin")
data[c("Treatment","Day","Replicate")] <- str_split_fixed(data$variable, "_", 3)

data$Treatment<-factor(data$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
data$Day <- as.numeric(data$Day)
data$order <- factor(data$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                          "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                          "Sphingobacteriales","Xanthomonadales","no support"))
#data$class <- factor(data$class, levels=c("Actinobacteria","Alphaproteobacteria","Bacilli","Betaproteobacteria","Flavobacteriia","Gammaproteobacteria","Sphingobacteriia","no support"))
#pal for organisms
palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")



pbarorg <- ggplot() +
  geom_bar(data=data, aes(y=value, x=variable, color=order, fill=order),stat="identity") +
  theme_classic(base_size=15) +
  theme(axis.text.x = element_blank()) +
  labs(y = "Percent Metagenomic Reads Mapped", x = "Samples Arranged by Time", fill = "Order") +
  guides(color= "none") +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg) +
  ggtitle("Relative Abundance") +
  facet_grid(. ~ Treatment, scales="free_x")

pbarorg




##add in mt data

data$bin_variable <- paste(data$bin, data$variable, sep = ".")
mt <- read.csv("CR_perentmtabundance.csv")
data <- merge(data, mt, by = "bin_variable")
data$mt_by_mg <- (data$percent_mt)/(data$value)


pbarmtmg <- ggplot() +
  geom_bar(data=data, aes(y=mt_by_mg, x=variable, color=order, fill=order),stat="identity") +
  theme_classic(base_size=15) +
  theme(axis.text.x = element_blank()) +
  labs(y = "MT/MG", x = "Samples Arranged by Time", fill = "Order") +
  guides(color= "none") +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg) +
  facet_grid(. ~ Treatment, scales="free_x")

pbarmtmg



#summary, add up percents per bin
summary <- data %>%
  group_by(Treatment, Replicate, Day, order) %>%
  summarise(sum=sum(percent))

#summary, get mean and sd between replicates
summary2 <- summary %>%
  group_by(Treatment, Day, order) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count=n())

summary2$Treatment<-factor(summary2$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
summary2$Day <- as.numeric(summary2$Day)
summary2$order <- factor(summary2$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                        "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                        "Sphingobacteriales","Xanthomonadales","no support"))
#summary2$class <- factor(summary2$class, levels=c("Actinobacteria","Alphaproteobacteria","Bacilli","Betaproteobacteria","Flavobacteriia","Gammaproteobacteria","Sphingobacteriia","no support"))
#pal for organisms
palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")


plineorg <- ggplot(summary2,aes(y=mean, x= Day, color=order)) + 
  geom_point(size=2.5) + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  theme_classic(base_size=15) + #make it pretty
  theme(#axis.text.x=element_blank(),
    legend.position="bottom") + 
  #theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
  ylab("Percent Reads Mapped") + 
  xlab("Time (Day)") +
  labs(fill='Cluster')+ #change the legend title
  scale_color_manual(values=palorg) +
  facet_grid( . ~ Treatment, scales = "free_x")
plineorg


#organism bar graph

#oops lost "variable"
summary$variable <- paste(summary$Treatment,summary$Day,summary$Replicate, sep="_")
summary$Treatment<-factor(summary$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
summary$order <- factor(summary$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                  "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                  "Sphingobacteriales","Xanthomonadales","no support"))

pbarorg <- ggplot() +
  geom_bar(data=summary, aes(y=sum, x=variable, color=order, fill=order),stat="identity") +
  theme_classic() +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg) +
  facet_grid(. ~ Treatment, scales="free_x")

pbarorg


#by cluster

data <- merge(org_table, cluster_dat, by = "bin")
data$Treatment<-factor(data$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
data$Day <- as.character(data$Day)
data$cluster <- as.character(data$cluster)


#summary, add up percents per bin
summary <- data %>%
  group_by(Treatment, Replicate, Day, cluster) %>%
  summarise(sum=sum(percent))

#summary, get mean and sd between replicates
summary2 <- summary %>%
  group_by(Treatment, Day, cluster) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count=n())

summary2$Treatment<-factor(summary2$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
summary2$cluster <- as.character(summary2$cluster)
summary2$Day <- as.numeric(summary2$Day)

summary2 <- summary2 %>%
  mutate(across('cluster', str_replace, '0', 'unclustered'))



palette <-c( "#56B4E9", "#999999","#E69F00","#0072B2", "#F0E442","#CC79A7")
palette2 <-c("#66a61e","#7570b3","#e7298a","#d95f02")
clustercol <- c("#bcd190","#9f9bca","#ee6aad","#f0a463","#f4dc97")


pline <- ggplot(summary2,aes(y=mean, x= Day, color=cluster)) + 
  geom_point(size=2.5) + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  theme_classic(base_size=15) + #make it pretty
  theme(#axis.text.x=element_blank(),
    legend.position="bottom") + 
  #theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
  ylab("Percent Reads Mapped") + 
  xlab("Time (Day)") +
  labs(fill='Cluster')+ #change the legend title
  scale_color_manual(values=clustercol) +
  facet_grid( . ~ Treatment, scales = "free_x")
pline



#########check in on unclustered bins

check <- merge(org_table, cluster_dat, by = "bin")
check <- subset(check, cluster == "0")



check <- check %>%
  group_by(Treatment, Replicate, Day, order) %>%
  summarise(sum=sum(percent))

check <- check %>%
  group_by(Treatment, Day, order) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count=n())
check$Treatment<-factor(check$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))


pline <- ggplot(check,aes(y=mean, x= Day, color=order)) + 
  geom_point(size=2.5) + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  theme_classic(base_size=15) + #make it pretty
  theme(#axis.text.x=element_blank(),
    legend.position="bottom") + 
  #theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
  ylab("Percent Reads Mapped") + 
  xlab("Time (Day)") +
  labs(fill='bin')+ #change the legend title
  ggtitle("Abundance of unclusterd bins") +
  scale_color_manual(values=palette) +
  facet_grid( . ~ Treatment, scales = "free_x")
pline


check2 <- merge(org_table, cluster_dat, by = "bin")
check2 <- subset(check2, cluster == "0")
check2 <- subset(check2, order == "Lactobacillales")
check2 <- subset(check2, Treatment == "MES")

hmm <- check2 %>%
  group_by(Day, Replicate) %>%
  summarise(sum=sum(percent))


######## count of clusters by class


cluster_dat$cluster <- as.character(cluster_dat$cluster)
cluster_dat <- cluster_dat %>%
  mutate(across('cluster', str_replace, '0', 'unclustered'))


summary2$order <- factor(summary2$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                        "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                        "Sphingobacteriales","Xanthomonadales","no support"))
cluster_dat$class <- factor(cluster_dat$class, levels=c("Actinobacteria","Alphaproteobacteria","Bacilli","Betaproteobacteria","Flavobacteriia","Gammaproteobacteria","Sphingobacteriia","no support"))
#pal for organisms
palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")


clustersummaryclass <- ggplot(cluster_dat, aes(x=cluster, color= class, fill = class)) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg)

clustersummaryclass


