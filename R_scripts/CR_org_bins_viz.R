library(ggplot2)
library(scales)
library(reshape2)
library(knitr)
library(RColorBrewer)
library(dplyr)
library(stringr)


#read in summary file
org_table <- read.csv("CR_org_total_table_summary", header=TRUE)
cluster_dat <- read.csv("CR_bin_tax_cluster.csv", header=TRUE)
cluster_dat2 <- cluster_dat[,c(1,9)]


data <- merge(org_table, cluster_dat, by = "bin")
data$Treatment<-factor(data$Treatment, levels = c("pH6","MES","G05","G10","pH8","pH10"))
data$Day <- as.character(data$Day)
data$cluster_total <- as.character(data$cluster_total)


#summary, add up percents per bin
summary <- data %>%
  group_by(Treatment, Replicate, Day, cluster_total) %>%
  summarise(sum=sum(percent))

#summary, get mean and sd between replicates
summary2 <- summary %>%
  group_by(Treatment, Day, cluster_total) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count=n())

summary2$Treatment<-factor(summary2$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
summary2$cluster_total <- as.character(summary2$cluster_total)
summary2$Day <- as.numeric(summary2$Day)

summary2 <- summary2 %>%
  mutate(across('cluster_total', str_replace, '0', 'unclustered'))



palette <-c( "#56B4E9", "#999999","#E69F00","#0072B2", "#F0E442","#CC79A7")
palette2 <-c("#66a61e","#7570b3","#e7298a","#d95f02")
clustercol <- c("#bcd190","#9f9bca","#ee6aad","#f0a463","#f4dc97")


pline <- ggplot(summary2,aes(y=mean, x= Day, color=cluster_total)) + 
  geom_point(size=2.5) + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  theme_classic(base_size=15) + #make it pretty
  theme(#axis.text.x=element_blank(),
    legend.position="bottom") + 
  #theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
  ylab("Relative Abundance (% Reads Mapped)") + 
  xlab("Time (Day)") +
  labs(color='Cluster')+ #change the legend title
  scale_color_manual(values=clustercol) +
  facet_grid( . ~ Treatment, scales = "free_x")
pline


##count bins by treatment

data2 <- data %>%
  group_by(Treatment, bin, Day, order) %>%
  summarise(count = n())

data2 <- data2 %>%
  group_by(bin,order) %>%
  summarise(count = n())

data2$order <- factor(data2$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                        "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                        "Sphingobacteriales","Xanthomonadales","no support"))

p_bin <- ggplot(data2, aes(x=order, color = bin, fill = order)) +
  geom_bar() +
  theme_classic(base_size=15) +
  theme(axis.text.x = element_text(angle= 45, hjust = 1)) +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=black) +
  labs(y = "Count", x = "Cluster", main = "Count of MAGs by Order", fill = "Order", main = "Count of Bins by Cluster") +
  guides(color = "none")
p_bin

######## count of clusters by order


cluster_dat$cluster_total <- as.character(cluster_dat$cluster_total)
cluster_dat <- cluster_dat %>%
  mutate(across('cluster_total', str_replace, '0', 'unclustered'))


cluster_dat$order <- factor(cluster_dat$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                        "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                        "Sphingobacteriales","Xanthomonadales","no support"))
cluster_dat$class <- factor(cluster_dat$class, levels=c("Actinobacteria","Alphaproteobacteria","Bacilli","Betaproteobacteria","Flavobacteriia","Gammaproteobacteria","Sphingobacteriia","no support"))
#pal for organisms
palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
black <- c("black","black","black","black","black","black","black","black","black","black",
           "black","black","black","black","black","black","black","black","black","black",
           "black","black","black","black","black","black","black","black","black","black",
           "black","black","black","black","black","black","black","black","black","black",
           "black","black","black","black","black","black","black","black","black","black",
           "black","black","black","black")

clustersummaryorder <- ggplot(cluster_dat, aes(x=cluster_total, color = bin, fill = order)) +
  geom_bar() +
  theme_classic(base_size=15) +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=black) +
  labs(y = "Count", x = "Cluster", main = "Count of MAGs by Order", fill = "Order", main = "Count of Bins by Cluster") +
  guides(color = "none")

clustersummaryorder


#########check in on unclustered bins

check <- merge(org_table, cluster_dat, by = "bin")
check <- subset(check, cluster_total == "0")



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
