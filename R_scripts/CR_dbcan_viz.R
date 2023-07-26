library(ggplot2)
library(reshape2)
library(knitr)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(ggridges)
library(ggpubr)

cazy_total <-read.csv("CR_dbcan_total_table.csv", header=TRUE)


cazy_total$Treatment<-factor(cazy_total$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
cazy_total$func <- factor(cazy_total$func, levels = c("cellulase", "betaglucanase", "amylase", "xylanase",
                                                      "ligninase", "GMC oxidoreductases", "other oxidases", "other"))

cazy_total$order <- factor(cazy_total$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                      "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                      "Sphingobacteriales","Xanthomonadales","no support"))


cazy_total_2 <- cazy_total[,c(7:18,20:23)]

AA <- c("GMC oxidoreductases","other oxidases","ligninase")
GH <- c("amylase","betaglucanase","cellulase","xylanase")
other <- c("other")

cazy_total_2 <- cazy_total_2 %>%
  mutate(func2 = case_when(
    func %in% AA ~ "AA",
    func %in% GH ~ "GH",
    func %in% other ~ "other"
  ))

cazy_total_2 <- subset(cazy_total_2, func != "other")

palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
palette <-c( "#56B4E9","#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")


#### boxplot

total_table <- cazy_total

data_trim <- total_table %>%
  group_by(Treatment,Day, Replicate) %>%
  summarise(sum=sum(percent),
            count=n())
data_trim$Day <- as.character(data_trim$Day)

#palette <-c( "#56B4E9", "#999999","#E69F00","#0072B2", "#F0E442","#CC79A7")

pbox <- ggplot(data_trim, aes(x=Day, y=sum, fill=Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values=palette) +
  #scale_color_manual(values=pal1) +
  theme_classic(base_size=20) +
  #ylim(0,19) +
  ylab("Percent reads mapped") + 
  xlab("Time (Day sampled)") +
  #theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(aspect.ratio=1)+
  ggtitle("CAZy Expression") +
  facet_wrap(  ~ Treatment , scales = "free_x")
pbox




#counts of org by function

cazy_trim <- cazy_total[,c(8,15,20)]
cazy_trim <-subset(cazy_trim, func != "other")
cazy_trim$func <- str_to_title(cazy_trim$func)
cazy_trim <- cazy_trim %>%
  group_by(func, order) %>%
  summarise(value=sum(value))
sums <- cazy_trim %>%
  group_by(func) %>%
  summarise(sum=sum(value))

cazy_trim <- merge(cazy_trim, sums, by = "func")

cazy_trim$percent <- ((cazy_trim$value)/(cazy_trim$sum))*100

sums$sum <- format(sums$sum, big.mark = ",", scientific=FALSE)
sums$percent <- 102

cazy_trim$func <- factor(cazy_trim$func, levels=c("Cellulase","Betaglucanase","Xylanase","Amylase","Ligninase","Gmc Oxidoreductases","Other Oxidases"))

pcountsbar <- ggplot() + 
  geom_bar(data=cazy_trim, aes(x=func, y = percent,fill = order, color=order), stat='identity') +
  geom_text(data=sums, aes(x=func, y= percent, label=sum)) +
  theme_classic(base_size=15) +
  theme(axis.text.x = element_text(angle= 25, hjust=1),
        aspect.ratio = 1) +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  labs (x = "CAZyme Function", y = "Percent Expression \n (Summed Across Treatment and Time)")+
  guides(color="none", fill=guide_legend("Order")) +
  ggtitle("Relative Contribution of Bacterial Orders to CAZyme Function")
pcountsbar




cazy_trim2 <- cazy_total[,c(8,15,20,21)]
cazy_trim2 <-subset(cazy_trim2, func != "other")
cazy_trim2$func <- str_to_title(cazy_trim2$func)
cazy_trim2 <- cazy_trim2 %>%
  group_by(func, order, Treatment) %>%
  summarise(value=sum(value))
sums2 <- cazy_trim2 %>%
  group_by(func, Treatment) %>%
  summarise(sum=sum(value))
cazy_trim2 <- merge(cazy_trim2, sums2, by = c("func"="func", "Treatment"="Treatment"))
cazy_trim2$percent <- ((cazy_trim2$value)/(cazy_trim2$sum))*100
sums2$sum <- format(sums2$sum, big.mark = ",", scientific=FALSE)
sums2$percent <- 103

cazy_trim2$func <- factor(cazy_trim2$func, levels=c("Cellulase","Betaglucanase","Xylanase","Amylase","Ligninase","Gmc Oxidoreductases","Other Oxidases"))


pcountsbar_treatment <- ggplot() + 
  geom_bar(data=cazy_trim2, aes(x=func, y = percent,fill = order, color=order), stat='identity') +
  geom_text(data=sums2, aes(x=func, y= percent, label=sum), position=position_dodge(width=0.9), size =3) +
  theme_classic(base_size=15) +
  theme(axis.text.x = element_text(angle= 25, hjust=1),
        aspect.ratio = 1) +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  labs (x = "CAZyme Function", y = "Percent Expression \n (Summed Across Time)")+
  guides(color="none", fill=guide_legend("Order")) +
  ggtitle("Relative Contribution of Bacterial Orders to CAZyme Function") +
  facet_grid(~Treatment, scales ="free")
pcountsbar_treatment



##focus on ligninase
cazy_trim3 <- subset(cazy_total, func == "ligninase")
cazy_trim3 <- cazy_trim3[,c(3,8:10,15,21:23)]
cazy_trim3 <- cazy_trim3 %>%
  group_by(CAZy, order, Treatment, Day, Replicate) %>%
  summarise(value=sum(value))
cazy_trim3 <- cazy_trim3 %>%
  group_by(CAZy, order, Treatment, Day) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            count=n())



pallig <- c("#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a")

plotter_ligcounts <- function(...){
  return(...%>%
           ggplot(aes(x=Day, y = mean, color = CAZy)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd)) +
  theme_classic(base_size=15) +
  scale_color_manual(values=pallig) +
  labs (x = "Day", y = "Number of Reads")+
  facet_grid(~Treatment, scales = "free")
  )
}

sphingo_lig <- subset(cazy_trim3, order == "Sphingobacteriales")
micrococ_lig <- subset(cazy_trim3, order == "Micrococcales")
burk_lig <- subset(cazy_trim3, order == "Burkholderiales")

psphingo <- plotter_ligcounts(sphingo_lig) + ggtitle("Sphingobacteriales Expression of AA1 and AA2 Contigs")
pmicrococ <- plotter_ligcounts(micrococ_lig) + ggtitle("Micrococcales Expression of AA1 and AA2 Contigs")
pburk <- plotter_ligcounts(burk_lig) + ggtitle("Burkholderiales Expression of AA1 and AA2 Contigs")

psphingo
pmicrococ
pburk 


ggarrange(psphingo + rremove("xlab"), 
          pburk + rremove("xlab"),
          pmicrococ,
       common.legend = TRUE, legend = "right", ncol= 1)




#####graph cazy function by time group all treatments #####

#sum up by func
data <- cazy_total_2 %>%
  group_by(func, Day, Treatment, Replicate) %>%
  summarise(sum=sum(percent))
#average across treatments and replicates
data <- data %>%
  group_by(func, Treatment, Day) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count=n())
data$Day <- as.character(data$Day)

plotter <- function(...){
  return(...%>%
           ggplot(aes(y=sum, x = Day, fill = Treatment)) +
           geom_boxplot() +
           #geom_point(size = 2.5) +
          # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
           theme_classic(base_size=15) +
           labs(fill='CAZy Function', y='Percent reads mapped (%)', x = 'Day')+
           scale_fill_manual(values=palette) +
           facet_grid(~ Treatment))
}


lignin <- subset(data, func =="ligninase")
plignin <- plotter(lignin) + ggtitle('Lignin-acting Contigs Expression')
plignin

gmcoxi <- subset(data, func =="GMC oxidoreductases")
pgmcoxi <- plotter(gmcoxi) + ggtitle('GMC oxidoreductas Contigs Expression')

otheroxi <- subset(data, func =="other oxidases")
potheroxi <- plotter(otheroxi) + ggtitle('Other Oxidase Contigs Expression')

cellu <- subset(data, func =="cellulase")
pcellu <- plotter(cellu) + ggtitle('Cellulase Contigs Expression')

betaglucanase <- subset(data, func =="betaglucanase")
pbetaglucanase <- plotter(betaglucanase) + ggtitle('Betaglucanase Contigs Expression')

amylase <- subset(data, func =="amylase")
pamylase <- plotter(amylase) + ggtitle('Amylase Contigs Expression')

xylanase <- subset(data, func =="xylanase")
pxylanase <- plotter(xylanase) + ggtitle('Xylanase Contigs Expression')

plignin
pgmcoxi 
potheroxi
pcellu
pbetaglucanase
pamylase
pxylanase


##ligninases have diff trend than other functional groups
#counts graphs


palorg2<- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")

plotter_counts <- function(...){
  return(...%>%
           ggplot(aes(x=variable, y = value, color = order, fill = order)) +
           geom_bar(stat="identity") +
           #scale_color_manual(values = palorg) +
           #scale_fill_manual(values = palorg) +
           theme_classic(base_size=15) +
           theme(axis.text.x = element_blank())+
           labs(x= "Samples, arranged by Day", y = "Read Counts") +
           guides(color = "none", fill = guide_legend("Order")) +
           facet_grid(~Treatment, scales = "free")
         )
}

lignin <- subset(cazy_total, func =="ligninase")
plignin <- plotter_counts(lignin) + ggtitle('Lignin-acting Contigs Expression') + scale_color_manual(values = palorg2) + scale_fill_manual(values = palorg2)
plignin

gmcoxi <- subset(cazy_total, func =="GMC oxidoreductases")
pgmcoxi <- plotter_counts(gmcoxi) + ggtitle('GMC oxidoreductas Contigs Expression') + scale_color_manual(values = palorg) + scale_fill_manual(values = palorg)

otheroxi <- subset(cazy_total, func =="other oxidases")
potheroxi <- plotter_counts(otheroxi) + ggtitle('Other Oxidase Contigs Expression') + scale_color_manual(values = palorg) + scale_fill_manual(values = palorg)

cellu <- subset(cazy_total, func =="cellulase")
pcellu <- plotter_counts(cellu) + ggtitle('Cellulase Contigs Expression') + scale_color_manual(values = palorg) + scale_fill_manual(values = palorg)

betaglucanase <- subset(cazy_total, func =="betaglucanase")
pbetaglucanase <- plotter_counts(betaglucanase) + ggtitle('Betaglucanase Contigs Expression') + scale_color_manual(values = palorg) + scale_fill_manual(values = palorg)

amylase <- subset(cazy_total, func =="amylase")
pamylase <- plotter_counts(amylase) + ggtitle('Amylase Contigs Expression') + scale_color_manual(values = palorg) + scale_fill_manual(values = palorg)

xylanase <- subset(cazy_total, func =="xylanase")
pxylanase <- plotter_counts(xylanase) + ggtitle('Xylanase Contigs Expression') + scale_color_manual(values = palorg) + scale_fill_manual(values = palorg)

pgmcoxi
potheroxi
pcellu #needs new color palette
pbetaglucanase
pamylase #needs new palette
pxylanase 


#### AA only ###
#get just AA
data2org <- subset(cazy_total, CAZy == "AA2")

top <- c("Sphingobacteriales")

data2org <- subset(data2org, order %in% top)
data2org$order <- factor(data2org$order, levels = c("Sphingobacteriales"))


#group by treatment, day ,rep
data2org <- data2org %>%
  group_by(order, Treatment, Day, Replicate) %>%
  summarise(sum=sum(percent))

data2org <- data2org %>%
  group_by(order, Treatment, Day) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count = n())


paa2top <- ggplot(data2org, aes(y=mean, x= Day, color = order)) +
  geom_point(size=2) +
  geom_errorbar(data=data2org, aes(ymax=mean+sd, ymin=mean-sd)) +
  theme_classic(base_size=15) +
  #scale_color_manual(values = c("#1f78b4")) +
  #guides(color="none") +
  labs(x = "Day", y = "Percent reads mapped (%)", title = "Expression of AA2 contigs classified to Sphingobacteriales") +
  facet_grid(~Treatment, scales = "free")
paa2top




palorg2 <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")

pligorg <- ggplot(data2org, aes(y=sum, x = variable, fill= order, color=order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=15) + 
  theme(#axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "bottom") +
  labs(fill='Order', y='Percent reads mapped (%)', title = "Expression of AA2 contigs by Order") +
  scale_fill_manual(values=palorg2) +
  scale_color_manual(values=palorg2) +
  facet_grid(~ Treatment, scales="free") +
  guides(color = "none", fill = guide_legend(ncol=10))
pligorg




unique(cazy_total$func)

#### GHonly ###

data4 <- subset(cazy_total_2, func2 == "GH")
data4 <- data4 %>%
  group_by(variable, func, Treatment, Day, Replicate) %>%
  summarise(sum=sum(percent))



data4 <- subset(data4, func == "amylase")

pgh <- ggboxplot(data4, x="Day", y="sum", fill = "Treatment", lwd=1) +
  theme_classic(base_size=22)+
  theme(aspect.ratio=1) +
  scale_fill_manual(values=palette) + 
  labs(y="Percent Reads", x="Day") + 
  ggtitle("amylase Expression") + 
  facet_grid( ~ Treatment) +
  stat_compare_means(method="t.test", aes(group=Day), label="p.signif")
pgh


#get just GH
data4org <- subset(cazy_total_2, func2 == "GH")
#group by treatment, day ,rep
data4org <- data4org %>%
  group_by(func, order, Treatment, Day, Replicate) %>%
  summarise(sum=sum(percent))
#average across all
data4org <- data4org %>%
  group_by(func, order, Treatment, Day) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count=n())

data4org <- subset(data4org, func == "cellulase")
org <- c("Enterobacterales",
        "Lactobacillales","Micrococcales",
        "Sphingobacteriales")

data4org <- subset(data4org, order %in% org)
pghorg <- ggplot(data4org, aes(y=mean, x = Day, fill= order, color=order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "bottom") +
  labs(fill='CAZy function', y='Percent reads mapped (%)') +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  facet_grid(~ Treatment) +
  guides(color = "none", fill = guide_legend(ncol=10))
pghorg


### proportion of amylase reads to a specific order ###



#within "variable" sum
#then mean% by "func"

#by order


table_trim <- cazy_total_2 %>%
  group_by(func,Treatment,counts,Day, Replicate, order) %>%
  summarise(sum=sum(percent))

table_trim <- table_trim %>%
  group_by(func,Treatment,order) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count=n())

#colorCount <-length(unique(table_trim$order))
table_trim$order <- factor(table_trim$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                    "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                    "Sphingobacteriales","Xanthomonadales","no support"))


#palette for functional graph
palfunc <-c("#5ab4ac","#d8b365","#8c510a","#fc8d59","#4575b4")
#pal for organisms
palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
#pal for clusters
clustercol <- c("#bcd190","#9f9bca","#ee6aad","#f0a463","#f4dc97")



porg <-ggplot(table_trim, aes(y=mean,x=Treatment,fill= order, color= order)) + 
  geom_bar(stat="identity") +
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "bottom") +
  labs(fill='Order', y='Average (across days sampled) Percent reads mapped (%)') +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  facet_grid( ~ func, scales = "free") +
  guides(color = "none", fill = guide_legend(ncol=10))

porg



#####order by day

cazy_trim <- subset(cazy_total, func == "cellulase")
table_trim <- cazy_trim %>%
  group_by(func,Treatment,counts,Day, Replicate, order) %>%
  summarise(sum=sum(percent))


table_trim <- table_trim %>%
  group_by(func,Treatment, Day, order) %>%
  summarise(mean = mean (sum),
            sd = sd(sum),
            count = n())


table_trim$order <- factor(table_trim$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                      "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                      "Sphingobacteriales","Xanthomonadales","no support"))

porderday <- ggplot(table_trim, aes(y=mean, x = Day, color = order)) +
  geom_point(size=2.5) + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  theme_classic(base_size=15) +
  scale_color_manual(values=palorg) +
  labs(color='Order', y='Percent reads mapped (%)') +
  ggtitle("Ligninase Expression") +
  facet_grid ( ~ Treatment, scales = "free")



porderday






######pl and ce

ce <- cazy_total %>%
  filter(grepl('CE', CAZy))

length(unique(ce$bin))

ce <- ce %>%
  group_by(variable, order, Treatment, Day, Replicate) %>%
  summarise(sum=sum(percent))


pce <- ggplot(ce, aes(y=sum, x = variable, fill= order, color=order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        legend.position = "bottom") +
  labs(fill='Order', y='Percent reads mapped (%)') +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  facet_grid(~ Treatment, scales = "free") +
  ggtitle("Carbohydrate Esterases") +
  guides(color = "none", fill = guide_legend(ncol=10))
pce

pl <- cazy_total %>%
  filter(grepl('PL', CAZy))

length(unique(pl$bin))

pl <- pl %>%
  group_by(variable, order, Treatment, Day, Replicate) %>%
  summarise(sum=sum(percent))

palorg3 <- c("#FFFFB3","#BEBADA","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
ppl <- ggplot(pl, aes(y=sum, x = variable, fill= order, color=order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        legend.position = "bottom") +
  labs(fill='Order', y='Percent reads mapped (%)') +
  scale_fill_manual(values=palorg3) +
  scale_color_manual(values=palorg3) +
  facet_grid(~ Treatment, scales = "free") +
  ggtitle("Polysaccharaide lyases") +
  guides(color = "none", fill = guide_legend(ncol=10))
ppl



gh <- cazy_total %>%
  filter(grepl('GH', CAZy))

length(unique(gh$bin))

gh <- gh %>%
  group_by(variable, order, Treatment, Day, Replicate) %>%
  summarise(sum=sum(percent))

pgh <- ggplot(gh, aes(y=sum, x = variable, fill= order, color=order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        legend.position = "bottom") +
  labs(fill='Order', y='Percent reads mapped (%)') +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  facet_grid(~ Treatment, scales = "free") +
  ggtitle("Glycoside hydrolases") +
  guides(color = "none", fill = guide_legend(ncol=10))
pgh

aa <- cazy_total %>%
  filter(grepl('AA', CAZy))

length(unique(aa$bin))

aa <- aa %>%
  group_by(variable, order, Treatment, Day, Replicate) %>%
  summarise(sum=sum(percent))

paa <- ggplot(aa, aes(y=sum, x = variable, fill= order, color=order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        legend.position = "bottom") +
  labs(fill='Order', y='Percent reads mapped (%)') +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  facet_grid(~ Treatment, scales = "free") +
  ggtitle("Auxiliary Activity") +
  guides(color = "none", fill = guide_legend(ncol=10))
paa


gt <- cazy_total %>%
  filter(grepl('GT', CAZy))

length(unique(gt$bin))

gt <- gt %>%
  group_by(variable, order, Treatment, Day, Replicate) %>%
  summarise(sum=sum(percent))

pgt <- ggplot(gt, aes(y=sum, x = variable, fill= order, color=order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        legend.position = "bottom") +
  labs(fill='Order', y='Percent reads mapped (%)') +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  facet_grid(~ Treatment, scales = "free") +
  ggtitle("Glycosyl transferases") +
  guides(color = "none", fill = guide_legend(ncol=10))
pgt

ggarrange(pgh, pgt, paa, pce, ppl, common.legend = TRUE, legend = "bottom")

#tax_matrix <- dcast(data=total_table, formula = bin~variable, fun.aggregate = sum, value.var="value")
#write.csv(tax_matrix, "CR_dbCAN_taxmatrix.csv", row.names = FALSE)


total_table <- cazy_total

data_trim <- total_table %>%
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
  ylab("Percent reads mapped") + 
  xlab("Time (Day sampled)") +
  #theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(aspect.ratio=1)+
  ggtitle("CAZy Expression") +
  facet_wrap(  ~ Treatment , scales = "free_x")
pbox


