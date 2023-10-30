library(ggplot2)
library(reshape2)
library(knitr)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(viridis)
library(cowplot)
library(ggpubr)
library(cowplot)

cazy_total <-read.csv("CR_dbcan_total_table.csv", header=TRUE)
unique(cazy_total$Substrate)

cazy_total$Treatment<-factor(cazy_total$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
cazy_total$Substrate <- factor(cazy_total$Substrate, levels = c("Cellulose", "Cellobiose", "Xylan",
                                                                "Starch", "Trehalose","Pectin","Lignin",
                                                                "Peptidoglycan","Chitin","Other"))

cazy_total$order <- factor(cazy_total$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                      "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                      "Sphingobacteriales","Xanthomonadales","no support"))

##read in cazy mg file
cazy_total_mg <-read.csv("CR_dbcan_total_table_mg.csv", header=TRUE)
cazy_total_mg$Treatment<-factor(cazy_total_mg$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
cazy_total_mg$Substrate <- factor(cazy_total_mg$Substrate, levels = c("Cellulose", "Cellobiose", "Xylan",
                                                          "Starch", "Trehalose","Pectin","Lignin",
                                                          "Peptidoglycan","Chitin","Other"))

palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
palette <-c( "#56B4E9","#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")

palorg <- c("#8DD3C7", #Actinomycetales
            "#FFFFB3", #Burkholderiales
            "#BEBADA", #Caulobacterales
            "#FB8072", #Corynebacteriales
            "#80B1D3", #Enterobacterales
            "#FDB462", #Flavobacteriales
            "#B3DE69", #Lactobacillales
            "#FCCDE5", #Micrococcales
            "#BC80BD", #Propionibacteriales
            "#CCEBC5", #Pseudomonadales
            "#e5c494", #Rhizobiales
            "#1f78b4", #Sphingobacteriales
            "#fb9a99", #Xanthomonadales
            "#D9D9D9" #no support
            )

palsubstrate <- c("#73A87C","#CE9486","#205D89","#CF784B","#FFEA59","#B7DBDB","#B5B867","#77674E","#292176")
palday <- c("#7DCEA0","#52BE80","#27AE60","#229954","#1E8449","#196F3D","#145A32")

### set labels so ggplot graphs only day and not entire "variable" string
## ok so turns out I have to set these manually (rip) because facet_grid() resets from beginning of vector for each facet
##and I have a different number replicates for days between treatments
## rip I can't believe I did this
## is anybody going to care about this
#no

mt_labels <- c("G05_1_A" = "1", "G05_1_B" = "1", "G05_1_C" = "1", "G05_1_D" = "1", "G05_1_E" = "1", "G05_1_F" = "1", 
               "G05_2_A" = "2", "G05_2_B" = "2", "G05_2_C" = "2", "G05_2_D" = "2", "G05_2_E" = "2", "G05_2_F" = "2", 
               "G05_3_A" = "3", "G05_3_B" = "3", "G05_3_C" = "3", "G05_3_D" = "3", "G05_3_E" = "3", "G05_3_F" = "3",
               "G05_4_B" = "4", "G05_4_C" = "4", "G05_4_D" = "4", "G05_4_E" = "4", "G05_4_F" = "4",
               "G05_5_A" = "5", "G05_5_B" = "5", "G05_5_C" = "5", "G05_5_D" = "5", "G05_5_E" = "5", "G05_5_F" = "5",
               "G05_6_A" = "6", "G05_6_B" = "6", "G05_6_C" = "6", "G05_6_D" = "6", "G05_6_E" = "6", "G05_6_F" = "6",
               "G05_7_A" = "7", "G05_7_B" = "7", "G05_7_C" = "7",
               "G10_1_A" = "1", "G10_1_B" = "1", "G10_1_C" = "1", "G10_1_D" = "1", "G10_1_E" = "1", "G10_1_F" = "1",
               "G10_2_A" = "2", "G10_2_B" = "2", "G10_2_C" = "2", "G10_2_D" = "2", "G10_2_E" = "2", "G10_2_F" = "2",
               "G10_3_A" = "3", "G10_3_B" = "3", "G10_3_C" = "3", "G10_3_D" = "3", "G10_3_E" = "3", "G10_3_F" = "3",
               "G10_4_A" = "4", "G10_4_B" = "4", "G10_4_C" = "4", "G10_4_D" = "4", "G10_4_E" = "4", "G10_4_F" = "4",
               "G10_5_A" = "5", "G10_5_B" = "5", "G10_5_C" = "5", "G10_5_D" = "5", "G10_5_E" = "5", "G10_5_F" = "5",
               "G10_6_A" = "6", "G10_6_B" = "6", "G10_6_C" = "6", "G10_6_D" = "6", "G10_6_E" = "6", "G10_6_F" = "7",
               "G10_7_B" = "7", "G10_7_C" = "7", "G10_7_E" = "7", "G10_7_F" = "7",
               "MES_1_C" = "1", "MES_1_E" = "1",
               "MES_2_A" = "2", "MES_2_B" = "2", "MES_2_C" = "2", "MES_2_D" = "2", "MES_2_E" = "2", "MES_2_F" = "2",
               "MES_3_A" = "3", "MES_3_B" = "3", "MES_3_C" = "3", "MES_3_D" = "3", "MES_3_F" = "3",
               "MES_4_A" = "4", "MES_4_B" = "4", "MES_4_C" = "4", "MES_4_F" = "4",
               "MES_5_A" = "5", "MES_5_B" = "5", "MES_5_C" = "5", "MES_5_D" = "5", "MES_5_E" = "5", "MES_5_F" = "5",
               "MES_6_A" = "6", "MES_6_B" = "6", "MES_6_C" = "6", "MES_6_D" = "6", "MES_6_E" = "6",
               "MES_7_A" = "7", "MES_7_B" = "7", "MES_7_C" = "7", "MES_7_D" = "7", "MES_7_E" = "7", "MES_7_F" = "7",
               "pH10_1_A" = "1", "pH10_1_C" = "1", "pH10_1_D" = "1", "pH10_1_E" = "1", "pH10_1_F" = "1",
               "pH10_2_A" = "2", "pH10_2_B" = "2", "pH10_2_C" = "2", "pH10_2_D" = "2", "pH10_2_E" = "2", "pH10_2_F" = "2",
               "pH10_3_C" = "3", "pH10_3_F" = "3",
               "pH10_4_F" = "4",
               "pH10_5_A" = "5", "pH10_5_B" = "5", "pH10_5_C" = "5", "pH10_5_D" = "5", "pH10_5_E" = "5", "pH10_5_F" = "5",
               "pH10_6_A" = "6", "pH10_6_B" = "6", "pH10_6_C" = "6", "pH10_6_E" = "6", "pH10_6_F" = "6",
               "pH10_7_B" = "7", "pH10_7_C" = "7", "pH10_7_D" = "7", "pH10_7_E" = "7", "pH10_7_F" = "7",
               "pH6_1_A" = "1", "pH6_1_B" = "1", "pH6_1_C" = "1", "pH6_1_D" = "1", "pH6_1_E" = "1", "pH6_1_F" = "1",
               "pH6_2_A" = "2", "pH6_2_B" = "2", "pH6_2_C" = "2", "pH6_2_D" = "2", "pH6_2_E" = "2", "pH6_2_F" = "2",
               "pH6_3_A" = "3", "pH6_3_B" = "3", "pH6_3_C" = "3", "pH6_3_D" = "3", "pH6_3_E" = "3", "pH6_3_F" = "3",
               "pH6_4_B" = "4", "pH6_4_C" = "4", "pH6_4_E" = "4", "pH6_4_F" = "4",
               "pH6_5_A" = "5", "pH6_5_B" = "5", "pH6_5_C" = "5", "pH6_5_D" = "5", "pH6_5_E" = "5", "pH6_5_F" = "5",
               "pH6_6_A" = "6", "pH6_6_B" = "6", "pH6_6_C" = "6", "pH6_6_D" = "6", "pH6_6_E" = "6", "pH6_6_F" = "6",
               "pH6_7_A" = "7", "pH6_7_B" = "7", "pH6_7_C" = "7", "pH6_7_D" = "7", "pH6_7_F" = "7",
               "pH8_1_B" = "1", "pH8_1_C" = "1", "pH8_1_D" = "1", "pH8_1_E" = "1", "pH8_1_F" = "1",
               "pH8_2_A" = "2", "pH8_2_B" = "2", "pH8_2_C" = "2", "pH8_2_D" = "2", "pH8_2_E" = "2", "pH8_2_F" = "2",
               "pH8_3_A" = "3", "pH8_3_D" = "3", "pH8_3_E" = "3", "pH8_3_F" = "3",
               "pH8_4_B" = "4", "pH8_4_C" = "4", "pH8_4_E" = "4", "pH8_4_F" = "4", 
               "pH8_5_A" = "5", "pH8_5_B" = "5", "pH8_5_C" = "5", "pH8_5_D" = "5", "pH8_5_E" = "5", "pH8_5_F" = "5",
               "pH8_6_A" = "6", "pH8_6_B" = "6", "pH8_6_C" = "6", "pH8_6_D" = "6", "pH8_6_E" = "6", "pH8_6_F" = "6",
               "pH8_7_A" = "7", "pH8_7_B" = "7", "pH8_7_C" = "7", "pH8_7_D" = "7", "pH8_7_F" = "7")
               
## omg i also need to do mg labels i'm gonna lose my mind
mg_labels <- c("G05_1_A" = "1", "G05_1_B" = "1", "G05_1_C" = "1", "G05_1_D" = "1", "G05_1_E" = "1", "G05_1_F" = "1",
               "G05_2_A" = "2", "G05_2_B" = "2", "G05_2_C" = "2", "G05_2_E" = "2", "G05_2_F" = "2",
               "G05_3_A" = "3", "G05_3_B" = "3", "G05_3_C" = "3", "G05_3_D" = "3", "G05_3_E" = "3", "G05_3_F" = "3",
               "G05_4_A" = "4", "G05_4_B" = "4", "G05_4_D" = "4", "G05_4_E" = "4", "G05_4_F" = "4",
               "G05_5_A" = "5", "G05_5_B" = "5", "G05_5_C" = "5", "G05_5_D" = "5", "G05_5_E" = "5", "G05_5_F" = "5",
               "G05_6_A" = "6", "G05_6_B" = "6", "G05_6_C" = "6", "G05_6_D" = "6", "G05_6_E" = "6", "G05_6_F" = "6",
               "G05_7_A" = "7", "G05_7_B" = "7", "G05_7_C" = "7",
               "G10_1_A" = "1", "G10_1_B" = "1", "G10_1_C" = "1", "G10_1_D" = "1", "G10_1_E" = "1", "G10_1_F" = "1",
               "G10_2_A" = "2", "G10_2_B" = "2", "G10_2_C" = "2", "G10_2_D" = "2", "G10_2_E" = "2", "G10_2_F" = "2",
               "G10_3_A" = "3", "G10_3_B" = "3", "G10_3_C" = "3", "G10_3_D" = "3", "G10_3_E" = "3", "G10_3_F" = "3",
               "G10_4_A" = "4", "G10_4_B" = "4", "G10_4_C" = "4", "G10_4_D" = "4", "G10_4_E" = "4", "G10_4_F" = "4",
               "G10_5_A" = "5", "G10_5_B" = "5", "G10_5_C" = "5", "G10_5_D" = "5", "G10_5_E" = "5", "G10_5_F" = "5",
               "G10_6_A" = "6", "G10_6_B" = "6", "G10_6_C" = "6", "G10_6_D" = "6", "G10_6_E" = "6", "G10_6_F" = "6",
               "G10_7_A" = "7", "G10_7_B" = "7", "G10_7_C" = "7", "G10_7_E" = "7", "G10_7_F" = "7",
               "MES_1_A" = "1", "MES_1_C" = "1", "MES_1_D" = "1", "MES_1_E" = "1", "MES_1_F" = "1",
               "MES_2_A" = "2", "MES_2_B" = "2", "MES_2_C" = "2", "MES_2_D" = "2", "MES_2_E" = "2", "MES_2_F" = "2",
               "MES_3_A" = "3", "MES_3_B" = "3", "MES_3_C" = "3", "MES_3_D" = "3", "MES_3_E" = "3", "MES_3_F" = "3",
               "MES_4_A" = "4", "MES_4_B" = "4", "MES_4_C" = "4", "MES_4_D" = "4", "MES_4_E" = "4", "MES_4_F" = "4",
               "MES_5_A" = "5", "MES_5_B" = "5", "MES_5_C" = "5", "MES_5_D" = "5", "MES_5_E" = "5", "MES_5_F" = "5",
               "MES_6_A" = "6", "MES_6_B" = "6", "MES_6_C" = "6", "MES_6_D" = "6", "MES_6_E" = "6", "MES_6_F" = "6",
               "MES_7_A" = "7", "MES_7_B" = "7", "MES_7_C" = "7", "MES_7_D" = "7", "MES_7_E" = "7", "MES_7_F" = "7",
               "pH10_1_A" = "1", "pH10_1_B" = "1", "pH10_1_C" = "1", "pH10_1_D" = "1", "pH10_1_E" = "1", "pH10_1_F" = "1",
               "pH10_2_A" = "2", "pH10_2_B" = "2", "pH10_2_C" = "2", "pH10_2_D" = "2", "pH10_2_E" = "2", "pH10_2_F" = "2",
               "pH10_3_A" = "3", "pH10_3_B" = "3", "pH10_3_C" = "3", "pH10_3_D" = "3", "pH10_3_F" = "3",
               "pH10_4_A" = "4", "pH10_4_B" = "4", "pH10_4_C" = "4", "pH10_4_D" = "4", "pH10_4_E" = "4", "pH10_4_F" = "4",
               "pH10_5_F" = "5",
               "pH10_6_A" = "6", "pH10_6_B" = "6", "pH10_6_C" = "6", "pH10_6_D" = "6", "pH10_6_E" = "6", "pH10_6_F" = "6",
               "pH10_7_A" = "7", "pH10_7_B" = "7", "pH10_7_C" = "7", "pH10_7_D" = "7", "pH10_7_E" = "7", "pH10_7_F" = "7",
               "pH6_1_A" = "1", "pH6_1_B" = "1", "pH6_1_C" = "1", "pH6_1_D" = "1", "pH6_1_E" = "1", "pH6_1_F" = "1",
               "pH6_2_A" = "2", "pH6_2_B" = "2", "pH6_2_C" = "2", "pH6_2_D" = "2", "pH6_2_E" = "2", "pH6_2_F" = "2",
               "pH6_3_A" = "3", "pH6_3_B" = "3", "pH6_3_C" = "3", "pH6_3_D" = "3", "pH6_3_E" = "3", "pH6_3_F" = "3",
               "pH6_4_A" = "4", "pH6_4_B" = "4", "pH6_4_C" = "4", "pH6_4_D" = "4", "pH6_4_E" = "4", "pH6_4_F" = "4",
               "pH6_5_A" = "5", "pH6_5_B" = "5", "pH6_5_C" = "5", "pH6_5_D" = "5", "pH6_5_E" = "5", "pH6_5_F" = "5",
               "pH6_6_A" = "6", "pH6_6_B" = "6", "pH6_6_D" = "6", "pH6_6_E" = "6", "pH6_6_F" = "6",
               "pH6_7_A" = "7", "pH6_7_B" = "7", "pH6_7_C" = "7", "pH6_7_D" = "7", "pH6_7_F" = "7",
               "pH8_1_A" = "1", "pH8_1_B" = "1", "pH8_1_C" = "1", "pH8_1_D" = "1", "pH8_1_E" = "1", "pH8_1_F" = "1",
               "pH8_2_A" = "2", "pH8_2_B" = "2", "pH8_2_C" = "2", "pH8_2_D" = "2", "pH8_2_E" = "2", "pH8_2_F" = "2",
               "pH8_3_A" = "3", "pH8_3_B" = "3", "pH8_3_C" = "3", "pH8_3_D" = "3", "pH8_3_E" = "3", "pH8_3_F" = "3",
               "pH8_4_A" = "4", "pH8_4_B" = "4", "pH8_4_C" = "4", "pH8_4_D" = "4", "pH8_4_E" = "4", "pH8_4_F" = "4",
               "pH8_5_A" = "5", "pH8_5_B" = "5", "pH8_5_C" = "5", "pH8_5_D" = "5", "pH8_5_E" = "5", "pH8_5_F" = "5",
               "pH8_6_A" = "6", "pH8_6_B" = "6", "pH8_6_C" = "6", "pH8_6_D" = "6", "pH8_6_E" = "6", "pH8_6_F" = "6",
               "pH8_7_A" = "7", "pH8_7_B" = "7", "pH8_7_C" = "7", "pH8_7_D" = "7", "pH8_7_F" = "1")

#### summarize data for boxplot


data_trim <- cazy_total %>%
  group_by(Treatment,Day, Replicate) %>%
  summarise(percent=sum(percent))
#write.csv(data_trim, "cazy_expression_summary.csv", row.names = FALSE)


####### BAR PLOTS OF RELATIVE CAZY CONTRIBUTION BY ORG ########

cazy_trim <- cazy_total[,c(8,15,24)]
cazy_trim <-subset(cazy_trim, Substrate != "Other")
cazy_trim <- cazy_trim %>%
  group_by(Substrate, order) %>%
  summarise(value=sum(value))
sums <- cazy_trim %>%
  group_by(Substrate) %>%
  summarise(sum=sum(value))

cazy_trim <- merge(cazy_trim, sums, by = "Substrate")

cazy_trim$percent <- ((cazy_trim$value)/(cazy_trim$sum))*100

sums$sum <- format(sums$sum, big.mark = ",", scientific=FALSE)
sums$percent <- 102

cazy_trim$Substrate <- factor(cazy_trim$Substrate, levels=c("Cellulose", "Cellobiose", "Xylan",
                                                            "Starch", "Trehalose","Pectin","Lignin",
                                                            "Peptidoglycan","Chitin","Other"))


pcountsbar <- ggplot() + 
  geom_bar(data=cazy_trim, aes(x=Substrate, y = percent,fill = order, color=order), stat='identity') +
  geom_text(data=sums, aes(x=Substrate, y= percent, label=sum)) +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(angle= 25, hjust=1),
        aspect.ratio = 1) +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  labs (x = "CAZyme Grouping", y = "Percent Expression \n (Summed Across Treatment and Time)")+
  guides(color="none", fill=guide_legend("Order")) +
  ggtitle("Relative Contribution of Bacterial Orders to CAZyme Grouping")
pcountsbar




cazy_trim2 <- cazy_total[,c(8,15,20,24)]
cazy_trim2 <-subset(cazy_trim2, Substrate != "Other")
cazy_trim2 <- cazy_trim2 %>%
  group_by(Substrate, order, Treatment) %>%
  summarise(value=sum(value))
sums2 <- cazy_trim2 %>%
  group_by(Substrate, Treatment) %>%
  summarise(sum=sum(value))
cazy_trim2 <- merge(cazy_trim2, sums2, by = c("Substrate"="Substrate", "Treatment"="Treatment"))
cazy_trim2$percent <- ((cazy_trim2$value)/(cazy_trim2$sum))*100
sums2$sum <- format(sums2$sum, big.mark = ",", scientific=FALSE)
sums2$percent <- 103

cazy_trim2$Substrate <- factor(cazy_trim2$Substrate, levels=c("Cellulose", "Cellobiose", "Xylan",
                                                              "Starch", "Trehalose","Pectin","Lignin",
                                                              "Peptidoglycan","Chitin","Other"))


pcountsbar_treatment <- ggplot() + 
  geom_bar(data=cazy_trim2, aes(x=Substrate, y = percent,fill = order, color=order), stat='identity') +
  #geom_text(data=sums2, aes(x=Substrate, y= percent, label=sum), position=position_dodge(width=0.9), size =3) +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(angle= 25, hjust=1, size = 12),
        #aspect.ratio = 1,
        legend.position = "bottom") +
  scale_fill_manual(values=palorg) +
  scale_color_manual(values=palorg) +
  labs (x = "CAZyme Grouping", y = "Percent Expression \n (Summed Across Time)")+
  guides(color="none", fill=guide_legend("Order")) +
  ggtitle("Relative Contribution of Bacterial Orders to CAZyme Grouping") +
  facet_wrap(.~Treatment, scales ="free")
pcountsbar_treatment

ggsave("pcountsbar.tiff", device = "tiff", dpi = 700)
ggsave("pcountsbar.png", device = "png", dpi = 700)
ggsave("pcountsbar.pdf", device = "pdf", dpi = 700)


###############REL ABUND+EXPRESSOF CAZY t #############

data_sub_time <- cazy_total[,c(7,10,20,24)]

data_sub_time <- subset(data_sub_time, Substrate != "Other")

data_sub_time <- data_sub_time %>%
  group_by(variable, Treatment, Substrate) %>%
  summarise(percent=sum(percent))


p_sub_time <- ggplot(data=data_sub_time, aes(x=variable, y= percent, color = Substrate, fill = Substrate)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = palsubstrate) +
  scale_fill_manual(values = palsubstrate) +
  labs (x = "Samples, Arranged by Day", y = "Relative Expression (%)")+
  guides(color="none", fill=guide_legend("CAZy Grouping")) +
  scale_x_discrete(labels = mt_labels) +
  facet_grid( ~ Treatment, scales = "free")
p_sub_time  

############ MG #############

data_sub_time_mg <- cazy_total_mg[,c(4,7,18,19)]

data_sub_time_mg <- subset(data_sub_time_mg, Substrate != "Other")

data_sub_time_mg <- data_sub_time_mg %>%
  group_by(variable, Treatment, Substrate) %>%
  summarise(percent=sum(percent))


p_sub_time_mg <- ggplot(data=data_sub_time_mg, aes(x=variable, y= percent, color = Substrate, fill = Substrate)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = palsubstrate) +
  scale_fill_manual(values = palsubstrate) +
  labs (x = "Samples, Arranged by Day", y = "Relative Abundnace (%)")+
  guides(color="none", fill=guide_legend("CAZy Grouping")) +
  scale_x_discrete(labels = mg_labels) +
  facet_grid( ~ Treatment, scales = "free")
p_sub_time_mg  

ggarrange(p_sub_time_mg , p_sub_time, ncol=1, common.legend = TRUE, legend = "bottom")

ggsave("cr_dbcan_relabun_express.tiff", device = "tiff", dpi = 700)
ggsave("cr_dbcan_relabun_express.png", device = "png", dpi = 700)
ggsave("cr_dbcan_relabun_express.pdf", device = "pdf", dpi = 700)

############################# Xylan ##################################

palce <- c("#89B151","#FFEA59","#B7DBDB","#D9B253","#6CA18F","#CE9486","#335B8E","#C9573C","#999999")
palgh <- c("#E69EA1","#DED9A2","#F6C83E","#4D5B28","#DB4372","#B77E60","#5785C1","#8CBEA3", "#a6cee3", "#999999")

##focus on xylan
cazy_trim4 <- subset(cazy_total, Substrate == "Xylan")
cazy_trim4$CAZy <- factor(cazy_trim4$CAZy, levels=c("CE1",  "CE2",  "CE3",  "CE4",  "CE5",  "CE6",  "CE7","CE12", "CE16", "GH10","GH5_22","GH30","GH30_2","GH30_3","GH30_7","GH30_8","GH67","GH115","GH120"))

cazy_trim4 <- cazy_trim4[,c(1,2,7:10,15,20:24)]

#separate into CE and GH
cazy_trim4$cazy_type <- substr(cazy_trim4$CAZy, 1,2)

cazy_trim5 <- cazy_trim4 %>%
  group_by(cazy_type, order, variable, Treatment, Day, Replicate) %>%
  summarise(value=sum(percent))

pxylan_counts <- ggplot() +
  geom_bar(data=cazy_trim5, aes(x=variable, y = value, color = order, fill = order), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg) +
  scale_x_discrete(labels = mt_labels) +
  labs(x = "Samples, Arranged by Day", y = "Relative Expression (%)") +
  guides(color = "none", fill = guide_legend("Order")) +
  facet_grid(cazy_type~Treatment, scales ="free")

pxylan_counts

#ggsave("xylan_expression_byorg.png")
#ggsave("xylan_expression_byorg.pdf")

### xylan ce

cazy_trim6 <- subset(cazy_trim4, cazy_type == "CE")

cazy_trim6 <- cazy_trim6 %>%
  group_by(variable, Day, CAZy, order, Treatment) %>%
  summarize(percent =sum(percent),
            value = sum(value),
            mean=mean(value),
            sd=sd(value))



pce_counts <- ggplot()+
  geom_bar(data=cazy_trim6, aes(x=variable, y = percent, color = CAZy, fill = CAZy), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_x_discrete(labels = mt_labels) +
  labs (x = "Samples, Arranged by Day", y = "Relative Expression (%)")+
  scale_color_manual(values=palce) +
  scale_fill_manual(values=palce)+
  ggtitle("Expression of Xylan CE CAZys") +
  facet_grid(~Treatment, scales = "free_x")

pce_counts



pce_org <- ggplot() +
  geom_bar(data=cazy_trim6, aes(x=CAZy, y = value, color = order, fill = order), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "bottom") +
  labs (x = "CAZyme Family", y = "Number of Reads")+
  guides(color = "none", fill = guide_legend("Order")) +
  ggtitle("Count of Xylan CE CAZy Transcripts") +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg)
pce_org


## xylan gh

cazy_trim7 <- subset(cazy_trim4, cazy_type == "GH")
cazy_trim7 <- cazy_trim7 %>%
  group_by(CAZy, variable, order, Treatment) %>%
  summarize(percent = sum(percent),
            value =sum(value),
            mean=mean(value),
            sd=sd(value))


cazy_trim7$CAZy <- factor(cazy_trim7$CAZy, levels=c( "GH10","GH5_22","GH30","GH30_2","GH30_3","GH30_7","GH30_8","GH67","GH115","GH120" ))


pgh_counts <- ggplot()+
  geom_bar(data=cazy_trim7, aes(x=variable, y = percent, color = CAZy, fill = CAZy), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_x_discrete(labels = mt_labels) +
  labs(x = "Samples, Arranged by Day", y = "Relative Expression (%)")+
  scale_color_manual(values=palgh) +
  scale_fill_manual(values=palgh) +
  ggtitle("Expression of Xylan GH CAZys") +
  facet_grid(~Treatment, scales = "free")

pgh_counts



pgh_org <- ggplot() +
  geom_bar(data=cazy_trim7, aes(x=CAZy, y = value, color = order, fill = order), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs (x = "CAZyme Family", y = "Number of Reads")+
  guides(color = "none", fill = guide_legend("Order")) +
  ggtitle("Count of Xylan GH CAZy Transcripts") +
  scale_color_manual(values=palorggh) +
  scale_fill_manual(values=palorggh)
pgh_org





p1 <- ggarrange(pce_counts, pgh_counts, ncol = 1, labels = c("A","C"), font.label = list(size = 25))
p1
p2 <- ggarrange(pce_org, pgh_org, ncol = 1, common.legend = TRUE, legend = "right",
                labels = c("B","D"), font.label = list(size = 25))
p2


xylan_total <- plot_grid(p1,p2, align = "v", ncol = 2, axis = "tb", rel_widths = c(20, 10))
xylan_total

ggsave("xylan_expression.tiff", device = "tiff", dpi = 700)
ggsave("xylan_expression.png", device = "png", dpi = 700)
ggsave("xylan_expression.pdf", device = "pdf", dpi = 700)


### stats, skip if graphing 
cazy_trim7 <- subset(cazy_trim4, cazy_type == "GH")
cazy_trim7 <- cazy_trim7 %>%
  group_by(CAZy, order) %>%
  summarize(value =sum(value),
            mean=mean(value),
            sd=sd(value))
cazy_trim7 <- subset(cazy_trim7, CAZy == "GH10")
cazy_trim8 <- cazy_trim7 %>%
  summarise(value=sum(value))
cazy_trim7$percent <- (cazy_trim7$value)/(108502) * 100

############################# Lignin ##################################

cazy_trim3 <- subset(cazy_total, Substrate == "Lignin")
length(unique(cazy_trim3$CAZy))

cazy_trim3$CAZy <- factor(cazy_trim3$CAZy, levels=c("AA1","AA1_1","AA1_2","AA1_3","AA2","AA3","AA3_1","AA3_2","AA3_3","AA3_4","AA4","AA5","AA5_1","AA5_2","AA6","AA7","AA12"))

cazy_trim3 <- cazy_trim3 %>%
  group_by(variable, CAZy, order, Treatment) %>%
  summarize(percent=sum(percent),
            value =sum(value),
            mean=mean(value),
            sd=sd(value))

pallignin1 <- c("#fa9fb5","#f768a1","#dd3497","#ae017e", #AA1, pink
                "#335B8E", #AA2, blue
                #### no AA3 ###
                "#984ea3", #AA4, purple
                "#fdae6b","#f16913","#8c2d04", #AA5, orange
                "#DFBA47", #AA6, yellow
                "#6C8645", #AA7, green
                "#e41a1c" #AA12, red
)

pallignin2 <- c("#DFBA47","#175149","#8CBEA3","#CB654F","#F4B5BD")


cazy_trim8 <- filter(cazy_trim3, !CAZy %in% c("AA3","AA3_1","AA3_2","AA3_3","AA3_4"))
cazy_trim9 <- filter(cazy_trim3, CAZy %in% c("AA3","AA3_1","AA3_2","AA3_3","AA3_4"))


plignin_counts1 <- ggplot()+
  geom_bar(data=cazy_trim8, aes(x=variable, y = percent, color = CAZy, fill = CAZy), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_x_discrete(labels = mt_labels) +
  labs (x = "Samples, Arranged by Day", y = "Relative Expression (%)")+
  scale_color_manual(values=pallignin1) +
  scale_fill_manual(values=pallignin1)+
  ggtitle("Expression of Lignin CAZys") +
  facet_grid(~Treatment, scales = "free")

plignin_counts1


pligni_org1 <- ggplot() +
  geom_bar(data=cazy_trim8, aes(x=CAZy, y = value, color = order, fill = order), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs (x = "CAZyme Family", y = "Number of Reads")+
  guides(color = "none", fill = guide_legend("Order")) +
  ggtitle("Count of Lignin CAZy Transcripts") +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg)
pligni_org1




plignin_counts2 <- ggplot()+
  geom_bar(data=cazy_trim9, aes(x=variable, y = percent, color = CAZy, fill = CAZy), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_x_discrete(labels = mt_labels) +
  labs (x = "Samples, Arranged by Day", y = "Relative Expression (%)")+
  scale_color_manual(values=pallignin2) +
  scale_fill_manual(values=pallignin2)+
  ggtitle("Expression of Lignin AA3 GH CAZys") +
  facet_grid(~Treatment, scales = "free")

plignin_counts2

pligni_org2 <- ggplot() +
  geom_bar(data=cazy_trim9, aes(x=CAZy, y = value, color = order, fill = order), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs (x = "CAZyme Family", y = "Number of Reads")+
  guides(color = "none", fill = guide_legend("Order")) +
  ggtitle("Count of Lignin AA3 CAZy Transcripts") +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg)
pligni_org2



p1 <- ggarrange(plignin_counts1,plignin_counts2, ncol = 1, labels = c("A","C"), font.label = list(size = 25))
p1
p2 <- ggarrange(pligni_org1, pligni_org2, ncol = 1, common.legend = TRUE, legend = "right",
                labels = c("B","D"), font.label = list(size = 25))
p2


lignin_total <- plot_grid(p1,p2, align = "h", ncol = 2, axis = "tb", rel_widths = c(20, 10))
lignin_total


ggsave("lignin_expression.tiff", device = "tiff", dpi = 700)
ggsave("lignin_expression.png", device = "png", dpi = 700)
ggsave("lignin_expression.pdf", device = "pdf", dpi = 700)




######### Cellulose I guess #######################


cazy_trim10 <- subset(cazy_total, Substrate == "Cellulose")
cazy_trim11 <- subset(cazy_total, Substrate == "Cellobiose")

length(unique(cazy_trim10$bin))
length(unique(cazy_trim11$bin))


unique(cazy_trim10$CAZy)
unique(cazy_trim11$CAZy)

cazy_trim10$CAZy <- factor(cazy_trim10$CAZy, levels=c("GH5","GH5_1","GH5_2","GH5_4","GH5_5","GH5_25","GH5_26","GH5_38","GH5_39","GH5_46","GH6","GH8","GH9","GH44","GH74","AA10"))

cazy_trim11$CAZy <- factor(cazy_trim11$CAZy, levels=c("GH1","GH3","GH116" ))


cazy_trim10 <- cazy_trim10 %>%
  group_by(variable, CAZy, order, Treatment) %>%
  summarize(percent = sum(percent),
            value =sum(value),
            mean=mean(value),
            sd=sd(value))

cazy_trim11 <- cazy_trim11 %>%
  group_by(variable, CAZy, order, Treatment) %>%
  summarize(percent = sum(percent),
            value =sum(value),
            mean=mean(value),
            sd=sd(value))

palcellu <- c("#B57EDC", "#A1CAF1","#6495ED","#318CE7","#007FFF","#1560BD","#0047AB","#08457E","#003366","#002147",
              "#AF4E24","#E54E21","#FDD262","#6C8645","#54D8B1","#9C964A")

palorgcellu <- c("#8DD3C7", #Actinomycetales
                 "#FFFFB3", #Burkholderiales
                 "#BEBADA", #Caulobacterales
                 "#80B1D3", #Enterobacterales
                 "#FDB462", #Flavobacteriales
                 "#B3DE69", #Lactobacillales
                 "#FCCDE5", #Micrococcales
                 "#BC80BD", #Propionibacteriales
                 "#e5c494", #Rhizobiales
                 "#1f78b4", #Sphingobacteriales
                 "#fb9a99", #Xanthomonadales
                 "#D9D9D9" #no support
)

pcellulose_counts <- ggplot()+
  geom_bar(data=cazy_trim10, aes(x=variable, y = percent, color = CAZy, fill = CAZy), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_x_discrete(labels = mt_labels) +
  labs (x = "Samples, Arranged by Day", y = "Relative Expression (%)")+
  ggtitle("Expression of Cellulose CAZys") +
  scale_color_manual(values=palcellu) +
  scale_fill_manual(values=palcellu)+
  facet_grid(~Treatment, scales = "free")

pcellulose_counts

length(unique(cazy_trim10$order))


pcellulose_org <- ggplot() +
  geom_bar(data=cazy_trim10, aes(x=CAZy, y = value, color = order, fill = order), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "non") +
  labs (x = "CAZyme Family", y = "Number of Reads")+
  guides(color = "none", fill = guide_legend("Order")) +
  ggtitle("Count of Cellulose CAZy Transcripts") +
  scale_color_manual(values=palorgcellu) +
  scale_fill_manual(values=palorgcellu)
pcellulose_org


pcellobio <- c("#3B9AB2","#EBCC2A","#F21A00")
pcellobiose_counts <- ggplot()+
  geom_bar(data=cazy_trim11, aes(x=variable, y = percent, color = CAZy, fill = CAZy), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_x_discrete(labels = mt_labels) +
  labs (x = "Samples, Arranged by Day", y = "Relative Expression (%)")+
  ggtitle("Expression of Cellobiose GH CAZys") +
  scale_color_manual(values=pcellobio) +
  scale_fill_manual(values=pcellobio)+
  facet_grid(~Treatment, scales = "free")

pcellobiose_counts

pcellobiose_org <- ggplot() +
  geom_bar(data=cazy_trim11, aes(x=CAZy, y = value, color = order, fill = order), stat='identity') +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs (x = "CAZyme Family", y = "Number of Reads")+
  guides(color = "none", fill = guide_legend("Order")) +
  ggtitle("Count of Cellobiose Transcripts") +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg)
pcellobiose_org

p1 <- ggarrange(pcellulose_counts,pcellobiose_counts, ncol = 1, labels = c("A","C"), font.label = list(size = 25))
p1
p2 <- ggarrange(pcellulose_org, pcellobiose_org, ncol = 1, common.legend = TRUE, legend = "right",
                labels = c("B","D"), font.label = list(size = 25))
p2


cellulose_total <- plot_grid(p1,p2, align = "h", ncol = 2, axis = "tb", rel_widths = c(20, 10))
cellulose_total

ggsave("cellulose_expression.tiff", device = "tiff", dpi = 700)
ggsave("cellulose_expression.png", device = "png", dpi = 700)
ggsave("cellulose_expression.pdf", device = "pdf", dpi = 700)


##########



#write matrix for diversity plots
total_table <- cazy_total

Substrate_matrix <-dcast(data=total_table, formula = CAZy~variable, fun.aggregate = sum, value.var="value")
tax_matrix <- dcast(data=total_table, formula = bin~variable, fun.aggregate = sum, value.var="value")


Substrate_matrix_mg <-dcast(data=cazy_total_mg, formula = CAZy~variable, fun.aggregate = sum, value.var="value")
#write.csv(tax_matrix, "CR_dbCAN_taxmatrix.csv", row.names = FALSE)
#write.csv(Substrate_matrix, "CR_dbCAN_Substratematrix_cazy.csv", row.names = FALSE)
#write.csv(Substrate_matrix_mg, "CR_dbCAN_funcmatrix_cazy_mg.csv", row.names = FALSE)


################# OLD, Ignore ##################################


################ #focus on lacto, corne, sphingo in G10 #######################

cazy_mg <- subset(cazy_total_mg, Treatment == "G10" & order %in% c("Corynebacteriales", "Lactobacillales","Sphingobacteriales")& Substrate != "Other")

cazy_mg <- cazy_mg %>%
  group_by(variable, Day, order, Substrate) %>%
  summarise(sum=sum(percent))



cazy_g10 <- subset(cazy_total, Treatment == "G10" & order %in% c("Corynebacteriales", "Lactobacillales","Sphingobacteriales")& Substrate != "Other")
cazy_g10 <- cazy_g10 %>%
  group_by(variable, Day, order, Substrate) %>%
  summarise(sum=sum(percent))

unique(cazy_g10$Substrate)

pg10 <- ggplot(data=cazy_g10 , aes(x=variable, y = sum, color = Substrate, fill = Substrate)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=20) +
  #ylim(0,1.1) +
  theme(axis.text.x = element_blank()) +
  scale_color_manual(values = palsubstrate)+
  scale_fill_manual(values = palsubstrate) +
  labs(y = "Percent Transcripts Mapped", x = "Samples, Arranged by Day")+
  facet_grid(~order, scales = "free") +
  ggtitle("Relative Expression of CAZymes in G10 Treatment")
pg10


pg10_mg <- ggplot(data=cazy_mg , aes(x=variable, y = sum, color = Substrate, fill = Substrate)) +
  geom_bar(stat="identity") +
  theme_classic(base_size=20) +
  #ylim(0,1.1) +
  theme(axis.text.x = element_blank()) +
  scale_color_manual(values = palsubstrate)+
  scale_fill_manual(values = palsubstrate) +
  labs(y = "Percent Reads Mapped", x = "Samples, Arranged by Day")+
  facet_grid(~order, scales = "free") +
  ggtitle("Relative Abundance of CAZymes in G10 Treatment")
pg10_mg


cazy_g10 <- ggarrange(pg10_mg, pg10, common.legend=TRUE, ncol=1, nrow=2, legend = "right")
cazy_g10
annotate_figure(cazy_g10, top=text_grob("Abundance and Expression of CAZy Annotations, Treatment G10", size = 20))


############### BOXPLOT OF MG VS MY XYLAN DYNAMICS ###################

palorgcorynshingo <- c("#FB8072","#B3DE69","#1f78b4")
palorgcorynshingo2 <- c("#B3DE69","#1f78b4")

xyl_ph6 <- subset(cazy_total, Treatment == "pH6" & Substrate == "Xylan")
xyl_ph6$cazy_type <- substr(xyl_ph6$CAZy, 1,2)


#group xylan reads
xyl_ph6 <- xyl_ph6 %>%
  group_by(Day, order, variable, Substrate, cazy_type, Replicate) %>%
  summarise(value=sum(value))

#sum up xylan reads per sample
xyl_ph6_sum <- xyl_ph6 %>%
  group_by(variable) %>%
  summarise(sum=sum(value))

#get percents
xyl_ph6 <- merge(xyl_ph6, xyl_ph6_sum, by = "variable")
xyl_ph6$percent <- (xyl_ph6$value/xyl_ph6$sum)*100


xyl_ph6_summary <- xyl_ph6 %>%
  group_by(Day, order, Substrate, cazy_type) %>%
  summarise(mean=mean(percent),
            sd=sd(percent),
            count=n())

xyl_ph6$Type <- "MT"
xyl_ph6 <- subset(xyl_ph6 , order %in% c("Corynebacteriales", "Lactobacillales", "Sphingobacteriales"))


#work on mg

xyl_ph6_mg <- subset(cazy_total_mg, Treatment == "pH6" & Substrate == "Xylan")
xyl_ph6_mg$cazy_type <- substr(xyl_ph6_mg$CAZy, 1,2)

#group xylan reads
xyl_ph6_mg <- xyl_ph6_mg %>%
  group_by(Day, order, variable, Substrate, Replicate, cazy_type) %>%
  summarise(value=sum(value))

#sum up xylan reads per sample
xyl_ph6_mg_sum <- xyl_ph6_mg %>%
  group_by(variable) %>%
  summarise(sum=sum(value))

#get percents
xyl_ph6_mg <- merge(xyl_ph6_mg, xyl_ph6_mg_sum, by = "variable")
xyl_ph6_mg$percent <- (xyl_ph6_mg$value/xyl_ph6_mg$sum)*100

xyl_ph6_mg_summary <- xyl_ph6_mg %>%
  group_by(Day, order, Substrate, cazy_type) %>%
  summarise(mean=mean(percent),
            sd=sd(percent),
            count=n())

xyl_ph6_mg$Type <- "MG"
xyl_ph6_mg <- subset(xyl_ph6_mg , order %in% c("Corynebacteriales", "Lactobacillales", "Sphingobacteriales"))

## merge xyl_ph6 data

xyl_ph6_total <- rbind(xyl_ph6, xyl_ph6_mg)
xyl_ph6_total$Day <- as.factor(xyl_ph6_total$Day)


xyl_ph6_cory <- subset(xyl_ph6_total, order == "Corynebacteriales")
xyl_ph6_sphing <- subset(xyl_ph6_total, order == "Sphingobacteriales")
xyl_ph6_lacto <- subset(xyl_ph6_total, order == "Lactobacillales")


xyl_ph6_sphing_summary <- xyl_ph6_sphing %>%
  group_by(Day, cazy_type, Type) %>%
  summarise(mean = mean(percent),
            sd= sd(percent))

xyl_ph6_lacto_summary <- xyl_ph6_lacto %>%
  group_by(Day, cazy_type, Type) %>%
  summarise(mean = mean(percent),
            sd= sd(percent))





cazy.labels <- c("Carbodyrate Esterases", "Glycoside Hydrolases")
names(cazy.labels) <- c("CE", "GH")


pph6xyl_cory <- ggplot(xyl_ph6_cory, aes(x = Day, y = percent, fill = Type)) +
  geom_boxplot() +
  theme_classic(base_size=20) +
  labs(y = "Relative Amount of Xylan CAZys (%)", x = "Day") +
  stat_compare_means(aes(group=Type), label = "p.signif", method = "t.test") +
  scale_fill_manual(values = c("#1E88E5", "#FFC107")) +
  guides(fill=guide_legend("Dataset")) +
  ggtitle("Xylan CAZys Mapped to Sphingobacteriales MAGs") +
  facet_grid(~cazy_type, labeller = labeller(cazy_type = cazy.labels))
pph6xyl_cory

pph6xyl_sph <- ggplot(xyl_ph6_sphing, aes(x = Day, y = percent, fill = Type)) +
  geom_boxplot() +
  theme_classic(base_size=20) +
  labs(y = "Relative Amount of Xylan CAZys (%)", x = "Day") +
  stat_compare_means(aes(group=Type), label = "p.signif", method = "t.test") +
  scale_fill_manual(values = c("#1E88E5", "#FFC107")) +
  guides(fill=guide_legend("Dataset")) +
  ggtitle("pH6: Xylan CAZys Mapped to Sphingobacteriales MAGs") +
  facet_grid(~cazy_type, labeller = labeller(cazy_type = cazy.labels))
pph6xyl_sph



pph6xyl_lacto <- ggplot(xyl_ph6_lacto, aes(x = Day, y = percent, fill = Type)) +
  geom_boxplot() +
  theme_classic(base_size=20) +
  labs(y = "Relative Amount of Xylan CAZys (%)", x = "Day") +
  stat_compare_means(aes(group=Type), label = "p.signif", method = "t.test") +
  scale_fill_manual(values = c("#1E88E5", "#FFC107")) +
  guides(fill=guide_legend("Dataset")) +
  ggtitle("pH6: Xylan CAZys Mapped to Lactobacillales MAGs") +
  facet_grid(~cazy_type, labeller = labeller(cazy_type = cazy.labels))
pph6xyl_lacto


ggarrange(pph6xyl_sph, pph6xyl_lacto, nrow= 2, common.legend = TRUE, legend = "bottom")

#ggsave("ph6_xylan.png")
#ggsave("ph6_xylan.pdf")

####### xylan expression bar by time

palce <- c("#89B151","#FFEA59","#B7DBDB","#D9B253","#6CA18F","#CE9486","#335B8E","#C9573C","#999999")
palgh <- c("#E69EA1","#DED9A2","#F6C83E","#4D5B28","#DB4372","#B77E60","#5785C1","#8CBEA3", "#a6cee3", "#999999")
palorggh <-c("#FFFFB3","#BEBADA","#80B1D3","#B3DE69","#FCCDE5","#BC80BD","#e5c494","#1f78b4","#D9D9D9")

d_xyl<- subset(cazy_total, Substrate == "Xylan")
d_xyl$cazy_type <- substr(d_xyl$CAZy, 1,2)
length(unique(d_xyl$CAZy))



#group xylan reads
d_xyl <- d_xyl %>%
  group_by(order, variable, CAZy, cazy_type, Treatment, Day, Replicate) %>%
  summarise(value=sum(value))

#sum up xylan reads per sample
d_xyl_sum <- d_xyl %>%
  group_by(variable) %>%
  summarise(sum=sum(value))

#get percents
d_xyl <- merge(d_xyl, d_xyl_sum, by = "variable")
d_xyl$percent <- (d_xyl$value/d_xyl$sum)*100


d_xyl_summary <- d_xyl %>%
  group_by(CAZy, cazy_type) %>%
  summarise(mean=mean(percent),
            sd=sd(percent),
            count=n())


d_xyl <- subset(d_xyl, Treatment == "G10")
d_xyl_ce <- subset(d_xyl, cazy_type == "CE")
d_xyl_gh <- subset(d_xyl, cazy_type == "GH")
d_xyl_gh$CAZy <- factor(d_xyl_gh$CAZy, levels=c( "GH10","GH5_22","GH30","GH30_2","GH30_3","GH30_7","GH30_8","GH67","GH115","GH120" ))

p_xyl_ce <- ggplot(d_xyl_ce, aes(x=variable, y = value, fill = CAZy, color = CAZy)) +
  geom_bar(stat="identity") +
  theme_classic(base_size= 15) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = palce) +
  scale_color_manual(values = palce) +
  labs (x = "Samples, Arranged by Day", y = "Number of Reads")+
  facet_grid(~Treatment, scales = "free")

p_xyl_ce

p_xyl_ce_org <- ggplot(d_xyl_ce, aes(x=variable, y = value, fill = order, color = order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size= 15) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = palorg) +
  scale_color_manual(values = palorg) +
  labs (x = "Samples, Arranged by Day", y = "Number of Reads")+
  guides(color="none", fill=guide_legend("Order")) +
  facet_grid(~Treatment, scales = "free")

p_xyl_ce_org


p_xyl_gh <- ggplot(d_xyl_gh, aes(x=variable, y = value, fill = CAZy, color = CAZy)) +
  geom_bar(stat="identity") +
  theme_classic(base_size= 15) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = palgh) +
  scale_color_manual(values = palgh) +
  labs (x = "Samples, Arranged by Day", y = "Number of Reads")+
  facet_grid(~Treatment, scales = "free")

p_xyl_gh

p_xyl_gh_org <- ggplot(d_xyl_gh, aes(x=variable, y = value, fill = order, color = order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size= 15) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = palorggh) +
  scale_color_manual(values = palorggh) +
  labs (x = "Samples, Arranged by Day", y = "Number of Reads")+
  guides(color="none", fill=guide_legend("Order")) +
  facet_grid(~Treatment, scales = "free")

p_xyl_gh_org

p_xyl_total <- ggarrange(p_xyl_ce, p_xyl_gh,
                         p_xyl_ce_org, p_xyl_gh_org, ncol = 2, nrow = 2,
                         labels = c("A","B","C","D"))
#p_xyl_total <- annotate(p_xyl_total, top=text_grob("Expression of Xylan-acting CAZys", size = 20))

p_xyl_total

