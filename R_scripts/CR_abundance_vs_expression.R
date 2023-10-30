library(ggplot2)
library(dplyr)
library(ggpubr)

#read in mg, read in mt
#make sure these are unfiltered for ribosomes
#pull just percent
#plot

mg <- read.csv("CR_org_total_table_mg.csv", header = TRUE)
mt <- read.csv("CR_org_total_table_mt_unfiltered.csv", header = TRUE)
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


#summarise before messing around with tables 

mg_summary <- mg %>%
  group_by(variable) %>%
  summarise(sum=sum(percent))
mg_summary <- mg_summary %>%
  summarise(mean=mean(sum),
            sd=sd(sum))

mg_summary_tax <- mg %>%
  group_by(order,Treatment, Day, Replicate) %>%
  summarise(percent=sum(percent))
mg_summary_tax <- mg_summary_tax %>%
  group_by(order, Treatment, Day) %>%
  summarise(mean = mean(percent),
            sd=sd(percent))


mt_summary <- mt %>%
  group_by(variable) %>%
  summarise(sum=sum(percent))
mt_summary <- mt_summary %>%
  summarise(mean=mean(sum),
            sd=sd(sum))

mt_summary_tax <- mt %>%
  group_by(order,Treatment, Day, Replicate) %>%
  summarise(percent=sum(percent))
mt_summary_tax <- mt_summary_tax %>%
  group_by(order,Treatment, Day) %>%
  summarise(mean = mean(percent),
            sd=sd(percent),
            count = n())

##########

#######

mt$Treatment<-factor(mt$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
mt$order <- factor(mt$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                          "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                          "Sphingobacteriales","Xanthomonadales","no support"))


mg$Treatment<-factor(mg$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
mg$order <- factor(mg$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                      "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                      "Sphingobacteriales","Xanthomonadales","no support"))
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


##bar graphs

pbar_mt<- ggplot(mt, aes(x=variable, y = percent, color = order, fill = order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(size = 10)) +
  scale_color_manual(values=palorg, guide = "none") +
  scale_fill_manual(values=palorg) +
  labs(y = "Percent Reads Mapped", x = "Samples Arranged by Time", fill = "Order") +
  ggtitle("Relative Expression") +
  scale_x_discrete(labels = mt_labels) +
  facet_grid (~Treatment, scales = "free")
pbar_mt

pbar_mg<- ggplot(mg, aes(x=variable, y = percent, color = order, fill = order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(size = 10)) +
  scale_color_manual(values=palorg, guide = "none") +
  scale_fill_manual(values=palorg) +
  labs(y = "Percent Reads Mapped", x = "Samples Arranged by Time", fill = "Order") +
  ggtitle("Relative Abundance") +
  scale_x_discrete(labels = mg_labels) +
  facet_grid (~Treatment, scales = "free")
pbar_mg

ggarrange(pbar_mg, pbar_mt, nrow=2,  labels = "AUTO", font.label = list(size=25),
          common.legend = TRUE, legend = "bottom")


ggsave("cr_relabun_express.tiff", device = "tiff", dpi = 700)
ggsave("cr_relabun_express.png", device = "png", dpi = 700)
ggsave("cr_relabun_express.pdf", device = "pdf", dpi = 700)


##expression vs abundance
mg$bin_variable <- paste(mg$bin, mg$variable, sep = "_")
mt$bin_variable <- paste(mt$bin, mt$variable, sep = "_")


mg <- mg[,c(1,2,5,10,15:18)]
mt <- mt[,c(5,18)]


colnames(mg)[3] <- "percent_mg"
colnames(mt)[1] <- "percent_mt"

data <- merge(mt, mg, by = "bin_variable")

data$Treatment<-factor(data$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
data$order <- factor(data$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                          "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                          "Sphingobacteriales","Xanthomonadales","no support"))

plot_facet <- ggplot() +
  geom_point(data=data, aes(x=percent_mg, y = percent_mt, color = order), size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept=0) +
  #coord_fixed() +
  theme_classic(base_size=15) +
  theme(legend.position="bottom") +
  scale_color_manual(values=palorg) +
  labs(y="Relative Expression (%)", x="Relative Abundance (%)", color = "Order") +
  ggtitle("Relative Expression vs. Abundance of MAGs") +
  facet_wrap(~Treatment, ncol = 3)
plot_facet

ggsave("cr_abun_vs_express.tiff", device = "tiff", dpi = 700)
ggsave("cr_abun_vs_express.png", device = "png", dpi = 700)
ggsave("cr_abun_vs_express.pdf", device = "pdf", dpi = 700)

##### numbers
mt <- read.csv("CR_org_total_table_mt_unfiltered.csv", header = TRUE)
mt_summary <- mt %>%
  group_by(Treatment, Day, Replicate, order) %>%
  summarise(value=sum(percent))

mt_summary <- subset(mt_summary, Treatment == "MES")

mt_summary <- mt_summary %>%
  group_by(order) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            count=n())
