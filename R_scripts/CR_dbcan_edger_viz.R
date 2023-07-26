library(ggplot2)
library(dplyr)
library(ggpubr)

data <- read.csv("CR_dbcan_edger_earlyvslate.csv")
data <- filter(data, PValue<0.05)
data <- na.omit(data)


data$Comparison <- factor(data$Comparison, levels = c("pH6","pH8","pH10","MES","G05","G10"))
data$func <- factor(data$func, levels = c("other","GMC oxidoreductases","amylase","other oxidases","betaglucanase","ligninase","cellulase","xylanase"))
data$cluster_total <- as.factor(data$cluster_total)

#define into large groups, AA and hydrolases

AA <- c("GMC oxidoreductases","other oxidases","ligninase")
GH <- c("amylase","betaglucanase","cellulase","xylanase")
other <- c("other")

data <- data %>%
  mutate(func2 = case_when(
    func %in% AA ~ "AA",
    func %in% GH ~ "GH",
    func %in% other ~ "other"
  ))


summary <- data %>%
  group_by(phylum) %>%
  summarise(count=n())

data$phylum <- factor(data$phylum, levels = c("Actinobacteria","Proteobacteria","Bacteroidetes","Firmicutes"))

data$class <- factor(data$class, levels =c ("Actinobacteria","Betaproteobacteria","Alphaproteobacteria","Sphingobacteriia","Bacilli","Gammaproteobacteria","Flavobacteriia","no support"))

data$order <- factor(data$order, levels = c("Lactobacillales","Micrococcales","Burkholderiales","Corynebacteriales","Sphingobacteriales",
                                            "Rhizobiales","no support","Pseudomonadales","Caulobacterales","Enterobacterales",
                                            "Flavobacteriales","Xanthomonadales","Propionibacteriales"))
data$func <- factor(data$func, levels = c("ligninase","GMC oxidoreductases","other oxidases","cellulase","betaglucanase","amylase","xylanase","other"))


palorg <- c("#B3DE69", #green
            "#FCCDE5", #pink
            "#FFFFB3", #yellow
            "#FB8072", #coral
            "#1f78b4", #blue
            "#e5c494", #tan
            "#D9D9D9",  #grey
            "#CCEBC5", #mint
            "#BEBADA", #lavender
            "#80B1D3", #cornflower
            "#FDB462", #orange
            "#fb9a99", #blush
            #"#8DD3C7", #teal
            "#BC80BD" #purple
            )
palette <-c( "#999999","#7fc97f","#beaed4","#fdc086","#ffff99")

#this one is more saturated, use for class
#palorg <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#80b1d3","#e6ab02","#a6761d","#666666")

data1 <- subset(data, func != "other")

p <- ggplot(data1 %>%
              arrange(order),
            aes(x=logCPM, y=logFC, color = order)) +
  geom_point(size=2) +
  scale_color_manual(values=palorg) +
  theme_classic(base_size = 15) +
  #coord_fixed() +
  ggtitle("Differential expression of CAZy contigs, Days 1-3 vs. Days 4-7") +
  facet_grid(Comparison ~ func)
p

########### Hydrolases ############

#define non GH
data <- data %>%
  mutate(is_GH = case_when(
    func2 == "GH" ~ "yes",
    func2 != "GH" ~ "no"
  ))


##makes two subsets, one of "no" and another of just GH contigs

data1 <- subset(data, is_GH == "no")
data2 <- subset(data, func2 == "GH")



pGH <- ggplot() +
  #geom_point(data = data1, aes(x = logCPM, y=logFC), color = "8C8C8C", size=2) +
  geom_point(data = data2 %>%
               arrange(order), aes(x=logCPM, y=logFC, color = order), size = 2) +
  scale_color_manual(values=palorg) +
  theme_classic(base_size = 15) +
  coord_fixed() +
  ggtitle("Differential Expression of GH contigs early vs late") +
  facet_grid( ~ Comparison)
pGH

############### Auxiliary Activity ####

#define non AA
data <- data %>%
  mutate(is_AA = case_when(
    func2 == "AA" ~ "yes",
    func2 != "AA" ~ "no"
  ))




data3 <- subset(data, is_AA == "no")
data4 <- subset(data, func2 == "AA")



pAA <- ggplot() +
  #geom_point(data = data3, aes(x = logCPM, y=logFC), color = "8C8C8C", size=2) +
  geom_point(data = data4 %>%
               arrange(order), aes(x=logCPM, y=logFC, color = order), size = 2) +
  scale_color_manual(values=palorg) +
  theme_classic(base_size = 15) +
  coord_fixed() +
  ggtitle("Differential Expression of AA contigs early vs late") +
  facet_grid(~ Comparison)
pAA

ggarrange(pGH, pAA, common.legend = TRUE)

d_lig <- subset(data, func == "ligninase")
d_gmc <- subset(data, func == "GMC oxidoreductases")
d_oxi <- subset(data, func == "other oxidases")
d_AA4 <- subset(data, CAZy == "AA4")
d_AA6 <- subset(data, CAZy == "AA6")
d_AA7 <- subset(data, CAZy == "AA7")
d_AA3 <- subset(data, CAZy == "AA3_1")





palorg2 <- c("#FCCDE5", #pink
            "#FFFFB3", #yellow
            "#FB8072", #coral
            "#1f78b4", #blue
            "#e5c494", #tan
            "#D9D9D9",  #grey
            "#CCEBC5", #mint
            "#BEBADA", #lavender
            "#80B1D3", #cornflower
            "#FDB462", #orange
            "#fb9a99", #blush
            #"#8DD3C7", #teal
            "#BC80BD" #purple
)

palette2 <- c("#419fde","#1f78b4","#1f78b4","#144c73","#144c73")
palette3 <- c("#cab2d6","#b2df8a","#fdbf6f","#fb9a99","#a6cee3")
p_lig <- ggplot() +
  geom_point(data = d_lig %>%
               arrange(order), aes(x=logCPM, y=logFC, color = CAZy), size = 5, alpha = 0.9) +
  scale_color_manual(values=palette3) +
  theme_classic(base_size = 15) +
  coord_fixed() +
  labs(y="log Fold Change", x = "Average log Counts per Million", color = "CAZy") +
  ggtitle("Differential Expression of Lignin-acting Contigs Days 1-3 vs Days 4-7") +
  facet_grid(~ Comparison)
p_lig


d_ligsphingo <- subset(d_lig, order == "Sphingobacteriales")

p_ligsph <- ggplot() +
  geom_point(data = d_ligsphingo %>%
               arrange(order), aes(x=logCPM, y=logFC, color = CAZy), size = 5, alpha = 0.9) +
  scale_color_manual(values=palette2) +
  theme_classic(base_size = 15) +
  coord_fixed() +
  labs(y="log Fold Change", x = "Average log Counts per Million", color = "Order") +
  ggtitle("Differential Expression of Sphingobacteriales Classified Lignin-acting Expression Days 1-3 vs Days 4-7") +
  facet_grid(~ Comparison)

p_ligsph


