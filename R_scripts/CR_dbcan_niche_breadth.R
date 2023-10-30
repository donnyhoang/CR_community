library(ggplot2)
library(cowplot)
library(reshape2)
library(knitr)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(ggpubr)
library(broom)
library(MicroNiche)

palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
palette <-c( "#56B4E9","#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")


tax <- read.csv("CR_bin_tax_cluster.csv", header = TRUE)

cazy_mt <-read.csv("CR_dbcan_total_table_30m.csv", header=TRUE)
length(unique(cazy_mt$variable))
cazy_mt <- subset(cazy_mt, Substrate != "Other")
cazy_mt$Treatment<-factor(cazy_mt$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))

cazy_mg <-read.csv("CR_dbcan_total_table_14m.csv", header=TRUE)
length(unique(cazy_mg$variable))
cazy_mg <- subset(cazy_mg, Substrate != "Other")

#I want a matrix treatment, where cazy grouping is "environment"

sampleinfo <- c("Cellulose", "Cellobiose", "Xylan",
                "Starch", "Trehalose","Pectin","Lignin",
                "Peptidoglycan","Chitin")

### Subset into separate dataframes
G05_mt <- subset(cazy_mt, Treatment == "G05")
G10_mt <- subset(cazy_mt, Treatment == "G10")
MES_mt <- subset(cazy_mt, Treatment == "MES")
pH6_mt <- subset(cazy_mt, Treatment == "pH6")
pH8_mt <- subset(cazy_mt, Treatment == "pH8")
pH10_mt <- subset(cazy_mt, Treatment == "pH10")


### Combine df, name them for easy tracking
mt_list <- list(G05_mt, G10_mt, MES_mt, pH6_mt, pH8_mt, pH10_mt)
names(mt_list) <- c("G05", "G10", "MES", "pH6", "pH8", "pH10")

# check
mt_list[3]

# cast df into matrix
cast_matrix <- function(x) {
  dcast(x, formula = bin ~ Substrate, fun.aggregate = sum, value.var = "value")
  
}

matrix_mt <- lapply(mt_list, cast_matrix)

#check
matrix_mt[1]


# calculate levin's niche breadth
get_Bn <- function(x) {
  levins.Bn(x, 9, sampleinfo)

}

levinsBn_mt <- lapply(matrix_mt, get_Bn)

#check
levinsBn_mt[1]


## Write name (Treatment) into new column
# actually, didn't need this
#levinsBn_list_mt <- Map(cbind, levinsBn_mt, Treatment = names(levinsBn_mt))

# Write into one new big df
total_bn_mt <- do.call(rbind, lapply(levinsBn_mt, data.frame))
total_bn_mt[c("Treatment", "bin")] <- str_split_fixed(rownames(total_bn_mt), "[.]", 2)
rownames(total_bn_mt) <- NULL
total_bn_mt$Type <- "MT"
total_bn_mt <- filter(total_bn_mt , P.adj < 0.05)



###### Repeat for MG


### Subset into separate dataframes
G05_mg <- subset(cazy_mg, Treatment == "G05")
G10_mg <- subset(cazy_mg, Treatment == "G10")
MES_mg <- subset(cazy_mg, Treatment == "MES")
pH6_mg <- subset(cazy_mg, Treatment == "pH6")
pH8_mg <- subset(cazy_mg, Treatment == "pH8")
pH10_mg <- subset(cazy_mg, Treatment == "pH10")


### Combine df, name them for easy tracking
mg_list <- list(G05_mg, G10_mg, MES_mg, pH6_mg, pH8_mg, pH10_mg)
names(mg_list) <- c("G05", "G10", "MES", "pH6", "pH8", "pH10")

# check
mg_list[3]

# cast df into matrix
matrix_mg <- lapply(mg_list, cast_matrix)

#check
matrix_mg[1]


# calculate levin's niche breadth
levinsBn_mg <- lapply(matrix_mg, get_Bn)

#check
levinsBn_mg[1]


# Write into one new big df
total_bn_mg <- do.call(rbind, lapply(levinsBn_mg, data.frame))
total_bn_mg[c("Treatment", "bin")] <- str_split_fixed(rownames(total_bn_mg), "[.]", 2)
rownames(total_bn_mg) <- NULL
total_bn_mg$Type <- "MG"
total_bn_mg <- filter(total_bn_mg , P.adj < 0.05)

## Merge MG and MT Bn calculations
total_bn <- rbind(total_bn_mg, total_bn_mt)
total_bn <- merge(total_bn, tax, by = "bin")

total_bn$order <- factor(total_bn$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                  "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                  "Sphingobacteriales","Xanthomonadales","no support"))

total_bn$Treatment<-factor(total_bn$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
palette <-c( "#56B4E9","#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")

#no comparison for enterobaterales MG vs MT, remove
#total_bn <- subset(total_bn, order != "Enterobacterales")

total_bn_summary <- total_bn %>%
  group_by(order, Type) %>%
  summarise(mean=mean(Bn),
            sd = sd(Bn),
            count=n())

#trim_bn <- subset(total_bn, Substrate == "Pectin" & order %in% c("Corynebacteriales", "Lactobacillales", "Sphingobacteriales"))

bnplot <- ggplot(total_bn, aes(x=order, y = Bn, fill = Type)) +
  geom_boxplot() +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size =13),
        legend.position = "bottom") +
  labs(y = "Levin's Niche Breadth Index", x = "Order", fill = "Dataset") +
  stat_compare_means(aes(group=Type), label = "p.signif", method = "t.test") +
  scale_fill_manual(values = c("#1E88E5", "#FFC107"), labels = c("Metagenomic", "Metatranscriptomic")) +
  ggtitle("Niche Breadth of Bacterial Orders")
bnplot 
ggsave("bnplot.tiff", device = "tiff", dpi = 700)
ggsave("bnplot.png", device = "png", dpi = 700)
ggsave("bnplot.pdf", device = "pdf", dpi = 700)

########### Levin's overlap


levins_overlap <-  function(x) {
  temp <- levins.overlap(as.data.frame(x))
  melt(as.data.frame(temp))
}

lo_mt <- lapply(matrix_mt, levins_overlap)
lo_total_mt <- do.call(rbind, lapply(lo_mt, data.frame))
lo_total_mt[c("Treatment", "unused")] <- str_split_fixed(rownames(lo_total_mt), "[.]", 2)
rownames(lo_total_mt) <- NULL
lo_total_mt <- lo_total_mt[,-c(5)]
colnames(lo_total_mt) <- c("bin", "biny", "value", "Treatment")
lo_total_mt <- subset(lo_total_mt, bin != biny)
lo_total_mt$Treatment<-factor(lo_total_mt$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
lo_total_mt <- merge(lo_total_mt, tax, by = "bin")

#order order
lo_total_mt$order <- factor(lo_total_mt$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                                        "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                        "Sphingobacteriales","Xanthomonadales","no support"))
#order bins
lo_total_mt$bin <- factor(lo_total_mt$bin, levels=c("MES_1_D_cluster_DBSCAN_round2_17",
                                                    "G05_2_C_cluster_DBSCAN_round3_1",
                                                    "G10_7_C_cluster_DBSCAN_round1_6",
                                                    "pH10_1_E_cluster_DBSCAN_round2_19",
                                                    "pH10_2_A_cluster_DBSCAN_round6_5",
                                                    "pH8_7_D_cluster_DBSCAN_round9_5",
                                                    "pH10_2_A_cluster_DBSCAN_round2_4",
                                                    "G05_3_E_cluster_DBSCAN_round7_3",
                                                    "pH6_2_E_cluster_DBSCAN_round5_7",
                                                    "pH10_1_E_cluster_DBSCAN_round1_4",
                                                    "pH10_1_E_cluster_DBSCAN_round2_2",
                                                    "G05_1_E_cluster_DBSCAN_round3_4",
                                                    "G05_5_C_cluster_DBSCAN_round2_2",
                                                    "G10_2_C_cluster_DBSCAN_round1_10",
                                                    "MES_1_E_cluster_DBSCAN_round6_3",
                                                    "pH10_6_C_cluster_DBSCAN_round2_10",
                                                    "pH6_1_D_cluster_DBSCAN_round2_7",
                                                    "pH8_6_D_cluster_DBSCAN_round6_12",
                                                    "G05_6_E_cluster_DBSCAN_round5_10",
                                                    "G05_6_E_cluster_DBSCAN_round7_6",
                                                    "G05_7_A_cluster_DBSCAN_round8_2",
                                                    "G10_3_A_cluster_DBSCAN_round5_0",
                                                    "G10_3_D_cluster_DBSCAN_round210_2",
                                                    "G10_6_B_cluster_DBSCAN_round179_0",
                                                    "G10_7_A_cluster_DBSCAN_round160_0",
                                                    "G10_7_C_cluster_DBSCAN_round7_14",
                                                    "G10_7_C_cluster_DBSCAN_round7_2",
                                                    "G10_7_E_cluster_DBSCAN_round81_0",
                                                    "MES_2_C_cluster_DBSCAN_round3_1",
                                                    "MES_7_C_cluster_DBSCAN_round4_2",
                                                    "MES_7_D_cluster_DBSCAN_round362_17",
                                                    "pH10_6_A_cluster_DBSCAN_round525_32",
                                                    "pH10_6_E_cluster_DBSCAN_round5_14",
                                                    "pH10_7_A_cluster_DBSCAN_round2_25",
                                                    "pH10_7_F_cluster_DBSCAN_round4_20",
                                                    "pH6_1_A_cluster_DBSCAN_round1_4",
                                                    "pH6_4_F_cluster_DBSCAN_round203_0",
                                                    "G05_7_C_cluster_DBSCAN_round4_18",
                                                    "pH8_2_E_cluster_DBSCAN_round1_4",
                                                    "G05_6_E_cluster_DBSCAN_round6_6",
                                                    "G10_6_E_cluster_DBSCAN_round5_1",
                                                    "pH10_4_D_cluster_DBSCAN_round1_1",
                                                    "pH10_7_B_cluster_DBSCAN_round161_0",
                                                    "pH10_7_C_cluster_DBSCAN_round5_6",
                                                    "pH10_7_C_cluster_DBSCAN_round7_1",
                                                    "pH10_7_C_cluster_DBSCAN_round7_4",
                                                    "pH10_7_E_cluster_DBSCAN_round3_20",
                                                    "pH8_7_C_cluster_DBSCAN_round6_1",
                                                    "G05_2_B_cluster_DBSCAN_round2_17",
                                                    "pH10_7_F_cluster_DBSCAN_round5_0",
                                                    "pH6_6_D_cluster_DBSCAN_round9_1",
                                                    "MES_7_E_cluster_DBSCAN_round2_10",
                                                    "MES_5_C_cluster_DBSCAN_round8_4",
                                                    "pH8_7_F_cluster_DBSCAN_round202_0"
))

#annotations
org <- tax[,c(1,6)]
org$order <- factor(org$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                        "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                        "Sphingobacteriales","Xanthomonadales","no support"))
org$bin <- factor(org$bin, levels=c("MES_1_D_cluster_DBSCAN_round2_17",
                                    "G05_2_C_cluster_DBSCAN_round3_1",
                                    "G10_7_C_cluster_DBSCAN_round1_6",
                                    "pH10_1_E_cluster_DBSCAN_round2_19",
                                    "pH10_2_A_cluster_DBSCAN_round6_5",
                                    "pH8_7_D_cluster_DBSCAN_round9_5",
                                    "pH10_2_A_cluster_DBSCAN_round2_4",
                                    "G05_3_E_cluster_DBSCAN_round7_3",
                                    "pH6_2_E_cluster_DBSCAN_round5_7",
                                    "pH10_1_E_cluster_DBSCAN_round1_4",
                                    "pH10_1_E_cluster_DBSCAN_round2_2",
                                    "G05_1_E_cluster_DBSCAN_round3_4",
                                    "G05_5_C_cluster_DBSCAN_round2_2",
                                    "G10_2_C_cluster_DBSCAN_round1_10",
                                    "MES_1_E_cluster_DBSCAN_round6_3",
                                    "pH10_6_C_cluster_DBSCAN_round2_10",
                                    "pH6_1_D_cluster_DBSCAN_round2_7",
                                    "pH8_6_D_cluster_DBSCAN_round6_12",
                                    "G05_6_E_cluster_DBSCAN_round5_10",
                                    "G05_6_E_cluster_DBSCAN_round7_6",
                                    "G05_7_A_cluster_DBSCAN_round8_2",
                                    "G10_3_A_cluster_DBSCAN_round5_0",
                                    "G10_3_D_cluster_DBSCAN_round210_2",
                                    "G10_6_B_cluster_DBSCAN_round179_0",
                                    "G10_7_A_cluster_DBSCAN_round160_0",
                                    "G10_7_C_cluster_DBSCAN_round7_14",
                                    "G10_7_C_cluster_DBSCAN_round7_2",
                                    "G10_7_E_cluster_DBSCAN_round81_0",
                                    "MES_2_C_cluster_DBSCAN_round3_1",
                                    "MES_7_C_cluster_DBSCAN_round4_2",
                                    "MES_7_D_cluster_DBSCAN_round362_17",
                                    "pH10_6_A_cluster_DBSCAN_round525_32",
                                    "pH10_6_E_cluster_DBSCAN_round5_14",
                                    "pH10_7_A_cluster_DBSCAN_round2_25",
                                    "pH10_7_F_cluster_DBSCAN_round4_20",
                                    "pH6_1_A_cluster_DBSCAN_round1_4",
                                    "pH6_4_F_cluster_DBSCAN_round203_0",
                                    "G05_7_C_cluster_DBSCAN_round4_18",
                                    "pH8_2_E_cluster_DBSCAN_round1_4",
                                    "G05_6_E_cluster_DBSCAN_round6_6",
                                    "G10_6_E_cluster_DBSCAN_round5_1",
                                    "pH10_4_D_cluster_DBSCAN_round1_1",
                                    "pH10_7_B_cluster_DBSCAN_round161_0",
                                    "pH10_7_C_cluster_DBSCAN_round5_6",
                                    "pH10_7_C_cluster_DBSCAN_round7_1",
                                    "pH10_7_C_cluster_DBSCAN_round7_4",
                                    "pH10_7_E_cluster_DBSCAN_round3_20",
                                    "pH8_7_C_cluster_DBSCAN_round6_1",
                                    "G05_2_B_cluster_DBSCAN_round2_17",
                                    "pH10_7_F_cluster_DBSCAN_round5_0",
                                    "pH6_6_D_cluster_DBSCAN_round9_1",
                                    "MES_7_E_cluster_DBSCAN_round2_10",
                                    "MES_5_C_cluster_DBSCAN_round8_4",
                                    "pH8_7_F_cluster_DBSCAN_round202_0"
))

annotatex <- ggplot(org) +
  geom_bar(aes(x=bin, y = 1, fill = order), stat = "identity", width=1) +
  scale_fill_manual(values = palorg) +
  theme_void(base_size = 15) +
  guides(fill=guide_legend("Order")) +
  theme(
    panel.spacing.x = unit(1, "mm"))
annotatex

annotatey <- ggplot(org) +
  geom_bar(aes(x=1, y = bin, fill = order), stat = "identity", width=1) +
  scale_fill_manual(values = palorg) +
  theme_void() +
  theme(legend.position = "none")
annotatey


lo_total_mt_summary <- lo_total_mt %>%
  group_by(bin, order) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            count=n())

#write.csv(lo_total_mt_summary, "CR_dbcan_lo_summary.csv", row.names = FALSE)




lo_total_mt$order <- factor(lo_total_mt$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales","Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                                        "Sphingobacteriales","Xanthomonadales","no support"))

p_lo_dot_mt <- ggplot(lo_total_mt, aes(x=bin, y = value, color = Treatment)) +
  geom_jitter(alpha = 0.7) +
  theme_classic(base_size=25) +
  theme(axis.text.x.bottom = element_blank(), 
        panel.spacing.x = unit(1, "mm"),
        axis.title.x = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
        #axis.text.x=element_text(angle = 45, hjust = 1)
  ) +
  scale_color_manual(values=palette) +
  labs(y = "Levin's Overlap") +
  ggtitle("Levin's Overlap of MAG CAZy Expression")
p_lo_dot_mt

legend <- plot_grid(get_legend(p_lo_dot_mt), get_legend(annotatex), ncol = 1)

p1 <- p_lo_dot_mt + theme(legend.position="none")
p2 <- annotatex + theme(legend.position="none")

dot <- plot_grid(p1,p2, align = "v", ncol = 1, axis = "tb", rel_heights = c(10, 0.3))
dot2 <- plot_grid(dot, legend, nrow = 1, rel_widths = c(10,1.5))

dot2

ggsave("lo_dot.tiff", device = "tiff", dpi = 700)
ggsave("lo_dot.png", device = "png", dpi = 700)
ggsave("lo_dot.pdf", device = "pdf", dpi = 700)





#pH6_2_E_cluster_DBSCAN_round5_7
#pH10_4_D_cluster_DBSCAN_round1_1

low_lo <- subset(lo_total_mt, bin == "pH10_4_D_cluster_DBSCAN_round1_1")

p_low_lo <- ggplot(low_lo, aes(x = Treatment, y = value, fill = Treatment)) +
  geom_boxplot() +
  theme_classic(base_size = 15) +
  #theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(y = "Levin's Overlap", x = "Treatment") +
  stat_compare_means(method = "anova") +
  scale_fill_manual(values = palette) +
  ggtitle("Levin's Overlap of Selected Rhizobiales MAG")

p_low_lo

ggsave("low_lo.tiff", device = "tiff", dpi = 700)
ggsave("low_lo.png", device = "png", dpi = 700)
ggsave("low_lo.pdf", device = "pdf", dpi = 700)


p_lo_hist_mt <- ggplot() +
  #geom_histogram(data=lo_total_mt, aes(x=value, fill = order, color = order), position="identity", alpha = 0.6) +
  geom_density(data=lo_total_mt, aes(x=value, fill = order, color = order), position="identity", alpha = 0.6) +
  theme_classic(base_size=15) +
  scale_color_manual(values=palorg) +
  scale_fill_manual(values=palorg) +
  labs(x="Levin's Niche Overlap", y = "Density") +
  guides(fill=guide_legend("Order"), color= "none") +
  ggtitle("Density Plot of Levin's Overlap")
p_lo_hist_mt
ggsave("pdensity.tiff", device = "tiff", dpi = 700)
ggsave("pdensity.png", device = "png", dpi = 700)
ggsave("pdensity.pdf", device = "pdf", dpi = 700)

### Redo, focusing on specific CAZY

#Subset into substrate+treatment
Cellulose_G05_mt <- subset(cazy_mt, Substrate == "Cellulose" & Treatment == "G05")
Cellulose_G10_mt <- subset(cazy_mt, Substrate == "Cellulose" & Treatment == "G10")
Cellulose_MES_mt <- subset(cazy_mt, Substrate == "Cellulose" & Treatment == "MES")
Cellulose_pH6_mt <- subset(cazy_mt, Substrate == "Cellulose" & Treatment == "pH6")
Cellulose_pH8_mt <- subset(cazy_mt, Substrate == "Cellulose" & Treatment == "pH8")
Cellulose_pH10_mt <- subset(cazy_mt, Substrate == "Cellulose" & Treatment == "pH10")

Xylan_G05_mt <- subset(cazy_mt, Substrate == "Xylan" & Treatment == "G05")
Xylan_G10_mt <- subset(cazy_mt, Substrate == "Xylan" & Treatment == "G10")
Xylan_MES_mt <- subset(cazy_mt, Substrate == "Xylan" & Treatment == "MES")
Xylan_pH6_mt <- subset(cazy_mt, Substrate == "Xylan" & Treatment == "pH6")
Xylan_pH8_mt <- subset(cazy_mt, Substrate == "Xylan" & Treatment == "pH8")
Xylan_pH10_mt <- subset(cazy_mt, Substrate == "Xylan" & Treatment == "pH10")

Lignin_G05_mt <- subset(cazy_mt, Substrate == "Lignin" & Treatment == "G05")
Lignin_G10_mt <- subset(cazy_mt, Substrate == "Lignin" & Treatment == "G10")
Lignin_MES_mt <- subset(cazy_mt, Substrate == "Lignin" & Treatment == "MES")
Lignin_pH6_mt <- subset(cazy_mt, Substrate == "Lignin" & Treatment == "pH6")
Lignin_pH8_mt <- subset(cazy_mt, Substrate == "Lignin" & Treatment == "pH8")
Lignin_pH10_mt <- subset(cazy_mt, Substrate == "Lignin" & Treatment == "pH10")

Cellobiose_G05_mt <- subset(cazy_mt, Substrate == "Cellobiose" & Treatment == "G05")
Cellobiose_G10_mt <- subset(cazy_mt, Substrate == "Cellobiose" & Treatment == "G10")
Cellobiose_MES_mt <- subset(cazy_mt, Substrate == "Cellobiose" & Treatment == "MES")
Cellobiose_pH6_mt <- subset(cazy_mt, Substrate == "Cellobiose" & Treatment == "pH6")
Cellobiose_pH8_mt <- subset(cazy_mt, Substrate == "Cellobiose" & Treatment == "pH8")
Cellobiose_pH10_mt <- subset(cazy_mt, Substrate == "Cellobiose" & Treatment == "pH10")

Starch_G05_mt <- subset(cazy_mt, Substrate == "Starch" & Treatment == "G05")
Starch_G10_mt <- subset(cazy_mt, Substrate == "Starch" & Treatment == "G10")
Starch_MES_mt <- subset(cazy_mt, Substrate == "Starch" & Treatment == "MES")
Starch_pH6_mt <- subset(cazy_mt, Substrate == "Starch" & Treatment == "pH6")
Starch_pH8_mt <- subset(cazy_mt, Substrate == "Starch" & Treatment == "pH8")
Starch_pH10_mt <- subset(cazy_mt, Substrate == "Starch" & Treatment == "pH10")

Trehalose_G05_mt <- subset(cazy_mt, Substrate == "Trehalose" & Treatment == "G05")
Trehalose_G10_mt <- subset(cazy_mt, Substrate == "Trehalose" & Treatment == "G10")
Trehalose_MES_mt <- subset(cazy_mt, Substrate == "Trehalose" & Treatment == "MES")
Trehalose_pH6_mt <- subset(cazy_mt, Substrate == "Trehalose" & Treatment == "pH6")
Trehalose_pH8_mt <- subset(cazy_mt, Substrate == "Trehalose" & Treatment == "pH8")
Trehalose_pH10_mt <- subset(cazy_mt, Substrate == "Trehalose" & Treatment == "pH10")

Pectin_G05_mt <- subset(cazy_mt, Substrate == "Pectin" & Treatment == "G05")
Pectin_G10_mt <- subset(cazy_mt, Substrate == "Pectin" & Treatment == "G10")
Pectin_MES_mt <- subset(cazy_mt, Substrate == "Pectin" & Treatment == "MES")
Pectin_pH6_mt <- subset(cazy_mt, Substrate == "Pectin" & Treatment == "pH6")
Pectin_pH8_mt <- subset(cazy_mt, Substrate == "Pectin" & Treatment == "pH8")
Pectin_pH10_mt <- subset(cazy_mt, Substrate == "Pectin" & Treatment == "pH10")

Peptidoglycan_G05_mt <- subset(cazy_mt, Substrate == "Peptidoglycan" & Treatment == "G05")
Peptidoglycan_G10_mt <- subset(cazy_mt, Substrate == "Peptidoglycan" & Treatment == "G10")
Peptidoglycan_MES_mt <- subset(cazy_mt, Substrate == "Peptidoglycan" & Treatment == "MES")
Peptidoglycan_pH6_mt <- subset(cazy_mt, Substrate == "Peptidoglycan" & Treatment == "pH6")
Peptidoglycan_pH8_mt <- subset(cazy_mt, Substrate == "Peptidoglycan" & Treatment == "pH8")
Peptidoglycan_pH10_mt <- subset(cazy_mt, Substrate == "Peptidoglycan" & Treatment == "pH10")

Chitin_G05_mt <- subset(cazy_mt, Substrate == "Chitin" & Treatment == "G05")
Chitin_G10_mt <- subset(cazy_mt, Substrate == "Chitin" & Treatment == "G10")
Chitin_MES_mt <- subset(cazy_mt, Substrate == "Chitin" & Treatment == "MES")
Chitin_pH6_mt <- subset(cazy_mt, Substrate == "Chitin" & Treatment == "pH6")
Chitin_pH8_mt <- subset(cazy_mt, Substrate == "Chitin" & Treatment == "pH8")
Chitin_pH10_mt <- subset(cazy_mt, Substrate == "Chitin" & Treatment == "pH10")



### Combine df, name them for easy tracking
substrate_mt_list <- list(Cellulose_G05_mt, Cellulose_G10_mt, Cellulose_MES_mt, Cellulose_pH6_mt, Cellulose_pH8_mt, Cellulose_pH10_mt,
                          Xylan_G05_mt, Xylan_G10_mt, Xylan_MES_mt, Xylan_pH6_mt, Xylan_pH8_mt, Xylan_pH10_mt,
                          Lignin_G05_mt, Lignin_G10_mt, Lignin_MES_mt, Lignin_pH6_mt, Lignin_pH8_mt, Lignin_pH10_mt,
                          Cellobiose_G05_mt, Cellobiose_G10_mt, Cellobiose_MES_mt, Cellobiose_pH6_mt, Cellobiose_pH8_mt, Cellobiose_pH10_mt,
                          Starch_G05_mt, Starch_G10_mt, Starch_MES_mt, Starch_pH6_mt, Starch_pH8_mt, Starch_pH10_mt,
                          Trehalose_G05_mt, Trehalose_G10_mt, Trehalose_MES_mt, Trehalose_pH6_mt, Trehalose_pH8_mt, Trehalose_pH10_mt,
                          Chitin_G05_mt, Chitin_G10_mt, Chitin_MES_mt, Chitin_pH6_mt, Chitin_pH8_mt, Chitin_pH10_mt,
                          Peptidoglycan_G05_mt, Peptidoglycan_G10_mt, Peptidoglycan_MES_mt, Peptidoglycan_pH6_mt, Peptidoglycan_pH8_mt, Peptidoglycan_pH10_mt,
                          Pectin_G05_mt, Pectin_G10_mt, Pectin_MES_mt, Pectin_pH6_mt, Pectin_pH8_mt, Pectin_pH10_mt)
names(substrate_mt_list) <- c("Cellulose_G05", "Cellulose_G10", "Cellulose_MES", "Cellulose_pH6", "Cellulose_pH8", "Cellulose_pH10",
                              "Xylan_G05", "Xylan_G10", "Xylan_MES", "Xylan_pH6", "Xylan_pH8", "Xylan_pH10",
                              "Lignin_G05", "Lignin_G10", "Lignin_MES", "Lignin_pH6", "Lignin_pH8", "Lignin_pH10",
                              "Cellobiose_G05", "Cellobiose_G10", "Cellobiose_MES", "Cellobiose_pH6", "Cellobiose_pH8", "Cellobiose_pH10",
                              "Starch_G05", "Starch_G10", "Starch_MES", "Starch_pH6", "Starch_pH8", "Starch_pH10",
                              "Trehalose_G05", "Trehalose_G10", "Trehalose_MES", "Trehalose_pH6", "Trehalose_pH8", "Trehalose_pH10",
                              "Chitin_G05", "Chitin_G10", "Chitin_MES", "Chitin_pH6", "Chitin_pH8", "Chitin_pH10",
                              "Peptidoglycan_G05", "Peptidoglycan_G10", "Peptidoglycan_MES", "Peptidoglycan_pH6", "Peptidoglycan_pH8", "Peptidoglycan_pH10",
                              "Pectin_G05", "Pectin_G10", "Pectin_MES", "Pectin_pH6", "Pectin_pH8", "Pectin_pH10")

# check
substrate_mt_list[3]

#cast into matrix

cast_matrix2 <- function(x) {
  dcast(x, formula = bin ~ CAZy, fun.aggregate = sum, value.var = "value")
  
}

substrate_matrix_mt <- lapply(substrate_mt_list, cast_matrix2)

#calculate levin's niche breadth

# calculate levin's niche breadth
get_Bn2 <- function(x) {
  R <- (length(colnames(x)) - 1)
  sampleInfo <- colnames(x)
  sampleInfo <- sampleInfo[-1]
  levins.Bn(x, R, sampleInfo)
  
}

levinsBn_mt <- lapply(substrate_matrix_mt, get_Bn2)

#check
levinsBn_mt[21]


## Write name (Treatment) into new column
# actually, didn't need this
#levinsBn_list_mt <- Map(cbind, levinsBn_mt, Treatment = names(levinsBn_mt))

# Write into one new big df
total_bn_mt <- do.call(rbind, lapply(levinsBn_mt, data.frame))
total_bn_mt[c("temp", "bin")] <- str_split_fixed(rownames(total_bn_mt), "[.]", 2)
total_bn_mt[c("Substrate", "Treatment")] <- str_split_fixed(total_bn_mt$temp, "_", 2)
rownames(total_bn_mt) <- NULL
total_bn_mt$Type <- "MT"
total_bn_mt <- filter(total_bn_mt , P.adj < 0.05)
total_bn_mt <- subset(total_bn_mt, Below.LOQ != "Y")
total_bn_mt <- total_bn_mt[,-c(5)]

### Repeat for metagenomic data

#Subset into substrate+treatment
Cellulose_G05_mg <- subset(cazy_mg, Substrate == "Cellulose" & Treatment == "G05")
Cellulose_G10_mg <- subset(cazy_mg, Substrate == "Cellulose" & Treatment == "G10")
Cellulose_MES_mg <- subset(cazy_mg, Substrate == "Cellulose" & Treatment == "MES")
Cellulose_pH6_mg <- subset(cazy_mg, Substrate == "Cellulose" & Treatment == "pH6")
Cellulose_pH8_mg <- subset(cazy_mg, Substrate == "Cellulose" & Treatment == "pH8")
Cellulose_pH10_mg <- subset(cazy_mg, Substrate == "Cellulose" & Treatment == "pH10")

Xylan_G05_mg <- subset(cazy_mg, Substrate == "Xylan" & Treatment == "G05")
Xylan_G10_mg <- subset(cazy_mg, Substrate == "Xylan" & Treatment == "G10")
Xylan_MES_mg <- subset(cazy_mg, Substrate == "Xylan" & Treatment == "MES")
Xylan_pH6_mg <- subset(cazy_mg, Substrate == "Xylan" & Treatment == "pH6")
Xylan_pH8_mg <- subset(cazy_mg, Substrate == "Xylan" & Treatment == "pH8")
Xylan_pH10_mg <- subset(cazy_mg, Substrate == "Xylan" & Treatment == "pH10")

Lignin_G05_mg <- subset(cazy_mg, Substrate == "Lignin" & Treatment == "G05")
Lignin_G10_mg <- subset(cazy_mg, Substrate == "Lignin" & Treatment == "G10")
Lignin_MES_mg <- subset(cazy_mg, Substrate == "Lignin" & Treatment == "MES")
Lignin_pH6_mg <- subset(cazy_mg, Substrate == "Lignin" & Treatment == "pH6")
Lignin_pH8_mg <- subset(cazy_mg, Substrate == "Lignin" & Treatment == "pH8")
Lignin_pH10_mg <- subset(cazy_mg, Substrate == "Lignin" & Treatment == "pH10")

Cellobiose_G05_mg <- subset(cazy_mg, Substrate == "Cellobiose" & Treatment == "G05")
Cellobiose_G10_mg <- subset(cazy_mg, Substrate == "Cellobiose" & Treatment == "G10")
Cellobiose_MES_mg <- subset(cazy_mg, Substrate == "Cellobiose" & Treatment == "MES")
Cellobiose_pH6_mg <- subset(cazy_mg, Substrate == "Cellobiose" & Treatment == "pH6")
Cellobiose_pH8_mg <- subset(cazy_mg, Substrate == "Cellobiose" & Treatment == "pH8")
Cellobiose_pH10_mg <- subset(cazy_mg, Substrate == "Cellobiose" & Treatment == "pH10")

Starch_G05_mg <- subset(cazy_mg, Substrate == "Starch" & Treatment == "G05")
Starch_G10_mg <- subset(cazy_mg, Substrate == "Starch" & Treatment == "G10")
Starch_MES_mg <- subset(cazy_mg, Substrate == "Starch" & Treatment == "MES")
Starch_pH6_mg <- subset(cazy_mg, Substrate == "Starch" & Treatment == "pH6")
Starch_pH8_mg <- subset(cazy_mg, Substrate == "Starch" & Treatment == "pH8")
Starch_pH10_mg <- subset(cazy_mg, Substrate == "Starch" & Treatment == "pH10")

Trehalose_G05_mg <- subset(cazy_mg, Substrate == "Trehalose" & Treatment == "G05")
Trehalose_G10_mg <- subset(cazy_mg, Substrate == "Trehalose" & Treatment == "G10")
Trehalose_MES_mg <- subset(cazy_mg, Substrate == "Trehalose" & Treatment == "MES")
Trehalose_pH6_mg <- subset(cazy_mg, Substrate == "Trehalose" & Treatment == "pH6")
Trehalose_pH8_mg <- subset(cazy_mg, Substrate == "Trehalose" & Treatment == "pH8")
Trehalose_pH10_mg <- subset(cazy_mg, Substrate == "Trehalose" & Treatment == "pH10")

Pectin_G05_mg <- subset(cazy_mg, Substrate == "Pectin" & Treatment == "G05")
Pectin_G10_mg <- subset(cazy_mg, Substrate == "Pectin" & Treatment == "G10")
Pectin_MES_mg <- subset(cazy_mg, Substrate == "Pectin" & Treatment == "MES")
Pectin_pH6_mg <- subset(cazy_mg, Substrate == "Pectin" & Treatment == "pH6")
Pectin_pH8_mg <- subset(cazy_mg, Substrate == "Pectin" & Treatment == "pH8")
Pectin_pH10_mg <- subset(cazy_mg, Substrate == "Pectin" & Treatment == "pH10")

Peptidoglycan_G05_mg <- subset(cazy_mg, Substrate == "Peptidoglycan" & Treatment == "G05")
Peptidoglycan_G10_mg <- subset(cazy_mg, Substrate == "Peptidoglycan" & Treatment == "G10")
Peptidoglycan_MES_mg <- subset(cazy_mg, Substrate == "Peptidoglycan" & Treatment == "MES")
Peptidoglycan_pH6_mg <- subset(cazy_mg, Substrate == "Peptidoglycan" & Treatment == "pH6")
Peptidoglycan_pH8_mg <- subset(cazy_mg, Substrate == "Peptidoglycan" & Treatment == "pH8")
Peptidoglycan_pH10_mg <- subset(cazy_mg, Substrate == "Peptidoglycan" & Treatment == "pH10")

Chitin_G05_mg <- subset(cazy_mg, Substrate == "Chitin" & Treatment == "G05")
Chitin_G10_mg <- subset(cazy_mg, Substrate == "Chitin" & Treatment == "G10")
Chitin_MES_mg <- subset(cazy_mg, Substrate == "Chitin" & Treatment == "MES")
Chitin_pH6_mg <- subset(cazy_mg, Substrate == "Chitin" & Treatment == "pH6")
Chitin_pH8_mg <- subset(cazy_mg, Substrate == "Chitin" & Treatment == "pH8")
Chitin_pH10_mg <- subset(cazy_mg, Substrate == "Chitin" & Treatment == "pH10")



### Combine df, name them for easy tracking
substrate_mg_list <- list(Cellulose_G05_mg, Cellulose_G10_mg, Cellulose_MES_mg, Cellulose_pH6_mg, Cellulose_pH8_mg, Cellulose_pH10_mg,
                          Xylan_G05_mg, Xylan_G10_mg, Xylan_MES_mg, Xylan_pH6_mg, Xylan_pH8_mg, Xylan_pH10_mg,
                          Lignin_G05_mg, Lignin_G10_mg, Lignin_MES_mg, Lignin_pH6_mg, Lignin_pH8_mg, Lignin_pH10_mg,
                          Cellobiose_G05_mg, Cellobiose_G10_mg, Cellobiose_MES_mg, Cellobiose_pH6_mg, Cellobiose_pH8_mg, Cellobiose_pH10_mg,
                          Starch_G05_mg, Starch_G10_mg, Starch_MES_mg, Starch_pH6_mg, Starch_pH8_mg, Starch_pH10_mg,
                          Trehalose_G05_mg, Trehalose_G10_mg, Trehalose_MES_mg, Trehalose_pH6_mg, Trehalose_pH8_mg, Trehalose_pH10_mg,
                          Chitin_G05_mg, Chitin_G10_mg, Chitin_MES_mg, Chitin_pH6_mg, Chitin_pH8_mg, Chitin_pH10_mg,
                          Peptidoglycan_G05_mg, Peptidoglycan_G10_mg, Peptidoglycan_MES_mg, Peptidoglycan_pH6_mg, Peptidoglycan_pH8_mg, Peptidoglycan_pH10_mg,
                          Pectin_G05_mg, Pectin_G10_mg, Pectin_MES_mg, Pectin_pH6_mg, Pectin_pH8_mg, Pectin_pH10_mg)
names(substrate_mg_list) <- c("Cellulose_G05", "Cellulose_G10", "Cellulose_MES", "Cellulose_pH6", "Cellulose_pH8", "Cellulose_pH10",
                              "Xylan_G05", "Xylan_G10", "Xylan_MES", "Xylan_pH6", "Xylan_pH8", "Xylan_pH10",
                              "Lignin_G05", "Lignin_G10", "Lignin_MES", "Lignin_pH6", "Lignin_pH8", "Lignin_pH10",
                              "Cellobiose_G05", "Cellobiose_G10", "Cellobiose_MES", "Cellobiose_pH6", "Cellobiose_pH8", "Cellobiose_pH10",
                              "Starch_G05", "Starch_G10", "Starch_MES", "Starch_pH6", "Starch_pH8", "Starch_pH10",
                              "Trehalose_G05", "Trehalose_G10", "Trehalose_MES", "Trehalose_pH6", "Trehalose_pH8", "Trehalose_pH10",
                              "Chitin_G05", "Chitin_G10", "Chitin_MES", "Chitin_pH6", "Chitin_pH8", "Chitin_pH10",
                              "Peptidoglycan_G05", "Peptidoglycan_G10", "Peptidoglycan_MES", "Peptidoglycan_pH6", "Peptidoglycan_pH8", "Peptidoglycan_pH10",
                              "Pectin_G05", "Pectin_G10", "Pectin_MES", "Pectin_pH6", "Pectin_pH8", "Pectin_pH10")

# check
substrate_mg_list[3]

#cast into matrix

substrate_matrix_mg <- lapply(substrate_mg_list, cast_matrix2)

#calculate levin's niche breadth

levinsBn_mg <- lapply(substrate_matrix_mg, get_Bn2)

#check
levinsBn_mg[21]


# Write into one new big df
total_bn_mg <- do.call(rbind, lapply(levinsBn_mg, data.frame))
total_bn_mg[c("temp", "bin")] <- str_split_fixed(rownames(total_bn_mg), "[.]", 2)
total_bn_mg[c("Substrate", "Treatment")] <- str_split_fixed(total_bn_mg$temp, "_", 2)
rownames(total_bn_mg) <- NULL
total_bn_mg$Type <- "MG"
total_bn_mg <- filter(total_bn_mg , P.adj < 0.05)
total_bn_mg <- subset(total_bn_mg, Below.LOQ != "Y")
total_bn_mg <- total_bn_mg[,-c(5)]




#calculate levin's overlap

#calc and write into large dataframe
substrate_lo_mt <- lapply(substrate_matrix_mt, levins_overlap)
substrate_lo_total_mt <- do.call(rbind, lapply(substrate_lo_mt, data.frame))

#fix up columns
substrate_lo_total_mt[c("temp", "unused")] <- str_split_fixed(rownames(substrate_lo_total_mt), "[.]", 2)
substrate_lo_total_mt[c("Substrate", "Treatment")] <- str_split_fixed(substrate_lo_total_mt$temp, "_", 2)
rownames(substrate_lo_total_mt) <- NULL
substrate_lo_total_mt <- substrate_lo_total_mt[,-c(4,5)]
colnames(substrate_lo_total_mt) <- c("bin", "biny", "value","Substrate","Treatment")

#looks like I have some bins below LOQ (marked with *), filter them out

substrate_lo_total_mt$LOQ <- str_sub(substrate_lo_total_mt$bin, -1)
substrate_lo_total_mt <- subset(substrate_lo_total_mt, LOQ != "*")
substrate_lo_total_mt <- substrate_lo_total_mt[,-c(6)]

substrate_lo_total_mt$LOQ <- str_sub(substrate_lo_total_mt$biny, -1)
substrate_lo_total_mt <- subset(substrate_lo_total_mt, LOQ != "*")
substrate_lo_total_mt <- substrate_lo_total_mt[,-c(6)]


substrate_lo_total_mt$Type <- "MT"



###############calculate levin's overlap for MG ########################

#calc and write into large dataframe
substrate_lo_mg <- lapply(substrate_matrix_mg, levins_overlap)
substrate_lo_total_mg <- do.call(rbind, lapply(substrate_lo_mg, data.frame))

#fix up columns
substrate_lo_total_mg[c("temp", "unused")] <- str_split_fixed(rownames(substrate_lo_total_mg), "[.]", 2)
substrate_lo_total_mg[c("Substrate", "Treatment")] <- str_split_fixed(substrate_lo_total_mg$temp, "_", 2)
rownames(substrate_lo_total_mg) <- NULL
substrate_lo_total_mg <- substrate_lo_total_mg[,-c(4,5)]
colnames(substrate_lo_total_mg) <- c("bin", "biny", "value","Substrate","Treatment")

#looks like I have some bins below LOQ (marked with *), filter them out

substrate_lo_total_mg$LOQ <- str_sub(substrate_lo_total_mg$bin, -1)
substrate_lo_total_mg <- subset(substrate_lo_total_mg, LOQ != "*")
substrate_lo_total_mg <- substrate_lo_total_mg[,-c(6)]


substrate_lo_total_mg$LOQ <- str_sub(substrate_lo_total_mg$biny, -1)
substrate_lo_total_mg <- subset(substrate_lo_total_mg, LOQ != "*")
substrate_lo_total_mg <- substrate_lo_total_mg[,-c(6)]


substrate_lo_total_mg$Type <- "MG"



########### Merge MG and MT################## 
substrate_lo <- rbind(substrate_lo_total_mg, substrate_lo_total_mt)

#factors
substrate_lo$Treatment<-factor(substrate_lo$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
substrate_lo$Substrate <- factor(substrate_lo$Substrate, levels = c("Cellulose", "Xylan", "Lignin"))

##merge
substrate_lo <- merge(substrate_lo, tax, by = "bin")


palorg3 <- c("#FB8072","#B3DE69","#1f78b4")

substrate_lo <- subset(substrate_lo, order %in% c("Corynebacteriales","Lactobacillales","Sphingobacteriales"))
substrate_lo$order <- factor(substrate_lo$order, levels = c("Corynebacteriales","Lactobacillales","Sphingobacteriales"))
substrate_lo <- subset(substrate_lo, Substrate %in% c("Cellulose", "Xylan", "Lignin"))


p_substrate_lo <- ggplot() +
  geom_density(data=substrate_lo, aes(x=value, fill = order, color = order), position="identity", alpha = 0.6) +
  theme_classic(base_size=15) +
  scale_color_manual(values=palorg3) +
  scale_fill_manual(values=palorg3) +
  labs(x="Levin's Niche Overlap", y = "Density") +
  guides(fill=guide_legend("Order"), color= "none") +
  facet_grid(Type ~Substrate) +
  ggtitle("Levin's Niche Overlap of Selected Bacterial Orders")
p_substrate_lo

ggsave("lo_substrate_density.tiff", device = "tiff", dpi = 700)
ggsave("lo_substrate_density.png", device = "png", dpi = 700)
ggsave("lo_substrate_density.pdf", device = "pdf", dpi = 700)


p_sub_lo_box <- ggplot(data = substrate_lo, aes (x = order, y = value, fill = Type)) +
  geom_boxplot() +
  theme_classic(base_size=15) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(y = "Levin's Niche Overlap", x = "Order", fill = "Dataset") +
  scale_fill_manual(values = c("#1E88E5", "#FFC107"), labels = c("Metagenomic", "Metatranscriptomic")) +
  facet_grid(~Substrate) +
  stat_compare_means(aes(group = Type), label = "p.signif", method = "t.test") +
  ggtitle("Niche Overlap of Selected Bacterial Orders")
p_sub_lo_box

ggsave("lo_substrate_box.tiff", device = "tiff", dpi = 700)
ggsave("lo_substrate_box.png", device = "png", dpi = 700)
ggsave("lo_substrate_box.pdf", device = "pdf", dpi = 700)




substrate_lo_summary <- substrate_lo %>%
  group_by (order, Substrate, Type) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            count=n())

substrate_lo_stat <- substrate_lo %>%
  group_by(order, Substrate) %>%
  summarise(ttest = t.test(value[Type == "MG"], value[Type == "MT"])$p.value)


############################## What is Sphingo doing with lignin? ################################

lignin_mt <- subset(cazy_mt, Substrate == "Lignin")

lignin_mt <- lignin_mt %>%
  group_by(bin, order, CAZy) %>%
  summarise(sum = sum(value))
lignin_mt$Type <- "MT"


lignin_mg <- subset(cazy_mg, Substrate == "Lignin")
lignin_mg <- lignin_mg %>%
  group_by(bin, order, CAZy) %>%
  summarise(sum = sum(value))
lignin_mg$Type <- "MG"


lignin_mt_trim <- subset(lignin_mt, !CAZy %in% c("AA3", "AA3_1", "AA3_2", "AA3_3", "AA3_4"))



plig_mt <- ggplot(lignin_mt_trim, aes(y = CAZy, x = bin, fill = sum)) +
  geom_tile() +
  theme(axis.text.x = element_blank())

plig_mt




legend <- plot_grid(get_legend(plig_mt), get_legend(annotatex), ncol = 1)

p1 <- plig_mt + theme(legend.position="none")
p2 <- annotatex + theme(legend.position="none")

lig <- plot_grid(p1,p2, align = "v", ncol = 1, axis = "tb", rel_heights = c(10, 0.3))
plot_grid(lig , legend, nrow = 1, rel_widths = c(10,1.5))


lignin_total$order <- factor(lignin_total$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                        "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                        "Sphingobacteriales","Xanthomonadales","no support"))
lignin_total$bin <- factor(lignin_total$bin, levels=c("MES_1_D_cluster_DBSCAN_round2_17",
                                    "G05_2_C_cluster_DBSCAN_round3_1",
                                    "G10_7_C_cluster_DBSCAN_round1_6",
                                    "pH10_1_E_cluster_DBSCAN_round2_19",
                                    "pH10_2_A_cluster_DBSCAN_round6_5",
                                    "pH8_7_D_cluster_DBSCAN_round9_5",
                                    "pH10_2_A_cluster_DBSCAN_round2_4",
                                    "G05_3_E_cluster_DBSCAN_round7_3",
                                    "pH6_2_E_cluster_DBSCAN_round5_7",
                                    "pH10_1_E_cluster_DBSCAN_round1_4",
                                    "pH10_1_E_cluster_DBSCAN_round2_2",
                                    "G05_1_E_cluster_DBSCAN_round3_4",
                                    "G05_5_C_cluster_DBSCAN_round2_2",
                                    "G10_2_C_cluster_DBSCAN_round1_10",
                                    "MES_1_E_cluster_DBSCAN_round6_3",
                                    "pH10_6_C_cluster_DBSCAN_round2_10",
                                    "pH6_1_D_cluster_DBSCAN_round2_7",
                                    "pH8_6_D_cluster_DBSCAN_round6_12",
                                    "G05_6_E_cluster_DBSCAN_round5_10",
                                    "G05_6_E_cluster_DBSCAN_round7_6",
                                    "G05_7_A_cluster_DBSCAN_round8_2",
                                    "G10_3_A_cluster_DBSCAN_round5_0",
                                    "G10_3_D_cluster_DBSCAN_round210_2",
                                    "G10_6_B_cluster_DBSCAN_round179_0",
                                    "G10_7_A_cluster_DBSCAN_round160_0",
                                    "G10_7_C_cluster_DBSCAN_round7_14",
                                    "G10_7_C_cluster_DBSCAN_round7_2",
                                    "G10_7_E_cluster_DBSCAN_round81_0",
                                    "MES_2_C_cluster_DBSCAN_round3_1",
                                    "MES_7_C_cluster_DBSCAN_round4_2",
                                    "MES_7_D_cluster_DBSCAN_round362_17",
                                    "pH10_6_A_cluster_DBSCAN_round525_32",
                                    "pH10_6_E_cluster_DBSCAN_round5_14",
                                    "pH10_7_A_cluster_DBSCAN_round2_25",
                                    "pH10_7_F_cluster_DBSCAN_round4_20",
                                    "pH6_1_A_cluster_DBSCAN_round1_4",
                                    "pH6_4_F_cluster_DBSCAN_round203_0",
                                    "G05_7_C_cluster_DBSCAN_round4_18",
                                    "pH8_2_E_cluster_DBSCAN_round1_4",
                                    "G05_6_E_cluster_DBSCAN_round6_6",
                                    "G10_6_E_cluster_DBSCAN_round5_1",
                                    "pH10_4_D_cluster_DBSCAN_round1_1",
                                    "pH10_7_B_cluster_DBSCAN_round161_0",
                                    "pH10_7_C_cluster_DBSCAN_round5_6",
                                    "pH10_7_C_cluster_DBSCAN_round7_1",
                                    "pH10_7_C_cluster_DBSCAN_round7_4",
                                    "pH10_7_E_cluster_DBSCAN_round3_20",
                                    "pH8_7_C_cluster_DBSCAN_round6_1",
                                    "G05_2_B_cluster_DBSCAN_round2_17",
                                    "pH10_7_F_cluster_DBSCAN_round5_0",
                                    "pH6_6_D_cluster_DBSCAN_round9_1",
                                    "MES_7_E_cluster_DBSCAN_round2_10",
                                    "MES_5_C_cluster_DBSCAN_round8_4",
                                    "pH8_7_F_cluster_DBSCAN_round202_0"
))
