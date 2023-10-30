library(ggplot2)

mag <- read.csv("CR_dbcan_MAG_annotation_total.csv", header = TRUE)

#count annotations, keep bin and CAZy separate
mag <- mag %>%
  group_by(bin,CAZy, Substrate, order) %>%
  summarise(count=n())

#sum up, group by Substrate and sum up for MAGs of same order
mag <- mag %>%
  group_by(bin, Substrate, order) %>%
  summarise(sum=sum(count))

mag_summary <- mag %>%
  group_by(Substrate, order) %>%
  summarise(mean=mean(sum),
            sd=sd(sum),
            count = n())



mag$Substrate <- factor(mag$Substrate, levels = c("Cellulose", "Cellobiose", "Xylan",
                                                  "Starch", "Trehalose","Pectin","Lignin",
                                                  "Peptidoglycan","Chitin","Other"))

mag$order <- factor(mag$order, levels=c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                        "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                        "Sphingobacteriales","Xanthomonadales","no support"))

#& order %in% c("Corynebacteriales", "Lactobacillales","Sphingobacteriales")

mag <- subset(mag, Substrate != "Other" & order %in% c("Corynebacteriales", "Lactobacillales","Sphingobacteriales"))



palorgcorynshingo <- c("#FB8072","#B3DE69","#1f78b4")

mag$bin <- factor(mag$bin, levels = c("G05_3_E_cluster_DBSCAN_round7_3", #Coryne
                                      "pH6_2_E_cluster_DBSCAN_round5_7",
                                      "G05_1_E_cluster_DBSCAN_round3_4", #Lacto
                                      "G05_5_C_cluster_DBSCAN_round2_2",
                                      "G10_2_C_cluster_DBSCAN_round1_10",
                                      "MES_1_E_cluster_DBSCAN_round6_3",
                                      "pH10_6_C_cluster_DBSCAN_round2_10",
                                      "pH6_1_D_cluster_DBSCAN_round2_7",
                                      "pH8_6_D_cluster_DBSCAN_round6_12",
                                      "G05_2_B_cluster_DBSCAN_round2_17", #Sphingo
                                      "pH10_7_F_cluster_DBSCAN_round5_0",
                                      "pH6_6_D_cluster_DBSCAN_round9_1"
))

mag_plot <- ggplot() +
  geom_point(data = mag, aes(x=Substrate, y = bin, color = order, size = sum), alpha = 0.9) +
  scale_size(range = c(1, 10), name = "Number of CAZy Annotations") +
  theme_classic(base_size=20) +
  theme(axis.text.x =element_text(angle = 25, hjust = 1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  #scale_color_manual(values = palorg) +
  scale_color_manual(values = palorgcorynshingo) +
  labs(y = "MAGs", color = "Order") +
  ggtitle("CAZy Annotation of MAGs")
mag_plot

ggsave("mag_bubble.tiff", device = "tiff", dpi = 700)
ggsave("mag_bubble.png", device = "png", dpi = 700)
ggsave("mag_bubble.pdf", device = "pdf", dpi = 700)
