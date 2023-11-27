library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggtext)

data <- read.csv("CR_dbcan_edger_time_treatments.csv")
data <- filter(data, PValue<0.05)
data <- na.omit(data)
data <- subset(data, Substrate != "Other")

data$Treatment <- factor(data$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))

data$order <- factor(data$order, levels = c("Actinomycetales","Burkholderiales","Caulobacterales","Corynebacteriales","Enterobacterales","Flavobacteriales",
                                            "Lactobacillales","Micrococcales","Propionibacteriales","Pseudomonadales","Rhizobiales",
                                            "Sphingobacteriales","Xanthomonadales","no support"))
data$Substrate <- factor(data$Substrate, levels = c("Cellulose", "Cellobiose", "Xylan",
                                                                "Starch", "Trehalose","Pectin","Lignin",
                                                                "Peptidoglycan","Chitin","Other"))

time.labels <- c("Day 1 to 2", "Day 2 to 3", "Day 3 to 4", "Day 4 to 5", "Day 5 to 6", "Day 6 to 7")
names(time.labels) <- c("1_to_2", "2_to_3", "3_to_4", "4_to_5", "5_to_6", "6_to_7")

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
palette <-c( "#56B4E9","#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")


palsubstrate <- c("#73A87C","#CE9486","#205D89","#CF784B","#FFEA59","#B7DBDB","#B5B867","#77674E","#292176")
palday <- c("#7DCEA0","#52BE80","#27AE60","#229954","#1E8449","#196F3D","#145A32")


##################### Start with Cellulose ################


data_Cellulose <- subset(data, Substrate == "Cellulose")
data_Cellulose <- subset(data_Cellulose, Time %in% c("1_to_2", "3_to_4", "5_to_6"))

#replace GH5 subfamily with GH5 only, since I don't see obvious difference when I 




data_Cellulose$CAZy <- factor(data_Cellulose$CAZy, levels=c("GH5","GH5_1","GH5_2","GH5_4","GH5_5","GH5_25","GH5_26","GH5_38","GH5_39","GH5_46","GH6","GH8","GH9","GH44","GH74","AA10"))

palcellu <- c("#B57EDC", "#A1CAF1","#6495ED","#318CE7","#007FFF","#1560BD","#0047AB","#003366","#002147",
              "#AF4E24","#E54E21","#FDD262","#6C8645","#54D8B1","#9C964A")

palorg_Cellulose <- c(
                      "#FFFFB3", #Burkholderiales
                      "#BEBADA", #Caulobacterales
                      "#80B1D3", #Enterobacterales
                      "#FDB462", #Flavobacteriales
                      "#B3DE69", #Lactobacillales
                      "#FCCDE5", #Micrococcales
                      "#BC80BD", #Propionibacteriales
                      "#e5c494", #Rhizobiales
                      "#1f78b4", #Sphingobacteriales
                      "#D9D9D9" #no support)
)

palorg_Cellulose2 <- c(
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
  "#D9D9D9" #no support)
)

p_Cellulose_order <- ggplot(data_Cellulose %>%
                    arrange(order),
                  aes(x=logCPM, y=logFC, color = order)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palorg_Cellulose) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Order", ncol = 3)) +
  ggtitle("Differential Expression of Cellulose CAZys, Colored by Order") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))

p_Cellulose_order

p_Cellulose_cazy <- ggplot(data_Cellulose %>%
                             arrange(order),
                           aes(x=logCPM, y=logFC, color = CAZy)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palcellu) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("CAZy")) +
  ggtitle("Differential Expression of Cellulose CAZys, Colored by CAZy") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))


p_Cellulose_cazy


legend_Cellulose <- plot_grid(get_legend(p_Cellulose_order), get_legend(p_Cellulose_cazy), ncol = 2)

p1 <- p_Cellulose_order + theme(legend.position="none")
p2 <- p_Cellulose_cazy + theme(legend.position="none")

total_Cellulose <- plot_grid(p1,p2, align = "v", ncol = 2, axis = "tb", rel_heights = c(10, 10))
plot_grid(total_Cellulose, legend_Cellulose, ncol = 1, rel_heights = c(15,3))

#ggsave("cellulose_de.png")
#ggsave("cellulose_de.pdf")


#################### Xylan ###########################

data_Xylan <- subset(data, Substrate == "Xylan")
data_Xylan$cazy_type <- substr(data_Xylan$CAZy, 1,2)

data_Xylan$CAZy <- factor(data_Xylan$CAZy, levels=c("CE1",  "CE2",  "CE3",  "CE4",  "CE5",  "CE6",  "CE7","CE12", "CE16", "GH10","GH5_22","GH30","GH30_2","GH30_3","GH30_7","GH30_8","GH67","GH115","GH120"))


data_Xylan  <- subset(data_Xylan , Time %in% c("1_to_2", "3_to_4", "5_to_6"))

data_Xylan_gh <- subset(data_Xylan , cazy_type == "GH")
data_Xylan_ce <- subset(data_Xylan , cazy_type == "CE")

palce <- c("#89B151","#FFEA59","#B7DBDB","#D9B253","#6CA18F","#CE9486","#335B8E","#C9573C","#999999")
palgh <- c("#E69EA1","#DED9A2","#F6C83E","#4D5B28","#DB4372","#B77E60","#5785C1","#8CBEA3", "#a6cee3", "#999999")

palorg_Xylan <- c(
  "#FFFFB3", #Burkholderiales
  "#BEBADA", #Caulobacterales
  "#80B1D3", #Enterobacterales
  "#B3DE69", #Lactobacillales
  "#FCCDE5", #Micrococcales
  "#e5c494", #Rhizobiales
  "#1f78b4" #Sphingobacteriales
)


palorg_Xylan2 <- c(
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


palorg_Xylan3 <- c(
  "#FFFFB3", #Burkholderiales
  "#BEBADA", #Caulobacterales
  "#FB8072", #Corynebacteriales
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


palorg_Xylan4 <- c(
  "#BEBADA", #Caulobacterales
  "#B3DE69", #Lactobacillales
  "#FCCDE5", #Micrococcales
  "#1f78b4" #Sphingobacteriales
)





############ CE ##################

p_Xylan_order_ce <- ggplot(data_Xylan_ce %>%
                          arrange(order),
                        aes(x=logCPM, y=logFC, color = order)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palorg_Xylan2) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Order", ncol = 3)) +
  ggtitle("Differential Expression of Xylan CE CAZys, Colored by Order") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))

p_Xylan_order_ce


p_Xylan_cazy_ce <- ggplot(data_Xylan_ce %>%
                         arrange(order),
                       aes(x=logCPM, y=logFC, color = CAZy)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palce) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("CAZy")) +
  ggtitle("Differential Expression of Xylan CE CAZys, Colored by CAZy") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))


p_Xylan_cazy_ce


legend_Xylan_ce <- plot_grid(get_legend(p_Xylan_order_ce), get_legend(p_Xylan_cazy_ce), ncol = 2)

p1 <- p_Xylan_order_ce + theme(legend.position="none")
p2 <- p_Xylan_cazy_ce + theme(legend.position="none")

total_Xylan_ce <- plot_grid(p1,p2, align = "v", ncol = 2, axis = "tb", rel_heights = c(10, 10))
plot_grid(total_Xylan_ce, legend_Xylan_ce, ncol = 1, rel_heights = c(15,3))

#ggsave("xylan_ce_de.png")
#ggsave("xylan_ce_de.pdf")


#################### GH ################



p_Xylan_order_gh <- ggplot(data_Xylan_gh %>%
                             arrange(order),
                           aes(x=logCPM, y=logFC, color = order)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palorg_Xylan) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Order", ncol = 3)) +
  ggtitle("Differential Expression of Xylan GH CAZys, Colored by Order") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))

p_Xylan_order_gh


p_Xylan_cazy_gh <- ggplot(data_Xylan_gh %>%
                            arrange(order),
                          aes(x=logCPM, y=logFC, color = CAZy)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palgh) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("CAZy")) +
  ggtitle("Differential Expression of Xylan GH CAZys, Colored by CAZy") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))


p_Xylan_cazy_gh

legend_Xylan_gh <- plot_grid(get_legend(p_Xylan_order_gh), get_legend(p_Xylan_cazy_gh), ncol = 2)

p1 <- p_Xylan_order_gh + theme(legend.position="none")
p2 <- p_Xylan_cazy_gh + theme(legend.position="none")

total_Xylan_gh <- plot_grid(p1,p2, align = "v", ncol = 2, axis = "tb", rel_heights = c(10, 10))
plot_grid(total_Xylan_gh, legend_Xylan_gh, ncol = 1, rel_heights = c(15,3))

#ggsave("xylan_gh_de.png")
#ggsave("xylan_gh_de.pdf")





#################### Lignin ###########################

data_Lignin <- subset(data, Substrate == "Lignin")
data_Lignin$CAZy <- factor(data_Lignin$CAZy, levels=c("AA1","AA1_1","AA1_2","AA1_3","AA2","AA3","AA3_1","AA3_2","AA3_3","AA3_4","AA4","AA5","AA5_1","AA5_2","AA6","AA7","AA12"))

data_Lignin <- subset(data_Lignin, Time %in% c("1_to_2", "3_to_4", "5_to_6"))
data_Lignin_aa3 <- subset(data_Lignin, CAZy %in% c("AA3", "AA3_1", "AA3_2", "AA3_3", "AA3_4"))
data_Lignin <- subset(data_Lignin, ! CAZy %in% c("AA3", "AA3_1", "AA3_2", "AA3_3", "AA3_4"))


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


palorg_Lignin <- c(
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


p_Lignin_order <- ggplot(data_Lignin %>%
                           arrange(order),
                         aes(x=logCPM, y=logFC, color = order)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palorg_Lignin) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  theme(axis.text = element_text(size = 10)) +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Order", ncol = 3)) +
  ggtitle("Differential Expression of Lignin CAZys, Colored by Order") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))

p_Lignin_order


p_Lignin_cazy <- ggplot(data_Lignin %>%
                          arrange(order),
                        aes(x=logCPM, y=logFC, color = CAZy)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=pallignin1 ) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("CAZy")) +
  ggtitle("Differential Expression of Lignin CAZys, Colored by CAZy") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))


p_Lignin_cazy


legend_Lignin <- plot_grid(get_legend(p_Lignin_order), get_legend(p_Lignin_cazy), ncol = 2)

p1 <- p_Lignin_order + theme(legend.position="none")
p2 <- p_Lignin_cazy + theme(legend.position="none")

total_Lignin <- plot_grid(p1,p2, align = "v", ncol = 2, axis = "tb", rel_heights = c(10, 10))
plot_grid(total_Lignin, legend_Lignin, ncol = 1, rel_heights = c(15,3))

#ggsave("lignin_de.png")
#ggsave("lignin_de.pdf")





### AA 3

p_Lignin_order_aa3 <- ggplot(data_Lignin_aa3 %>%
                           arrange(order),
                         aes(x=logCPM, y=logFC, color = order)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palorg_Lignin) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  theme(axis.text = element_text(size = 10)) +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Order", ncol = 3)) +
  ggtitle("Differential Expression of AA3 CAZys, Colored by Order") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))

p_Lignin_order_aa3


p_Lignin_cazy_aa3 <- ggplot(data_Lignin_aa3 %>%
                          arrange(order),
                        aes(x=logCPM, y=logFC, color = CAZy)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=pallignin2 ) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom") +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("CAZy")) +
  ggtitle("Differential Expression of AA3 CAZys, Colored by CAZy") +
  facet_grid(Treatment ~ Time, labeller = labeller(Time = time.labels))


p_Lignin_cazy_aa3

legend_Lignin_aa3 <- plot_grid(get_legend(p_Lignin_order_aa3), get_legend(p_Lignin_cazy_aa3), ncol = 2)

p1 <- p_Lignin_order_aa3 + theme(legend.position="none")
p2 <- p_Lignin_cazy_aa3 + theme(legend.position="none")

total_Lignin_aa3 <- plot_grid(p1,p2, align = "v", ncol = 2, axis = "tb", rel_heights = c(10, 10))
plot_grid(total_Lignin_aa3, legend_Lignin_aa3, ncol = 1, rel_heights = c(15,3))


#ggsave("lignin_aa3_de.png")
#ggsave("lignin_aa3_de.pdf")


########### MAIN FIGURE FOCUS ON XYLAN, LIGNIN DAY 1 TO 2 PH6 #######

pallignin3 <- c("#fa9fb5","#f768a1","#dd3497","#ae017e", #AA1, pink
                "#335B8E", #AA2, blue
                #### no AA3 ###
                "#984ea3", #AA4, purple
                "#fdae6b","#8c2d04", #AA5, orange
                "#DFBA47", #AA6, yellow
                "#6C8645", #AA7, green
                "#e41a1c" #AA12, red
)




data_Xylan_main <- subset(data_Xylan, Treatment %in% c("pH6") & Time == "1_to_2" & cazy_type == "CE")

data_Lignin_main <- subset(data_Lignin, Treatment %in% c("pH6") & Time == "1_to_2")

cazy_labels <- c("CE" = "Xylan Carbohydrate Esterase", "GH" = "Glycoside Hydrolase")



p_main_Xylan_order <- ggplot(data_Xylan_main %>%
                               arrange(order),
                             aes(x=logCPM, y=logFC, color = order)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palorg_Xylan3) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Order")) +
  ggtitle("Colored by Order") +
  facet_grid( ~ cazy_type, labeller = labeller(cazy_type = cazy_labels))
p_main_Xylan_order

p_main_Xylan_ce <- ggplot(data_Xylan_main %>%
                                  arrange(order),
                                aes(x=logCPM, y=logFC, color = CAZy)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palce) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Xylan CE CAZy")) +
  ggtitle("Colored by CAZy") +
  facet_grid( ~ cazy_type, labeller = labeller(cazy_type = cazy_labels))

p_main_Xylan_ce



p_main_Lignin_order <- ggplot(data_Lignin_main %>%
                                arrange(order),
                              aes(x=logCPM, y=logFC, color = order)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=palorg_Lignin) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Order")) +
  ggtitle("Colored by Order") +
  facet_grid( ~ Substrate)
p_main_Lignin_order

p_main_Lignin_cazy <- ggplot(data_Lignin_main %>%
                               arrange(order),
                             aes(x=logCPM, y=logFC, color = CAZy)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed") +
  scale_color_manual(values=pallignin3) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  labs(y="Log Fold Change", x = "Average Log CPM (Counts per Million)") +
  guides(color=guide_legend("Lignin CAZy")) +
  ggtitle("Colored by CAZy") +
  facet_grid( ~ Substrate)

p_main_Lignin_cazy


title <- ggdraw() +
  draw_label(
    "Differential Expression of Selected CAZys, Day 1 to 2, pH 6 Treatment",
    x = 0,
    hjust = 0
  )


#l1 <-plot_grid(get_legend(p_main_Lignin_order))
#l2 <- plot_grid(get_legend(p_main_Xylan_ce))
#l3 <- plot_grid(get_legend(p_main_Lignin_cazy))

#l1
#ggsave("l1.png", dpi = 700)
#l2
#ggsave("l2.png", dpi = 700)
#l3
#ggsave("l3.png", dpi = 700)


legend_main <- plot_grid(get_legend(p_main_Lignin_order), get_legend(p_main_Xylan_ce),get_legend(p_main_Lignin_cazy), ncol = 1)

legend_main


p1 <- p_main_Xylan_order + theme(legend.position="none")
p2 <- p_main_Lignin_order + theme(legend.position="none")
p3 <- p_main_Xylan_ce + theme(legend.position="none")
p4 <- p_main_Lignin_cazy + theme(legend.position="none")

main_org <- plot_grid(p1, p2, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("A", "B"))
main_cazy <- plot_grid(p3, p4, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("C", "D"))
main_de <- plot_grid(main_org,main_cazy, align = "v", ncol = 1, axis = "tb", rel_heights = c(10, 10))
main_de_legend <- plot_grid(main_de, legend_main, nrow = 1, rel_widths = c(10,3))

main_de_legend
title <- ggdraw() +
  draw_label(
    "Differential Expression of Selected CAZys, Day 1 to 2, pH 6 Treatment",
    size = 20,
    x = 0,
    hjust = 0
  )




plot_grid(title, main_de_legend, ncol = 1, rel_heights = c(0.1, 1))

#ggsave("cr_main_de.png")
#ggsave("cr_main_de.pdf")

############### quanitfy DE ########

xy_temp <- data_Xylan_main[,-c(20)]
lig_temp <- data_Lignin_main

de_dat <- rbind(xy_temp, lig_temp)

de_dat <- de_dat %>%
  mutate(change = case_when(
    logFC > 0 ~ "Positive",
    logFC < 0 ~ "Negative"
  ))

de_dat$change <- factor(de_dat$change, levels=c("Positive", "Negative"))

de_dat_sum <- de_dat %>%
  group_by(Substrate, change) %>%
  summarise(sum = n())

de_dat_sum$subchange <- paste(de_dat_sum$Substrate, de_dat_sum$change, sep = "_")
de_dat_sum <- de_dat_sum[,-c(1,2)]

de_dat <- de_dat %>%
  group_by(CAZy, Substrate, order, change) %>%
  summarise(count = n())

de_dat$subchange <- paste(de_dat$Substrate, de_dat$change, sep = "_")

de_dat <- merge(de_dat, de_dat_sum, by = "subchange")
de_dat$percent <- (de_dat$count/de_dat$sum)*100



p_de <- ggplot(de_dat, aes(x=change, y = percent, fill = order)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=palorg_Lignin) +
  labs(y="Percent (%)", x = "Fold Change Direction") +
  guides(color=guide_legend("Order")) +
  ggtitle("Proportion of DE CAZys from Day 1 to Day 2 in Treatment pH 6") +
  facet_grid(~Substrate)
p_de



de_dat_xy <- subset(de_dat, Substrate == "Xylan")
de_dat_lig <- subset(de_dat, Substrate == "Lignin")

de_dat_xy_summary <- de_dat_xy %>%
  group_by(CAZy, change) %>%
  summarise(sumper=sum(percent),
            sumnum=sum(count),
            sd=sd(percent))


p_de_xy <- ggplot(de_dat_xy, aes(x=change, y = count, fill = CAZy)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=palce) +
  labs(y="Count", x = "Fold Change Direction") +
  guides(color=guide_legend("Xylan CE CAZy")) +
  ggtitle("Count of Xylan CE CAZys") +
  facet_grid(~Substrate)
p_de_xy


p_de_bar_xy <- ggplot(de_dat_xy, aes(x=CAZy, y = count, fill = order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  scale_fill_manual(values=palorg_Xylan3) +
  labs(y="Count", x = "CAZy") +
  guides(fill=guide_legend("Order")) +
  ggtitle("Count of Xylan CE CAZys") +
  facet_grid(~change)
p_de_bar_xy



p_de_lig <- ggplot(de_dat_lig, aes(x=change, y = count, fill = CAZy)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=pallignin3) +
  labs(y="Count", x = "Fold Change Direction") +
  guides(color=guide_legend("Lignin CAZy")) +
  ggtitle("Count of Lignin CAZys") +
  facet_grid(~Substrate)
p_de_lig

p_de_bar_lig <- ggplot(de_dat_lig, aes(x=CAZy, y = count, fill = order)) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  scale_fill_manual(values=palorg_Lignin) +
  labs(y="Count", x = "CAZy") +
  guides(fill=guide_legend("Order")) +
  ggtitle("Count of Lignin CAZys") +
  facet_grid(~change)
p_de_bar_lig

p_de_xy
p_de_bar_xy
p_de_lig
p_de_bar_lig


legend_de <- plot_grid(get_legend(p_de_bar_lig), get_legend(p_de_xy),get_legend(p_de_lig), ncol = 1)

p1 <- p_de_xy + theme(legend.position="none")
p2 <- p_de_lig + theme(legend.position="none")
p3 <- p_de_bar_xy + theme(legend.position="none")
p4 <- p_de_bar_lig + theme(legend.position="none")

de_org <- plot_grid(p1, p2, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("A", "B"))
de_cazy <- plot_grid(p3, p4, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("C", "D"))
de_de <- plot_grid(de_org,de_cazy, align = "v", ncol = 1, axis = "tb", rel_heights = c(10, 10))
de_de_legend <- plot_grid(de_de, legend_de, nrow = 1, rel_widths = c(10,3))

#de_de_legend
title <- ggdraw() +
  draw_label(
    "Count of Differentially Expressed CAZys, Day 1 to 2, pH 6 Treatment",
    size = 20,
    x = 0,
    hjust = 0
  )




plot_grid(title, de_de_legend, ncol = 1, rel_heights = c(0.1, 1))

#ggsave("cr_de_de.png")
#ggsave("cr_de_de.pdf")


