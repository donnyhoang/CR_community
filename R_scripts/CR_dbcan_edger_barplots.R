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


######################################### de every thing #####################


de_dat <- data %>%
  mutate(change = case_when(
    logFC > 0 ~ "Positive",
    logFC < 0 ~ "Negative"
  ))

de_dat$change <- factor(de_dat$change, levels=c("Positive", "Negative"))

de_dat_sum <- de_dat %>%
  group_by(Substrate, change, Time, Treatment) %>%
  summarise(sum = n())

de_dat_sum$subchange <- paste(de_dat_sum$Treatment, de_dat_sum$Time, de_dat_sum$Substrate, de_dat_sum$change, sep = "_")

de_dat_sum_neg <- subset(de_dat_sum, change == "Negative")
de_dat_sum <- subset(de_dat_sum, change != "Negative")
de_dat_sum_neg$sum <- de_dat_sum_neg$sum*(-1)

de_dat_sum <- rbind(de_dat_sum, de_dat_sum_neg)
de_dat_sum <- de_dat_sum[,-c(1:4)]

de_dat <- de_dat %>%
  group_by(CAZy, Substrate, order, change, Treatment, Time, logFC) %>%
  summarise(count = n())

de_dat$subchange <- paste(de_dat$Treatment, de_dat$Time, de_dat$Substrate, de_dat$change, sep = "_")

de_dat <- merge(de_dat, de_dat_sum, by = "subchange")

de_dat_neg <- subset(de_dat, change == "Negative")
de_dat <- subset(de_dat, change != "Negative")
de_dat_neg$count <- de_dat_neg$count*(-1)
de_dat <- rbind(de_dat_neg, de_dat)
de_dat$abscount <- abs(de_dat$count)
de_dat$percent <- (de_dat$count/de_dat$sum)*100
de_dat$Time <- factor(de_dat$Time, levels = rev(unique(de_dat$Time)))



############# Cellulose

palcellu <- c("#B57EDC", "#A1CAF1","#6495ED","#318CE7","#007FFF","#1560BD","#0047AB","#003366","#002147",
              "#AF4E24","#E54E21","#FDD262","#6C8645","#54D8B1","#9C964A")
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



de_dat_Cellulose <- subset(de_dat, Substrate == "Cellulose")

de_dat_Cellulose$CAZy <- factor(de_dat_Cellulose$CAZy, levels=c("GH5","GH5_1","GH5_2","GH5_4","GH5_5","GH5_25","GH5_26","GH5_38","GH5_39","GH5_46","GH6","GH8","GH9","GH44","GH74","AA10"))



p_de_Cellulose_org <- ggplot(de_dat_Cellulose, aes(x = Time, y = count, fill = order)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed" ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = palorg_Cellulose2) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Order")) +
  #ggtitle("Counts of Differentially Expressed Cellulose CAZys") +
  facet_wrap(~Treatment)
p_de_Cellulose_org


p_de_Cellulose <- ggplot(de_dat_Cellulose, aes(x = Time, y = count, fill = CAZy)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed" ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = palcellu) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Cellulose CAZy")) +
  #ggtitle("Counts of Differentially Expressed Cellulose CAZys") +
  facet_wrap(~Treatment)
p_de_Cellulose

legend_Cellulose <- plot_grid(get_legend(p_de_Cellulose_org), get_legend(p_de_Cellulose), ncol = 2)

p1 <- p_de_Cellulose_org + theme(legend.position="none")
p2 <- p_de_Cellulose + theme(legend.position="none")

plots <- plot_grid(p1, p2, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("A", "B"))
plot_legend <- plot_grid(plots, legend_Cellulose, ncol=1, rel_heights = c(10,3)) 
title <- ggdraw() +
  draw_label(
    "Count of Differentially Expressed Cellulose CAZys",
    size = 20,
    x = 0.01,
    hjust = 0
  )

plot_grid(title, plot_legend, ncol = 1, rel_heights = c(0.1, 1))


#ggsave("de_counts_Cellulose.png")
#ggsave("de_counts_Cellulose.pdf")

############### Xylan

palce <- c("#89B151","#FFEA59","#B7DBDB","#D9B253","#6CA18F","#CE9486","#335B8E","#C9573C","#999999")
palgh <- c("#E69EA1","#DED9A2","#F6C83E","#4D5B28","#DB4372","#B77E60","#5785C1","#8CBEA3", "#a6cee3", "#999999")

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



palorg_Xylan <- c(
  "#FFFFB3", #Burkholderiales
  "#BEBADA", #Caulobacterales
  "#80B1D3", #Enterobacterales
  "#B3DE69", #Lactobacillales
  "#FCCDE5", #Micrococcales
  "#e5c494", #Rhizobiales
  "#1f78b4", #Sphingobacteriales
  "#D9D9D9" #no support
)

de_dat_Xylan <- subset(de_dat, Substrate == "Xylan")
de_dat_Xylan$CAZy <- factor(de_dat_Xylan$CAZy, levels=c("CE1",  "CE2",  "CE3",  "CE4",  "CE5",  "CE6",  "CE7","CE12", "CE16", "GH10","GH5_22","GH30","GH30_2","GH30_3","GH30_7","GH30_8","GH67","GH115","GH120"))
de_dat_Xylan$cazy_type <- substr(de_dat_Xylan$CAZy, 1,2)
de_dat_Xylan_gh <- subset(de_dat_Xylan , cazy_type == "GH")
de_dat_Xylan_ce <- subset(de_dat_Xylan , cazy_type == "CE")



p_de_Xylan_org_ce <- ggplot(de_dat_Xylan_ce, aes(x = Time, y = count, fill = order)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed" ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = palorg_Xylan2) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Order")) +
  ggtitle("Counts of Differentially Expressed Xylan CE CAZys")+
  facet_wrap(~Treatment)
p_de_Xylan_org_ce


p_de_Xylan_ce <- ggplot(de_dat_Xylan_ce, aes(x = Time, y = count, fill = CAZy)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed" ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = palce) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Xylan CAZy")) +
  ggtitle("Counts of Differentially Expressed CE Xylan CAZys")+
  facet_wrap(~Treatment)
p_de_Xylan_ce



p_de_Xylan_org_gh <- ggplot(de_dat_Xylan_gh, aes(x = Time, y = count, fill = order)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed" ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = palorg_Xylan) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Order")) +
  ggtitle("Counts of Differentially Expressed Xylan GH CAZys")+
  facet_wrap(~Treatment)
p_de_Xylan_org_gh


p_de_Xylan_gh <- ggplot(de_dat_Xylan_gh, aes(x = Time, y = count, fill = CAZy)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed" ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = palgh) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Xylan CAZy")) +
  ggtitle("Counts of Differentially Expressed Xylan GH CAZys")+
  facet_wrap(~Treatment)
p_de_Xylan_gh



legend_Xylan_ce <- plot_grid(get_legend(p_de_Xylan_org_ce), get_legend(p_de_Xylan_ce), ncol = 2)

p1 <- p_de_Xylan_org_ce + theme(legend.position="none")
p2 <- p_de_Xylan_ce + theme(legend.position="none")

plots <- plot_grid(p1, p2, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("A", "B"))
plot_legend <- plot_grid(plots, legend_Xylan_ce, ncol=1, rel_heights = c(10,3)) 
title <- ggdraw() +
  draw_label(
    "Count of Differentially Expressed Xylan CE CAZys",
    size = 20,
    x = 0.01,
    hjust = 0
  )


plot_grid(title, plot_legend, ncol = 1, rel_heights = c(0.1, 1))

#ggsave("de_counts_Xylan_ce.png")
#ggsave("de_counts_Xylan_ce.pdf")


legend_Xylan_gh <- plot_grid(get_legend(p_de_Xylan_org_gh), get_legend(p_de_Xylan_gh), ncol = 2)

p1 <- p_de_Xylan_org_gh + theme(legend.position="none")
p2 <- p_de_Xylan_gh + theme(legend.position="none")

plots <- plot_grid(p1, p2, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("A", "B"))
plot_legend <- plot_grid(plots, legend_Xylan_gh, ncol=1, rel_heights = c(10,3)) 
title <- ggdraw() +
  draw_label(
    "Count of Differentially Expressed Xylan GH CAZys",
    size = 20,
    x = 0.01,
    hjust = 0
  )

plot_grid(title, plot_legend, ncol = 1, rel_heights = c(0.1, 1))

#ggsave("de_counts_Xylan_gh.png")
#ggsave("de_counts_Xylan_gh.pdf")


################################ LIGNIN ########################################

de_dat_Lignin <- subset(de_dat, Substrate == "Lignin")
de_dat_Lignin$CAZy <- factor(de_dat_Lignin$CAZy, levels=c("AA1","AA1_1","AA1_2","AA1_3","AA2","AA3","AA3_1","AA3_2","AA3_3","AA3_4","AA4","AA5","AA5_1","AA5_2","AA6","AA7","AA12"))
de_dat_Lignin_aa3 <- subset(de_dat_Lignin, CAZy %in% c("AA3", "AA3_1", "AA3_2", "AA3_3", "AA3_4"))
de_dat_Lignin <- subset(de_dat_Lignin, ! CAZy %in% c("AA3", "AA3_1", "AA3_2", "AA3_3", "AA3_4"))


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

p_de_Lignin_org <- ggplot(de_dat_Lignin, aes(x = Time, y = count, fill = order)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed" ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = palorg_Lignin) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Order")) +
  ggtitle("Counts of Differentially Expressed Lignin CAZys")+
  facet_wrap(~Treatment)
p_de_Lignin_org


p_de_Lignin <- ggplot(de_dat_Lignin, aes(x = Time, y = count, fill = CAZy)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = pallignin1) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Lignin CAZy")) +
  ggtitle("Counts of Differentially Expressed Lignin CAZys")+
  facet_wrap(~Treatment)
p_de_Lignin

p_de_Lignin_org_aa3 <- ggplot(de_dat_Lignin_aa3, aes(x = Time, y = count, fill = order)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = palorg_Lignin) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("Order")) +
  ggtitle("Counts of Differentially Expressed AA3 CAZys")+
  facet_wrap(~Treatment)
p_de_Lignin_org_aa3


p_de_Lignin_aa3 <- ggplot(de_dat_Lignin_aa3, aes(x = Time, y = count, fill = CAZy)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom") +
  coord_flip() +
  scale_x_discrete(labels = time.labels) + 
  scale_fill_manual(values = pallignin2) +
  labs(y="Counts", x = "Time") +
  guides(fill=guide_legend("AA3 CAZy")) +
  ggtitle("Counts of Differentially Expressed AA3 CAZys")+
  facet_wrap(~Treatment)
p_de_Lignin_aa3



legend_Lignin <- plot_grid(get_legend(p_de_Lignin_org), get_legend(p_de_Lignin), ncol = 2)

p1 <- p_de_Lignin_org + theme(legend.position="none")
p2 <- p_de_Lignin + theme(legend.position="none")

plots <- plot_grid(p1, p2, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("A", "B"))
plot_legend <- plot_grid(plots, legend_Lignin, ncol=1, rel_heights = c(10,3)) 
title <- ggdraw() +
  draw_label(
    "Count of Differentially Expressed Lignin CAZys",
    size = 20,
    x = 0.01,
    hjust = 0
  )

plot_grid(title, plot_legend, ncol = 1, rel_heights = c(0.1, 1))


#ggsave("de_counts_Lignin.png")
#ggsave("de_counts_Lignin.pdf")



legend_Lignin_aa3 <- plot_grid(get_legend(p_de_Lignin_org_aa3), get_legend(p_de_Lignin_aa3), ncol = 2)

p1 <- p_de_Lignin_org_aa3 + theme(legend.position="none")
p2 <- p_de_Lignin_aa3 + theme(legend.position="none")

plots <- plot_grid(p1, p2, align = "h", ncol = 2, axis = "tb", rel_heights = c(10,10), labels = c("A", "B"))
plot_legend <- plot_grid(plots, legend_Lignin_aa3, ncol=1, rel_heights = c(10,3)) 
title <- ggdraw() +
  draw_label(
    "Count of Differentially Expressed AA3 CAZys",
    size = 20,
    x = 0.01,
    hjust = 0
  )

plot_grid(title, plot_legend, ncol = 1, rel_heights = c(0.1, 1))


#ggsave("de_counts_Lignin_aa3.png")
#ggsave("de_counts_Lignin_aa3.pdf")


################### MAIN FIGURE #############################


p_bar_Cellulose_cazy <- ggplot(de_dat_Cellulose, aes(x = change, y = abscount, fill = CAZy)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palcellu) +
  facet_grid(~Substrate)
p_bar_Cellulose_cazy

p_bar_Cellulose <- ggplot(de_dat_Cellulose, aes(x = CAZy, y = abscount, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palorg_Cellulose2) +
  facet_grid(~change)
p_bar_Cellulose

