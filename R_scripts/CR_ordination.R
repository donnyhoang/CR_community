library(BiodiversityR) #this loads vegan too
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyverse)
library(ade4)
library(usedist)
library(ggpubr)
library(purrr)
library(cowplot)
# https://github.com/kylebittinger/usedist
#reference https://www.rpubs.com/RGrieger/545184

#ordinate organism abundance, expression, and dbcan expression

ph <- read.csv("CR_pH_measurements.csv")
ph$variable <- paste(ph$Treatment, ph$Day, ph$Replicate, sep = "_")

palorg <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#e5c494","#1f78b4","#fb9a99","#D9D9D9")
palette <-c( "#56B4E9","#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")


############ MAG Abundance, MG ############
otu_mg <- as.data.frame(t(read.csv("CR_matrix_tax_mg.csv", header = TRUE, row.names = 1)))

otu_mg <- as.data.frame(t(otu_mg))

meta_mg <- as.data.frame(row.names(otu_mg))
colnames(meta_mg)[1] <- "variable"

meta_mg <- merge(meta_mg,ph , by="variable")

ordination.model_mg <- metaMDS(otu_mg, distance='bray', k=3)
data.bin.fit_mg <- envfit(ordination.model_mg, otu_mg, permutations = 999)
data.env.fit_mg <- envfit(ordination.model_mg, meta_mg, permutations = 999, na.rm = TRUE)

#write out NMDS coords
data_mg <- as.data.frame(scores(ordination.model_mg, display = "sites"))
data_mg <- cbind(data_mg, meta_mg)
data_mg$Treatment<-factor(data_mg$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))

#write out bin data
data_bin_mg <- as.data.frame(scores(data.bin.fit_mg, display = "vectors"))
data_bin_mg <- cbind(data_bin_mg, bin = rownames(data_bin_mg))
data_bin_mg <- cbind(data_bin_mg, pval = data.bin.fit_mg$vectors$pvals)
#I want order name for bins
tax <- read.csv("CR_bin_tax_cluster.csv", header = TRUE)
data_bin_mg <- merge(data_bin_mg, tax, by = "bin")
sig_data_bin_mg <-subset(data_bin_mg, pval <= 0.05)

#write out env data
data_env_mg <- as.data.frame(scores(data.env.fit_mg, display = "vectors"))
data_env_mg <- cbind(data_env_mg, env_variables = rownames(data_env_mg))
data_env_mg <- cbind(data_env_mg, pval = data.env.fit_mg$vectors$pvals)

##### plot

nmds_mg <- ggplot(data_mg, aes(x=NMDS1, y = NMDS2)) +
  geom_point(aes(x=NMDS1, y = NMDS2, color = Treatment),size = 3) +
  #coord_fixed() +
  theme_classic(base_size =15) +
  scale_color_manual(values=palette)
nmds_mg

#add bin vectors, label with order
nmds_mg +
  geom_segment(data=sig_data_bin_mg, aes (x = 0, xend=NMDS1, y = 0, yend=NMDS2), arrow = arrow(length=unit(0.25, "cm")), color = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data=sig_data_bin_mg, aes(x=NMDS1, y = NMDS2, label = order), cex = 3, direction = "both", segment.size=0.25) +
  ggtitle("Ordination of Sample MAG Abundance \n Vectors Labeled by MAG's Order Classification")

#add env vector

nmds_mg <- nmds_mg +
  geom_segment(data=data_env_mg, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow=arrow(length = unit(0.25, "cm")), color = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data=data_env_mg, aes(x=NMDS1, y = NMDS2, label = env_variables), cex = 5, direction = "both", segment.size=0.25) +
  ggtitle("Ordination of Sample MAG Abundance")

############ MAG expression, MT ############
otu_mt <- as.data.frame(t(read.csv("CR_matrix_tax_mt_unfiltered.csv", header = TRUE, row.names = 1)))

meta_mt <- as.data.frame(row.names(otu_mt))
colnames(meta_mt)[1] <- "variable"

meta_mt <- merge(meta_mt,ph , by="variable")

ordination.model <- metaMDS(otu_mt, distance='bray', k=3)
data.bin.fit <- envfit(ordination.model, otu_mt, permutations = 999)
data.env.fit <- envfit(ordination.model, meta_mt, permutations = 999, na.rm = TRUE)

#write out NMDS coords
data <- as.data.frame(scores(ordination.model, display = "sites"))
data <- cbind(data, meta_mt)
data$Treatment<-factor(data$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))

#write out bin data
data_bin <- as.data.frame(scores(data.bin.fit, display = "vectors"))
data_bin <- cbind(data_bin, bin = rownames(data_bin))
data_bin <- cbind(data_bin, pval = data.bin.fit$vectors$pvals)
#I want order name for bins
tax <- read.csv("CR_bin_tax_cluster.csv", header = TRUE)
data_bin <- merge(data_bin, tax, by = "bin")
sig_data_bin <-subset(data_bin, pval <= 0.05)

#write out env data
data_env <- as.data.frame(scores(data.env.fit, display = "vectors"))
data_env <- cbind(data_env, env_variables = rownames(data_env))
data_env <- cbind(data_env, pval = data.env.fit$vectors$pvals)

##### plot

nmds_mt <- ggplot(data, aes(x=NMDS1, y = NMDS2)) +
  geom_point(aes(x=NMDS1, y = NMDS2, color = Treatment),size = 3) +
  #coord_fixed() +
  theme_classic(base_size =15) +
  scale_color_manual(values=palette)
nmds_mt

#add bin vectors, label with order
nmds_mt +
  geom_segment(data=sig_data_bin, aes (x = 0, xend=NMDS1, y = 0, yend=NMDS2), arrow = arrow(length=unit(0.25, "cm")), color = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data=sig_data_bin, aes(x=NMDS1, y = NMDS2, label = order), cex = 3, direction = "both", segment.size=0.25) +
  ggtitle("Ordination of Sample MAG Expression \n Vectors Labeled by MAG's Order Classification")

#add env vector

nmds_mt <- nmds_mt +
  geom_segment(data=data_env, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow=arrow(length = unit(0.25, "cm")), color = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data=data_env, aes(x=NMDS1, y = NMDS2, label = env_variables), cex = 5, direction = "both", segment.size=0.25) +
  ggtitle("Ordination of Sample MAG Expression")

nmds_mt


nmds_mag <- ggarrange(nmds_mg, nmds_mt, ncol = 1, common.legend = TRUE)
nmds_mag

################## Repeat for CAZy Expression #######################

otu_cazy <- as.data.frame(t(read.csv("CR_dbCAN_funcmatrix_cazy.csv", header = TRUE, row.names = 1)))


meta_cazy  <- as.data.frame(row.names(otu_cazy))
colnames(meta_cazy)[1] <- "variable"

meta_cazy <- merge(meta_cazy,ph , by="variable")


ordination.model <- metaMDS(otu_cazy, distance='bray', k=3)
#data.bin.fit <- envfit(ordination.model, otu_cazy, permutations = 999)
data.env.fit <- envfit(ordination.model, meta_cazy, permutations = 999, na.rm = TRUE)

#write out NMDS coords
data <- as.data.frame(scores(ordination.model, display = "sites"))
data <- cbind(data, meta_cazy)
data$Treatment<-factor(data$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))


#write out env data
data_env <- as.data.frame(scores(data.env.fit, display = "vectors"))
data_env <- cbind(data_env, env_variables = rownames(data_env))
data_env <- cbind(data_env, pval = data.env.fit$vectors$pvals)

##### plot

nmds_cazy <- ggplot(data, aes(x=NMDS1, y = NMDS2)) +
  geom_point(aes(x=NMDS1, y = NMDS2, color = Treatment),size = 3) +
  #stat_ellipse(aes(color = Treatment)) +
  #coord_fixed() +
  theme_classic(base_size =15) +
  ggtitle("CAZyme Expression") +
  scale_color_manual(values=palette)
nmds_cazy


nmds_cazy +
  geom_segment(data=data_env, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow=arrow(length = unit(0.25, "cm")), color = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data=data_env, aes(x=NMDS1, y = NMDS2, label = env_variables), cex = 5, direction = "both", segment.size=0.25) +
  ggtitle("Ordination of CAZyme Expression")

### Repeat for CAZy Abundance #######

otu_cazy_mg <- as.data.frame(t(read.csv("CR_dbCAN_funcmatrix_cazy_mg.csv", header = TRUE, row.names = 1)))


meta_cazy_mg  <- as.data.frame(row.names(otu_cazy_mg))
colnames(meta_cazy_mg)[1] <- "variable"

meta_cazy_mg <- merge(meta_cazy_mg,ph , by="variable")


ordination.model_mg <- metaMDS(otu_cazy_mg, distance='bray', k=3)
data.env.fit_mg <- envfit(ordination.model_mg, meta_cazy_mg, permutations = 999, na.rm = TRUE)

#write out NMDS coords
data_mg <- as.data.frame(scores(ordination.model_mg, display = "sites"))
data_mg <- cbind(data_mg, meta_cazy_mg)
data_mg$Treatment<-factor(data_mg$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))
data_mg$Day <- as.factor(data_mg$Day)

#write out env data
data_env_mg <- as.data.frame(scores(data.env.fit_mg, display = "vectors"))
data_env_mg <- cbind(data_env_mg, env_variables = rownames(data_env_mg))
data_env_mg <- cbind(data_env_mg, pval = data.env.fit_mg$vectors$pvals)

##### plot

nmds_cazy_mg <- ggplot(data_mg, aes(x=NMDS1, y = NMDS2)) +
  geom_point(aes(x=NMDS1, y = NMDS2, color = Treatment),size = 3) +
  #stat_ellipse(aes(color = Treatment)) +
  #coord_fixed() +
  theme_classic(base_size =15) +
  scale_color_manual(values=palette)+
  ggtitle("CAZyme Abundance")
nmds_cazy_mg


nmds_cazy_mg +
  geom_segment(data=data_env_mg, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow=arrow(length = unit(0.25, "cm")), color = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data=data_env_mg, aes(x=NMDS1, y = NMDS2, label = env_variables), cex = 5, direction = "both", segment.size=0.25) +
  ggtitle("CAZyme Abundance")

nmds_cazy_total <- ggarrange(nmds_cazy_mg, nmds_cazy, common.legend = TRUE, legend = "right", labels = "AUTO", font.label = list(size=25))
nmds_cazy_total

p1 <- nmds_cazy_total

#ggsave("nmds_cazy_total.png", dpi = 600, width = 8, heigh = 6, units = "in")

########### stats ####

#MT stats
dist_mt <- vegdist(otu_mt)
attach(meta_mt)  

anoTreatment_mt <- anosim(dist_mt, Treatment)
anoReplicate_mt <- anosim(dist_mt, Replicate)
anoDay_mt <- anosim(dist_mt, Day)
#anopH_mt <- anosim(dist_mt, pH)

summary(anoTreatment_mt)
summary(anoReplicate_mt)
summary(anoDay_mt)

#cazy stats
dist_cazy <- vegdist(otu_cazy)
attach(meta_cazy)  


anoTreatment_cazy <- anosim(dist_cazy, Treatment)
anoReplicate_cazy <- anosim(dist_cazy, Replicate)
anoDay_cazy <- anosim(dist_cazy, Day)
#anopH_cazy <- anosim(dist_cazy, pH)

summary(anoTreatment_cazy)
summary(anoDay_cazy)
summary(anoReplicate_cazy)

#cazy stats mg
dist_cazy_mg <- vegdist(otu_cazy_mg)
attach(meta_cazy_mg)  


anoTreatment_cazy_mg <- anosim(dist_cazy_mg, Treatment)
anoReplicate_cazy_mg <- anosim(dist_cazy_mg, Replicate)
anoDay_cazy_mg <- anosim(dist_cazy_mg, Day)


summary(anoTreatment_cazy_mg)
summary(anoDay_cazy_mg)
summary(anoReplicate_cazy_mg)

Dataset <- c("CAZy Abundance","CAZy Abundance","CAZy Abundance","CAZy Expression","CAZy Expression","CAZy Expression")
Variables <- c("Treatment", "Day", "Replicate")
ANOSIM_R <- c("0.1329","0.229","-0.000874","0.14","0.238","0.002942")
ANOSIM_p <- c("0.001","0.001","0.491","0.001","0.001","0.343")

stat_table <- data.frame(Dataset, Variables, ANOSIM_R, ANOSIM_p)
colnames(stat_table) <- c("Dataset", "Variable", "ANOSIM", "ANOSIM p-value")

p_stat_table <- ggtexttable(stat_table, rows = NULL,
                            theme = ttheme("light"))
p_stat_table

p4 <-  p_stat_table

#https://stats.oarc.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
#compare matrix
#https://jkzorz.github.io/2019/07/08/mantel-test.html


#mantel test of distance matrices

#mt_cazy <- mantel(dist_mt, dist_cazy, method = "spearman", permutations = 9999, na.rm = TRUE)
#mt_cazy

#mantel stat is R: 0.7424 and pvalue is: 1e-04
# strong relationsgip between cazy dissimilarity matrix and mag dissimilarily matrix
#as samples because more dissimilar in mag, they also become more dissimilar in cazy


#make env dist for mt and cazy
#mt
mt_ph <- meta_mt$pH
mt_day <- meta_mt$Day

distph_mt <- dist(mt_ph, method = "euclidean")
distday_mt <- dist(mt_day, method = "euclidean")

#mag_ph_mt <- mantel(dist_mt, distph_mt, method = "spearman", permutations = 9999, na.rm = TRUE)
#mag_day_mt <- mantel(dist_mt, distday_mt, method = "spearman", permutations = 9999, na.rm = TRUE)

mag_ph_mt
mag_day_mt

#cazy

cazy_ph <- meta_cazy$pH
cazy_day <- meta_cazy$Day

distph_cazy <- dist(cazy_ph, method = "euclidean")
distday_cazy <- dist(cazy_day, method = "euclidean")

#mag_ph_cazy <- mantel(dist_cazy, distph_cazy, method = "spearman", permutations = 9999, na.rm = TRUE)
#mag_day_cazy <- mantel(dist_cazy, distday_cazy, method = "spearman", permutations = 9999, na.rm = TRUE)

mag_ph_cazy
mag_day_cazy



# compare cazy mt and mg dist from centroids
attach(meta_cazy_mg)
cazy_mg_t <- dist_to_centroids(dist_cazy_mg, Treatment)
cazy_mg_t$type <- "MG"
cazy_mg_d <- dist_to_centroids(dist_cazy_mg, Day)
cazy_mg_d$type <- "MG"
attach(meta_cazy)
cazy_mt_t <- dist_to_centroids(dist_cazy, Treatment)
cazy_mt_t$type <- "MT"
cazy_mt_d <- dist_to_centroids(dist_cazy, Day)
cazy_mt_d$type <- "MT"

centroid_treatment <- rbind(cazy_mg_t,cazy_mt_t)
centroid_treatment[c("Treatment", "Day", "Replicate")] <- str_split_fixed(centroid_treatment$Item, "_", 3)
centroid_day <- rbind(cazy_mg_d, cazy_mt_d)
centroid_day[c("Treatment", "Day", "Replicate")] <- str_split_fixed(centroid_day$Item, "_", 3)

centroid_treatment_summary <- centroid_treatment %>%
  group_by(type) %>%
  summarise(mean=mean(CentroidDistance),
            sd=sd(CentroidDistance),
            count = n())

centroid_day_summary <- centroid_day %>%
  group_by(type) %>%
  summarise(mean=mean(CentroidDistance),
            sd=sd(CentroidDistance),
            count = n())

centroid_treatment$CentroidGroup <-factor(centroid_treatment$CentroidGroup, levels = c("pH6","pH8","pH10","MES","G05","G10"))
centroid_day$Treatment<-factor(centroid_day$Treatment, levels = c("pH6","pH8","pH10","MES","G05","G10"))


palday <- c("#7DCEA0","#52BE80","#27AE60","#229954","#1E8449","#196F3D","#145A32")
p_cent_treat <- ggplot(data=centroid_treatment, aes(x=type, y = CentroidDistance)) +
  geom_jitter(data=centroid_treatment, aes(color = Day), width = 0.2, size = 2) +
  geom_boxplot(outlier.shape = NA, lwd = 1, alpha = 0) + #rm outliers bc they are plotted with jitter
  theme_classic(base_size=15) +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = palday) +
  #coord_fixed() +
  labs(y= "Distance from Centroid") + 
  stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.5) +
  ggtitle("Grouped by Treatment")
p_cent_treat


centroid_day$CentroidGroup <- as.factor(centroid_day$CentroidGroup)
palette <-c( "#56B4E9","#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")

p_cent_day <- ggplot(data=centroid_day, aes(x=type, y = CentroidDistance, color)) +
  geom_jitter(data=centroid_day, aes(color = Treatment), width = 0.2, size = 2) +
  geom_boxplot(outlier.shape = NA, lwd = 1, alpha = 0) +
  theme_classic(base_size=15) +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = palette) +
  #coord_fixed() +
  labs(y= "Distance from Centroid") + 
  stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.5) +
  ggtitle("Grouped by Day")
p_cent_day

p_cent <- ggarrange(p_cent_treat, p_cent_day + rremove("y.title"), legend = "bottom")
p2 <- p_cent


#ggsave("p_cent.png", dpi = 600, width = 8, heigh = 6, units = "in")



#compare cazy mt and mg dist between centroids
# need to input two lists, so group by treatment, then manually create distances, then combine results


start <- c("pH6","pH6","pH6","pH6","pH6","pH8","pH8","pH8","pH8","pH10","pH10","pH10","MES","MES","G05")
end <- c("pH8","pH10","MES","G05","G10","pH10","MES","G05","G10","MES","G05","G10","G05","G10","G10")

distance_between_mg <- map2_dfr(start, end, ~ {
  idx1 <- which(meta_cazy_mg$Treatment == .x)
  idx2 <- which(meta_cazy_mg$Treatment == .y)
  tibble(
    start_centroid = .x,
    end_centroid = .y,
    distance = dist_between_centroids(dist_cazy_mg, idx1, idx2)
  )
})

distance_between_mg$type <- "MG"


distance_between_mt <- map2_dfr(start, end, ~ {
  idx1 <- which(meta_cazy$Treatment == .x)
  idx2 <- which(meta_cazy$Treatment == .y)
  tibble(
    start_centroid = .x,
    end_centroid = .y,
    distance = dist_between_centroids(dist_cazy, idx1, idx2)
  )
})

distance_between_mt$type <- "MT"

distances_cent <- rbind(distance_between_mg, distance_between_mt)

start <- c("1","1","1","1","1","1","2","2","2","2","2","3","3","3","3","4","4","4","5","5","6")
end <- c("2","3","4","5","6","7","3","4","5","6","7","4","5","6","7","5","6","7","6","7","7")

meta_cazy_mg$Day <- as.factor(meta_cazy_mg$Day)
meta_cazy$Day <- as.factor(meta_cazy$Day)

distance_between_day_mg <- map2_dfr(start, end, ~ {
  idx1 <- which(meta_cazy_mg$Day == .x)
  idx2 <- which(meta_cazy_mg$Day == .y)
  tibble(
    start_centroid = .x,
    end_centroid = .y,
    distance = dist_between_centroids(dist_cazy_mg, idx1, idx2)
  )
})

distance_between_day_mg$type <- "MG"


distance_between_day_mt <- map2_dfr(start, end, ~ {
  idx1 <- which(meta_cazy$Day == .x)
  idx2 <- which(meta_cazy$Day == .y)
  tibble(
    start_centroid = .x,
    end_centroid = .y,
    distance = dist_between_centroids(dist_cazy, idx1, idx2)
  )
})

distance_between_day_mt$type <- "MT"

distances_cent_day <- rbind(distance_between_day_mg, distance_between_day_mt)



distances_cent_summary <- distances_cent %>%
  group_by(type) %>%
  summarise(mean=mean(distance),
            sd=sd(distance),
            count = n())
distances_cent_day_summary <- distances_cent_day %>%
  group_by(type) %>%
  summarise(mean=mean(distance),
            sd=sd(distance),
            count = n())


# viz ##


p_centd_treat <- ggplot(data=distances_cent, aes(x=type, y = distance)) +
  geom_jitter(data=distances_cent, width = 0.2, size = 2) +
  geom_boxplot(outlier.shape = NA, lwd = 1, alpha = 0) + #rm outliers bc they are plotted with jitter
  theme_classic(base_size=15) +
  theme(axis.title.x = element_blank()) +
  labs(y= "Distance Between Centroid") + 
  stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.5) +
  ggtitle("Grouped by Treatment")
p_centd_treat



p_centd_day <- ggplot(data=distances_cent_day, aes(x=type, y = distance)) +
  geom_jitter(data=distances_cent_day, width = 0.2, size = 2) +
  geom_boxplot(outlier.shape = NA, lwd = 1, alpha = 0) +
  theme_classic(base_size=15) +
  theme(axis.title.x = element_blank()) +
  labs(y= "Distance Between Centroid") + 
  stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.5) +
  ggtitle("Grouped by Day")
p_centd_day

p_cent2 <- ggarrange(p_centd_treat, p_centd_day + rremove("y.title"), legend = "bottom")
p_cent2

p3 <- p_cent2

#ggsave("p_cent2.png", dpi = 600, width = 8, heigh = 6, units = "in")

################# stick everything into one plot

fig <- ggdraw() +
  draw_plot(p1, x = 0, y = 0.5, width = 0.6, height = 0.5) +
  draw_plot(p4, x = 0.6, y = 0.5, width = 0.4, height = 0.5) +
  draw_plot(p2, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(p3, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("","C","D","E"), size = 25,
                  x=c(0,0.6,0,0.5),
                  y=c(0.5,1,0.5,0.5))
fig

ggsave("cr_ordination.tiff", device = "tiff", dpi = 700)
ggsave("cr_ordination.png", device = "png", dpi = 700)
ggsave("cr_ordination.pdf", device = "pdf", dpi = 700)
