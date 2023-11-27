#graphs of fake data illustrating the idea of realized vs fundamental niche

library(ggplot2)

set.seed(3)

x <-runif(500)
y <- 4 * x ^ 2 + rnorm(length(x), sd = 4)

data <- data.frame(x, y)

p <- ggplot(data, aes(x,y)) +
  theme_classic(base_size = 20) +
  stat_ellipse(level = 0.5) +
  stat_ellipse(level = 0.99)
p


data2 <- data
data2$x <- (data2$x) + 1
data2$y <- (data2$y) * .5
data2$group <- "B"
data$group <- "A"
df <- rbind(data, data2)


#blue "#1E88E5"
#yellow "#FFC107"

p1 <- ggplot(data, aes(x,y)) +
  ylim(-10, 13) +
  xlim(-0.5, 2.5) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  stat_ellipse(geom = "polygon", fill = "#1E88E5", alpha = 0.25, level = 0.99) +
  stat_ellipse(geom = "polygon", fill = "#FFC107", alpha = 0.25, level = 0.5) +
  labs(x = "Envrionmental Condition", y = "Resource")
p1

ggsave("niche_fun_v_realized.png", device = "png", dpi = 700)

p2 <- ggplot(data, aes(x,y)) +
  ylim(-10, 13) +
  xlim(-0.5, 2.5) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  stat_ellipse(geom = "polygon", fill = "#1E88E5", alpha = 0.25, level = 0.99)+
  labs(x = "Envrionmental Condition", y = "Resource")
p2

ggsave("niche.png", device = "png", dpi = 700)

p3 <- ggplot(df, aes(x,y, fill = group)) +
  ylim(-10, 13) +
  xlim(-0.5, 2.5) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#1E88E5", "#FFC107")) +
  stat_ellipse(geom = "polygon", alpha = 0.25, level = 0.99)+
  labs(x = "Envrionmental Condition", y = "Resource")
p3

ggsave("niche_fun_v_realized2.png", device = "png", dpi = 700)



###make an illustrative plot for functional redundancy

x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
y <- c(8,8,8,8,8,8,8,8,8,8,1,2,3,4,5,6,7,8,9,10)
group <- c("A","A","A","A","A","A","A","A","A","A",
           "B","B","B","B","B","B","B","B","B","B")

data <- data.frame(x,y,group)


p_fr <- ggplot(data, aes(x,y, color = group)) +
  geom_line(linetype = "dashed", lwd = 2) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("#1E88E5", "#FFC107")) +
  labs(x = "Taxonomic Diversity", y = "Functional Diversity") +
  ggtitle("Functional Redundancy?")
p_fr
ggsave("fr_plot.png", device = "png", dpi = 700)

# niche breadth vs overlap
library(reshape2)
set.seed(3)


d1 <- as.data.frame(rnorm(1000, 35, 20))
d1$x <- row.names(d1)
colnames(d1)[1] <- "y"
d1$group <- "A"
d2 <- as.data.frame(rnorm(1000, 65, 20))
d2$x <- row.names(d2)
d2$group <- "B"
colnames(d2)[1] <- "y"

dtotal <- rbind(d1,d2)



p1 <- ggplot(d1, aes(x=y, color = group, fill = group)) +
  ylim(0, 0.023)+
  xlim(-25, 135) +
  geom_density(alpha = 0.4) +
  theme_classic(base_size = 20) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("#F6C83E")) +
  scale_fill_manual(values = c("#F6C83E")) +
  labs(x = "Resource Gradient", y = "Resource Use")
p1 

ggsave("niche_breadth.png", device = "png", dpi = 700)


p3 <- ggplot(dtotal, aes(x=y, color = group, fill = group)) +
  ylim(0, 0.023)+
  xlim(-25, 135) +
  geom_density(alpha = 0.4) +
  theme_classic(base_size = 20) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("#F6C83E", "#E69EA1")) +
  scale_fill_manual(values = c("#F6C83E", "#E69EA1")) +
  labs(x = "Resource Gradient", y = "Resource Use")
p3
ggsave("niche_overlap.png", device = "png", dpi = 700)

