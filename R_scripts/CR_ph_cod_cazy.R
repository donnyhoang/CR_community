library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(reshape2)

data <- read.csv("CR_pH_measurements.csv", header=TRUE)
shapeorder <- c(65, 66, 67, 68, 69, 70, 4)
data <- subset(data, Treatment != "innoc")
data2 <- subset(data, Replicate == "X")
data <- subset(data, Replicate != "X")

data$Treatment <- factor(data$Treatment, levels = c("pH6", "pH8", "pH10", "MES", "G05", "G10"))
data2$Treatment <- factor(data2$Treatment, levels = c("pH6", "pH8", "pH10", "MES", "G05", "G10"))

#remove NA
data <- na.omit(data)

data <- data %>%
  group_by(Day, Treatment) %>%
  summarise(mean=mean(pH),
            sd=sd(pH),
            count=n())

palette <-c( "#56B4E9","#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")

ph_plot <- ggplot(data = data, aes(x=Day,y=mean, color = Treatment))+
  geom_line() +
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
  geom_point(data=data2, aes(x=Day, y = pH, color = Treatment), shape = "X", size = 4) +
  geom_line(data=data2, aes(x=Day, y = pH, color = Treatment)) +
  scale_color_manual(values =palette) +
  labs(y = "pH") +
  #theme(aspect.ratio=1) +
  theme_classic(base_size=15) +
  facet_grid(~Treatment) +
  ggtitle("pH Measurements")

ph_plot

############### COD ##########################

#load in data
data <- read.csv("COD_total_no_uninnoc.csv")

#take in only relevant columns
data <- data[,c(1,2,3,4,9,10,11)]
head(data)

#change Day to factor
data$Day <- as.factor(data$Day)

#melt for easy manipulation and graphing
data <- melt(data,id=c("Sample","Treatment","Replicate","Day"))
head(data)

#calc mean and sd
data<-data %>%
  group_by(Treatment, Day, Replicate) %>%
  summarise(count=n(),
            mean=mean(value),
            sd=sd(value))
head(data)

#get ready to graph and stats
palette <-c( "#56B4E9", "#999999","#E69F00","#0072B2", "#F0E442","#CC79A7") 
#palette2 <-c("#e8f7ff","#e2e2e2","#ffda88","#9edcff","#fffbbf","#f8d6e9") 
#data$Treatment<-factor(data$Treatment, levels = c("pH 6","pH 6 + MES","pH 6 + 5% glucose","pH 6 + 10% glucose","pH 8","pH 10"))

data$Treatment <- factor(data$Treatment, levels = c("pH6", "pH8", "pH10", "MES", "G05", "G10"))
palette <-c( "#56B4E9", "#F0E442","#CC79A7", "#999999","#E69F00","#0072B2")


############## by treatment

p <-  ggplot(data=data, aes(x=Day,y=mean,color=Treatment)) +
  geom_point(size=3)+
  geom_errorbar(data=data, aes(ymin=mean-sd,ymax=mean+sd),width=0.2,position=position_dodge(0.05)) +
  theme_classic(base_size=15) +
  #theme(aspect.ratio=1) +
  #theme(legend.position = "bottom") +
  scale_colour_manual(values=palette) +
  facet_grid(~ Treatment) +
  labs(title="Remaining COD", y ="Remaining COD (% mg O2/L)")
p



#################### CAZy expressiom
cazy <- read.csv("cazy_expression_summary.csv")

cazy_summary <- cazy %>%
  group_by(Treatment,Day) %>%
  summarise(mean=mean(percent),
            sd=sd(percent),
            count=n())

cazy$Day <- as.character(cazy$Day)
cazy$Treatment <- factor(cazy$Treatment, levels = c("pH6", "pH8", "pH10", "MES", "G05", "G10"))
#palette <-c( "#56B4E9", "#999999","#E69F00","#0072B2", "#F0E442","#CC79A7")

pbox <- ggplot(cazy, aes(x=Day, y=percent, fill=Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values=palette) +
  #scale_color_manual(values=pal1) +
  theme_classic(base_size=15) +
  #ylim(0,19) +
  ylab("Percent reads mapped") + 
  xlab("Day") +
  #theme(axis.text.x=element_text(angle=45,hjust=1)) +
  #theme(aspect.ratio=1)+
  ggtitle("CAZy Expression") +
  facet_grid(  ~ Treatment , scales = "free_x")
pbox

#pull out numbers to quickly reference in paper
d1 <- subset(data_trim, Day =="1")
d7 <- subset(data_trim, Day =="7")

p1 <- ph_plot
p2 <- p 
p3 <- pbox 

fig<- plot_grid(p1,p2,p3, nrow = 3, labels = "AUTO", label_size = 25,
                 rel_widths = c(10,10,10))
fig


ggsave("ph_cod_cazy.tiff", device = "tiff", dpi = 700)
ggsave("ph_cod_cazy.png", device = "png", dpi = 700)
ggsave("ph_cod_cazy.pdf", device = "pdf", dpi = 700)


ph_cod <- plot_grid(p1, p2, nrow = 2, rel_widths = c(10,10))
ph_cod
ggsave("ph_cod.png", device = "png", dpi = 700)
