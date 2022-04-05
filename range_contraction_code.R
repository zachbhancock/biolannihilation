####Range contraction R code####
setwd("~/Desktop/Range Contractions/Revised")
library(ggplot2)
library(ggpubr)
library(ggbreak)

#code for making random contraction patterns#
N <- 88
landscape <- matrix("g", nrow=20, ncol=20)
landscape2 <- landscape
landscape2[sample(1:length(landscape),N,replace=FALSE)] <- "d"
landscape3 <- landscape2
placestoreplace <- which(landscape3=="g")
idstoreplace <- sample(placestoreplace,N,replace=FALSE)
landscape3 <- replace(landscape3,idstoreplace,"d")
landscape4 <- landscape3
placestoreplace.2 <- which(landscape4=="g")
idstoreplace.2 <- sample(placestoreplace.2, N, replace=FALSE)
landscape4 <- replace(landscape4, idstoreplace.2, "d")
landscape5 <- landscape4
placestoreplace.3 <- which(landscape5=="g")
idstoreplace.3 <- sample(placestoreplace.3, N, replace=FALSE)
landscape5 <- replace(landscape5, idstoreplace.3, "d")
write.table(landscape2, "random_frag_1.txt", row.names=FALSE, col.names=FALSE, sep=",")
write.table(landscape3, "random_frag_2.txt", row.names=FALSE, col.names=FALSE, sep=",")
write.table(landscape4, "random_frag_3.txt", row.names=FALSE, col.names=FALSE, sep=",")
write.table(landscape5, "random_frag_4.txt", row.names=FALSE, col.names=FALSE, sep=",")

#analyzing heterozygosity spatially

cbPalette <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

time800 <- read.table("random_1_800.txt", header=TRUE, sep=",")
time650 <- read.table("random_1_650.txt", header=TRUE, sep=",")
time550 <- read.table("random_1_550.txt", header=TRUE, sep=",")
time450 <- read.table("random_1_450.txt", header=TRUE, sep=",")
time350 <- read.table("random_1_350.txt", header=TRUE, sep=",")
time0 <- read.table("random_1_0.txt", header=TRUE, sep=",")
random_1 <- rbind(time0, time350, time450, time550, time650, time800)
ggplot(random_1, aes(x=x, y=y, color=relative, size=het)) + scale_color_viridis() + geom_point(alpha=0.8) + facet_wrap(~time)

#comparisons between patterns
before.random <- read.table("random_1_800.txt", header=TRUE, sep=",")
mid.random <- read.table("random_1_550.txt", header=TRUE, sep=",")
end.random <- read.table("random_1_350.txt", header=TRUE, sep=",")
after.random <- read.table("random_1_0.txt", header=TRUE, sep=",")
random <- rbind(before.random, mid.random, end.random, after.random)
random$pattern <- "random"

before.shrinkage <- read.table("shrinkage_800.txt", header=TRUE, sep=",")
mid.shrinkage <- read.table("shrinkage_550.txt", header=TRUE, sep=",")
end.shrinkage <- read.table("shrinkage_350.txt", header=TRUE, sep=",")
after.shrinkage <- read.table("shrinkage_0.txt", header=TRUE, sep=",")
shrinkage <- rbind(before.shrinkage, mid.shrinkage, end.shrinkage, after.shrinkage)
shrinkage$time <- factor(shrinkage$time, levels=c("800", "550", "350", "0"))
shrinkage$pattern <- "shrinkage"

before.amputation <- read.table("amputation_800.txt", header=TRUE, sep=",")
mid.amputation <- read.table("amputation_550.txt", header=TRUE, sep=",")
end.amputation <- read.table("amputation_350.txt", header=TRUE, sep=",")
after.amputation <- read.table("amputation_0.txt", header=TRUE, sep=",")
amputation <- rbind(before.amputation, mid.amputation, end.amputation, after.amputation)
amputation$time <- factor(amputation$time, levels=c("800", "550", "350", "0"))
amputation$pattern <- "amputation"

before.frag <- read.table("frag_800.txt", header=TRUE, sep=",")
mid.frag <- read.table("frag_550.txt", header=TRUE, sep=",")
end.frag <- read.table("frag_350.txt", header=TRUE, sep=",")
after.frag <- read.table("frag_0.txt", header=TRUE, sep=",")
fragmentation <- rbind(before.frag, mid.frag, end.frag, after.frag)
fragmentation$time <- factor(fragmentation$time, levels=c("800", "550", "350", "0"))
fragmentation$pattern <- "fragmentation"

patterns <- rbind(shrinkage, amputation, fragmentation, random)
patterns$time <- factor(patterns$time, levels = c(800, 550, 350, 0))

#ANOVA and Tukey HSD tests to determine groups
amp.model <- aov(het~time, data=amputation)
tukey.amp.model <- TukeyHSD(amp.model, "time")
frag.model <- aov(het~time, data=fragmentation)
tukey.frag.model <- TukeyHSD(frag.model, "time")
shrink.model <- aov(het~time, data=shrinkage)
tukey.shrink.model <- TukeyHSD(shrink.model, "time")

ind.het <- ggplot(patterns, aes(x=time, y = relative, fill = pattern)) + theme_bw() + scale_fill_manual(values=cbPalette) + geom_boxplot() + xlab("Timesteps") + ylab("Individual \u03C0 / Average \u03C0") + theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.title=element_blank(), legend.text=element_text(size=10), legend.position="top") + geom_hline(yintercept=1.0, linetype="dashed")

ggplot(patterns, aes(x=time, y=relative, fill=pattern)) + geom_boxplot()

#demographic analysis

random.demo <- read.table("random_1_results.txt", header=TRUE, sep=",")
random.demo$pattern <- "random"
shrinkage.demo <- read.table("shrinkage_results.txt", header=TRUE, sep=",")
shrinkage.demo$pattern <- "shrinkage"
amputation.demo <- read.table("amputation_results.txt", header=TRUE, sep=",")
amputation.demo$pattern <- "amputation"
frag.demo <- read.table("fragmentation_results.txt", header=TRUE, sep=",")
frag.demo$pattern <- "fragmentation"
demo <- rbind(shrinkage.demo, amputation.demo, frag.demo)

psize <- ggplot(demo, aes(x = generation, y = pop_size, color = pattern)) + theme_bw() + scale_color_manual(values=cbPalette) + geom_point(size=1, alpha=0.8) + geom_line() + geom_vline(xintercept=20100, linetype="dashed") + geom_vline(xintercept=20200, linetype="dashed") + geom_vline(xintercept=20300, linetype="dashed") + geom_vline(xintercept=20400, linetype="dashed") + xlab("Timesteps") + ylab("Population size") + theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.text=element_text(size=12), legend.title=element_blank()) + guides(colour = guide_legend(override.aes = list(size=2)))

mean_off <- ggplot(demo, aes(x = generation, y = mean_off, color = pattern)) + theme_bw() + scale_color_manual(values=cbPalette) + geom_point(size=1, alpha=0.8) + geom_vline(xintercept=20100, linetype="dashed") + geom_vline(xintercept=20200, linetype="dashed") + geom_vline(xintercept=20300, linetype="dashed") + geom_vline(xintercept=20400, linetype="dashed") + xlab("Timesteps") + ylab("Mean Offspring") + theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.text=element_text(size=12), legend.title=element_blank()) + guides(colour = guide_legend(override.aes = list(size=2))) + geom_smooth()

relat <- ggplot(demo, aes(x = generation, y = relatedness, color = pattern)) + theme_bw() + scale_color_manual(values=cbPalette) + geom_point(size=1, alpha=0.4) + geom_vline(xintercept=20100, linetype="dashed") + geom_vline(xintercept=20200, linetype="dashed") + geom_vline(xintercept=20300, linetype="dashed") + geom_vline(xintercept=20400, linetype="dashed") + xlab("Timesteps") + ylab("Relatedness") + theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.text=element_text(size=12), legend.title=element_blank()) + geom_smooth(method="lm") + guides(colour = guide_legend(override.aes = list(size=2)))

meanage <- ggplot(demo, aes(x = generation, y = mean_age, color = pattern)) + theme_bw() + scale_color_manual(values=cbPalette) + geom_point(size=1, alpha=0.4) + geom_vline(xintercept=20100, linetype="dashed") + geom_vline(xintercept=20200, linetype="dashed") + geom_vline(xintercept=20300, linetype="dashed") + geom_vline(xintercept=20400, linetype="dashed") + xlab("Timesteps") + ylab("Mean age") + theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.text=element_text(size=12), legend.title=element_blank()) + geom_smooth(method="lm") + guides(colour = guide_legend(override.aes = list(size=2)))

maxage <- ggplot(demo, aes(x = pop_size, y = max_age, color = pattern)) + theme_bw() + scale_color_manual(values=cbPalette) + geom_point(size=1, alpha=0.8) + xlab("Population size") + ylab("Max age") + theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.text=element_text(size=12), legend.title=element_blank()) + guides(colour = guide_legend(override.aes = list(size=2)))

panels <- ggarrange(psize, relat, meanage, mean_off, labels=c("A", "B", "C", "D"), ncol=2, nrow=2, common.legend=TRUE, legend="top")
panels

shrink.landscape <- ggplot(after.shrinkage, aes(x=x, y=y, size=het, color=relative)) + theme_bw() + geom_point() + scale_x_continuous(limit=c(0,20)) + scale_y_continuous(limit=c(0,20)) + labs(title="shrinkage") + theme(plot.title=element_text(hjust=0.5, size=14), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), axis.text=element_text(size=10))

amp.landscape <- ggplot(after.amputation, aes(x=x, y=y, size=het, color=relative)) + theme_bw() + geom_point() + scale_x_continuous(limit=c(0,20)) + scale_y_continuous(limit=c(0,20)) + labs(title="amputation") + theme(plot.title=element_text(hjust=0.5, size=14), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), axis.text=element_text(size=10))

frag.landscape <- ggplot(after.frag, aes(x=x, y=y, size=het, color=relative)) + theme_bw() + geom_point() + scale_x_continuous(limit=c(0,20)) + scale_y_continuous(limit=c(0,20)) + labs(title="fragmentation") + theme(plot.title=element_text(hjust=0.5, size=14), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), axis.text=element_text(size=10))

ggarrange(amp.landscape, frag.landscape, shrink.landscape, ind.het, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)

#TMRCA

"#E69F00", "#009E73" 

amp.length <- read.table("amputation_length.txt", header=TRUE, sep=",")
amp.pre <- as.data.frame(amp.length$length.pre)
names(amp.pre)[1] <- "length"
amp.pre <- na.omit(amp.pre)
amp.during <- as.data.frame(amp.length$length.during)
names(amp.during)[1] <- "length"
amp.during <- na.omit(amp.during)
amp.post <- as.data.frame(amp.length$length.post)
names(amp.post)[1] <- "length"
amp.post <- na.omit(amp.post)
amp.both <- data.frame(variable = c(rep("pre-contraction", length(amp.pre$length)), rep("mid-contraction", length(amp.during$length)),
                              rep("post-contraction", length(amp.post$length))), 
                 value=c(amp.pre$length, amp.during$length, amp.post$length))
amp.both$variable <- factor(amp.both$variable, levels = c("pre-contraction", "mid-contraction", "post-contraction"))
amp.tmrca <- ggplot(amp.both, aes(x=value, fill=variable)) + scale_fill_manual(values = c("#E69F00", "#009E73", "#0072B2")) + geom_histogram(position="dodge") + xlim(0, 300000) + theme(legend.title = element_blank(), axis.text = element_text(size=10), legend.text = element_text(size=12), axis.title = element_text(size=12)) + ggtitle("amputation") + xlab("haplotype length") + ylab("count")


+ annotate(geom="text", x = 300000, y=6000, label="mean length = 42445", color="#E69F00") + annotate(geom="text", x = 300000, y=5750, label="# of haplotypes = 23559", color="#E69F00") + annotate(geom="text", x = 300000, y=5500, label="mean length = 34684", color="#009E73") + annotate(geom="text", x = 300000, y=5250, label="# of haplotypes = 28831", color="#009E73")

ggplot(amp.both, aes(x=variable, y=value)) + geom_boxplot()

frag.length <- read.table("frag_length.txt", header=TRUE, sep=",")
frag.pre <- as.data.frame(frag.length$length.pre)
names(frag.pre)[1] <- "length"
frag.pre <- na.omit(frag.pre)
frag.post <- as.data.frame(frag.length$length.post)
names(frag.post)[1] <- "length"
frag.post <- na.omit(frag.post)
frag.both <- data.frame(variable = c(rep("pre-contraction", length(frag.pre$length)), 
                                     rep("post-contraction", length(frag.post$length))), 
                        value=c(frag.pre$length, frag.post$length))
frag.tmrca <- ggplot(frag.both, aes(x=value, fill=variable)) + scale_fill_manual(values = c("#E69F00", "#009E73")) + geom_histogram() + theme(legend.title = element_blank(), axis.text = element_text(size=10), legend.text = element_text(size=12), axis.title = element_text(size=12)) + ggtitle("fragmentation") + ylab("count") + xlab("haplotype length") + xlim(0, 400000) + annotate(geom="text", x = 300000, y=10000, label="mean length = 60336", color="#E69F00") + annotate(geom="text", x = 300000, y=9500, label="# of haplotypes = 16570", color="#E69F00") + annotate(geom="text", x = 300000, y=9000, label="mean length = 34284", color="#009E73") + annotate(geom="text", x = 300000, y=8500, label="# of haplotypes = 29168", color="#009E73")


