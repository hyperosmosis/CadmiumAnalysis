setwd("~/Documents/Cui/Cadmium/")

library(tidyverse)
library(reshape2)
library(ggpubr)
library(qiime2R)
library(grid)
library(gridExtra)
library(RColorBrewer)

# Figure 2, Shannon Diversity

## Prepare data

shannon = read.csv("shannon_clean.csv", header = T, stringsAsFactors = F)
colnames(shannon) = c("sample.id", 1, 1112, 2223, 3334, 4445, 5556, 6667, 7778,
                      8889, 10000, "Genotype", "Chemical", "Sex")

group_mean <- shannon %>%
  group_by(Genotype, Chemical, Sex) %>%
  summarise_if(is.numeric, mean, na.rm=T)

group_sd <- shannon %>%
  group_by(Genotype, Chemical, Sex) %>%
  summarise_if(is.numeric, function(x) sd(x)/sqrt(length(x)))

mean_m = melt(group_mean, 
              id.vars=c('Genotype', 'Chemical', 'Sex'), 
              variable.name='Iteration')

sd_m = melt(group_sd,
            id.vars=c('Genotype', 'Chemical', 'Sex'),
            variable.name='Iteration')

shannon_m = data.frame(mean_m, sd = sd_m[,5])
shannon_m$Cadmium[shannon_m$Chemical=="ctrl"] = "Control"
shannon_m$Cadmium[shannon_m$Chemical=="lowCd"] = "Low"
shannon_m$Cadmium[shannon_m$Chemical=="highCd"] = "High"
shannon_m$GenoSex = paste(shannon_m$Genotype, shannon_m$Sex, sep=" ")
shannon_m$Cadmium <- factor(shannon_m$Cadmium, levels = c("Control", "Low", "High"))

## Make line plots

linePlot <- function(genotype, sex, title){
  ggplot(shannon_m[shannon_m$Genotype== genotype & shannon_m$Sex==sex,],
         aes(x=Iteration, y = value, group = Cadmium))+
    xlab("Sequences") + ylab("Shannon Index") + 
    ggtitle(title) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(x=Iteration, ymin=value-sd, ymax=value+sd, color = Cadmium),width=0.25)+
    geom_line(aes(color = Cadmium)) + geom_point(aes(color = Cadmium)) + 
    #scale_color_discrete(breaks=c("Control","Low","High")) + 
    scale_color_manual(values = c("#9999a0", "#1b2aed", "#f6131e"))
}

p1 = linePlot("E3", "M", "Alpha Diversity: ApoE3-KI Male")
p2 = linePlot("E3", "F", "Alpha Diversity: ApoE3-KI Female")
p3 = linePlot("E4", "M", "Alpha Diversity: ApoE4-KI Male")
p4 = linePlot("E4", "F", "Alpha Diversity: ApoE4-KI Female")

png("AlphaDiversity.png", height = 6, width = 8, units = "in", res = 350)
ggarrange(p1, p3, p2, p4, ncol=2, nrow=2, labels = c("a", "b", "c", "d"),
          common.legend = TRUE, legend="bottom")
dev.off()

png("AlphaDiversity_Basal.png", height = 3.5 , width = 6, units = "in", res = 250)
ggplot(shannon_m[shannon_m$Cadmium == "Control",],
       aes(x=Iteration, y = value, group = GenoSex))+
  xlab("Sequences") + ylab("Shannon Index") + 
  ggtitle("Alpha Diversity: Vehicle groups only") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  geom_errorbar(aes(x=Iteration, ymin=value-sd, ymax=value+sd, color = GenoSex), width=0.25)+
  geom_line(aes(color = GenoSex)) + geom_point(aes(color = GenoSex)) + 
  labs(color='Group') + scale_color_manual(values = c("#9999a0", "#1b2aed", "#f6131e", "chartreuse3"))
dev.off()

# Figure 3: Beta Diversity, 3 Beta-Diversity plot

pco = read_qza("weighted_unifrac_pcoa_results.qza")
pco_bray = read_qza("bray_curtis_pcoa_results.qza")


betaPlot = function(qza, group, title, color){
  pcoa_dat = qza$data$Vectors[,c(1:3)]
  pcoa_dat$Group = substr(pcoa_dat$SampleID,1,nchar(as.character(pcoa_dat$SampleID))-1)
  pcoa_dat$Sex = substr(pcoa_dat$Group, nchar(as.character(pcoa_dat$Group)),
                        nchar(as.character(pcoa_dat$Group)))
  pcoa_dat$Chemical = substr(pcoa_dat$Group, 3, nchar(as.character(pcoa_dat$Group))-1)
  pcoa_dat$Genotype = substr(pcoa_dat$Group, 1, 2)
  
  pcoa_dat$Chemical <- revalue(pcoa_dat$Chemical, c("ctrl" = "Control",
                                                    "lowCd" = "Low",
                                                    "highCd" = "High"))
  pcoa_dat$Chemical <- factor(pcoa_dat$Chemical, levels = c("Control", "Low", "High"))
  ggplot(pcoa_dat, aes(x=PC1, y=PC2, color=pcoa_dat[,group])) +
    geom_point() +  
    xlab(paste("PC1: ", round(100*pco$data$ProportionExplained[1]), "%")) +
    ylab(paste("PC2: ", round(100*pco$data$ProportionExplained[2]), "%")) +
    labs(color=colnames(pcoa_dat)[group]) + theme_bw() + 
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=color)
}

sex_weighted = betaPlot(pco, 5, "Weighted Unifrac Grouped by Sex")
chem_weighted = betaPlot(pco, 6, "Weighted Unifrac Grouped by Chemical")
geno_weighted = betaPlot(pco, 7, "Weighted Unifrac Grouped by Genotype")

sex_bray = betaPlot(pco_bray, 5, "Separated by Sex", c("#9999a0", "#1b2aed"))
chem_bray = betaPlot(pco_bray, 6, "Separated by Chemical", c("#9999a0", "#1b2aed", "#f6131e"))
geno_bray = betaPlot(pco_bray, 7, "Separated by Genotype", c("#9999a0", "#1b2aed"))

pdf("BetaDiversity.pdf", height = 8, width = 4)
plot_grid(sex_bray, geno_bray, chem_bray, nrow = 3, align = "v",
          labels = "AUTO")
dev.off()

# Figure 4: Compositional stacked bar plot

## Level 3:

level3 = read.csv("level-3.csv", header = T, stringsAsFactors = F)
level6 = read.csv("level-6.csv", header = T, stringsAsFactors = F)
level7 = read.csv("level-7.csv", header = T, stringsAsFactors = F)

#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#col = sample(color, 120)

stackedBar <- function(data, proportion = T){
  average <- data %>%
    group_by(Genotype, Chemical, Sex) %>%
    summarise_if(is.numeric, mean, na.rm=T)
  
  average$Sum = rowSums(average[ ,sapply(average, is.numeric)])
  
  average_m = melt(average[,c(1:(dim(average)[2]-3), dim(average)[2])], 
                 id.vars=c('Genotype', 'Chemical', 'Sex', 'Sum'), 
                 variable.name='Bacteria')
  
  average_m$Bacteria = sapply(average_m$Bacteria, function(x) strsplit(gsub("\\D__", "",x), ";")[[1]])
  average_m$Bacteria = gsub("Bacteria.", "", average_m$Bacteria)
  average_m$Bacteria = gsub(".", " ", average_m$Bacteria, fixed=TRUE)
  average_m$Bacteria[average_m$Bacteria==" "] = "Bacteria"
  
  average_m$Group = paste(average_m$Genotype, average_m$Chemical, average_m$Sex, sep = "")
  
  average_m$Group = factor(average_m$Group, 
                           levels = c("E3ctrlM","E3lowCdM","E3highCdM",
                                      "E3ctrlF","E3lowCdF","E3highCdF",
                                      "E4ctrlM","E4lowCdM","E4highCdM",
                                      "E4ctrlF","E4lowCdF","E4highCdF"))
  stand_dev <- data %>%
    group_by(Genotype, Chemical, Sex) %>%
    summarise_if(is.numeric, sd, na.rm=T)
  
  stand_dev_m = melt(stand_dev[,c(1:25)],
               id.vars=c('Genotype', 'Chemical', 'Sex'),
               variable.name='Bacteria')
  
  if(proportion == T){
    average_prop= average_m %>%
      group_by(Genotype, Chemical, Sex) %>%
      mutate(Proportion = (value/Sum)*100)
    
    ggplot(average_prop, aes(x = Group, y = Proportion, fill=Bacteria)) +
      geom_bar(stat='identity') + 
      scale_fill_manual(values = col) + 
      guides(fill = guide_legend(nrow = 22)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 60, hjust = 0.75, size = 12)) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.ticks.x=element_blank(),
            panel.border=element_blank(),
            axis.title.x=element_blank(),
            axis.line=element_blank())
  }
  
  else{
    ggplot(average_m, aes(x = Group, y = value, fill=Bacteria)) +
      geom_bar(stat='identity') + 
      scale_fill_manual(values = col) + 
      guides(fill = guide_legend(nrow = 22)) +
      theme_bw()
  }
}

col <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow3",
         "darkorange4","brown")

png("StackedBarLevel3.png", width = 10, height = 7, units = "in", res = 250)
stackedBar(level3, proportion = T) # 22 different groups
dev.off()

png("StackedBarLevel3_NoProp.png", width = 11.5, height = 6.5, units = "in", res = 250)
stackedBar(level3, proportion = F) # 22 different groups
dev.off()

stackedBar(level6, proportion = T)
stackedBar(level7, proportion = T)


# Figure 5, bar plots

library(multcomp)

E3F = read.table("ANCOM/sig_ancom_E3F_l7.txt", header=T,sep="\t")
E3M = read.table("ANCOM/sig_ancom_E3M_l7.txt", header = T, sep="\t")
E4F = read.table("ANCOM/sig_ancom_E4F_l7.txt", header = T, sep="\t")
E4M = read.table("ANCOM/sig_ancom_E4M_l7.txt", header = T, sep="\t")

bacPlot <- function(data){
  title = gsub(".", " ", colnames(data)[2], fixed = T)
  colnames(data)[2] = "Value"
  
  dat = data %>%
    group_by(Group) %>%
    summarise(Average = mean(Value), se = sd(Value)/sqrt(n()), max = max(Value))
  
  dat$Group = factor(dat$Group, levels = c("ctrl","lowCd","highCd"))
  
  dat = dat %>%
    arrange(Group)

  p <- ggplot(data = dat, aes(x=Group, y = Average, fill = Group)) + 
    geom_bar(stat="identity") + 
    geom_errorbar(aes(ymin=Average-se, ymax=Average+se), width=.1) +
    guides(fill=FALSE) + ggtitle(paste(strwrap(title,width = 40 ),
                                       collapse = "\n")) + 
    theme(legend.position="none") + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
          axis.title.y = element_blank()) +
    scale_fill_manual(values=c("#9999a0", "#1b2aed", "#f6131e"))

  colnames(data)[2] = "Average"
  
  my_comparisons = list(c("ctrl", "lowCd"), c("ctrl", "highCd"))
  
  p + stat_compare_means(data = data, comparisons = my_comparisons,
                          label = "p.signif", method= "t.test")

}

# Function for box plot with points

makeBox <- function(data){
  plot_list <- lapply(2:(ncol(data)), function(i) {
    
    name <- colnames(data[i])
    title <- gsub(".", " ", colnames(data)[i], fixed = T)
    data$Group <- factor(data$Group, levels = c("ctrl", "lowCd", "highCd"))
    ggplot(data[,c(1, i)], aes(x=Group, y=data[,i], fill=Group)) +
      geom_boxplot() +
      scale_fill_manual(values=c("#9999a0", "#1b2aed", "#f6131e")) +
      geom_jitter(color="black", size=2.5, alpha=0.7, width = 0.2) + theme_classic() +
      theme(legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
            axis.title.y = element_blank()) +
      ggtitle(title)
    
  })
  return(plot_list)
}

E4F = E4F[complete.cases(E4F),]
E4M = E4M[complete.cases(E4M),]

pdf("E3F_sig.pdf", width = 4, height = 4)
ggarrange(plotlist = makeBox(E3F))
dev.off()

pdf("E3M_sig.pdf", width = 8, height = 4)
ggarrange(plotlist = makeBox(E3M))
dev.off()

pdf("E4F_sig.pdf", width = 8, height = 4)
ggarrange(plotlist = makeBox(E4F), ncol = 2)
dev.off()

pdf("E4M_sig.pdf", width = 12, height = 8)
ggarrange(plotlist = makeBox(E4M), ncol = 3, nrow = 2)
dev.off()


# Cadmium RNA-seq figures

genes = c("Cyp2b13", "Cyp2c69", "Cyp4a12a", "Cyp4a12b",
          "Sult2a1", "Sult2a5", "Sult2a6")

counts = as.data.frame(t(counts[counts$Gene_ID %in% genes, c(1, 6:41)]))
colnames(counts) = genes
counts = counts[-1,]

counts$Group = substr(rownames(counts),1, nchar(as.character(rownames(counts)))-2)
counts[,1:7] = sapply(counts[,1:7],as.numeric)

mean = counts %>%
  group_by(Group) %>%
  summarise_all("mean")

se = counts %>%
  group_by(Group) %>%
  summarise_if(is.numeric, function(x) sd(x)/sqrt(length(x)))


## Function for barplot

ll = list()
for (i in 2:8){
  dat = data.frame(Group = mean[,1], mean = mean[,i], se = se[,i])
  colnames(dat) = c("Group", "Mean", "SE")
  dat$Genotype = substr(dat$Group, 3, 7)
  dat$Sex = substr(dat$Group, 1, 1)
  
  order = c("F_ApoE3_CN", "F_ApoE3_LowCd", "F_ApoE3_HighCd",
            "F_ApoE4_CN", "F_ApoE4_LowCd", "F_ApoE4_HighCd",
            "M_ApoE3_CN", "M_ApoE3_LowCd", "M_ApoE3_HighCd",
            "M_ApoE4_CN", "M_ApoE4_LowCd", "M_ApoE4_HighCd")
  name = colnames(mean)[i]
  
  dat = dat %>%
    slice(match(order, Group))
  
  p1 = ggplot(dat[dat$Sex == "M",], aes(x = Group, y = Mean, fill = Genotype)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(x=Group, ymin= Mean - SE, ymax= Mean + SE),width=0.25)+ 
    ggtitle(paste(name, "Male", sep = ": ")) + theme_bw() +
    ylab("Count") + xlab("Treatment") + 
    scale_x_discrete(labels = c('Control','Low Cd','High Cd',
                                'Control','Low Cd','High Cd')) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 = ggplot(dat[dat$Sex == "F",], aes(x = Group, y = Mean, fill = Genotype)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(x=Group, ymin= Mean - SE, ymax= Mean + SE),width=0.25)+ 
    ggtitle(paste(name, "Female", sep = ": ")) + theme_bw() + 
    ylab("Count") + xlab("Treatment") + 
    scale_x_discrete(labels = c('Control','Low Cd','High Cd',
                                'Control','Low Cd','High Cd')) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  ll[[i-1]] = list(p1, p2)
}

png("../Cadmium/Genes.png", height = 18, width = 8, units = "in", res = 250)
ggarrange(plotlist=unlist(ll, recursive = F), ncol = 2, nrow = 7, common.legend = TRUE, legend="bottom")
dev.off()

## SCFA

scfa = read.csv("SCFA_Cadmium.csv", header = T, stringsAsFactors = F,
                check.names = F)
scfa$Chemical <- revalue(scfa$Chemical, c("Control" = "Control",
                                                  "Low Cd" = "Low",
                                                  "High Cd" = "High"))
scfa$Chemical <- factor(scfa$Chemical, levels = c("Control", "Low", "High"))

scfaPlot <- function(dat, scfa){
  dat = dat[dat$SCFA == scfa,]
  
  ggplot(dat, aes(x = Group, y = Value, fill = Chemical)) + 
    geom_bar(stat="identity",position="dodge") + 
    geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=0.2, 
                  position=position_dodge(.9)) +
    ggtitle(scfa) +
    ylab("Value at μg/mL") + theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values=c("#599ad3", "#f9a65a", "#cc7058"))
    
}

p1 = scfaPlot(scfa, "Lactate")
p2 = scfaPlot(scfa, "Butyric acid")
p3 = scfaPlot(scfa, "2-Methylbutyric acid")

png("SCFA.png", height = 4.5, width = 12, units = "in", res = 250)
ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, 
          labels = "AUTO", legend="bottom")
dev.off()

pdf("Lactate.pdf", width = 5, height = 4)
p1
dev.off()

library(grid)

scfa <- read.csv("SCFA_Cadmium.csv", header = T)
ApoE3M <- scfa[scfa$Group == "ApoE3M",]
ApoE3F <- scfa[scfa$Group == "ApoE3F",]
ApoE4M <- scfa[scfa$Group == "ApoE4M",]
ApoE4F <- scfa[scfa$Group == "ApoE4F",]

scfaPlot <- function(data){
  ll <- list()
  
  for (i in 4:ncol(data)){
    dat = data[,c(3,i)]
    title = gsub(".", " ", colnames(dat)[2], fixed = T)
    colnames(dat)[2] = "Value"
    dat$Group = factor(dat$Chemical, levels = c("Control","Low","High"))
    
    p <- ggplot(data = dat, aes(x=Group, y = Value, fill = Group)) + 
      geom_boxplot() + ggtitle(title) + 
      scale_fill_manual(values=c("#9999a0", "#1b2aed", "#f6131e")) +
      geom_jitter(color="black", size=2.5, alpha=0.7, width = 0.2) + theme_classic() +
      theme(legend.position="none",
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
            axis.title.y = element_blank())
    
    my_comparisons = list(c("Control", "Low"), c("Control", "High"))
    
    ll[[i - 3]] <- p + stat_compare_means(data = dat, comparisons = my_comparisons,
                           label = "p.signif", method= "t.test")
  }
  return(ll)
  
}

E3M <- scfaPlot(ApoE3M)[-5]
E3F <- scfaPlot(ApoE3F)[-5]
E4M <- scfaPlot(ApoE4M)[-5]
E4F <- scfaPlot(ApoE4F)[-5]

p1 <- do.call("grid.arrange", c(E3M, ncol=5, top = "ApoE3 Males"))
p2 <- do.call("grid.arrange", c(E4M, ncol=5, top = "ApoE4 Males"))
p3 <- do.call("grid.arrange", c(E3F, ncol=5, top = "ApoE3 Females"))
p4 <- do.call("grid.arrange", c(E4F, ncol=5, top = "ApoE4 Females"))

pdf("SCFA_all_M.pdf", height = 12, width = 12)
ggarrange(p1, p2, nrow = 2, labels = "AUTO")
dev.off()

pdf("SCFA_all_F.pdf", height = 12, width = 12)
ggarrange(p3, p4, nrow = 2, labels = "AUTO")
dev.off()

# Quick stats

t.test(Lactate~Chemical, data = ApoE4M[ApoE4M$Chemical%in%c("Control", "High"),c(3,11)])
t.test(ApoE4F[ApoE4F$Chemical%in%c("Control", "Low"),c(5)]~Chemical, data = ApoE4F[ApoE4F$Chemical%in%c("Control", "Low"),c(3,5)])

# Boxplots for supplementary figures
stool <- read.csv("StoolBacteria.csv", header = T)
stool$SexGeno <- substr(stool$sample, 1, 3)
stool$SexGeno <- factor(stool$SexGeno, levels = c("ME3", "ME4", "FE3", "FE4"))
stool$Cd <- substr(stool$sample, 4, 5)
stool$Cd <- factor(stool$Cd, levels = c("CN", "LO", "HI"))

Universal <- ggplot(data = stool, aes(x=SexGeno, y = unc.ddcq, fill = Cd)) + 
  geom_boxplot() + ggtitle("Universal Bacteria") + 
  ylab("ddCQ/5 ng DNA") + 
  scale_fill_manual(values=c("#9999a0", "#1b2aed", "#f6131e")) +
  geom_jitter(color="black", size=1, alpha=0.7, width = 0.2) + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12))

amucin <- ggplot(data = stool, aes(x=SexGeno, y = akk_ddcq, fill = Cd)) + 
  geom_boxplot() + ggtitle("A. muciniphila") + 
  ylab("ddCQ/5 ng DNA") + 
  scale_fill_manual(values=c("#9999a0", "#1b2aed", "#f6131e")) +
  geom_jitter(color="black", size=1, alpha=0.7, width = 0.2) + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12))

pdf("Figures/Supplemental Figures/UniversalBacteria.pdf", height = 4, width = 8)
ggarrange(Universal, amucin, common.legend = T)
dev.off()

cyp <- read.csv("CypGenes.csv", header = T, stringsAsFactors = F)
cyp$SexGeno <- substr(cyp$Samples, 1, 3)
cyp$SexGeno <- factor(cyp$SexGeno, levels = c("ME3", "ME4", "FE3", "FE4"))
cyp$Cd <- substr(cyp$Samples, 4, 5)
cyp$Cd <- factor(cyp$Cd, levels = c("CN", "LO", "HI"))

makeBox <- function(data, y, title){
  ggplot(data = data, aes(x=SexGeno, y = y, fill = Cd)) + 
    geom_boxplot() + ggtitle(title) + 
    ylab("mRNA (% of Gadph)") + 
    scale_fill_manual(values=c("#9999a0", "#1b2aed", "#f6131e")) +
    geom_jitter(color="black", size=1, alpha=0.7, width = 0.2) + theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 12)) 
}

cyp_list <- list()
for (i in 2:7){
  dat <- cyp[, c(i, 8:9)]
  title <- colnames(cyp)[i]
  cyp_list[[i-1]] <- makeBox(dat, dat[,], title)
}

cyp_list <- lapply(2:7, function(i)
       makeBox(cyp[, c(i, 8:9)], cyp[,i], title))

pdf("Figures/Supplemental Figures/Cyp_qPCR.pdf", height = 9, width = 7)
ggarrange(plotlist = cyp_list, nrow = 3, ncol = 2, common.legend = T)
dev.off()
