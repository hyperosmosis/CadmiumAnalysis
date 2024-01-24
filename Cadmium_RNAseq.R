setwd("~/Documents/Cui/Cadmium-RNAseq/")
counts = read.csv("liver_cadmium.csv", header = T, stringsAsFactors = F)
counts.1 = read.table("liver_test", header = T, sep = "\t", stringsAsFactors = F)
## Clean up counts file

counts$Chr = sapply(strsplit(counts$Chr,split=";") , "[[", 1)
counts$Start = sapply(strsplit(counts$Start,split=";"), function(x) min(x))
counts$End = sapply(strsplit(counts$End,split=";"), function(x) max(x))
counts$Strand = sapply(strsplit(counts$Strand,split=";") , "[[", 1)
rownames(counts) = counts$Gene_ID

treatments = rep(c("M_ApoE3_CN", "M_ApoE3_LowCd", "M_ApoE3_HighCd", 
                   "M_ApoE4_CN", "M_ApoE4_LowCd", "M_ApoE4_HighCd",
                   "F_ApoE3_CN", "F_ApoE3_LowCd", "F_ApoE3_HighCd", 
                   "F_ApoE4_CN", "F_ApoE4_LowCd", "F_ApoE4_HighCd",
                   "ileum_F_ApoE3_CN", "ileum_F_ApoE3_HighCd",
                   "ileum_F_ApoE4_CN", "ileum_F_ApoE4_HighCd"), each = 3)
colnames(counts) = c("Gene_ID", "Start", "End", "Strand", "Length",
                     paste(treatments, rep(c(1:3), 16), sep= "_"))

## Load up DESeq2

library(DESeq2)

cad_design = data.frame(row.names = colnames(counts[,6:53]), 
                        genotype = c(rep(rep(c("ApoE3", "ApoE4"), each = 9), 2), 
                                      rep(c("ApoE3", "ApoE4"), each = 6)),
                        chemical = c(rep(rep(c("none", "low", "high"), each = 3), 4),
                                     rep(rep(c("none", "high"), each = 3), 2)),
                        sex = c(rep(c("M", "F"), each = 18), rep("F", 12)),
                        libType = rep("paired-end",48))
cad_design = cad_design[c(2:48),]

genotype = as.factor(cad_design$genotype)
chemical = as.factor(cad_design$chemical)
sex = as.factor(cad_design$sex)

## Combine factors
cad_design$geno_chem = paste(genotype,chemical, sep = "_")
cad_design$geno_chem_sex = paste(genotype,chemical, sex,sep = "_")

geno_chem = as.factor(cad_design$geno_chem)
geno_chem_sex = as.factor(cad_design$geno_chem_sex)

genotype <- relevel(genotype, ref = "ApoE3")
chemical <- relevel(chemical, ref = "none")
sex <- relevel(sex, ref = "M")
geno_chem <- relevel(geno_chem, ref = "ApoE3_none")
geno_chem_sex <- relevel(geno_chem_sex, ref = "ApoE3_none_M")

coldata <- data.frame(row.names=colnames(counts[,7:53]), genotype, chemical, sex, geno_chem, geno_chem_sex)

## Use DESeq2

dds <- DESeqDataSetFromMatrix(countData=counts[,7:53], colData=coldata, design=~genotype + chemical + sex)
dds <- DESeqDataSetFromMatrix(countData=counts[,7:53], colData=coldata, design=~0 + geno_chem_sex)

dds <- estimateSizeFactors(dds)
filter <- rowSums( counts(dds, normalized=TRUE) >= 4 ) >= 2 
dds <- dds[filter,]
dds <- DESeq(dds)

# Contrasts!

diffCondition <- function(factor, x,y){ # Factor must be in quotes!
dds <- results(dds, contrast=c(factor,x,y))
  dds = as.data.frame(dds)
  dds = dds[complete.cases(dds),]
  sig = dds[dds$padj <= 0.05,]
  unsig = dds[dds$padj > 0.05,]
  sig = sig[order(sig$padj),]
  return(sig)
}

ApoE3_M_Low = diffCondition("geno_chem_sex", "ApoE3_low_M","ApoE3_none_M")
ApoE3_M_High = diffCondition("geno_chem_sex", "ApoE3_high_M","ApoE3_none_M")
ApoE4_M_Low = diffCondition("geno_chem_sex", "ApoE4_low_M","ApoE4_none_M")
ApoE4_M_High = diffCondition("geno_chem_sex", "ApoE4_high_M","ApoE4_none_M")

ApoE3_F_Low = diffCondition("geno_chem_sex", "ApoE3_low_F","ApoE3_none_F")
ApoE3_F_High = diffCondition("geno_chem_sex", "ApoE3_high_F","ApoE3_none_F")
ApoE4_F_Low = diffCondition("geno_chem_sex", "ApoE4_low_F","ApoE4_none_F")
ApoE4_F_High = diffCondition("geno_chem_sex", "ApoE4_high_F","ApoE4_none_F")

save(ApoE3_M_Low, ApoE3_M_High, ApoE4_M_Low, ApoE4_M_High,
     ApoE3_F_Low, ApoE3_F_High, ApoE4_F_Low, ApoE4_F_High, file = "ApoE_comparisons.Rdata")

genes = c(rownames(rbind(ApoE3_M_Low, ApoE3_M_High, ApoE4_M_Low, ApoE4_M_High,
                         ApoE3_F_Low, ApoE3_F_High, ApoE4_F_Low, ApoE4_F_High)))

write.table(genes, file = "genes.txt", sep = "\t", row.names = F)


## Heatmaps of genes of interest

library(made4)
library(gplots)
library(dplyr)

genes = read.csv("GenesofInterest.csv", header = F)
transporters = counts[counts$Gene_ID %in% genes[genes[,2]=="transporters",][,1],][,6:41]
inflam = counts[counts$Gene_ID %in% genes[genes[,2]=="inflammation",][,1],][,6:41]
drug = counts[counts$Gene_ID %in% genes[genes[,2]=="drug metabolizing enzymes",][,1],][,6:41]

heatplot(transporters, dend = "row")
heatplot(inflam, dend = "row")
heatplot(drug, dend = "row")

save(transporters, inflam, drug, file = "GenesofInterest.Rdata")

# Find mean of groups

helpMean <- function(data){
  dat = as.data.frame(t(data))
  dat$Group = substr(rownames(dat),1,nchar(rownames(dat))-2)
  
  new_dat <- dat %>%
    group_by(Group) %>%
    summarize_if(is.numeric, mean)
  
  new_dat = as.data.frame(t(new_dat))
  colnames(new_dat) = unique(dat$Group)
  final <- mutate_all(new_dat[-1,], function(x) as.numeric(as.character(x)))
  rownames(final) = colnames(dat)[1:length(colnames(dat))-1]
  
  return(final)
}

m_inflam <- helpMean(inflam)
m_transport <- helpMean(transporters)
m_drug <- helpMean(drug)

pdf("m_inflam.pdf", height = 3.5, width = 8)
heatplot(m_inflam, dend = "row")
dev.off()

pdf("m_transporters.pdf", height = 6, width = 8)
heatplot(m_transport, dend = "row")
dev.off()

pdf("m_drug.pdf", height = 6, width = 8)
heatplot(m_drug, dend = "row")
dev.off()


load("ApoE_comparisons.Rdata")

getSig <- function(data){
  p.val <- matrix(FALSE, nrow = nrow(data), ncol = ncol(data))
  rownames(p.val) <- rownames(data)
  
  p.val[,2] <- rownames(data)%in%rownames(ApoE3_M_Low)
  p.val[,3] <- rownames(data)%in%rownames(ApoE3_M_High)
  p.val[,5] <- rownames(data)%in%rownames(ApoE4_M_Low)
  p.val[,6] <- rownames(data)%in%rownames(ApoE4_M_High)
  p.val[,8] <- rownames(data)%in%rownames(ApoE3_F_Low)
  p.val[,9] <- rownames(data)%in%rownames(ApoE3_F_High)
  p.val[,11] <- rownames(data)%in%rownames(ApoE4_F_Low)
  p.val[,12] <- rownames(data)%in%rownames(ApoE4_F_High)
  
  p.val <- matrix(ifelse(as.vector(p.val) == TRUE, 
                            "*", ""), 
                     ncol = ncol(data), nrow = nrow(data))
  return(p.val)
}

t(apply(m_inflam, 1, getZ))
