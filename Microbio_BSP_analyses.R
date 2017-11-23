# Body size phenotypes comprehensively assess cardiometabolic risk and refine the association between obesity and gut microbiot
# de la Cuesta-Zuluaga et al., 2017
# International Journal of Obesity
# doi: 10.1038/ijo.2017.281
# 
# BMI: body mass index
# CHS: Cardiometabolic Health Status
# BSP: Body Size Phenotype
#
# clear workspace
rm(list = ls())
set.seed(5600)

# Libraries
library(GUniFrac) # UniFracs
library(ade4) # s.class
library(phytools) # loads tree - read.newick
library(ggplot2) # graphs
library(car) # Anova
library(cowplot) # plot_grid
library(ggdendro) # Dendrogram
library(BiodiversityR) # Alpha diversity
library(reshape2) # melt
library(plyr) # modify dataframes
library(doBy) # summary tables
library(asbio) # Confidence interval for median
library(qvalue) # FDR correction

# Functions ----
# Quasipoisson GLM for the association between microbiota and some trait (e.g., BMI, CHS, BSP) ----
#
# Implements a quasipoisson generalized linear model to evaluate the association between OTU abundance an a given factor
microbiota_glm = function(otu_abs_freq, otu_rel_freq, variable, qvalue_cutoff, median_otu_abund, taxonomy, trans = NA){
  # Loads OTU tables and selects OTUs according to median_otu_abund
  median_abund = apply(otu_rel_freq, MARGIN = 1, FUN = median)
  abundant_otus = t(otu_abs_freq)[median_abund >= (median_otu_abund/100), ]
  # Transforms the factor to log, asin_sqrt or none
  if(trans == "asin"){
    y.trans = y.trans = asin(sqrt(variable/100))
  } else if (trans == "log") {
    y.trans = log(variable)
  } else {
    y.trans = variable
  }
  # Runs glm on each OTU
  lm_meta = apply(abundant_otus, MARGIN = 1, FUN = function(x) glm_otu = glm(x ~ y.trans, family = quasipoisson, maxit = 100))
  #
  # Creates an output table with taxonomy, coefficients and p- and q-values
  coef_lm = lapply(X = lm_meta, FUN = function(x) x$coefficients[2])
  anova_lm_meta = lapply(X = lm_meta, FUN = function(x) Anova(x, type = "II"))
  otus_meta = sapply(anova_lm_meta, FUN = function(x) x$`Pr(>Chisq)`[1])
  fdr_corrected = qvalue(otus_meta)
  sig_otus_meta = subset(anova_lm_meta, fdr_corrected$qvalues <= qvalue_cutoff)
  sig_otus_coef = subset(coef_lm, fdr_corrected$qvalues <= qvalue_cutoff)
  result_lm = sapply(X = sig_otus_meta, FUN = function(x) x$`Pr(>Chisq)`[1])
  result_lm = data.frame(Tax = taxonomy[row.names(taxonomy) %in% names(result_lm),2], Coef = unlist(sig_otus_coef, use.names = F), p = result_lm, q = fdr_corrected$qvalues[fdr_corrected$qvalues <= qvalue_cutoff])
  result_lm
}

# These two functions are used to summarize metadata into a table, using summaryBy
m_sd = function(x){
  m = round(mean(x, na.rm = T), 1)
  s = round(sd(x, na.rm = T), 1)
  paste(m, "?", s)
}

me_iqr = function(x){
  me = round(median(x, na.rm = T)*100, 1)
  iqr = round(IQR(x, na.rm = T)*100, 1)
  paste(me, "?", iqr)
}

prop_vari = function(x, y){
  table_a =  table(y, x)
  table_b = t(round(prop.table(table_a, margin = 1)*100, 1))
  table_b
}

# Load data ----
setwd(dir = ".../bsp")

#* Problems with phytools in linux.
# Source code from read.newick is in the file with same name
# Loads read.newick
#source("read.newick.R")

# Metadata
microbio.meta = read.table(file = "microbio_selected.meta", header = T, sep = "\t", dec = ",", row.names = 1)

# OTU table
microbio.otus = read.table(file = "microbio_selected.otus", header = T, sep = "\t", row.names = 1)

# Phylogenetic tree
microbio.tree = read.newick(file = "microbio_selected.tre")

# Taxonomy
microbio.taxonomy = read.table("microbio_selected.taxonomy", sep = "\t", row.names = 1, header = T)

# Co-abundance groups of microbes (CAGs)
microbio.cags = read.table("microbio_selected.cags", sep = "\t", row.names = 1, header = T)

# Phylum abundance
microbio.phyla = t(read.table("phyla.txt", sep = "\t", row.names = 1, header = T))

# CHS and BSP data
classification_table = read.table("bsp_classification.txt", sep = "\t", row.names = 1, header = T)

# OTU count rarefaction
# By default, minimum number rowSums(microbio.rare)
microbio.rare = Rarefy(microbio.otus)$otu.tab.rff

# Relative frecuencies
# Obtained from the non-rarefied absolute frecuencies
microbio.relative = t(microbio.otus/rowSums(microbio.otus))

# Removes replicate samples from the dataset
replicate_samples = c("MI_093_H12", "MI_008_H2", "MI_130_H2", "MI_198_H2", "MI_458_H2")
pos_replicates = c(9, 95, 132, 201, 445)
microbio.meta = microbio.meta[-pos_replicates,]
microbio.otus = microbio.otus[-pos_replicates,]
microbio.rare = microbio.rare[-pos_replicates,]
microbio.relative = microbio.relative[,-pos_replicates]
microbio.phyla = microbio.phyla[-pos_replicates,]

# Appending CHS and BSP to the metadata table
# Only include samples with microbiota data
microbio.meta$bsp = classification_table[rownames(microbio.otus), "bsp_class"]
microbio.meta$chs = classification_table[rownames(microbio.otus), "chs_class"]

# Level of factors
microbio.meta$bmi_class = factor(microbio.meta$bmi_class, levels = c("Normoweight", "Overweight", "Obese"))

microbio.meta$chs = factor(microbio.meta$chs, levels = c("Healthy", "Abnormal"))

microbio.meta$bsp = factor(microbio.meta$bsp, levels = c("Normoweight-Healthy", "Normoweight-Abnormal", "Overweight-Healthy", "Overweight-Abnormal", "Obese-Healthy", "Obese-Abnormal"))
count(microbio.meta$bsp) # Figure 1

# Clean subset ----
# This subset excludes smokers and medicament users
# To run the code with this subset, uncomment the following lines

# microbio.clean = rownames(microbio.meta[microbio.meta$smoker == "No" & microbio.meta$medicament == "No",])
# 
# microbio.otus = microbio.otus[microbio.clean,]
# microbio.rare = microbio.rare[microbio.clean,]
# microbio.relative = microbio.relative[,microbio.clean]
# microbio.phyla = microbio.phyla[microbio.clean,]
# microbio.meta = microbio.meta[microbio.clean,]
# classification_table = classification_table[microbio.clean,]
# microbio.cags = microbio.cags[microbio.clean,]


# Descriptive statistics of biochemical and anthropometric variables (Table 1) ----
continuous_variables = summaryBy(age + bmi + waist + body_fat + systolic_bp + diastolic_bp + cholesterol + HDL + LDL + triglycerides + hsCRP + glucose + insulin + glycosylated_hg + HOMA_IR + adiponectin + calories + mets ~ bsp, data = microbio.meta, FUN = m_sd)

prop_vari(microbio.meta$sex, microbio.meta$bsp)
prop_vari(microbio.meta$smoker, microbio.meta$bsp)
prop_vari(microbio.meta$medicament, microbio.meta$bsp)

# Description of abnormalities
# Observed abnormalities
abnormalities_table = classification_table[rownames(microbio.otus), 3:8]

sex_abnormalities = data.frame(Abnormalities = rowSums(abnormalities_table), Sex = microbio.meta$sex)

obs_abnormalities_plot = ggplot(sex_abnormalities, aes(x = Abnormalities)) + geom_bar(aes(fill = microbio.meta$sex), position = "dodge") + theme_bw() + labs(x = "Observed abnormalities", y = "Number of individuals") + scale_fill_discrete(name = "Sex")


# Frecuency of each abnormality
prop_abnormalities_total = colSums(abnormalities_table)/sum(colSums(abnormalities_table))
prop_abnormalities_total = data.frame(Abnormality = c ("Blood pressure", "Triglycerides", "HDL", "Glucose", "Inflammation", "Insulin resistance"), Proportion = prop_abnormalities_total, Sex = rep("Total", 6))
rownames(prop_abnormalities_total) = NULL

prop_abnormalities_male = colSums(abnormalities_table[microbio.meta$sex == "Male",])/sum(colSums(abnormalities_table[microbio.meta$sex == "Male",]))
prop_abnormalities_male = data.frame(Abnormality = c ("Blood pressure", "Triglycerides", "HDL", "Glucose", "Inflammation", "Insulin resistance"), Proportion = prop_abnormalities_male, Sex = rep("Male", 6))
rownames(prop_abnormalities_male) = NULL

prop_abnormalities_female = colSums(abnormalities_table[microbio.meta$sex == "Female",])/sum(colSums(abnormalities_table[microbio.meta$sex == "Female",]))
prop_abnormalities_female = data.frame(Abnormality = c ("Blood pressure", "Triglycerides", "HDL", "Glucose", "Inflammation", "Insulin resistance"), Proportion = prop_abnormalities_female, Sex = rep("Female", 6))
rownames(prop_abnormalities_female) = NULL

full_prop_abnormalities = rbind(prop_abnormalities_total, prop_abnormalities_male, prop_abnormalities_female)

freq_abnormalities_plot = ggplot(full_prop_abnormalities, aes(x = Sex, y = Proportion)) + geom_col(aes(fill = Abnormality)) + theme_bw() + labs(x = "Sex", y = "Relative frequency") 

# Figure 2
plot_grid(obs_abnormalities_plot, freq_abnormalities_plot, nrow = 1, ncol = 2, labels = c("A", "B"), rel_widths = c(1.7,1))


# Alpha diversity analysis -----
richness_all = diversityresult(x = microbio.rare, index = "richness", method = "each site")
shannon_all = diversityresult(x = microbio.rare, index = "Shannon", method = "each site")

alpha_data = data.frame(richness = richness_all$richness, shannon = shannon_all$Shannon, cont_BMI = microbio.meta$bmi, cat_BMI = microbio.meta$bmi_class, CHS = microbio.meta$chs, BSP = microbio.meta$bsp)

# Richness
aov_richness_BMI = aov(richness_all$richness ~ microbio.meta$bmi_class)
summary(aov_richness_BMI)

aov_richness_CHS = aov(richness_all$richness ~ microbio.meta$chs)
summary(aov_richness_CHS)

aov_richness_BSP = aov(richness_all$richness ~ microbio.meta$bsp)
summary(aov_richness_BSP)

lm_richness_bmi = lm(richness_all$richness ~ microbio.meta$bmi)
Anova(lm_richness_bmi)

lm_richness_chs = lm(richness_all$richness ~ microbio.meta$chs)
Anova(lm_richness_chs)

lm_richness_interaction = lm(richness_all$richness ~ microbio.meta$bmi:microbio.meta$chs)
Anova(lm_richness_interaction)
summary(lm_richness_interaction)

lm_richness_bsp = lm(richness_all$richness ~ microbio.meta$bsp)
Anova(lm_richness_bsp)

# Shannon diversity index
aov_shannon_BMI = aov(shannon_all$Shannon ~ microbio.meta$bmi_class)
summary(aov_shannon_BMI)

aov_shannon_CHS = aov(shannon_all$Shannon ~ microbio.meta$chs)
summary(aov_shannon_CHS)

aov_shannon_BSP = aov(shannon_all$Shannon ~ microbio.meta$bsp)
summary(aov_shannon_BSP)

lm_shannon_bmi = lm(shannon_all$Shannon ~ microbio.meta$bmi)
Anova(lm_shannon_bmi)

lm_shannon_chs = lm(shannon_all$Shannon ~ microbio.meta$chs)
Anova(lm_shannon_chs)

lm_shannon_interaction = lm(shannon_all$Shannon ~ microbio.meta$bmi:microbio.meta$chs)
Anova(lm_shannon_interaction)

lm_shannon_bsp = lm(shannon_all$Shannon ~ microbio.meta$bsp)
Anova(lm_shannon_bsp)


# Plots
# Richness
# Raw
lm_richness_BMI = ggplot(alpha_data, aes(x=cont_BMI, y=richness)) + geom_vline(xintercept = c(18.5, 25, 30), linetype="dashed") + geom_point() + geom_smooth(method='lm', se=F, color="grey50") + theme_bw() + labs(x = "BMI", y = "Species richness")

# BMI:CHS
lm_richness_BMIxCHS = ggplot(alpha_data, aes(x=cont_BMI, y=richness, color=CHS, shape=CHS)) + geom_vline(xintercept = c(18.5, 25, 30), linetype="dashed") + geom_point() + scale_color_grey(start=0.5, end=0.2) + geom_smooth(method='lm', se=F, fullrange=TRUE) + theme_bw() + labs(x = "BMI", y = "Species richness") + theme(legend.justification=c(1,0), legend.position=c(1,0))

# shannon
# Raw
lm_shannon_BMI = lm_shannon_BMI = ggplot(alpha_data, aes(x=cont_BMI, y=shannon)) + geom_vline(xintercept = c(18.5, 25, 30), linetype="dashed") + geom_point() + geom_smooth(method='lm', se=F, color="grey50") + theme_bw() + labs(x = "BMI", y = "Shannon Index")
# BMI:CHS
lm_shannon_BMIxCHS = ggplot(alpha_data, aes(x=cont_BMI, y=shannon, color=CHS, shape=CHS)) + geom_vline(xintercept = c(18.5, 25, 30), linetype="dashed") + geom_point() + scale_color_grey(start=0.5, end=0.2) + geom_smooth(method='lm', se=F, fullrange=TRUE) + theme_bw() + labs(x = "BMI", y = "Shannon Index") + theme(legend.justification=c(1,0), legend.position=c(1,0))

# Figure 3
plot_grid(lm_richness_BMI, lm_richness_BMIxCHS, lm_shannon_BMI, lm_shannon_BMIxCHS, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

# Summary tables (Table S4)
tab_alfa_bmi = summaryBy(richness_all$richness + shannon_all$Shannon ~ bmi_class, data = microbio.meta, FUN = m_sd)
tab_alfa_chs = summaryBy(richness_all$richness + shannon_all$Shannon ~ chs, data = microbio.meta, FUN = m_sd)
tab_alfa_bsp = summaryBy(richness_all$richness + shannon_all$Shannon ~ bsp, data = microbio.meta, FUN = m_sd)


# Beta-diversity analyses ----
# Computation of UniFrac distances

unifracs <- GUniFrac(microbio.rare, microbio.tree, alpha=c(0, 0.5, 1))$unifracs

dw <- unifracs[, , "d_1"]   # Weighted UniFrac
du <- unifracs[, , "d_UW"]  # Unweighted UniFrac    

# Distance matrices
dw.dist = as.dist(dw)
du.dist = as.dist(du)

# PERMANOVA
# weighted
adonis(dw.dist ~ microbio.meta$bmi_class)
adonis(dw.dist ~ microbio.meta$chs)
adonis(dw.dist ~ microbio.meta$bsp)

cmd_dw = cmdscale(dw.dist, k = 3, eig = T)
pcoa_dw = as.data.frame(cmd_dw$points)
var_pco1_dw = round(cmd_dw$eig[1]/sum(cmd_dw$eig), 4)* 100
var_pco2_dw = round(cmd_dw$eig[2]/sum(cmd_dw$eig), 4)* 100
var_pco3_dw = round(cmd_dw$eig[3]/sum(cmd_dw$eig), 4)* 100

graph.bmi.w = ggplot(cbind(pcoa_dw, microbio.meta$bmi_class), aes(V1, V2)) +   geom_point(aes(color = microbio.meta$bmi_class)) + stat_ellipse(type = "t", mapping = aes(color = microbio.meta$bmi_class), level = 0.75) + labs(x = paste("PCo1",var_pco1_dw, "%"), y = paste("PCo2",var_pco2_dw, "%"), colour = "BMI") + theme(legend.key.size = unit(0.5, "cm"))

graph.chs.w = ggplot(cbind(pcoa_dw, microbio.meta$chs), aes(V1, V2)) +   geom_point(aes(color = microbio.meta$chs)) + stat_ellipse(type = "t", mapping = aes(color = microbio.meta$chs), level = 0.75) + labs(x = paste("PCo1",var_pco1_dw, "%"), y = paste("PCo2",var_pco2_dw, "%"), colour = "CHS") + theme(legend.key.size = unit(0.5, "cm"))

graph.feno.w = ggplot(cbind(pcoa_dw, microbio.meta$bsp), aes(V1, V2)) +  geom_point(aes(color = microbio.meta$bsp)) + stat_ellipse(type = "t", mapping = aes(color = microbio.meta$bsp), level = 0.75) + labs(x = paste("PCo1",var_pco1_dw, "%"), y = paste("PCo2",var_pco2_dw, "%"), colour = "BSP") +   theme(legend.key.size = unit(0.5, "cm"))

# unweighted
adonis(du.dist ~ microbio.meta$bmi_class)
adonis(du.dist ~ microbio.meta$chs)
adonis(du.dist ~ microbio.meta$bsp)

cmd_du = cmdscale(du.dist, k = 3, eig = T)
pcoa_du = as.data.frame(cmd_du$points)
var_pco1_du = round(cmd_du$eig[1]/sum(cmd_du$eig), 4)* 100
var_pco2_du = round(cmd_du$eig[2]/sum(cmd_du$eig), 4)* 100
var_pco3_du = round(cmd_du$eig[3]/sum(cmd_du$eig), 4)* 100

graph.bmi.u = ggplot(cbind(pcoa_du, microbio.meta$bmi_class), aes(V1, V2)) +   geom_point(aes(color = microbio.meta$bmi_class)) + stat_ellipse(type = "t", mapping = aes(color = microbio.meta$bmi_class), level = 0.75) + labs(x = paste("PCo1",var_pco1_du, "%"), y = paste("PCo2",var_pco2_du, "%"), colour = "BMI") + theme(legend.key.size = unit(0.5, "cm"))

graph.chs.u = ggplot(cbind(pcoa_du, microbio.meta$chs), aes(V1, V2)) +   geom_point(aes(color = microbio.meta$chs)) + stat_ellipse(type = "t", mapping = aes(color = microbio.meta$chs), level = 0.75) + labs(x = paste("PCo1",var_pco1_du, "%"), y = paste("PCo2",var_pco2_du, "%"), colour = "CHS") + theme(legend.key.size = unit(0.5, "cm"))

graph.feno.u = ggplot(cbind(pcoa_du, microbio.meta$bsp), aes(V1, V2)) +  geom_point(aes(color = microbio.meta$bsp)) + stat_ellipse(type = "t", mapping = aes(color = microbio.meta$bsp), level = 0.75) + labs(x = paste("PCo1",var_pco1_du, "%"), y = paste("PCo2",var_pco2_du, "%"), colour = "BSP") +   theme(legend.key.size = unit(0.5, "cm"))

plot_grid(graph.bmi.w, graph.chs.w, graph.feno.w, graph.bmi.u, graph.chs.u, graph.feno.u, nrow = 2, ncol = 3, labels = c("A", "C", "E", "B", "D", "F"))

# Bacteroidetes-Firmicutes analysis ----
# Relative abundances and abundance ratio
bacteroidetes = microbio.phyla[,"k__Bacteria;p__Bacteroidetes"]
rownames(bacteroidetes) = NULL

firmicutes = microbio.phyla[,"k__Bacteria;p__Firmicutes"]
rownames(firmicutes) = NULL

bf_ratio = bacteroidetes/firmicutes

# Tests
# BMI
kruskal.test(firmicutes ~ microbio.meta$bmi_class)
kruskal.test(bacteroidetes ~ microbio.meta$bmi_class)
kruskal.test(bf_ratio ~ microbio.meta$bmi_class)

aggregate(firmicutes ~ microbio.meta$bmi_class, FUN = me_iqr)
aggregate(bacteroidetes ~ microbio.meta$bmi_class, FUN = me_iqr)
aggregate(bf_ratio ~ microbio.meta$bmi_class, FUN = me_iqr)

# CHS
kruskal.test(firmicutes ~ microbio.meta$chs)
kruskal.test(bacteroidetes ~ microbio.meta$chs)
kruskal.test(bf_ratio ~ microbio.meta$chs)

aggregate(firmicutes ~ microbio.meta$chs, FUN = me_iqr)
aggregate(bacteroidetes ~ microbio.meta$chs, FUN = me_iqr)
aggregate(bf_ratio ~ microbio.meta$chs, FUN = me_iqr)

# BSP
kruskal.test(firmicutes ~ microbio.meta$bsp)
kruskal.test(bacteroidetes ~ microbio.meta$bsp)
kruskal.test(bf_ratio ~ microbio.meta$bsp)

aggregate(firmicutes ~ microbio.meta$bsp, FUN = me_iqr)
aggregate(bacteroidetes ~ microbio.meta$bsp, FUN = me_iqr)
aggregate(bf_ratio ~ microbio.meta$bsp, FUN = me_iqr)

# Co-abundance group (CAG) analysis ----
cags_bsp = cbind(microbio.cags, chs = microbio.meta$chs, bmi_class = microbio.meta$bmi_class, bsp = microbio.meta$bsp)
cags_melt_bsp = melt(cags_bsp, measure.vars = 1:5, id.vars = 6:8)
names(cags_melt_bsp)[names(cags_melt_bsp) == 'variable'] = 'cag'

# Kruskal test (Figure 4)
# BMI
aggregate(microbio.cags$Prevotella ~ microbio.meta$bmi_class, FUN = median)
aggregate(microbio.cags$Lachnospiraceae ~ microbio.meta$bmi_class, FUN = median)
aggregate(microbio.cags$Pathobiont ~ microbio.meta$bmi_class, FUN = median)
aggregate(microbio.cags$Akkermansia ~ microbio.meta$bmi_class, FUN = median)
aggregate(microbio.cags$Ruminococcaceae ~ microbio.meta$bmi_class, FUN = median)

kruskal.test(microbio.cags$Prevotella ~ microbio.meta$bmi_class)
kruskal.test(microbio.cags$Lachnospiraceae ~ microbio.meta$bmi_class)
kruskal.test(microbio.cags$Pathobiont ~ microbio.meta$bmi_class)
kruskal.test(microbio.cags$Akkermansia ~ microbio.meta$bmi_class)
kruskal.test(microbio.cags$Ruminococcaceae ~ microbio.meta$bmi_class)

box_bmi_cag = ggplot(cags_melt_bsp, aes(x = bmi_class, y = value)) + geom_boxplot(notch = T) + facet_grid(. ~ cag) + theme_grey()  + theme(axis.text.x=element_text(angle=45, hjust=1)) + labs(x = "BMI", y = "Relative abundance")

# CHS
aggregate(microbio.cags$Prevotella ~ microbio.meta$chs, FUN = median)
aggregate(microbio.cags$Lachnospiraceae ~ microbio.meta$chs, FUN = median)
aggregate(microbio.cags$Pathobiont ~ microbio.meta$chs, FUN = median)
aggregate(microbio.cags$Akkermansia ~ microbio.meta$chs, FUN = median)
aggregate(microbio.cags$Ruminococcaceae ~ microbio.meta$chs, FUN = median)

kruskal.test(microbio.cags$Prevotella ~ microbio.meta$chs)
kruskal.test(microbio.cags$Lachnospiraceae ~ microbio.meta$chs)
kruskal.test(microbio.cags$Pathobiont ~ microbio.meta$chs)
kruskal.test(microbio.cags$Akkermansia ~ microbio.meta$chs)
kruskal.test(microbio.cags$Ruminococcaceae ~ microbio.meta$chs)

box_chs_cag = ggplot(cags_melt_bsp, aes(x = chs, y = value)) + geom_boxplot(notch = T) + facet_grid(. ~ cag) + theme_grey() + theme(axis.text.x=element_text(angle=45, hjust=1)) + labs(x = "CHS", y = "Relative abundance")

# BSP
aggregate(microbio.cags$Prevotella ~ microbio.meta$bsp, FUN = median)
aggregate(microbio.cags$Lachnospiraceae ~ microbio.meta$bsp, FUN = median)
aggregate(microbio.cags$Pathobiont~ microbio.meta$bsp, FUN = median)
aggregate(microbio.cags$Akkermansia ~ microbio.meta$bsp, FUN = median)
aggregate(microbio.cags$Ruminococcaceae ~ microbio.meta$bsp, FUN = median)

kruskal.test(microbio.cags$Prevotella ~ microbio.meta$bsp)
kruskal.test(microbio.cags$Lachnospiraceae ~ microbio.meta$bsp)
kruskal.test(microbio.cags$Pathobiont ~ microbio.meta$bsp)
kruskal.test(microbio.cags$Akkermansia ~ microbio.meta$bsp)
kruskal.test(microbio.cags$Ruminococcaceae ~ microbio.meta$bsp)

box_bsp_cag = ggplot(cags_melt_bsp, aes(x = bsp, y = value)) + geom_boxplot(notch = T) + facet_grid(. ~ cag) + theme_grey() + theme(axis.text.x=element_text(angle=45, hjust=1)) + labs(x = "BSP", y = "Relative abundance")

# Figure 4
plot_grid(box_bmi_cag, box_chs_cag, box_bsp_cag, nrow = 3, ncol = 1, labels = c("A", "B", "C"))


# OTU-level analysis
# Tables for the Lefse analysis
table_lefse_bmi = data.frame(BMI = microbio.meta$bmi_class)
table_lefse_bmi = rbind(t(table_lefse_bmi), t(microbio.rare), deparse.level = 2)
#write.table(table_lefse_bmi, "lefse_abnormalities_bmi.txt", sep = "\t", quote = F)
lefse_OTUs_bmi = c("Otu00011", "Otu00331", "Otu00026", "Otu00021", "Otu00025", "Otu00247", "Otu00065", "Otu00100", "Otu00102", "Otu00209", "Otu00062", "Otu00422", "Otu00050", "Otu00352", "Otu00098", "Otu00142", "Otu00037", "Otu00169", "Otu00396", "Otu00084", "Otu00092", "Otu00043", "Otu00111", "Otu00067", "Otu00186", "Otu00208", "Otu00387", "Otu00147", "Otu00230", "Otu00505", "Otu00622", "Otu00151", "Otu00540", "Otu00944", "Otu00271", "Otu00080", "Otu00438", "Otu00211", "Otu00362", "Otu00439", "Otu00166", "Otu00210", "Otu00149", "Otu00289", "Otu00550", "Otu00004", "Otu00003", "Otu00029", "Otu00024", "Otu00059", "Otu00053", "Otu00075", "Otu00179", "Otu00334", "Otu00162", "Otu00600")

lefse_OTUs_nonobese = c("Otu00011", "Otu00012", "Otu00331", "Otu00021", "Otu00026", "Otu00025", "Otu00102", "Otu00050", "Otu00062", "Otu00422", "Otu00352", "Otu00396", "Otu00084", "Otu00111", "Otu00208", "Otu00043", "Otu00387", "Otu00211", "Otu00540", "Otu00271", "Otu00438", "Otu00166", "Otu00439", "Otu00080", "Otu00505", "Otu00874", "Otu00289", "Otu00999", "Otu00362", "Otu00374", "Otu00158", "Otu00004", "Otu00003", "Otu00029", "Otu00059", "Otu00304", "Otu00053", "Otu00037", "Otu00142", "Otu00179", "Otu00083", "Otu00075", "Otu00162", "Otu00334", "Otu00141", "Otu00270", "Otu00908")

writeClipboard(as.character(microbio.taxonomy[lefse_OTUs_bmi,2]))
writeClipboard(as.character(microbio.taxonomy[lefse_OTUs_nonobese,2]))

table_lefse_chs = data.frame(CHS = microbio.meta$chs)
table_lefse_chs = rbind(t(table_lefse_chs), t(microbio.rare), deparse.level = 2)
#write.table(table_lefse_chs, "lefse_abnormalities_chs.txt", sep = "\t", quote = F)
lefse_OTUs_chs = c("Otu00032", "Otu00602", "Otu00104", "Otu00118", "Otu00275", "Otu00011", "Otu00007", "Otu00012", "Otu00331", "Otu00209", "Otu00026", "Otu00063", "Otu00100", "Otu00098", "Otu00060", "Otu00142", "Otu00103", "Otu00352", "Otu00422", "Otu00093", "Otu00210", "Otu00740", "Otu00341", "Otu00151", "Otu00265", "Otu00325", "Otu00149", "Otu00147", "Otu01246", "Otu03946", "Otu00639", "Otu02197", "Otu00388", "Otu00200", "Otu00999", "Otu00368", "Otu00446", "Otu00874", "Otu01353", "Otu00707", "Otu01003", "Otu02141", "Otu01041")
writeClipboard(as.character(microbio.taxonomy[lefse_OTUs_chs,2]))

table_lefse_bsp = data.frame(BSP = microbio.meta$bsp)
table_lefse_bsp = rbind(t(table_lefse_bsp), t(microbio.rare), deparse.level = 2)
#write.table(table_lefse_bsp, "lefse_abnormalities_bsp.txt", sep = "\t", quote = F)
lefse_OTUs_extremes = c("Otu00011", "Otu00012", "Otu00331", "Otu00028", "Otu00021", "Otu00026", "Otu00063", "Otu00209", "Otu00062", "Otu00102", "Otu00352", "Otu00318", "Otu00098", "Otu00422", "Otu02095", "Otu00111", "Otu00141", "Otu00484", "Otu00151", "Otu00208", "Otu01047", "Otu00438", "Otu00999", "Otu00388", "Otu00845", "Otu00636", "Otu02276", "Otu00210", "Otu00691", "Otu00540", "Otu00917", "Otu00874", "Otu00557", "Otu01276", "Otu00166", "Otu00609", "Otu00084", "Otu01056", "Otu00218", "Otu02232", "Otu00211", "Otu00486", "Otu00613", "Otu00505", "Otu00374", "Otu00830", "Otu00650", "Otu00500", "Otu00347", "Otu00487", "Otu00158", "Otu00694", "Otu00368", "Otu00550", "Otu00271", "Otu00605", "Otu00443", "Otu00760", "Otu01020", "Otu00004", "Otu00003", "Otu00020", "Otu00029", "Otu00059", "Otu00037", "Otu00035", "Otu00162", "Otu00334", "Otu00083", "Otu00118", "Otu00139", "Otu01195", "Otu00348")
writeClipboard(as.character(microbio.taxonomy[lefse_OTUs_extremes,2]))


# Biomarkers of CAGs
aggregate(microbio.relative["Otu00001",] ~ microbio.meta$bmi_class, FUN = me_iqr)
aggregate(microbio.relative["Otu00002",] ~ microbio.meta$bmi_class, FUN = me_iqr)

aggregate(microbio.relative["Otu00008",] ~ microbio.meta$bmi_class, FUN = me_iqr)
aggregate(microbio.relative["Otu00003",] ~ microbio.meta$bmi_class, FUN = me_iqr)

aggregate(microbio.relative["Otu00004",] ~ microbio.meta$bmi_class, FUN = me_iqr)
aggregate(microbio.relative["Otu00009",] ~ microbio.meta$bmi_class, FUN = me_iqr)

aggregate(microbio.relative["Otu00006",] ~ microbio.meta$bmi_class, FUN = me_iqr)
aggregate(microbio.relative["Otu00050",] ~ microbio.meta$bmi_class, FUN = me_iqr)

aggregate(microbio.relative["Otu00011",] ~ microbio.meta$bmi_class, FUN = me_iqr)
aggregate(microbio.relative["Otu00012",] ~ microbio.meta$bmi_class, FUN = me_iqr)

aggregate(microbio.relative["Otu00001",] ~ microbio.meta$chs, FUN = me_iqr)
aggregate(microbio.relative["Otu00002",] ~ microbio.meta$chs, FUN = me_iqr)

aggregate(microbio.relative["Otu00008",] ~ microbio.meta$chs, FUN = me_iqr)
aggregate(microbio.relative["Otu00003",] ~ microbio.meta$chs, FUN = me_iqr)

aggregate(microbio.relative["Otu00004",] ~ microbio.meta$chs, FUN = me_iqr)
aggregate(microbio.relative["Otu00009",] ~ microbio.meta$chs, FUN = me_iqr)

aggregate(microbio.relative["Otu00006",] ~ microbio.meta$chs, FUN = me_iqr)
aggregate(microbio.relative["Otu00050",] ~ microbio.meta$chs, FUN = me_iqr)

aggregate(microbio.relative["Otu00011",] ~ microbio.meta$chs, FUN = me_iqr)
aggregate(microbio.relative["Otu00012",] ~ microbio.meta$chs, FUN = me_iqr)

aggregate(microbio.relative["Otu00001",] ~ microbio.meta$bsp, FUN = me_iqr)
aggregate(microbio.relative["Otu00002",] ~ microbio.meta$bsp, FUN = me_iqr)

aggregate(microbio.relative["Otu00008",] ~ microbio.meta$bsp, FUN = me_iqr)
aggregate(microbio.relative["Otu00003",] ~ microbio.meta$bsp, FUN = me_iqr)

aggregate(microbio.relative["Otu00004",] ~ microbio.meta$bsp, FUN = me_iqr)
aggregate(microbio.relative["Otu00009",] ~ microbio.meta$bsp, FUN = me_iqr)

aggregate(microbio.relative["Otu00006",] ~ microbio.meta$bsp, FUN = me_iqr)
aggregate(microbio.relative["Otu00050",] ~ microbio.meta$bsp, FUN = me_iqr)

aggregate(microbio.relative["Otu00011",] ~ microbio.meta$bsp, FUN = me_iqr)
aggregate(microbio.relative["Otu00012",] ~ microbio.meta$bsp, FUN = me_iqr)


#GLM with BMI and abnormalities

bmi_glm = microbiota_glm(otu_abs_freq = microbio.rare, otu_rel_freq = microbio.relative, taxonomy = microbio.taxonomy, variable = log(microbio.meta$bmi), qvalue_cutoff = 0.1, median_otu_abund = 0.001, trans = "NA")

bmi_glm$rho = apply(microbio.rare[,rownames(bmi_glm)], MARGIN = 2, FUN = function(x) cor(x, log(microbio.meta$bmi), method = "s"))

