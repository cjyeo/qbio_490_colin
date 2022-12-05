"---
title: Final Project (R Section)
authors: Chloe Liu, Kailin Liu, Colin Yeo
---"
###SETTING UP
#install and load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
library(BiocManager)

if (!require("TCGAbiolinks", quietly = TRUE)) 
  BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

if (!require("SummarizedExperiment", quietly = TRUE)) 
  BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

if (!require("maftools", quietly = TRUE)) 
  BiocManager::install("maftools")
library(maftools)

if (!require(ggplot2)){ 
  install.packages("ggplot2")
}
library(ggplot2)

if (!require(survival)){ 
  install.packages("survival")
}
library(survival)

if (!require(survminer)){ 
  install.packages("survminer")
}
library(survminer)

if (!require(maftools)){ 
  install.packages("maftools")
}
library(maftools)

# make outputs folder
#first create a final_project_group4 folder through terminal, then run the following
dir.create("/Users/kathykliu/Desktop/qbio_490_kailin/final_project_group4/outputs")
setwd("/Users/kathykliu/Desktop/qbio_490_kailin/final_project_group4/outputs")


#--------------------------------------------------

# set up rna_query and rna_se
data_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")
GDCdownload(data_query)
data <- GDCprepare(data_query)


# setting up and saving clinical data frames
clinical_query <- GDCquery(project="TCGA-BRCA",data.category="Clinical",file.type="xml")
GDCdownload(clinical_query)

clinical <- GDCprepare_clinic(clinical_query,clinical.info="patient")
#change the bcr_patient_barcode column name, that way the MAF package can read our clinical file
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"


#query in the MAF files 
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query)

maf <- GDCprepare(maf_query)

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)
#--------------------------------------------------
#CREATING oncoplots to see patterns in gene mutations

#the following line will save the oncoplot to the outputs folder
jpeg("oncoplot.jpg")
#creating an oncoplot that shows the top 10 mutated genes
oncoplot(maf = maf_object,
         top = 10) 

#line below will save the oncoplot to the outputs folder
jpeg("oncoplot_specific.jpg")
# creating an oncoplot that only shows the mutations for TP53, CATA3, CDH1, and MAP3K1
oncoplot_specific<- oncoplot(maf = maf_object, genes = c("TP53","GATA3","CDH1","MAP3K1","PIK3CA"))

oncoplot_specific
dev.off()

#--------------------------------------------------
#CREATING somatic interactions plots to see relationships between mutated genes

#the following line will save the somatic interactions plot to the outputs folder
jpeg("somatic_interactions_plot.jpg")
#creating a somatic interactions plot
somaticInteractions(maf = maf_object, top = 25, pvalue = c(0.05, 0.1))

# paired fisher exact tests on our genes of interest
#      gene1  gene2       pValue oddsRatio  00 11  01  10         pAdj
# 1:   TP53  GATA3 6.109106e-14 0.1253992 518  9 116 321 1.414145e-12
# 2:   CDH1   TP53 1.172597e-12 0.1525908 517 11 319 117 2.527150e-11
# 3:   TP53 MAP3K1 7.195723e-07 0.2157246 561  9  73 321 1.450751e-05

#--------------------------------------------------
#SETTING UP FOR survival plots to see relationships between gene mutation and survival

#setting up data frames with only one desired gene
map_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "MAP3K1"]
gata_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "GATA3"]
cdh_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "CDH1"]
tp_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "TP53"]
pik_data <- maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol == "PIK3CA"]

#--------------------------------------------------
#CREATING survival plot using MAP3K1 data

#masking the MAP3K1 data frame for only patients that have MAP3K1 mutations or MAP3K1 mutations and TP53 mutations
map_mask <- ifelse(maf_object@clinical.data$Tumor_Sample_Barcode%in%levels(factor(map_data)),T,F)
map_clinical <- maf_object@clinical.data[map_mask,]
maptp <- intersect(map_data, tp_data)
map_clinical$mutation <- ifelse(map_clinical$Tumor_Sample_Barcode%in%maptp, "MAP3K1 and TP53", "MAP3K1")

#creating the time and event for the survival plot
map_clinical$survival_time <- ifelse(is.na(map_clinical$days_to_death), map_clinical$days_to_last_followup, map_clinical$days_to_death)
map_clinical$death_event <- ifelse(!is.na(map_clinical$vital_status),TRUE,FALSE)

#line below will save the survival plot to the outputs folder
jpeg("map_survival_plot.jpg")
surv_object_map <- Surv(time = as.numeric(map_clinical$survival_time), event = as.numeric(map_clinical$death_event))

# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
map_fit <- surv_fit( surv_object_map ~ map_clinical$mutation,
                     data = map_clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_map = ggsurvplot(map_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_map = survplot_map$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_map
dev.off()

#--------------------------------------------------
#CREATING survival plot using GATA3 data

#masking the GATA3 data frame for only patients that have GATA3 mutations or GATA3 mutations and TP53 mutations
gata_mask <- ifelse(maf_object@clinical.data$Tumor_Sample_Barcode%in%levels(factor(gata_data)),T,F)
gata_clinical <- maf_object@clinical.data[gata_mask,]
gatatp <- intersect(gata_data, tp_data)
gata_clinical$mutation <- ifelse(gata_clinical$Tumor_Sample_Barcode%in%gatatp, "GATA3 and TP53", "GATA3")

gata_clinical$survival_time <- ifelse(is.na(gata_clinical$days_to_death), gata_clinical$days_to_last_followup, gata_clinical$days_to_death)
gata_clinical$death_event <- ifelse(!is.na(gata_clinical$vital_status),TRUE,FALSE)

surv_object_gata <- Surv(time = as.numeric(gata_clinical$survival_time), event = as.logical(gata_clinical$death_event))

#line below will save the survival plot to the outputs folder
jpeg("gata_survival_plot.jpg")
# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
gata_fit <- surv_fit( surv_object_gata ~ gata_clinical$mutation,
                      data = gata_clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_gata = ggsurvplot(gata_fit, 
                           pval=TRUE, 
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                           legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_gata = survplot_gata$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))


KM_plot_gata
dev.off()

#--------------------------------------------------
#CREATING survival plot using CDH1 data

#masking the CDH1 data frame for only patients that have CDH1 mutations or MAP3K1 mutations and TP53 mutations
cdh_mask <- ifelse(maf_object@clinical.data$Tumor_Sample_Barcode%in%levels(factor(cdh_data)),T,F)
cdh_clinical <- maf_object@clinical.data[cdh_mask,]
cdhtp <- intersect(cdh_data, tp_data)
cdh_clinical$mutation <- ifelse(cdh_clinical$Tumor_Sample_Barcode%in%cdhtp, "CDH1 and TP53", "CDH1")

cdh_clinical$survival_time <- ifelse(is.na(cdh_clinical$days_to_death), cdh_clinical$days_to_last_followup, cdh_clinical$days_to_death)
cdh_clinical$death_event <- ifelse(!is.na(cdh_clinical$vital_status),TRUE,FALSE)

surv_object_cdh <- Surv(time = as.numeric(cdh_clinical$survival_time), event = as.logical(cdh_clinical$death_event))

#line below will save the survival plot to the outputs folder
jpeg("cdh_survival_plot.jpg")
# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
cdh_fit <- surv_fit( surv_object_cdh ~ cdh_clinical$mutation,
                     data = cdh_clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_cdh = ggsurvplot(cdh_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_cdh = survplot_cdh$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))


KM_plot_cdh
dev.off()

#--------------------------------------------------
#CREATING survival plot using PIK3CA data

#masking the MAP3K1 data frame for only patients that have MAP3K1 mutations or MAP3K1 mutations and TP53 mutations
pik_mask <- ifelse(maf_object@clinical.data$Tumor_Sample_Barcode%in%levels(factor(pik_data)),T,F)
pik_clinical <- maf_object@clinical.data[pik_mask,]
piktp <- intersect(pik_data, tp_data)
pik_clinical$mutation <- ifelse(pik_clinical$Tumor_Sample_Barcode%in%piktp, "PIK3CA and TP53", "PIK3CA")

#creating the time and event for the survival plot
pik_clinical$survival_time <- ifelse(is.na(pik_clinical$days_to_death), pik_clinical$days_to_last_followup, pik_clinical$days_to_death)
pik_clinical$death_event <- ifelse(!is.na(pik_clinical$vital_status),TRUE,FALSE)

#line below will save the survival plot to the outputs folder
jpeg("pik_survival_plot.jpg")
surv_object_pik <- Surv(time = as.numeric(pik_clinical$survival_time), event = as.numeric(pik_clinical$death_event))

# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
pik_fit <- surv_fit( surv_object_pik ~ pik_clinical$mutation,
                     data = pik_clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_pik = ggsurvplot(pik_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_pik = survplot_pik$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_pik
dev.off()


#--------------------------------------------------
#CREATING boxplots to more clearly see the difference in survival time

#MAP3K1
#line below will save the boxplot to the outputs folder
jpeg("survival_map_mutation_boxplot.jpg")
#creating the box plot with age as the y axis and menopause status as the x-axis
survival_map_mutation_boxplot <- boxplot(formula=as.numeric(map_clinical$survival_time) ~ map_clinical$mutation,
                                         data=clinical,
                                         xlab ="Mutation",
                                         ylab = "Survival Time (days)",
                                         main= "MAP TP53 Mutation vs Survival Time",
                                         cex.axis=0.5)
dev.off()

#GATA3
#line below will save the boxplot to the outputs folder
jpeg("survival_gata_mutation_boxplot.jpg")
#creating the box plot with mutation status as the y axis and survival time as the x-axis
survival_gata_mutation_boxplot <- boxplot(formula=as.numeric(gata_clinical$survival_time) ~ gata_clinical$mutation,
                                          data=clinical,
                                          xlab ="Mutation",
                                          ylab = "Survival Time (days)",
                                          main= "GATA TP53 Mutation vs Survival Time",
                                          cex.axis=0.5)
dev.off()

#CDH1
#line below will save the boxplot to the outputs folder
jpeg("survival_cdh_mutation_boxplot .jpg")
#creating the box plot with age as the y axis and menopause status as the x-axis
survival_cdh_mutation_boxplot <- boxplot(formula=as.numeric(cdh_clinical$survival_time) ~ cdh_clinical$mutation,
                                         data=clinical,
                                         xlab ="Mutation",
                                         ylab = "Survival Time (days)",
                                         main= "CDH TP53 Mutation vs Survival Time",
                                         cex.axis=0.5)
dev.off()

#PIK3CA
#line below will save the boxplot to the outputs folder
jpeg("survival_pik_mutation_boxplot.jpg")
#creating the box plot with age as the y axis and menopause status as the x-axis
survival_pik_mutation_boxplot <- boxplot(formula=as.numeric(pik_clinical$survival_time) ~ pik_clinical$mutation,
                                         data=clinical,
                                         xlab ="Mutation",
                                         ylab = "Survival Time (days)",
                                         main= "PIK3CA TP53 Mutation vs Survival Time",
                                         cex.axis=0.5)
dev.off()



