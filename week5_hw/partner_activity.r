#download and load dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install(version = "3.15")

install.packages("dplyr")
library(dplyr)

library(BiocManager)

if (!require("TCGAbiolinks", quietly = TRUE)) 
  BiocManager::install("TCGAbiolinks") # notice this is different from the typical "install.packages" command. that's because BiocManager is necessary to install and manage packages from the Bioconductor project

library(TCGAbiolinks)
  
library(ggplot2)


#set working directory
getwd()
setwd("/Users/colinyeo/Desktop/School/qbio_490_colin/analysis_data/")
knitr::opts_knit$set(root.dir = normalizePath("../analysis_data")) 

#read in data
clinical_query <- GDCquery(project="TCGA-BRCA", data.category="Clinical", file.type="xml")
clinical <- GDCprepare_clinic(clinical_query, clinical.info="patient")
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")

str(clinical) #view structure of clinical df
str(clinical.drug) #view structure of clinical.drug df
str(clinical.rad) #view structure of clinical.rad df

#view first few lines of each df
head(clinical)
head(clinical.drug)
head(clinical.rad)
clinical$age_at_initial_pathologic_diagnosis
clinical.drug$drug_name

#merge clinical and clinical.drug df's using patient ID as common column; needed to visualize using ggplot
total_clinical <- merge(clinical, clinical.drug, by="bcr_patient_barcode")

unique(total_clinical$drug_name) #there are over 400 drugs, which is too much to plot!
counts <- count(total_clinical,vars=  drug_name) #count the number of times each drug is used
counts
ordered_counts <- counts[order(counts$n),]
ordered_counts #order counts to obtain top 9 drugs administered
top_drugs = c("Cytoxan", "Tamoxifen", "Cyclophosphamide", "Arimidex", "Doxorubicin", "Adriamycin", "Taxoatere", "Paclitaxel", "Taxol", "Anastrozole")
top_frequency = ordered_counts[ordered_counts$vars %in% top_drugs,]
top_frequency

modified_clinical = droplevels(total_clinical[total_clinical$drug_name %in% top_drugs,]) #filters df to include only T9 drugs
head(modified_clinical)

#create plot of age at initial diagnosis vs. name of drug administered 
par(mar= c(8.5,8.5,7,7)) #sets margins of plot to ensure all labels fit
boxplot(age_at_initial_pathologic_diagnosis~drug_name, data= modified_clinical, las= 2, main="Age of Initial Pathologic Diagnosis and Drug Administered", xlab = "", ylab= "Age at Initial Pathologic Diagnosis")
mtext("Drug", side=1, line=4.5) #adds x label that doesn't overlap with drug names

#Survival Analysis for Age at Initial Pathologic Diagnosis
#install dependencies
if (!require("survival", quietly = TRUE)) 
  survival::install("survival")

library(survival)

if (!require("maftools", quietly = TRUE)) 
  install("maftools")

library(maftools)

if (!require("survminer", quietly = TRUE)) 
  install("survminer")

library(survminer)

#Conduct survivorship analysis
clinical$survival_time <- ifelse(is.na(clinical$days_to_death),clinical$days_to_last_followup,clinical$days_to_death)
clinical$survival_time


clinical$death_event <- ifelse(clinical$vital_status == "Dead", TRUE, FALSE)
clinical$death_event

surv_object_age <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)
age_fit <- surv_fit(surv_object_age~clinical$age_category, data = clinical)
survplot_age = ggsurvplot(age_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")
KM_plot_age = survplot_age$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_age

#Survival Analysis for Drug Administered
#Conduct survivorship analysis
modified_clinical$survival_time <- ifelse(is.na(modified_clinical$days_to_death),modified_clinical$days_to_last_followup,modified_clinical$days_to_death)
modified_clinical$survival_time


modified_clinical$death_event <- ifelse(modified_clinical$vital_status == "Dead", TRUE, FALSE)
modified_clinical$death_event

surv_object_age <- Surv(time = modified_clinical$survival_time,
                        event = modified_clinical$death_event)
drug_fit <- surv_fit(surv_object_age ~ modified_clinical$drug_name,
                    data = modified_clinical)
survplot_drug = ggsurvplot(drug_fit,
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")
KM_plot_drug = survplot_drug$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))

KM_plot_drug

#write total_clinical, modified_clinical df to csv 
  write.csv(total_clinical,"/Users/colinyeo/Desktop/School/qbio_490_colin/total_clinical.csv")
  write.csv(modified_clinical, "/Users/colinyeo/Desktop/School/qbio_490_colin/week5_hw/modified_clinical.csv")

