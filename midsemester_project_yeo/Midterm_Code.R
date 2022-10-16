# Author: Colin Yeo
# Title: Midterm Project

# Start by downloading/lodading libraries 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install(version = "3.15")

if (!require("TCGAbiolinks", quietly = TRUE)) 
  BiocManager::install("TCGAbiolinks")

if (!require("SummarizedExperiment", quietly = TRUE)) 
  BiocManager::install("SummarizedExperiment")

if (!require("maftools", quietly = TRUE)) 
  BiocManager::install("maftools")
  
if (!require("ggplot2", quietly = TRUE)) 
  install.packages("ggplot2")
  
if (!require("survival", quietly = TRUE)) 
  install.packages("survival")
  
if (!require("survminer", quietly = TRUE)) 
  install.packages("survminer")
  
library("BiocManager")
library("TCGAbiolinks")
library("SummarizedExperiment")
library("maftools")
library("ggplot2")
library("survival")
library("survminer")
  
# Create outputs folder
dir.create("/Users/colinyeo/desktop/School/qbio_490_colin/midsemester_project_yeo/outputs")  
  
# Query all the relevant data
knitr::opts_knit$set(root.dir = normalizePath("/Users/colinyeo/desktop/School/qbio_490_colin/midsemester_project_yeo/outputs")) 

# Get clinical data
clinical_query <- GDCquery(project="TCGA-BRCA", data.category="Clinical", file.type="xml")
GDCdownload(clinical_query)

clinical.drug <- GDCprepare_clinic(clinical_query, clinical.info="drug")
clinical <- GDCprepare_clinic(clinical_query, clinical.info="patient")

# Get MAF data
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

maf <- GDCprepare(maf_query)
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
maf_object <- read.maf(maf = maf, clinicalData = clinical, isTCGA = TRUE)
  
# 1. Kaplan-Meier Plot for type of procedure done vs. survival
# Need to find all the unique procedures done and remove any N/A's, spaces, or "other"
unique_procedures <- c("Lumpectomy", "Simple Mastectomy", "Modified Radical Mastectomy")

# Create mask for the three valid procedures & create a new dataframe with just these three procedures
procedure_mask <- ifelse(clinical$breast_carcinoma_surgical_procedure_name %in% unique_procedures, T, F)
clinical_kp <- clinical[procedure_mask, ]

# Now create survival time/death event columns for the plot
clinical_kp$survival_time <- ifelse(is.na(clinical_kp$days_to_death), clinical_kp$days_to_last_followup, clinical_kp$days_to_death)
clinical_kp$death_event <- ifelse(clinical_kp$vital_status == "Dead", TRUE, FALSE)

# Finally, create the KM plot... only saves properly if full file path is listed for some reason
jpeg("/Users/colinyeo/desktop/School/qbio_490_colin/midsemester_project_yeo/outputs/KMplot.jpg")
surv_object_age <- Surv(time = clinical_kp$survival_time, event = clinical_kp$death_event)
age_fit <- surv_fit(surv_object_age ~ clinical_kp$breast_carcinoma_surgical_procedure_name, data = clinical_kp)
survplot_age = ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right", title = "Survival plot based on surgical procedure type", legend.labs = unique_procedures)
KM_plot_age = survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16),legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_age
dev.off()

# 2. Bar Graph comparing the average survival of individuals with each treatment
jpeg("/Users/colinyeo/desktop/School/qbio_490_colin/midsemester_project_yeo/outputs/bargraph.jpg")

# Calculate the average survival times for each of the three procedures
lumpectomy_mask <- ifelse(clinical_kp$breast_carcinoma_surgical_procedure_name == "Lumpectomy", T, F)
lumpectomy_total <- clinical_kp$survival_time[lumpectomy_mask]
lumpectomy_average <- sum(lumpectomy_total) / length(lumpectomy_total)
radical_mask <- ifelse(clinical_kp$breast_carcinoma_surgical_procedure_name == "Modified Radical Mastectomy", T, F)
radical_total <- clinical_kp$survival_time[radical_mask]
radical_average <- sum(radical_total) / length(radical_total)
simple_mask <- ifelse(clinical_kp$breast_carcinoma_surgical_procedure_name == "Simple Mastectomy", T, F)
simple_total <- clinical_kp$survival_time[simple_mask]
simple_average <- sum(simple_total) / length(simple_total)

# Put each of these survival times in one vector and make a new data frame
avg_surv_times <- c(lumpectomy_average, simple_average, radical_average)
bar_df <- data.frame(unique_procedures, avg_surv_times)

# Graph the Bar Graph
bargraph = ggplot(data = bar_df, aes(x = unique_procedures, y = avg_surv_times)) + geom_bar(stat="identity")
bargraph = bargraph + labs(title = "Average Survival Time for each Type of Procedure", x = "Type of Procedure", y = "Average Survival Time (days)")
bargraph
dev.off()

# 3. Bar graph visualizing types of treatment with their drug counterparts
# These are the top four drugs used to treat breast cancer according to the clinical.drug data
top_drugs <- c("Cytoxan", "Tamoxifen", "Cyclophosphamide", "Arimidex")

# Rename the col names of clinical.drug, then merge them
colnames(clinical.drug)[1] <- "Tumor_Sample_Barcode"
drug_clinical_total <- merge(clinical_kp, clinical.drug, by="Tumor_Sample_Barcode")

# Mask out rows that don't have one of the top drugs in them then overwrite merged df
top_drug_mask <- ifelse(drug_clinical_total$drug_name %in% top_drugs, T, F)
drug_clinical_total <- drug_clinical_total[top_drug_mask, ]

# Create a new summary data frame to use in the barplot
type_of_procedure <- c("Lumpectomy", "Lumpectomy", "Lumpectomy", "Lumpectomy", "Modified Radical Mastectomy", "Modified Radical Mastectomy", "Modified Radical Mastectomy", "Modified Radical Mastectomy", "Simple Mastectomy", "Simple Mastectomy", "Simple Mastectomy", "Simple Mastectomy")
drug_used <- c("Cytoxan", "Tamoxifen", "Cyclophosphamide", "Arimidex", "Cytoxan", "Tamoxifen", "Cyclophosphamide", "Arimidex", "Cytoxan", "Tamoxifen", "Cyclophosphamide", "Arimidex")

# Have to count the number of appearances each drug makes under each category
# Recreate procedure masks for the new merged df
lumpectomy_mask <- ifelse(drug_clinical_total$breast_carcinoma_surgical_procedure_name == "Lumpectomy", T, F)
radical_mask <- ifelse(drug_clinical_total$breast_carcinoma_surgical_procedure_name == "Modified Radical Mastectomy", T, F)
simple_mask <- ifelse(drug_clinical_total$breast_carcinoma_surgical_procedure_name == "Simple Mastectomy", T, F)

# Now, count total number of patients with each type of procedure + type of drug
num_patients <- c(as.numeric(sum(lumpectomy_mask & drug_clinical_total$drug_name == top_drugs[1])), 
                  as.numeric(sum(lumpectomy_mask & drug_clinical_total$drug_name == top_drugs[2])),
                  as.numeric(sum(lumpectomy_mask & drug_clinical_total$drug_name == top_drugs[3])),
                  as.numeric(sum(lumpectomy_mask & drug_clinical_total$drug_name == top_drugs[4])),
                  as.numeric(sum(radical_mask & drug_clinical_total$drug_name == top_drugs[1])),
                  as.numeric(sum(radical_mask & drug_clinical_total$drug_name == top_drugs[2])),
                  as.numeric(sum(radical_mask & drug_clinical_total$drug_name == top_drugs[3])),
                  as.numeric(sum(radical_mask & drug_clinical_total$drug_name == top_drugs[4])),
                  as.numeric(sum(simple_mask & drug_clinical_total$drug_name == top_drugs[1])),
                  as.numeric(sum(simple_mask & drug_clinical_total$drug_name == top_drugs[2])),
                  as.numeric(sum(simple_mask & drug_clinical_total$drug_name == top_drugs[3])),
                  as.numeric(sum(simple_mask & drug_clinical_total$drug_name == top_drugs[4])))
stacked_df <- data.frame(factor(type_of_procedure), factor(drug_used), factor(num_patients))

# Create the stacked bar graph
jpeg("/Users/colinyeo/desktop/School/qbio_490_colin/midsemester_project_yeo/outputs/stackedbar.jpg")
stackedbar = ggplot(stacked_df, aes(x = type_of_procedure, y = num_patients, fill = drug_used)) + geom_bar(position="stack", stat="identity")
stackedbar = stackedbar + labs(title = "Drug Usage For Each Type of Procedure", x = "Type of Procedure", y = "Number of patients", fill= "Drug Used")
stackedbar
dev.off()

# 4. MAF Plot visualizing potential clustering between mutations and treatment types 
jpeg("/Users/colinyeo/desktop/School/qbio_490_colin/midsemester_project_yeo/outputs/oncoplot.jpg")
oncoplot(maf = maf_object,top = 10, clinicalFeatures = "breast_carcinoma_surgical_procedure_name")
dev.off()
