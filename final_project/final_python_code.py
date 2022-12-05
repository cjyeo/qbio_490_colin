# set working director to analysis data
import os
os.chdir('/Users/chloelyc/Desktop/final_project/output')

# import necessary library for analysis
import cptac
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
import seaborn as sns

# download the Brca dataset
cptac.download(dataset="Brca")
brca = cptac.Brca()

# extract the clinical, transcriptomic and protein data
clinical_data = brca.get_clinical()
transcriptomic_data = brca.get_transcriptomics()
protein_data = brca.get_proteomics()
protein_data.columns = protein_data.columns.get_level_values(0)

# Identify the genes (RNA, protein) shared between the two data sets
shared_rna_prot = np.intersect1d(transcriptomic_data.columns, protein_data.columns)

# Create the two data frames
rna_shared = transcriptomic_data.loc[:,shared_rna_prot]
prot_shared = protein_data.loc[:,shared_rna_prot]


# Construct the heat map
fig = plt.figure(figsize=(10,8), dpi=120)

ncomparisons = 5
gene_names = (["TP53","MAP3K1","GATA3", "CDH1", "PIK3CA"])
protein_data.loc[:,gene_names]

corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)

for g1 in gene_names:
    for g2 in gene_names:
        # calculate the spearman correlations between protein and RNA
        corr, pval = stats.spearmanr(rna_shared[g1], prot_shared[g2], nan_policy="omit")
        corr_df.loc[g1,g2] = corr

plot = sns.heatmap(
    corr_df,
    cmap='mako',
)
plot.set_xlabel('Protein', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)
plt.savefig("./heatmap.jpg")
plt.show()
