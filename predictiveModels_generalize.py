import gc
gc.collect()

import os
import pandas as pd
import scipy
import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt 
plt.rc("font", size=14)
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
import seaborn as sns
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.impute import KNNImputer


sns.set(style="white")
sns.set(style="whitegrid", color_codes=True)
pd.set_option('display.max_columns', None)
import openpyxl


#AUC data plot in wide format, y_scores needs to be a list
#Of the variables in the dataframe
def auc_plot(data, y_true, y_scores, figsize=(6,6), save_plot=False, leg_size=None):
    fin_dat = []
    fig, ax = plt.subplots(figsize=figsize)
    #Equality line
    plt.plot([0, 1], [0, 1], color='grey', lw=1, linestyle='-')
    plt.xlim([-0.03, 1.03])
    plt.ylim([-0.03, 1.03])
    tick_num = [i/10 for i in range(11)]
    plt.xticks(tick_num)
    plt.yticks(tick_num)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    ax.set_aspect(aspect='equal')
    for s in y_scores:
        fpr, tpr, thr = roc_curve(data[y_true], data[s])
        sub_dat = pd.DataFrame( zip(fpr, tpr, thr), 
                                columns=['FPR','TPR','Thresh'] )
        sub_dat['Var'] = s
        auc_val = roc_auc_score(data[y_true], data[s])
        sub_dat['AUC'] = auc_val
        fin_dat.append(sub_dat)
        ax.plot(fpr, tpr, lw=2, label=f'{s} (AUC={auc_val:0.2f})')
    if len(y_scores) == 1:
        plt.title(f'{y_scores[0]} (AUC={auc_val:0.2f})')
    else:
        plt.title('ROC Curves & AUC Scores')
        ax.legend(loc='lower right', fontsize=leg_size)
    if save_plot:
        plt.savefig(save_plot, dpi=500, bbox_inches='tight')
    plt.show()
    fin_df = pd.concat(fin_dat, axis=0)
    return fin_df


############### data preparation######################

dir0="H:/tmp/ML";
os.chdir(dir0);
workbook = openpyxl.load_workbook("Source Data.xlsx");
worksheets = workbook.sheetnames

ws_CPTAC=workbook["protein-gene_mapping_CPTAC"]
gene_mapping = pd.DataFrame(ws_CPTAC.values);
d_gene=gene_mapping.iloc[1:gene_mapping.shape[0], :3];
d_gene.columns=gene_mapping.iloc[0, :3]

geneIDs_CPTAC = d_gene["To"];



ws=workbook["expression_164proteins"]

df=pd.DataFrame(ws.values);
dff=df.iloc[1:df.shape[0], :];
dff.columns=df.iloc[0, :]

dff_gene = dff[["geneID", "Gene_Symbol"]]
dff = dff.set_index("Gene_Symbol");

## replace NA for cohort and Luminal-TN in the column geneID
dff.loc["cohort", "geneID"]="cohort";
dff.loc["Luminal-TN", "geneID"]="Luminal-TN";

geneIDs_inside= dff["geneID"];

### get the common IDs
geneIDs=list(set(geneIDs_CPTAC).intersection(set(geneIDs_inside)))

geneIDs=["cohort", "Luminal-TN"] + geneIDs;


#### get the genes with common geneIDs
dff0=dff[dff["geneID"].isin(geneIDs)]
dff0=dff0.iloc[:, 2:]

dff_t=dff0.T;
dff_t=dff_t[dff_t["Luminal-TN"].isin(["Luminal", "TN"])];

dff_t['subtype'] = 1*(dff_t['Luminal-TN'] == "TN")


### get CPTAC data
exp_training = dff_t.loc[dff_t["cohort"] == "Training", dff_t.columns[2:]]
exp_testing = dff_t.loc[dff_t["cohort"] == "Testing",dff_t.columns[2:]]

exp_training.shape
exp_testing.shape

### prepare CPTAC expression data
dd_gene = d_gene[d_gene["To"].isin(geneIDs)]

ws_CPTAC_exp = workbook["Global-Proteome-G3_CPTAC"]
exp_CPTAC0 = pd.DataFrame(ws_CPTAC_exp.values);

column_list=[0] + list(range(12, exp_CPTAC0.shape[1]))
exp_CPTAC = exp_CPTAC0.iloc[1:, column_list]
exp_CPTAC.columns=exp_CPTAC0.iloc[0, column_list]

### join
exp_CPTAC = dd_gene.merge(exp_CPTAC, left_on="From", right_on="accession_number");

CPTAC_samples=list(exp_CPTAC.columns[4:])
exp_CPTAC =exp_CPTAC[["From", "To"] + CPTAC_samples];
exp_CPTAC_gene = dff_gene.merge(exp_CPTAC, left_on = "geneID", right_on="To");

exp_CPTAC_gene=exp_CPTAC_gene[["Gene_Symbol"] + CPTAC_samples];
exp_CPTAC_gene = exp_CPTAC_gene.set_index("Gene_Symbol")

exp_CPTAC_gene[CPTAC_samples]=exp_CPTAC_gene[CPTAC_samples].astype(float);
### impute the data
imputer = KNNImputer(n_neighbors=10);
df00=pd.DataFrame(imputer.fit_transform(exp_CPTAC_gene));
df00.index=exp_CPTAC_gene.index;
df00.columns=exp_CPTAC_gene.columns;
### aggregate by Genes
exp_CPTAC_genes = df00.groupby("Gene_Symbol").median()
exp_CPTAC_genes.shape
### aggregate by duplicates
### first grep TCGA samples and convert the names
namess=exp_CPTAC_genes.columns;
bool_TCGA = namess.str.contains('TCGA') 
exp_CPTAC_genes = exp_CPTAC_genes[namess[bool_TCGA]];
namess=list(exp_CPTAC_genes.columns);
new_namess=["TCGA-" + name0[0:7] for name0 in namess];
exp_CPTAC_genes.columns = new_namess;

exp_CPTAC_t=exp_CPTAC_genes.T;
exp_CPTAC_t.reset_index(inplace=True);
exp_CPTAC_t= exp_CPTAC_t.rename(columns={"index":"samples"});

ws_clinical = workbook['TCGA_clinical_subtypes']

dd_clinical = pd.DataFrame(ws_clinical.values);

dd_clinicals=dd_clinical.iloc[1:, :];
dd_clinicals.columns=dd_clinical.iloc[0, :];
dd_clinicals = dd_clinicals[["samples", "Luminal-TN"]];

dd_clinicals = dd_clinicals[dd_clinicals["Luminal-TN"].isin(["Luminal", "TN"])]

dd_merge_CPTAC=dd_clinicals.merge(exp_CPTAC_t, left_on="samples", right_on="samples");
dd_merge_CPTAC['subtype'] = 1*(dd_merge_CPTAC['Luminal-TN'] == "TN")

dd_merge_CPTAC = dd_merge_CPTAC[list(exp_training.columns)];



#Independant variables
ind_vars = list(dd_merge_CPTAC.columns[1:])

# Dependent variable
y_var = 'subtype'


##### Next step: fit the model

### stuff the base model object in a dictionary, pipe in the same training data, 
### and fit the models. Then I can add in the predicted probabilities from each model into the test dataset.

###############################################
# Training three different models, Logit,
# Random Forest, and XGBoost

final_models = {}
final_models['XGB'] = XGBClassifier(n_estimators=100, max_depth=5)
final_models['RF'] = RandomForestClassifier(n_estimators=1000, max_depth=10, min_samples_split=50)
final_models['Logit'] = LogisticRegression(penalty='none', solver='newton-cg')


# Iterating over each model and fitting on train
for nm, mod in final_models.items():
    mod.fit(exp_training[ind_vars], dd_training[y_var])  ### extract individual independent variables, and dependent variables


# Adding predicted probabilities back into test dataset
for nm, mod in final_models.items():
    # Predicted probs out of sample
    dd_rest[nm] =  mod.predict_proba(dd_rest[ind_vars])[:,1]
    #recid_train[nm]=mod.predict_proba(recid_train[ind_vars][:, 1]
    
pred_prob_cols = list(final_models.keys()) #variable names

# binary_plots.cal_data_wide_group(pred_prob_cols, y_var, 'cohort', dd_all, bins=20, plot=True, save_plot='Cal4_all.png')


######### AUC plot

###metric to see the utility of a particular binary prediction model is the AUC stat. This has one interpretation in terms of the concordance stat, 
### an AUC of 0.7 means if you randomly picked a 0 case and a 1 case, the 1 case would have a higher value 70% of the time. 
#### So AUC is all about how well your prediction discriminates between the two classes.
### it tells how much the model is capable of distingushing between 2 classes
auc_plot(dd_training, y_var, ['Logit'], save_plot='AUC1_all_ASNS.png')
auc_plot_wide_group(dd_rest, y_var, pred_prob_cols, 'cohort', size=4, leg_size='x-small', save_plot='AUC4_all_ASNS.png')

# Multiple columns to show different models
# pred_prob_cols = list(final_models.keys()) #variable names
# binary_plots.auc_plot(dd_all, y_var, pred_prob_cols, save_plot='AUC2_all.png')

# binary_plots.cal_data_wide_group(pred_prob_cols, y_var, 'cohort', dd_all,  bins=4, plot=True, save_plot='Cal4_all.png')