######### This version is not for generalized use. The generalized codes can be requested 
#
__author__="Guisong Wang"
__copyright__ ="Copyright (C) 2022 Guisong Wang"
__license__ ="Public Domain"
__version__ = "1.0"


import gc
gc.collect()
import sys
print("User Current Version:-", sys.version)

import os
import re
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
from sklearn.model_selection import GridSearchCV


sns.set(style="white")
sns.set(style="whitegrid", color_codes=True)
pd.set_option('display.max_columns', None)
import openpyxl


#y_scores needs to be a list
#Of the variables in the dataframe
def auc_plot(data, y_true, y_scores, plt_title, figsize=(6,6), save_plot=False, leg_size=None):
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
        plt.title(plt_title)
        ax.legend(loc='lower right', fontsize=leg_size)
    if save_plot:
        plt.savefig(save_plot, dpi=500, bbox_inches='tight')
    plt.show()
    fin_df = pd.concat(fin_dat, axis=0)
    return fin_df


############### data preparation######################

dir0="H:/tmp/ML";
os.chdir(dir0);

##### find the common genes in in-hourse cohort and CPTAC cohort
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

############  get the common gene IDs
geneIDs=list(set(geneIDs_CPTAC).intersection(set(geneIDs_inside)))
geneIDs=["cohort", "Luminal-TN"] + geneIDs;

#### get the genes with common geneIDs
dff0=dff[dff["geneID"].isin(geneIDs)]
dff0=dff0.iloc[:, 2:]

dff_t=dff0.T;
dff_t=dff_t[dff_t["Luminal-TN"].isin(["Luminal", "TN"])];
dff_t['subtype'] = 1*(dff_t['Luminal-TN'] == "TN")

exp_training = dff_t.loc[dff_t["cohort"] == "Training", dff_t.columns[2:]]
exp_testing = dff_t.loc[dff_t["cohort"] == "Testing",dff_t.columns[2:]]
exp_training.shape
exp_testing.shape

################### prepare CPTAC expression data
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
ind_vars = list(dd_merge_CPTAC.columns[dd_merge_CPTAC.columns != "subtype"])

# Dependent variable
y_var = 'subtype'

exp_training[ind_vars] = exp_training[ind_vars].astype(float)
exp_testing[ind_vars] = exp_testing[ind_vars].astype(float)
dd_merge_CPTAC[ind_vars]=dd_merge_CPTAC[ind_vars].astype(float)

assert pd.notnull(exp_training).all().all()
assert pd.notnull(exp_testing).all().all()
assert pd.notnull(dd_merge_CPTAC).all().all()


######################### Build models ###############################

################## select genes by RANDOM FOREST####################################
### tuning for random forest using GridSearchCV
model_rf = RandomForestClassifier(random_state = 1, max_features = 70, n_jobs = 4)  ## 
params = {
    "max_depth": [5, 10, 15, 20],
    "min_samples_leaf": [2, 3, 4, 5, 9, 10],
    "n_estimators": [50, 100, 150, 200, 300]
}

grid_search = GridSearchCV(estimator = model_rf, param_grid = params, cv =5, n_jobs =4, scoring ="roc_auc")  ### stratified k fold

grid_search.fit(exp_training[ind_vars], exp_training[y_var])
grid_search.best_score_
model_rf_best = grid_search.best_estimator_
model_rf_best

### visualization of the tree
# from sklearn.tree import plot_tree
# plt.figure(figsize=(80,40))
# plot_tree(model_rf_best.estimators_[10], feature_names = ind_vars, class_names=['1', "0"],filled=True);

#### check the importance of features
model_rf_best.feature_importances_

imp_rf = pd.DataFrame({
    "Genes": ind_vars,
    "Imp": model_rf_best.feature_importances_
})
imp_rf.sort_values(by="Imp", ascending=False)

#### select features for RF
features_RF = list(imp_rf.loc[imp_rf["Imp"]>0.01, "Genes"])  ### 74 genes were selected

feat_importances = pd.Series(model_rf_best.feature_importances_, index=ind_vars)
feature_plot_rf = feat_importances.nlargest(len(features_RF)).plot(kind='barh')

feature_plot_rf.figure.savefig('importance_genes_RF_imp001.png', format="png", bbox_inches ="tight", dpi=300)

############################### select genes by eXtreme Gradient Booting, XGBoost applies a better regularization technique

model_xgb = XGBClassifier(
    objective= 'binary:logistic',
    n_jobs=4,
    random_state=1
)

params = {
    "max_depth": [5, 10],
    "n_estimators": [50, 100, 150, 200, 300],
    "eta": [0.1, 0.01, 0.05],  ### learning rate
    "gamma": [0.5, 1, 2, 5],  ### min split loss to make a further partition on a leaf node of the tree
    "min_child_weight": [1, 5],
    "subsample": [0.8, 1.0],  ### subsample ratio of the training instances
    "colsample_bytree": [0.8, 1.0] ### the subsampling ratio of columns
}

grid_search = GridSearchCV(
    estimator=model_xgb,
    param_grid=params,
    scoring = 'roc_auc',
    n_jobs = 4,
    cv = 5,
    verbose=True
)

grid_search = GridSearchCV(estimator = model_xgb, param_grid = params, cv =5, n_jobs =4, scoring ="roc_auc")  ### stratified k fold

grid_search.fit(exp_training[ind_vars], exp_training[y_var])
grid_search.best_score_
model_xgb_best = grid_search.best_estimator_
model_xgb_best

model_xgb_best.feature_importances_

imp_xgb = pd.DataFrame({
    "Genes": ind_vars,
    "Imp": model_xgb_best.feature_importances_
})
imp_xgb.sort_values(by="Imp", ascending=False)

features_xgb = list(imp_xgb.loc[imp_xgb["Imp"]>0.01, "Genes"])  ### 74 genes were selected

######## select features for XGB
feature_importances = pd.Series(model_xgb_best.feature_importances_, index=ind_vars)
feature_plot=feature_importances.nlargest(len(features_xgb)).plot(kind='barh')

feature_plot.figure.savefig('importance_genes_XGB_imp001.png', dpi=300, bbox_inches ="tight", format="png")

###### Build the model
final_models = {}
final_models['Logit_saga_l1'] = LogisticRegression(penalty='l1', solver='saga')  ### L1 penalty for feature selection
final_models['Logit_linear_l1'] = LogisticRegression(penalty='l1', solver='liblinear')

# Iterating over each model and fitting on train
for nm, mod in final_models.items():
    mod.fit(exp_training[ind_vars], exp_training[y_var])  ### extract individual independent variables, and dependent variables

features_saga_l1 = [final_models['Logit_saga_l1'].feature_names_in_[index] for index, value in enumerate(final_models['Logit_saga_l1'].coef_[0]) if (value!=0)]
features_linear_l1 = [final_models['Logit_linear_l1'].feature_names_in_[index] for index, value in enumerate(final_models['Logit_linear_l1'].coef_[0]) if (value!=0)]

final_models = {}
final_models['RF'] = model_rf_best
final_models['XGB'] = model_xgb_best
final_models['Logit_saga_l1'] = LogisticRegression(penalty='l1', solver='saga')
final_models['Logit_linear_l1'] = LogisticRegression(penalty='l1', solver='liblinear')


final_models['RF'].fit(exp_training[features_RF], exp_training[y_var])
final_models['XGB'].fit(exp_training[features_xgb], exp_training[y_var])
final_models['Logit_saga_l1'].fit(exp_training[features_saga_l1], exp_training[y_var])
final_models['Logit_linear_l1'].fit(exp_training[features_linear_l1], exp_training[y_var])

datasets=[exp_training, exp_testing, dd_merge_CPTAC];
db_names=["Training", "Testing", "CPTAC"]

############ Apply model to different cohorts and Performance Evaluation by ROC 
for i in range(len(datasets)):
    db=datasets[i];
    db_name = db_names[i];
    for nm, mod in final_models.items():
        # Predicted probs out of sample
        if re.match(nm, "RF"):
            db[nm] =  mod.predict_proba(db[features_RF])[:,1]
        if re.match(nm, "XGB"):
            db[nm] =  mod.predict_proba(db[features_xgb])[:,1]
        if re.match(nm, "Logit_saga_l1"):
            db[nm] =  mod.predict_proba(db[features_saga_l1])[:,1]
        if re.match(nm, "Logit_linear_l1"):
            db[nm] =  mod.predict_proba(db[features_linear_l1])[:,1]
    #recid_train[nm]=mod.predict_proba(recid_train[ind_vars][:, 1]
    filename="AUC_"+db_name+"_by_python.png";
    auc_plot(db, y_var, ["RF", "XGB",   "Logit_saga_l1",  "Logit_linear_l1"], db_name,  save_plot=filename)
