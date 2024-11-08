import os 
import scanpy as sc
import pandas as pd
import numpy as np
import celltypist
from celltypist import models
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
# from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef        
import matplotlib.pyplot as plt
import seaborn as sns


# os.chdir("") 
main_dir = 'HTAN/analysis/'
data_dir = '00.data/scRNA-seq/raw/'        
data_raw = sc.read(os.path.join(main_dir, data_dir, 'sce.raw.h5ad'), sparse=True, cleanup=False, dtype='float32')
output_dir = os.path.join(main_dir, '01.scRNA')    
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"{output_dir} has been built")
else:
    print(f"{output_dir} already existed")


max_gene = data_raw.shape[1]
number_genes = [300, 500, 800, 1000, 1200, 1500, 2000]  
mcc_scores_dict = {ng: [] for ng in number_genes}  
for number_gene in number_genes:
    adata = data_raw.copy()
    labels = adata.obs['CAFtype']
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    mcc_scores = []
    for train_index, test_index in kf.split(adata.X):
        train_set = adata[train_index]
        test_set = adata[test_index]
        print(train_set.shape[0])
        print(test_set.shape[0])
        sc.pp.normalize_total(train_set, target_sum = 1e4)
        sc.pp.log1p(train_set)
        train_set.X.expm1().sum(axis = 1)[:10]
        L2_train_model = celltypist.train(train_set, labels = 'CAFtype', n_jobs = 30, feature_selection = True,top_genes = number_gene) 
        sc.pp.normalize_total(test_set, target_sum = 1e4)
        sc.pp.log1p(test_set)
        # test_set.X.expm1().sum(axis = 1)[:10]
        predictions = celltypist.annotate(test_set, model = L2_train_model)
        test_set = predictions.to_adata(insert_prob = True)
        y_true = test_set.obs['CAFtype']
        y_pred = test_set.obs['predicted_labels']
        mcc = matthews_corrcoef(y_true, y_pred)
        print(f"Number of genes: {number_gene}, MCC: {mcc:.2f}")
        mcc_scores_dict[number_gene].append(mcc)


mcc_scores_df = pd.DataFrame.from_dict(mcc_scores_dict, orient='index').T
print(mcc_scores_df)
mcc_means = mcc_scores_df.mean()
print(mcc_means)
mcc_scores_df.to_csv(os.path.join(output_dir, "mcc_scores.csv"), index=False)