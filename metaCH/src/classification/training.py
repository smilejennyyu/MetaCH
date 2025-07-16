import numpy as np
import pandas as pd
import pickle
import joblib
import random
import os
from metaCH.src.utils import load_config
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc,roc_curve
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from flaml import AutoML
from collections import Counter

SEED = 1

# Set Python, NumPy, and OS seeds
random.seed(SEED)
np.random.seed(SEED)
os.environ["PYTHONHASHSEED"] = str(SEED)

# Control multi-threading randomness
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"



def seq_classifier_train_save(features_metk_seq, save_path):
    
    config = load_config()
    
    vectors_s, features_s, metadata_s = features_metk_seq

    vectors_s = pd.DataFrame(vectors_s, columns=features_s)
    dataset_s_classif_input = pd.concat([metadata_s, vectors_s], axis=1)


    features_seqclassif = dataset_s_classif_input.columns[
        dataset_s_classif_input.columns.str.contains(r'^(dg_|dbNSFP_|ge_)')
    ].tolist() + ['CANCER_TYPE']

    # ------------------
    # Annotate dataset for classifier
    dataset_s_classif_input['labels_CH_vs_ALL'] = ['Tumor' if i == 'Tumor' else "Tumor" if i=='CH putative driver' else 'Blood' for i in dataset_s_classif_input.cbp_driver_annotation]
    dataset_s_classif_input['labels_CHPD_vs_ALL'] = ['Tumor' if i == 'Tumor' else "Tumor" if i=='CH putative non-driver' else 'Blood' for i in dataset_s_classif_input.cbp_driver_annotation]
    dataset_s_classif_input['labels_CH_vs_Tumor'] = ['Tumor' if i == 'Tumor' else 'Blood' for i in dataset_s_classif_input.cbp_driver_annotation]

    prod_CH = AutoML()
    automl_settings = {
        "time_budget": 1200,  
        "metric": 'log_loss',
        "task": 'classification',
        "n_jobs": 1,  
        "n_splits": 5,
        "split_type": "stratified",
        "seed": SEED,  
        "estimator_list": ["lgbm", "xgboost", "rf"]  
    }
    prod_CH.fit(dataset_s_classif_input[features_seqclassif], dataset_s_classif_input.labels_CH_vs_ALL, **automl_settings)
    pickle.dump([prod_CH, features_seqclassif], open(save_path.format(n=2), 'wb'))


    prod_CHPD = AutoML()
    automl_settings = {
        "time_budget": 1200,  # in seconds
        "metric": 'log_loss',
        "task": 'classification',
        'n_jobs':1,
        "n_splits" :5,
        "split_type": "stratified",
    #     "class_weight": "balanced",
        "seed": SEED, 
        "estimator_list": ["lgbm", "xgboost", "rf"] 
    }
    prod_CHPD.fit(dataset_s_classif_input[features_seqclassif], dataset_s_classif_input.labels_CHPD_vs_ALL, **automl_settings)
    pickle.dump([prod_CHPD, features_seqclassif], open(save_path.format(n=1), 'wb'))


    
    return



def cfDNA_classifier_train_save(metk_f_list, save_path):
    config = load_config()
    
    vectors, features_, metadata = metk_f_list

    vectors = pd.DataFrame(vectors, columns=features_)
    dataset = pd.concat([metadata, vectors], axis=1)

    
    print(Counter(dataset.CANCER_TYPE))
    # ------- generate and merge patient level features and add it to the gene level embeddings
#     sample_id = 'Sample_ID'
    sample_id = 'Patient_ID'
    

    dataset_patient_level = dataset.groupby([sample_id]).mean().reset_index()[[sample_id] + list(dataset.columns[dataset.columns.str.match('dg_')]) ]
    dataset_patient_level.columns = ["pt_{}".format(i) for i in dataset_patient_level.columns]
    master_cfdna_dataset = pd.merge(dataset, dataset_patient_level, left_on=sample_id, right_on='pt_{}'.format(sample_id), how='left')

    dataset_patient_level_ge = dataset.groupby([sample_id]).mean().reset_index()[ [sample_id] + list(dataset.columns[dataset.columns.str.match('ge_')]) ]
    dataset_patient_level_ge.columns = ["pt_{}".format(i) for i in dataset_patient_level_ge.columns]
    master_cfdna_dataset = pd.merge(master_cfdna_dataset, dataset_patient_level_ge, left_on=sample_id, right_on='pt_{}'.format(sample_id), how='left')


    single_features = master_cfdna_dataset.columns[
        master_cfdna_dataset.columns.str.contains(r'^(dg_|dbNSFP_|ge_|pt_dg|pt_ge_)') 
    ].tolist() + ['CANCER_TYPE','VAF']        
    
    automl_settings = {
        "time_budget": 1200,  # in seconds
        "metric": 'log_loss',
        "task": 'classification',
        "n_jobs": 1,  
        "n_splits": 5,
        "split_type": "stratified",
        "seed": SEED,  
        "estimator_list": ["lgbm", "xgboost", "rf"]  
    }

    classifier = AutoML()
    classifier.fit(master_cfdna_dataset[single_features], master_cfdna_dataset.label, **automl_settings)

    pickle.dump([classifier, master_cfdna_dataset, single_features], open(
        save_path, 'wb' 
    ))

    
    return



def cfDNA_classifier_train(metk_f_list):
    config = load_config()
#     model_info = config['models'][version]
    

    vectors, features_, metadata = metk_f_list

    vectors = pd.DataFrame(vectors, columns=features_)
    dataset = pd.concat([metadata, vectors], axis=1)

    
    print(Counter(dataset.CANCER_TYPE))
    # ------- generate and merge patient level features and add it to the gene level embeddings
#     sample_id = 'Sample_ID'
    sample_id = 'Patient_ID'
    

    dataset_patient_level = dataset.groupby([sample_id]).mean().reset_index()[[sample_id] + list(dataset.columns[dataset.columns.str.match('dg_')]) ]
    dataset_patient_level.columns = ["pt_{}".format(i) for i in dataset_patient_level.columns]
    master_cfdna_dataset = pd.merge(dataset, dataset_patient_level, left_on=sample_id, right_on='pt_{}'.format(sample_id), how='left')

    dataset_patient_level_ge = dataset.groupby([sample_id]).mean().reset_index()[ [sample_id] + list(dataset.columns[dataset.columns.str.match('ge_')]) ]
    dataset_patient_level_ge.columns = ["pt_{}".format(i) for i in dataset_patient_level_ge.columns]
    master_cfdna_dataset = pd.merge(master_cfdna_dataset, dataset_patient_level_ge, left_on=sample_id, right_on='pt_{}'.format(sample_id), how='left')

    print(master_cfdna_dataset.label)
    # ------- train cfDNA classifier

    single_features = master_cfdna_dataset.columns[
        master_cfdna_dataset.columns.str.contains(r'^(dg_|dbNSFP_|ge_|pt_dg|pt_ge_)') 
    ].tolist() + ['CANCER_TYPE','VAF']        
    
    n_splits = 5
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    master_cfdna_dataset['CH_cfdna'] = np.nan

    fold_roc_aucs = []

#     X,y=master_cfdna_dataset[single_features], master_cfdna_dataset.label
    X = master_cfdna_dataset[single_features]
    label_map = {'Tumor': 0, 'Blood': 1}  # <- blood is positive class, mapped to 1
    y = master_cfdna_dataset.label.map(label_map)

    categorical_features = X.select_dtypes(include=['object']).columns.tolist()

    preprocessor = ColumnTransformer(
        transformers=[
            ('cat', OneHotEncoder(handle_unknown='ignore'), categorical_features)
        ],
        remainder='passthrough'  # Leave other columns unchanged
    )

    automl_settings = {
        "time_budget": 100,  # in seconds
        "metric": 'roc_auc',
        "task": 'classification',        
        "n_jobs": 1,  
        "n_splits": 5,
        "split_type": "stratified",
        "seed": SEED
    }


    # Step 4: Use StratifiedKFold to evaluate the best estimator with sklearn's pipeline
    n_splits = 5
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    # Initialize an empty list to store ROC AUC scores for each fold
    fold_roc_aucs = []
    fold_pr_aucs = []


    for fold, (train_index, val_index) in enumerate(skf.split(X, y)):
        print(f"Fold {fold + 1}/{n_splits}")

        # Split the data into training and validation sets for this fold
        X_train, X_val = X.iloc[train_index], X.iloc[val_index]
        y_train, y_val = y.iloc[train_index], y.iloc[val_index]

        # Fit the pipeline on the training data
        classifier_cv = AutoML()
        classifier_cv.fit(X_train, y_train, **automl_settings)    
        best_estimator = classifier_cv.model.estimator

        # Create a pipeline with the preprocessor and the best estimator
        pipeline = Pipeline(steps=[('preprocessor', preprocessor), ('model', best_estimator)])    
        pipeline.fit(X_train, y_train)

        # Get validation fold predictions (probabilities)
        val_pred_proba = pipeline.predict_proba(X_val)

        val_pred_proba_positive = val_pred_proba[:, 1]
    
        # Calculate ROC AUC score
        fpr, tpr, _ = roc_curve(y_val, val_pred_proba_positive, pos_label=1)

        master_cfdna_dataset.loc[val_index, 'CH_cfdna'] = val_pred_proba_positive
        master_cfdna_dataset.loc[val_index, 'fold'] = fold

        # Calculate ROC AUC score for this fold with 'Blood' as the positive class
        fold_roc_auc = roc_auc_score(y_val, val_pred_proba_positive)
        fold_roc_aucs.append(fold_roc_auc)
        print(f"Fold {fold + 1} ROC AUC: {fold_roc_auc:.4f}")

        # Calculate PR AUC for this fold with 'Blood' as the positive class
        precision, recall, _ = precision_recall_curve(y_val, val_pred_proba_positive, pos_label=1)
        fold_pr_auc = auc(recall, precision)
        fold_pr_aucs.append(fold_pr_auc)
        print(f"Fold {fold + 1} PR AUC: {fold_pr_auc:.4f}")


    #     plt.plot(fpr, tpr, alpha=0.3, label=f'Fold AUC = {roc_auc:.2f}')

    # Calculate the mean and standard deviation for ROC AUC and PR AUC across folds
    mean_cv_roc_auc = np.mean(fold_roc_aucs)
    std_cv_roc_auc = np.std(fold_roc_aucs)
    mean_cv_pr_auc = np.mean(fold_pr_aucs)
    std_cv_pr_auc = np.std(fold_pr_aucs)

    # Print the results
    print(f"Mean ROC AUC Score: {mean_cv_roc_auc:.4f} ± {std_cv_roc_auc:.4f}")
    print(f"Mean PR AUC Score: {mean_cv_pr_auc:.4f} ± {std_cv_pr_auc:.4f}")

    return master_cfdna_dataset


def meta_classifier_train_save(master_cfdna_dataset, metaclass_features, save_path):
    
    # Initialize Logistic Regression model with class balancing
    log_reg = LogisticRegression(class_weight='balanced', max_iter=1000)

    # Prepare empty lists to store fold results
    roc_auc_scores = []
    fold_pr_aucs=[]
    fold_f1=[]
    fold_best_precision=[]
    fold_best_recall=[]
    fold_best_thrsh=[]

    X = master_cfdna_dataset[metaclass_features]  
    y = master_cfdna_dataset['label'] 

    tprs = []
    recalls=[]
    mean_fpr = np.linspace(0, 1, 100)
    mean_recall = np.linspace(0, 1, 100)


    for fold in range(5):
        X_train = X[master_cfdna_dataset['fold'] != fold]
        y_train = y[master_cfdna_dataset['fold'] != fold]
        X_val = X[master_cfdna_dataset['fold'] == fold]
        y_val = y[master_cfdna_dataset['fold'] == fold]
        
        log_reg.fit(X_train, y_train)

        y_pred_proba = log_reg.predict_proba(X_val)[:, 0]  # Probability of the positive class

        # Calculate ROC AUC score
        fpr, tpr, _ = roc_curve(y_val, y_pred_proba, pos_label='Blood')
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0    

        roc_auc = roc_auc_score(y_val, 1-y_pred_proba)
        roc_auc_scores.append(roc_auc)

        precision, recall, thresholds = precision_recall_curve(y_val, y_pred_proba,pos_label='Blood')
        recalls.append(np.interp(mean_recall, recall[::-1], precision[::-1]))  # Interpolation

        fold_pr_auc = auc(recall, precision)    
        fold_pr_aucs.append(fold_pr_auc)

        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
        fold_f1.append(f1_scores)
        best_idx = f1_scores.argmax()
        best_threshold = thresholds[best_idx]
        best_precision = precision[best_idx]
        fold_best_precision.append(best_precision)
        best_recall = recall[best_idx]
        fold_best_recall.append(best_recall)
        fold_best_thrsh.append(best_threshold)

    # Calculate the average/std ROC AUC score across all folds
    mean_cv_roc_auc = np.mean(roc_auc_scores)
    std_cv_roc_auc = np.std(roc_auc_scores)
    mean_cv_pr_auc = np.mean(fold_pr_aucs)
    std_cv_pr_auc = np.std(fold_pr_aucs)

    print(f"Mean ROC AUC Score: {mean_cv_roc_auc:.2f} ± {std_cv_roc_auc:.3f}")
    print(f"Mean PR AUC Score: {mean_cv_pr_auc:.2f} ± {std_cv_pr_auc:.3f}")

    final_model = LogisticRegression(class_weight='balanced', max_iter=1000)

    # Fit the model on the entire dataset
    final_model.fit(X, y)


    pickle.dump([final_model, metaclass_features], open(
        save_path, 'wb' 
    ))
    
    return
    
    
    
if __name__ == "__main__":
    pass
