import yaml
import sys
import numpy as np
import pandas as pd
from sklearn import metrics

def load_config(config_path='config/config.yaml'):
    
    config_path = __file__.replace('/src/utils.py', '/config/config.yaml')

    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config



def pr_curve(labels, scores, pos_label, name):
    precision, recall, _ = metrics.precision_recall_curve(labels, scores, pos_label=pos_label)
    auc_val = metrics.auc(recall, precision)
    display = metrics.PrecisionRecallDisplay(precision, recall, estimator_name=f"{name} (auPR: {auc_val:.2f})")
    return display, auc_val

def roc_curve(labels, scores, pos_label, name):
    fpr, tpr, _ = metrics.roc_curve(labels, scores, pos_label=pos_label)
    auc_val = metrics.auc(fpr, tpr)
    display = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=auc_val, estimator_name=name)
    return display, auc_val

def add_patient_level_features(df, sample_col, pattern):
    selected = df.columns[df.columns.str.match(pattern)]
    agg = df.groupby(sample_col)[selected].mean().reset_index()
    agg.columns = [f"pt_{c}" for c in agg.columns]
    return agg


def seq_alter_to_endpos_vartype(variant_id):
    parts = variant_id.split('>')
    ref_allele = parts[0]
    tumor_seq_allele = parts[1]
    if (len(ref_allele) > len(tumor_seq_allele)) and (ref_allele[0]==tumor_seq_allele[0]):
        var_type='DEL';
        ref_allele = ref_allele[1:]  
        tumor_seq_allele = "-"
        end_position = len(ref_allele)
        
    elif len(ref_allele) < len(tumor_seq_allele) and (ref_allele[0]==tumor_seq_allele[0]):
        var_type='INS';
        ref_allele = "-"
        tumor_seq_allele = tumor_seq_allele[1:]
        end_position = 0
    else:
        if len(ref_allele) ==1:
            var_type='SNP';
            end_position=0
        else:
            var_type='ONV';
            end_position = len(ref_allele) - 1

    return pd.Series([var_type, end_position])