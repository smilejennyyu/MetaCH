import pandas as pd
import pickle
import joblib
from metaCH.src.utils import load_config, add_patient_level_features

def inference(metk_f_list, version, dbNSFP=False):
    config = load_config()
    model_info = config['models'][version]
    #------- load presaved extracted features by starspace
    vectors, features_, metadata = metk_f_list

    # metadata = pd.read_csv('./metk/metadata.txt', sep='\t', low_memory=False)
    vectors = pd.DataFrame(vectors, columns=features_)
    dataset = pd.concat([metadata, vectors], axis=1)

    
    # ------- generate and merge patient level features and add it to the gene level embeddings
    sample_id = 'Sample_ID'


    for patt in ['dg_', 'ge_']:
        pt_features = add_patient_level_features(dataset, sample_id, patt)
        dataset = pd.merge(dataset, pt_features, left_on=sample_id, right_on=f"pt_{sample_id}", how='left')

    # ------- Load and apply Stage 2 classifers
    classifier, _, ctdna_clf_features = pickle.load(open(
        model_info['classifier_ctdna'], 'rb' 
    ))

    [prod_Seq1, ch_seq1_features] = pickle.load(open(model_info['classifier_Seq1'], 'rb'))
    [prod_Seq2, ch_seq2_features] = pickle.load(open(model_info['classifier_Seq2'], 'rb'))

    if not dbNSFP:
        for feature in ctdna_clf_features:
            if "dbNSFP" in feature:
                dataset[feature] = 0

    dataset['CH_cfdna'] =  classifier.predict_proba(dataset[ctdna_clf_features])[:, 0]
    dataset['CH_seq1'] =prod_Seq1.predict_proba(dataset[ch_seq1_features])[:, 0]
    dataset['CH_seq2'] =prod_Seq2.predict_proba(dataset[ch_seq2_features])[:, 0]

    #------- generate and Meta-Classifier scores
    meta_clf, meta_feats = pickle.load(open(model_info['metaClassifier'], "rb"))
    dataset['ch_score_meta'] = meta_clf.predict_proba(dataset[meta_feats])[:, 0]
    
    results = dataset[config['parsing']['inference_required_columns']].copy()
    results['ch_score_meta']=dataset['ch_score_meta']

    return results

if __name__ == "__main__":
    pass

