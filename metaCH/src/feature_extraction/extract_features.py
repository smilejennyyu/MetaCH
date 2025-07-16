import os
import pandas as pd
import pickle
from tempfile import TemporaryDirectory
from metaCH.src.utils import load_config
from metk.extractor import FeatureExtractor
from metk.dataset import read
from metk.dataset import Dataset

def extract_features(dataset_table,version):
    
    config = load_config()
    
    model_info = config['models'][version]['metk']
    model_name=model_info.split('/')[-1]
    model_path='/'.join(model_info.split('/')[0:-1])
    
    if "<" in model_path:
        raise ValueError(f"Config parameter: metk , the path to the embedding model, has not been set. Please update it in config.yaml.")
        
    print(model_info)   
    
    dataset_table = dataset_table[dataset_table.Variant_Type.isin(['SNV', 'SNP','DEL','INS'])].reset_index(drop=True)
        
    dataset_table['Variant_Type']=dataset_table['Variant_Type'].replace({'SNV':'SNP'})
    
    variant_types_str='INS,DEL,SNP'
    
    dataset_table['table_unique_id_']='uid_'+dataset_table.index.astype(str)
    
    table = Dataset(dataset_table)
    
    ref_g=dataset_table['Reference_Genome'].values[0]
    with TemporaryDirectory() as metk_temp:
        
        metk = FeatureExtractor(
            reference_genome = ref_g,
            db = model_path,
            identifier = 'table_unique_id_',
            variant_types = variant_types_str,
            output_path = metk_temp,
            mutation_model = model_name,
            run_deepgesture=True,
            run_snpeff=True,
            run_dbnsfp=True,
            run_deepgesture_gene=True,
        )

        metk.extract_features(table)

        metk_f_list = pickle.load(
            open(metk_temp+'/mutation_features.pk', 'rb')
        )
        
        restricted_substr_ls = ['Polyphen2', 'CADD', 'VEST', 'REVEL', 'GenoCanyon','ClinPred']

        features = metk_f_list[1]
        
        filtered_features = [f for f in features if not any(sub in f for sub in restricted_substr_ls)]

        metk_f_list[1] = filtered_features
        return metk_f_list

if __name__ == "__main__":
    pass
