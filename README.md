# MetaCH
An artificial intelligence-based model for prediction of clonal hematopoiesis variants in cell-free DNA samples. 

Load the libraries
```python
import pandas as pd
from metach import chip_classifier
```

Required fields in the input data
```python
schema = {
    "Chromosome": 'Chromosome',
    "Start_Position": 'Start_Position',
    "End_Position": 'End_Position',
    "Reference_Allele": 'Reference_Allele',
    "Tumor_Allele_1": 'Tumor_Allele',
    "Tumor_Allele_2": 'Tumor_Allele',
    "Variant_Type": 'Mutation_Type',
    "Gene_Name": 'Hugo_Symbol',
    "Sample_ID": 'Sample_ID',
    'VAF': 'VAF',
    'CANCER_TYPE': 'cancer_type',
}
```

Format table for CHIP classifier
```python

table = clf.format_table(
    input_file='./datasets/Chabon_2020/Chabon_Nat2020.tsv'
)
```

Define parameters for CHIP classifier
```python

clf = chip_classifier(
    schema=schema,
    sep='\t',
    reference_genome='GRCh37',
    variant_types='Small insertion,Small deletion,Single base substitution',
    output_path='../analysis/',
    metk_db='./metk/',
    metk_model='dgv2.cbioportal.128.e500.bin:MIX_128',
    chip_model_dir='../data/'
)

```

Predict CHIP score 
```
embeddings = clf.get_embeddings(table)
predictions = clf.predict(embeddings)

```

for details please see the <code>./demo/</code>