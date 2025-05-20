# MetaCH
An artificial intelligence-based model for prediction of clonal hematopoiesis variants in cell-free DNA samples. 

## Requirements
We strongly recommend to use the docker container we prepared for running the MetaCH as this container has all the requirements installed including: 
* StarSpace
* SNPEff
* SNPSift
* FlaML

## Usage
### Jupyter Lab
To interactively run the MetaCH classifier, Initialize a jupyter lab inside the container and follow the next steps: 
#### Load the libraries
```python
import pandas as pd
from metach import chip_classifier
```

#### Required fields in the input data
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

#### Format table for CHIP classifier
```python

table = clf.format_table(
    input_file='./datasets/Chabon_2020/Chabon_Nat2020.tsv'
)
```

#### Define parameters for CHIP classifier
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

#### Predict CHIP score 
```
embeddings = clf.get_embeddings(table)
predictions = clf.predict(embeddings)

```
## Licensing
This software includes third-party components that are governed by separate license terms

While this tool is distributed under the MIT License, some included components are licensed under more restrictive terms, including non-commercial licenses.

>**⚠️ IMPORTANT:** Users are responsible for reviewing and complying with the licenses of all third-party components used by MetaCH and <a href='https://github.com/gaarangoa/METk/tree/main?tab=readme-ov-file#licensing'> METk</a>.
