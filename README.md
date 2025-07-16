# MetaCH
### An artificial intelligence-based model for prediction of clonal hematopoiesis variants in cell-free DNA samples. 
---

**MetaCH** is a machine learning framework for classifying the origin of cfDNA variants, distinguishing clonal hematopoiesis (CH) from tumor-derived mutations. This repository provides:

✅ Ready-to-use pre-trained classifiers  
✅ Scripts to generate variant-level predictions from user cfDNA data  
✅ Code to reproduce the results and figures from our publication  
✅ Integration with the [Mutation Enrichment Toolkit (METk)](https://github.com/gaarangoa/METk) for feature extraction for generation of general-purpose gene/variant embeddings

---

## 📦 Repository Structure

```
.                              
├── example_inference.ipynb    # ✅ Main notebook to test inference on user data using pre-trained models
├── example_input.csv          # Sample variant input table used in example_inference.ipynb
├── framework_training.ipynb   # Notebook to train classifiers whithin the multi-stage framework (e.g. from scratch or on new data)
├── generate_paper_figs.ipynb  # Notebook to reproduce publication figures by application of the framework on external validation datasets
├── metaCH                     
│   ├── config                 
│   │   ├── config.yaml        # Default settings (e.g. model paths, required data columns)
│   │   └── __init__.py        
│   ├── __init__.py            
│   └── src                    
│       ├── classification     
│       │   ├── inference.py   
│       │   ├── training.py    # Code to train models used in framework_training.ipynb
│       │   └── __init__.py
│       ├── feature_extraction # Feature extraction wrapper around METk
│       │   ├── extract_features.py # Wraps METk API for use in our pipeline
│       │   └── __init__.py
│       ├── utils.py           
│       └── __init__.py
├── models                     # Pre-trained classifier files
│   ├── cfDNA_classifier*.pk   # cfDNA classifier
│   ├── metaClassifier*.pk     # Meta-classifier
│   ├── seq1_classifier*.pk    # Sequence-based1 classifier (CH-O versus Others)
│   └── seq2_classifier*.pk    # Sequence-based2 classifier (CH-NO versus Others)
├── results                    
│   └── baselines
│       └── ssgan              
│           ├── chabon_ssgan_preds.csv
│           ├── chin_ssgan_preds.csv
│           ├── leal_ssgan_preds.csv
│           └── zhang_ssgan_preds.csv
├── README.md                  
└── setup.py                   
```

---

## 🚀 Quick Start

### 1. Setup Environment with METk Docker

We strongly recommend using the official METk Docker image for compatibility.

📎 Full instructions here: [METk Setup Guide](https://github.com/gaarangoa/METk#setup)

**Quick Example:**

```bash
docker run -it --rm \
    -p 8888:8888 \
    -v /path/to/metk/:/METk/ \
    -v /path/to/data/:/data/ \
    gaarangoa/chip_classifier:version-4.0.0 \
    jupyter notebook --NotebookApp.default_url=/lab/ --ip=0.0.0.0 --port=8888 --allow-root
```

Then inside the container:

```bash
pip install --user git+https://github.com/gaarangoa/METk.git
pip install --user git+https://github.com/gaarangoa/MetaCH.git
```

---

### 2. edit 🛠 `config.yaml` file

- **`path_info`**  
  File system paths to required datasets (e.g., benchmarking, Razavi cfDNA, cBioPortal/MSK).  
  **👉 Must be updated by the user before running any pipeline.**

- **`parsing`**  
  Lists the expected columns in input variant tables for training and inference.  
  **👉 Does not require changes unless user aims retraining/inference using a subset of columns**

- **`models`**  
  Paths to pre-trained model weights (METk embedding bin file and classifiers).  
  **👉 Update if models are stored elsewhere or you are using custom models.**


---

## ✅ Generating Predictions on Your Own cfDNA Data

We provide a worked example for applying the model to your own data in the notebook:

**`example_inference.ipynb`** – run this notebook to:
- Load an example variant table (e.g., `example_input.csv`)
- Extract features via METk
- Apply the framework in the inference mode by applying the pre-trained MetaCH classifiers
- Set the threshold for categorical Blood/Tumor predictions (default is 0.5)
- View and save prediction outputs

**Input:** Variant table in CSV format  
**Output:** Variant-level predictions indicating CH vs. tumor-derived origin

> ⚠️ Make sure you're running this inside the METk Docker container with both METk and MetaCH installed.

---

## ✅ Reproducing Results from the Paper

Use the following notebooks to replicate the analyses:

```bash
jupyter notebook framework_training.ipynb
jupyter notebook generate_paper_figs.ipynb
```

All pre-trained models required for reproduction are available in the `/models` directory except METk pretrained model for feature extraction which can be downloaded from [here](https://github.com/gaarangoa/METk?tab=readme-ov-file#download-metk-embeddings).

---

## 🛠 Available Pre-trained Models

| Model Name            | Description                                                     | 📚 Training Dataset               | 📁 path/link            |
| --------------------- | --------------------------------------------------------------- | --------------------------------- | ---------------------------- |
| `dgv2.cbioportal.128.e500.bin`        | Embedding model for genes/variants using co-occurrence patterns | TCGA + cBioPortal                 |  [METk pretrained model ](https://github.com/gaarangoa/METk?tab=readme-ov-file#download-metk-embeddings)   |
| `cfDNA_classifier.pk` | Primary classifier for distinguishing CH vs. tumor              | Razavi et al. cfDNA dataset       | `models/cfDNA_classifier.pk` |
| `metaClassifier.pk`   | Stacked model combining multiple feature types                  | Same as `cfDNA_classifier.pk`     | `models/metaClassifier.pk`   |
| `seq1_classifier.pk`  | Sequence-based classifier for CH-O vs. rest                     | `msk_impact_2017` + `msk_ch_2020` | `models/seq1_classifier.pk`  |
| `seq2_classifier.pk`  | Sequence-based classifier for CH-NO vs. rest                    | `msk_impact_2017` + `msk_ch_2020` | `models/seq2_classifier.pk`  |


---
| `metkmodel.pk`        | Embedding model for genes/variants using co-occurrence patterns | TCGA + cBioPortal                 |      |

## 📊 METk Utilization

 Feature extraction to generate Gene embeddings, Variant embeddings and Functional prediction scores relies on [METk](https://github.com/gaarangoa/METk).

METk supports both SNVs and INDELs.  
We strongly recommend using the [Docker setup](https://github.com/gaarangoa/METk#setup) for consistent execution.

---

## 📂 Data Availability

- **METk Training data:**
  - TCGA: [https://portal.gdc.cancer.gov](https://portal.gdc.cancer.gov)
  - cBioPortal studies: [https://www.cbioportal.org](https://www.cbioportal.org)
- **cfDNA classifier training dataset:** Available from Razavi et al.
- **sequence-based classifier dataset:** `msk_impact_2017`, `msk_ch_2020`
- **External validation sets:** Chin et al., Chabon et al., Leal et al., Zhang et al.
- **Comparison models:** Fairchild et al., SSGAN

Some datasets may require access requests.

---

## Licensing
This software includes third-party components that are governed by separate license terms

While this tool is distributed under the MIT License, some included components are licensed under more restrictive terms.

>**⚠️ IMPORTANT:** Users are responsible for reviewing and complying with the licenses of all third-party components used by MetaCH and <a href='https://github.com/gaarangoa/METk/tree/main?tab=readme-ov-file#licensing'> METk</a>.


---

## 📚 Citation

If you use MetaCH, please cite:

> **An artificial intelligence-based model for prediction of clonal hematopoiesis variants in cell-free DNA samples**  
> *Arango-Argoty, G., Haghighi, M., Sun, G.J., Choe, E.Y., Markovets, A., Barrett, J.C., Lai, Z. and Jacob, E.*  
> npj Precision Oncology 9.1 (2025): 147., DOI: [Link](https://www.nature.com/articles/s41698-025-00921-w.pdf)

---

## 📬 Contact

Please open a GitHub issue for needed support.

---

*This repository accompanies the publication of MetaCH and provides all necessary tools for variant origin prediction and result reproduction.*



