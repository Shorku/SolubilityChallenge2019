## About
This repository provides a model for 
"[Neural Mulliken Analysis: Molecular Graphs from Density Matrices for QSPR on Raw Quantum-Chemical Data](https://doi.org/10.1021/acs.jctc.5c00425)"
paper (also available as a 
[preprint](https://doi.org/10.26434/chemrxiv-2024-k2k3l-v3) which doesn't cover 
the model provided here). 

A notebook with minimal inference example is also available (see 
[rhnet2.ipynb](https://github.com/Shorku/SolubilityChallenge2019/blob/main/rhnet2.ipynb)).
For general RhNet2 model implementation see the 
[corresponding repository](https://github.com/Shorku/rhnet2).

## Contents
- [Repo structure](#repo-structure)
- [Model](#model)
- [Quantum Chemical Data](#quantum-chemical-data)

## Repo structure
```plaintext
+-- rhnet2.ipynb               # Data preprocessing, model loading and inference example
+-- data_utils.py              # Data preprocessing utilities
+-- RhNet2SC2_results.csv      # Predictions for Solubility Challenge (2019) test sets 
+-- model                      # RhNet2SC2 model
+-- data_example               # Examples of DFT calculations of the test sets compounds
|   +-- set2                   # Solubility Challenge (2019) SET2 compounds
|   |   +-- clofazimine.zip    # ORCA calculation for Clofazimine
```

## Model
The model provided here was trained and tested using test sets from
the Second Solubility Challenge (SC-2, 
[link](https://doi.org/10.1021/acs.jcim.9b00345),
[link](https://doi.org/10.1021/acs.jcim.0c00701)). 
Evaluation results are available in 
[RhNet2SC2_results.csv](https://github.com/Shorku/SolubilityChallenge2019/blob/main/RhNet2SC2_results.csv).
For model and fitting details see the corresponding
[paper](https://doi.org/10.1021/acs.jctc.5c00425)

## Quantum Chemical Data

🚧 🛠️ 🚧