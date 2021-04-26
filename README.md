# *Stylophora pistillata* single cell atlas

This repository contains the code for the analyses in the paper 

[A stony coral cell atlas illuminates the molecular and cellular basis of coral symbiosis, calcification, and immunity](https://doi.org/10.1016/j.cell.2021.04.005)    
Levy et al. Cell 2021

## Getting around

The files are structured in the following way.

The root directory contains Rmarkdown files with [metacell](https://tanaylab.github.io/metacell/) analysis of the six scRNA-seq datasets included in the paper:  
the three Stylophora life stages (coral, polyp, larva) generated in this study, as well as previously published *Nematostella vectensis* (Sebe-Pedros at al. Cell 2018), *Xenia sp.* (Hu et al. Nature 2020) and *Hydra vulgaris* (Siebert et al. Science 2019) datasets.  
There are two more Rmarkdown files with the cross-stages and cross-species comparative analysis.
```
├── 01_Metacell_coral.Rmd
├── 02_Metacell_polyp.Rmd 
├── 03_Metacell_larva.Rmd
├── 04_Cross_stages_comparison.Rmd
├── 05_Metacell_Nematostella.Rmd
├── 06_Metacell_Xenia.Rmd
├── 07_Metacell_Hydra.Rmd
├── 08_Cross_species_comparison.Rmd 
```

The output from each Rmarkdown file is saved in a separate directory. The output from metacell analyses is in the respective `clustering_*` directories, each of which contains a single-cell database (`scdb`) folder, and optionally other files ndependent of metacell clustering.(e.g. blacklisted genes list).  
```
├── clustering_coral
│   ├── scdb
│   ├── coral_bl_extended.txt
│   └── symbio_signal.rds
├── clustering_hydra
│   └── scdb
├── clustering_larva
│   └── scdb
├── clustering_nematostella
│   └── scdb
├── clustering_polyp
│   └── scdb
├── clustering_xenia
│   ├── scdb
│   └── xenia_bl.txt
``` 

Cross-species and cross-stages output is saved in the directories of the same name. Each directory contains subdirectories with specific analyses.  
```
├── cross_species
│   ├── broad_cell_type
│   └── cell_type
├── cross_stages
│   ├── genelists
│   ├── go
│   └── tree
```

Annotation files used in the analyses are organised in separate per-species directories (Stylophora annotations are in the `annotation` folder, others as indicated).
```
├── annotation
├── annotation_hydra
├── annotation_nematostella
├── annotation_xenia
```

Additional analysis and plotting functions that work on the output of metacell clustering are organized in the dedicated directory.

```
├── metacell_downstream_functions
│   ├── Cross_species_functions.R
│   ├── Downstream_functions.R
│   ├── Export_functions.R
│   ├── Modified_functions.R
│   └── Tree_functions.R
```

Finally, stand-alone functions and scripts that are submitted to HPC are in the `scripts` directory.

