# `untargeted-workflow.Rmd`

(This is just the `.rmd` version of `LC-MS` untargeted workflow, with many of the notes shifted into this readme.

## Intro

This R-script is suitable for non-targeted profiling experiments. Molecular features are extracted from MS data files, retention time aligned, grouped and filtered to generate a data matrix ready for downstream statistical analysis.

It is recommended to create a project directory (`"./MyProject"` or whatever you'd like to call the project) and copy this script there.

MS data files must be in open file formats such as:
  - `CDF` (tested, but doesn't support MSn for future releases)
  - `mzData` (tested, but larger file size relative to mzXML)
  - `mzML` (not-tested, but should work)
  - `mzXML` (tested and preferable)

Export proprietary formats directly from vendor software or use [ProteoWizard](http://proteowizard.sourceforge.net/index.shtml).


## Package Requirements
List of required packages below (handily presented as an install script)

```
install.packages('RANN')
install.packages('snow')
install.packages('pander')
install.packages('plyr')
install.packages("dplyr")
install.packages("ggpubr")
install.packages('gtools')
install.packages('sm')
install.packages('plotly', dependencies= TRUE)
install.packages('htmltools', dependencies= TRUE)
install.packages('tidyr', dependencies= TRUE)
install.packages('MSeasy', dependencies= TRUE)
install.packages("berryFunctions")
install.packages("magrittr")
install.packages("fda")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("xcms", version = "3.8")
BiocManager::install("CAMERA", version = "3.8")
BiocManager::install("modeest", version = "3.8")
```

## Project Directory Setup

To automatically determine sample classes, group files in folders according to their class e.g. `"./MyProject/MSfiles/Control"`, `"./MyProject/MSfiles/Treated"`.

Outputs intended for downstream analysis are stored in the working directory alongside this script. A subdirectory `/QC` will be created for various quality control outputs which are generally only required for troubleshooting and development. If `"FastMode"` is enabled some of these outputs won't be generated but the final (BP) matrix will be the same.

The first time a user runs this script software dependencies are installed
which will take some minutes.

The main sections of the script are listed below. It is intended that user
input is generally only required in the **Project Settings** and
**Instrument Parameters** sections.

* **Project Settings** - Project specific input
* **Acquisition Parameters** - Parameters relating to the acquisition conditions
* **Prepare Environment** - Installed required software etc
* **Create dataset** - Extract data from MS files
* **RT Alignment** - Aligns the extracted molecular features
* **Reconstruct Spectra** - Construct deconvoluted spectra
* **Export Matrix** - Filter peaks and export basepeak matrix

To clean things up and start afresh you may wish to run:

```
# Clear the console pane. Objects in the environment pane unaffected.
cat("\014")
# Clear the environment pane (the objects). Console pane unaffected.
rm(list = ls(all.names = TRUE))
# Invoke garbage collection to reclaim RAM
gc()
```

## 1. Project Settings
Prior to running this script create a project directory and specify it's
location. If the project directory is named "MStractorDemo" the MS files will
be loaded automatically. Otherwise, place MS files in a sub-directory named
"MSfiles" and group files according to their class into sub-directories of
MSfiles.
e.g.
```
./MyProject
      |--MSfiles
      |     |--Control
      |     |     |--MSfileC1.mzXML
      |     |     |--MSfileC2.mzXML
      |     |--Treatment
      |     |     |--MSfileT1.mzXML
      |     |     |--MSfileT2.mzXML
      |-Demo.R (this script)
```

Specify the path to the project folder. If left as "" a folder named
`MStractorDemo` will be created in the user's home directory.e.g.:
* `C:/users/%UserName%` on Windows
* `/home/$HOME` on Linux
