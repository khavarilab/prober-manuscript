# prober-manuscript

[![DOI](https://zenodo.org/badge/331781807.svg)](https://zenodo.org/badge/latestdoi/331781807)

Analysis code for the manuscript:

**PROBER Identifies Proteins Associated with Sequence-Specific DNA in Living Cells**<br>
Mondal, *et al*. 2022.

## Install conda environment

Use conda or mamba to avoid installing required packages independently.

```
mamba env create -f conda/env.yaml -n prober-analysis
conda activate prober-analysis
```

## Download ENCODE data used in study

```
source analysis/0-download_encode_data.sh
```

## Run analysis scripts

For example:

```
Rscript --vanilla analysis/1-PROBER-MS_contrasts.R
```

## Inspect results

```
ls results/
```

