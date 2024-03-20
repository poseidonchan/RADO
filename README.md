## RADO: Robust and Accurate DOublet detection

![scRADO](https://img.shields.io/badge/scRADO-v1.3-blue)![PyPI - Downloads](https://img.shields.io/pypi/dm/scRADO)![GitHub](https://img.shields.io/github/license/poseidonchan/RADO)

![Figure1](https://github.com/poseidonchan/RADO/blob/main/figures/Figure1.png)

### Installation
```bash
conda create -n RADO_env python=3.7
conda activate RADO_env
pip install umap-learn==0.5.3 #(to be compatible with python3.7)
pip install scRADO==1.3

```

### Usage

**For scRNA-seq data**

```python
from RADO import DoubletDetection
# adata (.H5AD file) is commmon data form in single-cell data analysis
adata = DoubletDetection(adata)
# filter out doublet
adata = adata[adata.obs['RADO_doublet_call']==0,]
```

Also see the [tutorial](https://github.com/poseidonchan/RADO/blob/main/tutorial.ipynb) please, for any other questions, raise issues please!

**For scATAC-seq data**

```python
from RADO import DoubletDetection
# Assume the adata.X is the peak matrix
adata = DoubletDetection(adata, atac_data=True)
# filter out doublet
adata = adata[adata.obs['RADO_doublet_call']==0,]
```

It will return an adata with predicted doublet score and doublet for each droplet in the dataset. The prediction can be found in adata.obs['RADO_doublet_score'] and adata.obs['RADO_doublet_call']. For doublet calling, 0 represents singlet and 1 represents doublet.

### Performance
**18 scRNA-seq datasets**

![Figure2](https://github.com/poseidonchan/RADO/blob/main/figures/Figure2.png)

**2 scATAC-seq datasets**

![Figure2](https://github.com/poseidonchan/RADO/blob/main/figures/Figure4.png)

### Datasets

The 16 scRNA-seq datasets were collected from the benchmarking paper of [Xi and Li](https://doi.org/10.1016/j.cels.2020.11.008). Datasets were transformed into **H5AD** format using [*sceasy*](https://github.com/cellgeni/sceasy). Processing script is [convertH5AD.R](https://github.com/poseidonchan/RADO/blob/main/convertH5AD.R).

The 2 DOGMA-seq datasets are from [Xu *et al*](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02698-8). Datasets with original singlet or doublet annotation need to be requested from the original authors.