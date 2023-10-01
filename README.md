## RADO: Robust and Accurate DOublet detection

![Figure1](https://github.com/poseidonchan/RADO/blob/main/figures/Figure1.png)

### Installation
```bash
conda create -n RADO_env python=3.7
# optional: conda env create -f RADO_env.yml
conda activate RADO_env
pip install scRADO==1.2
```

### Usage
```python
from RADO import RADO
adata = RADO(adata)
```

It will return an adata with predicted doublet score and doublet for each droplet in the dataset. The prediction can be found in adata.obs['RADO_doublet_score'] and adata.obs['RADO_doublet_call']. For doublet calling, 0 represents singlet and 1 represents doublet.

### Performance
**16 RNA-seq datasets**

| Method(AUPRC)  | cline-ch | HEK-HMEC-MULTI | hm-12k | hm-6k  | HMEC-orig-MULTI | HMEC-rep-MULTI | J293t-dm | mkidney-ch | nuc-MULTI | pbmc-1A-dm | pbmc-1B-dm | pbmc-1C-dm | pbmc-2ctrl-dm | pbmc-2stim-dm | pbmc-ch | pdx-MULTI | Average |
| -------------- | -------- | -------------- | ------ | ------ | --------------- | -------------- | -------- | ---------- | --------- | ---------- | ---------- | ---------- | ------------- | ------------- | ------- | --------- | ------- |
| DoubletsFinder | 0.3825   | 0.4765         | 0.9890 | 0.9908 | 0.3958          | 0.6035         | 0.2303   | 0.4605     | 0.4428    | 0.4869     | 0.2289     | 0.5368     | 0.6014        | 0.6446        | 0.6066  | 0.3937    | 0.5294  |
| RADO           | 0.4107   | 0.4974         | 0.9783 | 0.9978 | 0.4634          | 0.6104         | 0.2398   | 0.6236     | 0.4638    | 0.4961     | 0.3890     | 0.5684     | 0.6913        | 0.7021        | 0.6496  | 0.4493    | 0.5769  |
| scDblFinder    | 0.4242   | 0.5023         | 0.9533 | 0.9703 | 0.4954          | 0.6045         | 0.1557   | 0.5989     | 0.4480    | 0.5161     | 0.4119     | 0.5878     | 0.6550        | 0.6538        | 0.6535  | 0.3981    | 0.5643  |
| Scrublet       | 0.3789   | 0.4614         | 0.9108 | 0.9691 | 0.3981          | 0.4906         | 0.2553   | 0.5482     | 0.3580    | 0.2449     | 0.2023     | 0.3079     | 0.5635        | 0.5441        | 0.5260  | 0.2528    | 0.4632  |
| solo           |          |                |        |        |                 |                |          |            |           |            |            |            |               |               |         |           |         |
| vaeda          | 0.4067   | 0.4989         | 0.9483 | 0.9683 | 0.4829          | 0.6035         | 0.0978   | 0.5810     | 0.4398    | 0.0759     | 0.3799     | 0.5276     | 0.6733        | 0.6548        | 0.6072  | 0.4178    | 0.5227  |

### Datasets
The 16 RNA-seq datasets was colleted from the benchmarking paper of [Xi and Li](https://doi.org/10.1016/j.cels.2020.11.008). Datasets was transformed into **H5AD** format using [*sceasy*](https://github.com/cellgeni/sceasy). Processing script is convertH5AD.R.