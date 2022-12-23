## SimpleDoubletsDetection

### Procedures
- Simulate artificial doublets
	
	- Artifical doublets are the **average** gene expression profiles of two **randomly** selected real droplets in the dataset.
- Compute the simulated doublets proportion in the neighborhood (KNN score)
  - the neighborhood size is $\sqrt{2n}$
- Dimension Reduction
  - Use PCA to reduce the dimension, the first 10 PCs are used as the features
- Logisitic Regression

  - The training set is the combination of simulated and partial real datasets. Then the left samples of the real dataset is predicted, default is 5 fold validation. All samples are in the real datasets are labeled as singlets.
  - KNN score is concatenated to the PCA features and used in the training & test procedure.
  - The Logistic Regression is implemented with **L1** penalty using *saga* as the optimizer (faster on large dataset).
### Usage
Users can simply copy the function from SDD.py and call it:

```python
pred = SDD(adata)
```

Prediction is a numpy array object where 1 represents doublets and 0 represents singlets.

### Performance

To reproduce the performance of SimpleDoubletsDetection, run the following script

```bash
python reproduce.py --iter 1 --kfold 5 --seed 1234 --dir ./data/RNAseq/
```
```bash
usage: reproduce.py [-h] [--iter ITERATION] [--kfold KFOLD] [--seed SEED]
                    [--dir DIR]

optional arguments:
  -h, --help        show this help message and exit
  --iter ITERATION  iterations to perform PU-learning (default: 1)
  --kfold KFOLD     iterations to perform PU-learning (default: 5)
  --seed SEED       pin random seed to reproduce results (default: 1234)
  --dir DIR         specify the directory that contains all the datasets
                    (default: None)
```



| Method                  | cline-ch | HEK-HMEC-MULTI | hm-12k | hm-6k  | HMEC-orig-MULTI | HMEC-rep-MULTI | J293t-dm | mkidney-ch | nuc-MULTI | pbmc-1A-dm | pbmc-1B-dm | pbmc-1C-dm | pbmc-2ctrl-dm | pbmc-2stim-dm | pbmc-ch | pdx-MULTI | Average |
| ----------------------- | -------- | -------------- | ------ | ------ | --------------- | -------------- | -------- | ---------- | --------- | ---------- | ---------- | ---------- | ------------- | ------------- | ------- | --------- | ------- |
| scDblFinder             | 0.4242   | 0.5023         | 0.9533 | 0.9703 | 0.4954          | 0.6045         | 0.1557   | 0.5989     | 0.4480    | 0.5161     | 0.4119     | 0.5878     | 0.6550        | 0.6538        | 0.6535  | 0.3981    | 0.5643  |
| vaeda                   | 0.4067   | 0.4989         | 0.9483 | 0.9683 | 0.4829          | 0.6035         | 0.0978   | 0.5810     | 0.4398    | 0.0759     | 0.3799     | 0.5276     | 0.6733        | 0.6548        | 0.6072  | 0.4178    | 0.5227  |
| DoubletsFinder          | 0.3825   | 0.4765         | 0.9890 | 0.9908 | 0.3958          | 0.6035         | 0.2303   | 0.4605     | 0.4428    | 0.4869     | 0.2289     | 0.5368     | 0.6014        | 0.6446        | 0.6066  | 0.3937    | 0.5294  |
| Scrublet                | 0.3789   | 0.4614         | 0.9108 | 0.9691 | 0.3981          | 0.4906         | 0.2553   | 0.5482     | 0.3580    | 0.2449     | 0.2023     | 0.3079     | 0.5635        | 0.5441        | 0.5260  | 0.2528    | 0.4632  |
| SImpleDoubletsDetection | 0.4141   | 0.5274         | 0.9170 | 0.9997 | 0.4588          | 0.6150         | 0.1436   | 0.6303     | 0.4613    | 0.5024     | 0.3754     | 0.5555     | 0.6580        | 0.6509        | 0.6509  | 0.4401    | 0.5625  |

According to the benchmarking results done by other researchers such as [vaeda](https://doi.org/10.1093/bioinformatics/btac720) and [scDblFinder](https://f1000research.com/articles/10-979), SimpleDoubletsDetection achieves the ***comparable*** **average AUPRC** compared to state-of-the-art method  [scDblFinder](https://f1000research.com/articles/10-979) across 16 benchmarking datasets. 

### Datasets
All the 16 datasets was colleted from the benchmarking paper of [Xi and Li](https://doi.org/10.1016/j.cels.2020.11.008). Datasets was transformed into **H5AD** format using [*sceasy*](https://github.com/cellgeni/sceasy). Processing script is convertH5AD.R.