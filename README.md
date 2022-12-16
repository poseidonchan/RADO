## SimpleDoubletsDetection

### Procedures
- Simulate artificial doublets
	
	- Artifical doublets are the **average** gene expression profiles of two **randomly** selected real droplets in the dataset.
	
- Dimension Reduction

  - Use PCA to reduce the obtain the first 50 PCs of the **real** dataset and **then transform** the concatenation of **simulated and real** datasets.

- Logisitic Regression

  - The training set is the combination of simulated and real datasets. The sample size is 2 times the sample size of the real dataset. All samples are in the real datasets are labeled as singlets in the first iteration.
  - The Logistic Regression is implemented using stochasitic optimization with **Adam** as the optimizer (faster on large dataset).

- (Optional) Iterative training

  - In each iteration, samples with the predicted probability larger than 0.9 are considered as doublets in the next iteration training set.
### Usage
Users can simply copy the function from SDD.py and call it:

```python
pred = SDD(adata)
```

Prediction is a numpy array object where 1 represents doublets and 0 represents singlets.

### Performance

To reproduce the performance, run the following script

```bash
python performance_reproduce.py
```


| Random seed | J293t-dm | hm-6k  | hm-12k | cline-ch | HEK-HMEC-MULTI | HMEC-orig-MULTI | HMEC-rep-MULTI | mkidney-ch | nuc-MULTI | pbmc-1A-dm | pbmc-1B-dm | pbmc-1C-dm | pbmc-2ctrl-dm | pbmc-2stim-dm | pdx-MULTI | pbmc-ch | Average AUPRC |
| ----------- | -------- | ------ | ------ | -------- | -------------- | --------------- | -------------- | ---------- | --------- | ---------- | ---------- | ---------- | ------------- | ------------- | --------- | ------- | ------------- |
| 1234        | 0.6215   | 0.6525 | 0.7672 | 0.4861   | 0.4968         | 0.5216          | 0.6136         | 0.5577     | 0.5799    | 0.6309     | 0.5748     | 0.6307     | 0.6760        | 0.6692        | 0.5791    | 0.6530  | 0.6069        |
| 1235        | 0.4633   | 0.6535 | 0.7327 | 0.5062   | 0.5312         | 0.5569          | 0.6200         | 0.5543     | 0.5941    | 0.6266     | 0.5494     | 0.5518     | 0.6891        | 0.6829        | 0.5741    | 0.6751  | 0.5976        |
| 1236        | 0.5647   | 0.6603 | 0.7527 | 0.4899   | 0.4945         | 0.5250          | 0.6047         | 0.5644     | 0.5871    | 0.6237     | 0.5673     | 0.6197     | 0.6788        | 0.6786        | 0.5651    | 0.6536  | 0.6019        |

According to the benchmarking results done by other researchers such as [vaeda](https://doi.org/10.1093/bioinformatics/btac720) and [scDblFinder](https://f1000research.com/articles/10-979), SimpleDoubletsDetection achieves the **best average AUPRC** across 16 benchmarking datasets.

### Datasets
All the 16 datasets was colleted from the benchmarking paper of [Xi and Li](https://doi.org/10.1016/j.cels.2020.11.008). Datasets was transformed into **H5AD** format using [*sceasy*](https://github.com/cellgeni/sceasy). Processing script is convertH5AD.R.