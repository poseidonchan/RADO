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
python reproduce.py --iteration 2 --seed 1234 --dir ./data/RNAseq/
```
```bash
usage: reproduce.py [-h] [--iteration ITERATION] [--seed SEED]
                                [--dir DIR]

optional arguments:
  -h, --help            show this help message and exit
  --iteration ITERATION
                        iteration number of iterative training procedure,
                        usually 2 is enough (default: 2)
  --seed SEED           pin random seed to reproduce results (default: 1234)
  --dir DIR             specify the directory that contains all the datasets
                        (default: None)
```

**Without iterative training:**

| Random seed | J293t-dm | hm-6k  | hm-12k | cline-ch | HEK-HMEC-MULTI | HMEC-orig-MULTI | HMEC-rep-MULTI | mkidney-ch | nuc-MULTI | pbmc-1A-dm | pbmc-1B-dm | pbmc-1C-dm | pbmc-2ctrl-dm | pbmc-2stim-dm | pdx-MULTI | pbmc-ch | Average AUPRC |
| ----------- | -------- | ------ | ------ | -------- | -------------- | --------------- | -------------- | ---------- | --------- | ---------- | ---------- | ---------- | ------------- | ------------- | --------- | ------- | ------------- |
| 1234        | 0.6215   | 0.6525 | 0.7672 | 0.4861   | 0.4968         | 0.5216          | 0.6136         | 0.5577     | 0.5799    | 0.6309     | 0.5748     | 0.6307     | 0.6760        | 0.6692        | 0.5791    | 0.6530  | 0.6069        |
| 1235        | 0.4633   | 0.6535 | 0.7327 | 0.5062   | 0.5312         | 0.5569          | 0.6200         | 0.5543     | 0.5941    | 0.6266     | 0.5494     | 0.5518     | 0.6891        | 0.6829        | 0.5741    | 0.6751  | 0.5976        |
| 1236        | 0.5647   | 0.6603 | 0.7527 | 0.4899   | 0.4945         | 0.5250          | 0.6047         | 0.5644     | 0.5871    | 0.6237     | 0.5673     | 0.6197     | 0.6788        | 0.6786        | 0.5651    | 0.6536  | 0.6019        |
| 1237        | 0.4907   | 0.6671 | 0.7523 | 0.4780   | 0.4998         | 0.5064          | 0.6049         | 0.5611     | 0.5770    | 0.6004     | 0.5825     | 0.6184     | 0.6718        | 0.6792        | 0.5740    | 0.6605  | 0.5953        |
| 1238        | 0.4117   | 0.6510 | 0.7483 | 0.4903   | 0.5610         | 0.5315          | 0.6079         | 0.5648     | 0.6108    | 0.5192     | 0.5821     | 0.6221     | 0.6908        | 0.6756        | 0.5633    | 0.6620  | 0.5933        |

**With iterative training** (5 iteration):

| Random seed | J293t-dm   | hm-6k  | hm-12k | cline-ch   | HEK-HMEC-MULTI | HMEC-orig-MULTI | HMEC-rep-MULTI | mkidney-ch | nuc-MULTI  | pbmc-1A-dm | pbmc-1B-dm | pbmc-1C-dm | pbmc-2ctrl-dm | pbmc-2stim-dm | pdx-MULTI  | pbmc-ch    | Average AUPRC |
| ----------- | ---------- | ------ | ------ | ---------- | -------------- | --------------- | -------------- | ---------- | ---------- | ---------- | ---------- | ---------- | ------------- | ------------- | ---------- | ---------- | ------------- |
| 1234        | **0.6111** | 0.6459 | 0.7390 | **0.4833** | 0.4964         | **0.5390**      | **0.6158**     | 0.5528     | **0.5879** | **0.6204** | **0.5749** | **0.6204** | 0.6775        | 0.6741        | **0.5731** | 0.6535     | 0.6041        |
| 1235        | **0.4642** | 0.6563 | 0.7196 | **0.4734** | **0.5390**     | **0.5506**      | **0.6218**     | 0.5912     | **0.5790** | **0.6759** | **0.5650** | **0.6357** | 0.6874        | **0.6947**    | **0.5757** | **0.6627** | 0.6058        |
| 1236        | **0.5647** | 0.6463 | 0.7349 | **0.4931** | **0.5260**     | **0.5197**      | 0.6086         | 0.5720     | **0.5888** | **0.6210** | **0.5831** | **0.6207** | 0.6808        | **0.6804**    | **0.5681** | 0.6598     | 0.6042        |
| 1237        | **0.4907** | 0.6459 | 0.7409 | **0.4742** | **0.5198**     | **0.5289**      | **0.6162**     | 0.5748     | **0.5861** | **0.5991** | **0.5869** | **0.6196** | 0.6835        | **0.6866**    | **0.5647** | **0.6642** | 0.5989        |
| 1238        | **0.4588** | 0.6512 | 0.7334 | **0.4769** | **0.5426**     | **0.5381**      | 0.6067         | 0.5547     | **0.6043** | **0.6205** | **0.5832** | **0.6370** | 0.6873        | **0.6812**    | **0.5644** | **0.6627** | 0.6002        |

According to the benchmarking results done by other researchers such as [vaeda](https://doi.org/10.1093/bioinformatics/btac720) and [scDblFinder](https://f1000research.com/articles/10-979), SimpleDoubletsDetection achieves the ***comparable*** **average AUPRC** compared to state-of-the-art method  [scDblFinder](https://f1000research.com/articles/10-979) across 16 benchmarking datasets. Notably, SimpleDoubletsDetection can achieves ***10-13 best of 16*** datasets.

### Running time
All these results are produced on a MacBook Air with M2 chip. To predict all 16 datasets, SimpleDoubletsDetection **only cost 237 seconds** (without iterative training) to **288 seconds** (with 5 iteration), which is almost the fastest method.

### Datasets
All the 16 datasets was colleted from the benchmarking paper of [Xi and Li](https://doi.org/10.1016/j.cels.2020.11.008). Datasets was transformed into **H5AD** format using [*sceasy*](https://github.com/cellgeni/sceasy). Processing script is convertH5AD.R.