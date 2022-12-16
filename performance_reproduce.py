import anndata
import numpy as np
import scanpy as sc
from sklearn.neural_network import MLPClassifier
from sklearn.decomposition import PCA
from sklearn.metrics import auc, average_precision_score, precision_recall_curve
import time
def test_logistic_clf(seed):
    auprc_list = []
    for dataset in ["J293t-dm","hm-6k","hm-12k","cline-ch","HEK-HMEC-MULTI","HMEC-orig-MULTI",
                "HMEC-rep-MULTI","mkidney-ch","nuc-MULTI","pbmc-1A-dm","pbmc-1B-dm","pbmc-1C-dm",
                "pbmc-2ctrl-dm","pbmc-2stim-dm","pdx-MULTI","pbmc-ch"]:
        adata = sc.read("./data/RNAseq/"+dataset+".h5ad")
        adata.X = adata.X.todense()
        idx_1 = np.arange(len(adata))
        idx_2 = np.arange(len(adata))
        np.random.seed(seed)
        np.random.shuffle(idx_1)
        np.random.seed(seed+1)
        np.random.shuffle(idx_2)
        sim_x = np.array(adata[idx_2].X + adata[idx_1].X)/2
        adata_sim = sc.AnnData(sim_x,var=adata.var)
        adata_sim.obs['label'] = 'simulated'
        adata_new = anndata.concat([adata_sim,adata]).copy()
        adata_new.layers['counts'] = adata_new.X
        doublets_in_real = np.zeros(adata.shape[0]).astype('bool')
        sim_y = np.zeros((sim_x.shape[0]))
        sim_y[:] = 1
        real_x = adata.X
        sim_x = np.log(sim_x+1)
        real_x = np.log(real_x+1)
        pca = PCA(n_components=50)
        real_x = pca.fit_transform(np.array(real_x))
        sim_x = pca.transform(np.array(sim_x))
        train_x = np.concatenate((sim_x,real_x)).copy()
        test_x = train_x[-adata.shape[0]:,]
        for i in range(1):
            real_y = np.zeros((real_x.shape[0]))
            real_y[doublets_in_real] = 1
            train_y = np.concatenate((sim_y,real_y)).copy()

            nn = MLPClassifier(hidden_layer_sizes=(),
                               batch_size=128,
                               activation='logistic',
                               solver='adam',
                               alpha=1e-4,
                               learning_rate_init=1e-4,
                               learning_rate='adaptive',
                               early_stopping=True,
                               tol=0,
                               random_state=seed,max_iter=1000).fit(train_x, train_y)

            prob_y = nn.predict_proba(test_x)
            if i ==0:
                doublets_in_real = doublets_in_real | (prob_y[:,1]>0.9)
            else:
                doublets_in_real = doublets_in_real & (prob_y[:, 1] > 0.9)
            pred_y = nn.predict(test_x)
            label_y = np.array(adata.obs['label'] == 'doublet').astype('int')
            precision, recall, thresholds = precision_recall_curve(pred_y, label_y)
            auprc = auc(recall, precision)
            print(dataset+" AUPRC:",auprc)
        auprc_list.append(auprc)
    print("Average AUPRC across 16 benchmarking datasets:",np.array(auprc_list).mean())

    return np.array(auprc_list).mean()

start = time.time()
test_logistic_clf(1236)
end = time.time()

print((end-start)/16)

