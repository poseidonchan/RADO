import time
import anndata
import warnings
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import scale
from sklearn.model_selection import KFold
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from sklearn.metrics import auc, precision_recall_curve,


def main():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--iter", dest="iteration", default=1, type=int,
        help="iterations to perform PU-learning"
    )

    parser.add_argument(
        "--kfold", dest="kfold", default=5, type=int,
        help="iterations to perform PU-learning"
    )

    parser.add_argument(
        "--seed", dest="seed", default=1234, type=int,
        help="pin random seed to reproduce results"
    )
    parser.add_argument(
        "--dir",
        dest="dir",
        default=None,
        help="specify the directory that contains all the datasets",
    )
    args = parser.parse_args()
    auprc_list = []
    for dataset in ["cline-ch", "HEK-HMEC-MULTI", "hm-12k", "hm-6k", "HMEC-orig-MULTI",
                    "HMEC-rep-MULTI", "J293t-dm", "mkidney-ch", "nuc-MULTI", "pbmc-1A-dm", "pbmc-1B-dm", "pbmc-1C-dm",
                    "pbmc-2ctrl-dm", "pbmc-2stim-dm", "pbmc-ch", "pdx-MULTI"]:

        adata = sc.read(args.dir + dataset + ".h5ad")
        sim_x = SimulateDoublets(adata, flavor='sum', cluster_aware=False, seed=args.seed)

        sim_y = np.zeros((sim_x.shape[0]))
        sim_y[:] = 1
        real_x = adata.X
        """compute KNN score"""
        knn_score = naiveDoubletsScore(adata, sim_x, n_pcs=30)
        sim_knn_score = knn_score[:len(sim_x), ]
        real_knn_score = knn_score[len(sim_x):, ]
        """Normalize library size"""
        sim_library_size = sim_x.sum(axis=1).reshape(-1, 1)
        real_library_size = real_x.sum(axis=1).reshape(-1, 1)
        sim_x = 10000 * sim_x / sim_library_size
        real_x = 10000 * real_x / real_library_size
        """Log transformation"""
        if isinstance(real_x, np.ndarray):
            pass
        else:
            real_x = np.array(real_x.todense())
        sim_x = np.log2(sim_x + 1)
        real_x = np.log2(real_x + 1)
        """Scale"""
        train_x = np.concatenate([real_x, sim_x])
        train_x = scale(train_x.T).T
        """PCA"""
        pca = PCA(n_components=10)
        train_x = pca.fit_transform(train_x)
        sim_x = train_x[adata.shape[0]:, :]
        real_x = train_x[:adata.shape[0], :]
        # real_x = pca.fit_transform(real_x)
        # sim_x = pca.transform(np.array(sim_x))
        pu_pred = np.zeros((adata.shape[0], args.iteration))

        # display plot
        for k in range(args.iteration):
            kf = KFold(n_splits=args.kfold, shuffle=True)
            pu_pred_k = np.zeros((adata.shape[0], args.kfold))
            for i, (train_idx, test_idx) in enumerate(kf.split(real_x)):
                train_x = np.concatenate((sim_x, real_x[train_idx])).copy()
                train_knn_score = np.concatenate((sim_knn_score, real_knn_score[train_idx]))
                # train_library_size = np.concatenate((sim_library_size,real_library_size[train_idx]))
                train_x = np.concatenate((train_x, train_knn_score), axis=1)
                real_y = np.zeros((real_x.shape[0]))
                train_y = np.concatenate((sim_y, real_y[train_idx])).copy()
                test_x = real_x[test_idx]
                test_knn_score = real_knn_score[test_idx]
                # test_library_size = real_library_size[test_idx]
                test_x = np.concatenate((test_x, test_knn_score), axis=1)

                clf = LogisticRegression(penalty="l1", max_iter=1000,
                                         solver='saga', C=0.1,
                                         random_state=args.seed,
                                         n_jobs=-1).fit(train_x, train_y)

                pred_y = clf.predict_proba(test_x)[:, 1]
                pu_pred_k[test_idx, i] = pred_y

            pu_pred[:, k] = pu_pred_k.sum(axis=1)

        prob_pred_y = np.mean(pu_pred, axis=1)
        label_y = np.array(adata.obs['label'] == 'doublet').astype('int')
        precision, recall, thresholds = precision_recall_curve(label_y, prob_pred_y)
        auprc = auc(recall, precision)
        print(dataset + " AUPRC:", auprc)
        auprc_list.append(auprc)

    print("Average AUPRC across 16 benchmarking datasets:", np.array(auprc_list).mean())
    return np.array(auprc_list).mean()


def SimulateDoublets(adata: anndata.AnnData = None,
                     flavor: str = "average",
                     cluster_aware: bool = False,
                     fraction: float = 1.0,
                     seed: int = 0,
                     ) -> np.ndarray:
    if cluster_aware:
        adata = adata.copy()
        if isinstance(adata.X, np.ndarray):
            pass
        else:
            adata.X = np.array(adata.X.todense())
        adata.raw = adata
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        sc.tl.leiden(adata, resolution=0.3)
        idx_1 = np.arange(len(adata))
        idx_2 = np.arange(len(adata))
        np.random.seed(seed)
        np.random.shuffle(idx_1)
        np.random.seed(seed + 1)
        np.random.shuffle(idx_2)
        for i in range(adata.shape[0]):
            while adata.obs['leiden'][idx_1[i]] == adata.obs['leiden'][idx_2[i]]:
                idx_2[i] = (idx_2[i] + 17) % adata.shape[0]
        if flavor == 'average':
            sim_x = (adata.raw.X[idx_1] + adata.raw.X[idx_2]) / 2
            sim_x = sim_x[np.random.choice(len(adata), int(len(adata) * fraction))]
        elif flavor == 'sum':
            sim_x = (adata.raw.X[idx_1] + adata.raw.X[idx_2])
            sim_x = sim_x[np.random.choice(len(adata), int(len(adata) * fraction))]
        else:
            raise NotImplementedError
    else:
        if isinstance(adata.X, np.ndarray):
            pass
        else:
            adata.X = np.array(adata.X.todense())
        idx_1 = np.arange(len(adata))
        idx_2 = np.arange(len(adata))
        np.random.seed(seed)
        np.random.shuffle(idx_1)
        np.random.seed(seed + 1)
        np.random.shuffle(idx_2)
        if flavor == 'average':
            sim_x = (adata.X[idx_1] + adata.X[idx_2]) / 2
            sim_x = sim_x[np.random.choice(len(adata), int(len(adata) * fraction))]
        elif flavor == 'sum':
            sim_x = (adata.X[idx_1] + adata.X[idx_2])
            sim_x = sim_x[np.random.choice(len(adata), int(len(adata) * fraction))]
        else:
            raise NotImplementedError

    return np.array(sim_x)


def naiveDoubletsScore(adata, sim_x, n_pcs=32):
    adata_sim = sc.AnnData(X=sim_x, var=adata.var)
    adata_sim.obs['label'] = "simulated"
    adata_new = anndata.concat([adata_sim, adata])
    n_neighbors = int(np.sqrt(adata_new.shape[0]))
    if isinstance(adata_new.X, np.ndarray):
        pass
    else:
        adata_new.X = np.array(adata_new.X.todense())

    sc.pp.normalize_total(adata_new, target_sum=1e4)
    sc.pp.log1p(adata_new)
    sc.tl.pca(adata_new)
    sc.pp.neighbors(adata_new, n_neighbors=n_neighbors, n_pcs=n_pcs, metric="euclidean")
    simulated_score = []
    simulated = np.array(adata_new.obs['label'] == 'simulated')

    adata_new.obsp['distances'] = adata_new.obsp['distances'].todense()
    neighbors = (adata_new.obsp['distances'] > 0)
    for i in range(len(adata_new)):
        simulated_score.append((neighbors[i, :] & simulated).sum() / n_neighbors)

    return np.array(simulated_score).reshape(-1, 1)


if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print("Total time consumption is", end - start, "seconds")
