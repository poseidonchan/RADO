import anndata
import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn.decomposition import PCA

def SDD(adata: anndata.AnnData,
        npc_components: int = 50,
        randon_seed: int = 0,
        n_iter: int = 1,
        confident_threshold: float = 0.9,
        ):
    adata = adata.copy()
    if isinstance(adata.X, np.ndarray):
        pass
    else:
        adata.X = adata.X.todense()
    idx_1 = np.arange(len(adata))
    idx_2 = np.arange(len(adata))
    np.random.seed(randon_seed)
    np.random.shuffle(idx_1)
    np.random.seed(randon_seed + 1)
    np.random.shuffle(idx_2)
    sim_x = np.array(adata[idx_2].X + adata[idx_1].X) / 2
    doublets_in_real = np.zeros(adata.shape[0]).astype('bool')
    sim_y = np.zeros((sim_x.shape[0]))
    sim_y[:] = 1
    real_x = adata.X
    sim_x = np.log(sim_x + 1)
    real_x = np.log(real_x + 1)
    pca = PCA(n_components=npc_components)
    real_x = pca.fit_transform(np.array(real_x))
    sim_x = pca.transform(np.array(sim_x))
    train_x = np.concatenate((sim_x, real_x)).copy()
    test_x = train_x[-adata.shape[0]:, ]
    for i in range(n_iter):
        real_y = np.zeros((real_x.shape[0]))
        real_y[doublets_in_real] = 1
        train_y = np.concatenate((sim_y, real_y)).copy()

        nn = MLPClassifier(hidden_layer_sizes=(),
                           batch_size=128,
                           solver='adam',
                           alpha=1e-4,
                           learning_rate_init=1e-4,
                           learning_rate='adaptive',
                           early_stopping=True,
                           tol=0,
                           random_state=randon_seed, max_iter=1000).fit(train_x, train_y)

        prob_y = nn.predict_proba(test_x)
        if i == 0:
            doublets_in_real = doublets_in_real | (prob_y[:, 1] > confident_threshold)
        else:
            doublets_in_real = doublets_in_real & (prob_y[:, 1] > confident_threshold)
        pred_y = nn.predict(test_x)

    return pred_y


