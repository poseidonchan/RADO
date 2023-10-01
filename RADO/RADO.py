import copy
import torch
import random
import anndata
import numpy as np
import torch.nn as nn
import scanpy as sc
from sklearn.model_selection import KFold
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score
from torch.utils.data.sampler import Sampler
from sklearn.decomposition import PCA
from torch.optim import Adam
from tqdm import tqdm


import warnings
warnings.filterwarnings("ignore")


def reproducibility(seed=1):
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True


device = "cuda:0" if torch.cuda.is_available() else "cpu"


def RADO(adata: anndata.AnnData = None,
         seed: int = 1234,
         kfold: int = 5,
         atac_data: bool = False
         ) -> anndata.AnnData:
    """
    :param adata: the original.py scRNA-seq data, without any other preprocessing.
    :param seed: random seed of the model
    :param kfold: use k-1 folds as training set to train the model and predict doublets on the rest 1 fold.
    :return: anndata with annotated doublets score in adata.obs['pred']
    """

    """deal with atac data"""
    if atac_data:
        print("Selecting hihgly variable peaks for ATAC-seq data...")
        sc.pp.highly_variable_genes(adata, n_top_genes=int(adata.shape[1]*0.25), flavor='seurat_v3')
        adata = adata[:, adata.var.highly_variable].copy()
        print("Done!")
    """simulate doublets"""
    print("Simulating doublets...")
    simulation_ratio = (kfold - 1) / kfold
    sim_x = SimulateDoublets(adata, seed=seed, simulation_ratio=simulation_ratio)
    print("Done!")
    """compute KNN score"""
    print("Computing KNN doublet score...")
    knn_score = naiveDoubletsScore(adata, sim_x, n_pcs=30, atac_data=atac_data)
    sim_knn_score = knn_score[:len(sim_x), ]
    real_knn_score = knn_score[len(sim_x):, ]
    adata.obs['RADO_KNN_doublet_score'] = real_knn_score
    adata.obs['RADO_KNN_doublet_call'] = np.where(real_knn_score > 0.5, 1, 0)
    print("Done!")
    """generate label of the simulated data"""
    sim_y = np.zeros((sim_x.shape[0], 1))
    sim_y[:, 0] = 1  # all simulated data are doublets
    real_x = adata.X
    print(real_x.shape, sim_x.shape)
    """normalize library size"""
    print("Normalizing library size...")
    sim_library_size = sim_x.sum(axis=1).reshape(-1, 1)
    real_library_size = real_x.sum(axis=1).reshape(-1, 1)
    sim_x = 1000 * sim_x / sim_library_size
    real_x = 1000 * real_x / real_library_size
    print("Done!")
    """log transformation"""
    print("Doing log transformation...")
    if isinstance(real_x, np.ndarray):
        pass
    else:
        real_x = np.array(real_x.todense())
    sim_x = np.log2(sim_x + 1)
    real_x = np.log2(real_x + 1)
    print("Done!")
    """PCA"""
    print("Doing PCA...")
    pca = PCA(n_components=10)
    train_x = np.concatenate([real_x, sim_x])
    train_x = pca.fit_transform(np.array(train_x))
    sim_x = train_x[adata.shape[0]:, :]
    real_x = train_x[:adata.shape[0], :]
    print("Done!")
    """the predicted doublets score"""
    pred_y = np.zeros(adata.shape[0])
    """Kfold"""
    kf = KFold(n_splits=kfold, shuffle=True, random_state=seed)
    
    print("Training and predicting...")
    for i, (train_idx, test_idx) in enumerate(tqdm(kf.split(adata.X), total=kf.get_n_splits())):
        train_x = np.concatenate((sim_x, real_x[train_idx])).copy()
        train_knn_score = np.concatenate((sim_knn_score, real_knn_score[train_idx]))
        real_y = np.zeros((real_x.shape[0], 1))
        real_y[:, 0] = 0
        train_y = np.concatenate((sim_y, real_y[train_idx])).copy()
        test_x = real_x[test_idx]
        test_knn_score = real_knn_score[test_idx]
        """define the logistic regression model"""
        clf = DoubletsClassifier(train_x, train_y, train_knn_score,
                                 lr=1e-4, seed=seed)
        clf.fit()

        pred = clf.predict_proba(test_x, test_knn_score)[:, 0]

        pred_y[test_idx] = pred
    print("Done!")
    adata.obs['RADO_doublet_score'] = pred_y
    adata.obs['RADO_doublet_call'] = np.where(pred_y > 0.5, 1, 0)
    return adata


def SimulateDoublets(adata: anndata.AnnData = None,
                     simulation_ratio: float = 0.8,
                     seed: int = 0,
                     ) -> np.ndarray:
    """
    simulate doublets by averaing two random droplets' expression value
    :param adata: anndata with raw gene expression in adata.X, shape as (n, m)
    :param simulation_ratio: simulated doublets number is simulation_ratio*m
    :param seed: random seed
    :return: return a simulated expression matrix of doublets, shape as (n, (1+ratio)*m)
    """
    adata = adata.copy()
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
    sim_x = (adata.X[idx_1] + adata.X[idx_2]) / 2
    sim_x = sim_x[np.random.choice(len(adata), int(len(adata) * simulation_ratio))]
    
    return np.array(sim_x)


def naiveDoubletsScore(adata: anndata.AnnData,
                       sim_x: np.ndarray,
                       n_pcs: int = 30,
                       atac_data: bool = False,
                      ) -> np.ndarray:
    """
    this function inspired by DoubletsFinder (https://doi.org/10.1016/j.cels.2019.03.003), the heterogeneous doublets
    will be easily detected from the low dimension manifold.
    :param adata: input anndata
    :param sim_x: simulated doublets data
    :param n_pcs: use how many PCs to calculated neighborhood graph
    :return: doublets score based on KNN, the score can be viewed as probability, ranges from 0-1.
    """
    adata_sim = sc.AnnData(X=sim_x, var=adata.var)
    adata_sim.obs['label'] = "simulated"
    adata_new = anndata.concat([adata_sim, adata])
    n_neighbors = int(np.sqrt(adata_new.shape[0]))
    if isinstance(adata_new.X, np.ndarray):
        pass
    else:
        adata_new.X = np.array(adata_new.X.todense())
    
    if atac_data:
        sc.pp.normalize_total(adata_new, target_sum=1e4)
        sc.pp.log1p(adata_new)
        sc.tl.pca(adata_new)
        sc.pp.neighbors(adata_new, n_neighbors=n_neighbors, n_pcs=50, metric="euclidean")
    else:
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


class DoubletsClassifier():
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 add_feature: np.ndarray = None,
                 lr: float = 1e-4,
                 seed: int = 1234,
                 batch_size: int = 256,
                 margin: float = 1,
                 sampling_rate: float = 0.5,
                 validation_size: float = 0.2):

        super().__init__()

        self.x = torch.from_numpy(x).float()
        self.y = torch.from_numpy(y).float().to(device)
        self.add_feature = torch.from_numpy(add_feature).float()
        self.x = torch.concat([self.x, self.add_feature], dim=1).to(device)

        self.x_train, self.x_val, self.y_train, self.y_val = train_test_split(self.x, self.y, test_size=validation_size,
                                                                              random_state=42)
        self.add_feature_train = self.x_train[:, -1].clone().reshape(-1, 1)
        self.add_feature_val = self.x_val[:, -1].clone().reshape(-1, 1)
        self.x_train = self.x_train[:, :-1]
        self.x_val = self.x_val[:, :-1]

        self.bestmodel = None
        self.lastepoch = None

        self.validation_size = validation_size
        reproducibility(seed)
        self.clf = LogisticRegression(n_in=x.shape[1], n_out=y.shape[1],
                                      n_cat=1 if add_feature is not None else 0).to(device)
        self.lr = lr
        self.batch_size = batch_size
        self.margin = margin
        self.totalPNum = int(self.y_train.sum())  # total positive (labeled as 1) sample number
        self.posNum = int(sampling_rate * self.batch_size)

    def fit(self, iterations=100000):
        trainset = SimpleDataset(self.x_train, self.add_feature_train, self.y_train)
        valset = SimpleDataset(self.x_val, self.add_feature_val, self.y_val)

        epochs = int(iterations * self.batch_size / len(self.x))

        val_loader = DataLoader(valset, drop_last=True,
                                batch_size=len(self.x_val), shuffle=True)

        optimizer = Adam(self.clf.parameters(), lr=self.lr, weight_decay=0)

        train_loss = STOAPLOSS(margin=self.margin,
                               totalPNum=self.totalPNum, posNum=self.posNum,
                               data_length=len(trainset),
                               device=device)

        labels = self.y_train.cpu().long().numpy().reshape(-1)
        labels = list(labels)

        train_loader = DataLoader(trainset, batch_size=self.batch_size,
                                  sampler=AUPRCSampler(labels, self.batch_size, posNum=self.posNum))

        self.clf, train_loss, val_loss = self._training_stage(self.clf, [train_loader, val_loader], loss_fn=train_loss,
                                                              optimizer=optimizer, epochs=epochs)

    def _training_stage(self, clf, loaders, loss_fn=None, optimizer=None, epochs=128):
        clf_ = copy.deepcopy(clf)

        train_loss = []
        val_loss = []

        for i in (range(epochs)):
            train_epoch_loss = 0
            val_epoch_loss = 0
            clf.train()
            for k, (index, x, x_cat, y) in enumerate(loaders[0]):
                optimizer.zero_grad()
                state_dict = copy.deepcopy(clf.state_dict())
                y_pred = clf(x, x_cat)
                with torch.no_grad():
                    y_pred_ = clf_(x, x_cat)
                train_batch_loss = loss_fn(f_ps=y_pred[0:self.posNum], f_ns=y_pred[self.posNum:],
                                           f_ps_=y_pred_[0:self.posNum], f_ns_=y_pred_[self.posNum:],
                                           index_s=index[0:self.posNum])
                train_batch_loss.backward()
                optimizer.step()
                train_epoch_loss += train_batch_loss

                clf_.load_state_dict(state_dict)

            train_loss.append(train_epoch_loss.cpu().detach().numpy())
            torch.save(clf.state_dict(), "./model/model_epoch=" + str(i) + ".pth")

            for k, (index, x, x_cat, y) in enumerate(loaders[1]):
                y_pred = clf(x, x_cat)
                ap = average_precision_score(y.cpu().detach().numpy(), y_pred.cpu().detach().numpy())
                val_epoch_loss += ap

            val_loss.append(val_epoch_loss)

            if (len(val_loss) - np.argmax(val_loss) <= 10) or (val_loss[-1] < 0.7):
                pass
            else:
                # print("early stopping at epoch " + str(np.argmax(val_loss)))
                self.bestmodel = str(np.argmax(val_loss))
                break

            self.lastepoch = str(i)
        self.bestmodel = str(np.argmax(val_loss))
        # print(self.bestmodel)
        return clf, train_loss, val_loss

    def predict_proba(self, test_x, test_x_knn):
        if self.bestmodel is not None:
            # print("load best model")
            self.clf.load_state_dict(torch.load("./model/model_epoch=" + self.bestmodel + ".pth"))
        else:
            self.clf.load_state_dict(torch.load("./model/model_epoch=" + self.lastepoch + ".pth"))
        self.clf.eval()
        test_x = torch.from_numpy(test_x).float().to(device)
        test_x_knn = torch.from_numpy(test_x_knn).float().to(device)
        pred = self.clf(test_x, test_x_knn)
        return pred.cpu().detach().numpy()


class SimpleDataset(Dataset):
    def __init__(self,
                 X: np.ndarray,
                 X_add: np.ndarray,
                 Y: np.ndarray = None,
                 ):
        self.X = X
        self.Y = Y
        self.targets = Y.cpu().numpy()
        self.X_add = X_add

        self.pos_indices = np.flatnonzero(self.targets == 1)
        self.pos_index_map = {}
        for i, idx in enumerate(self.pos_indices):
            self.pos_index_map[idx] = i

    def __len__(self):
        return len(self.X)

    def __getitem__(self, index):
        if self.Y is not None:
            x = self.X[index]
            y = self.Y[index]
            x_add = self.X_add[index]

            index = self.pos_index_map[index] if index in self.pos_indices else -1

            return index, x, x_add, y
        else:
            x = self.X[index].float()
            x_add = self.X_add[index]
            return x, x_add


class LogisticRegression(nn.Module):
    def __init__(self,
                 n_in: int,
                 n_out: int,
                 n_cat: int = None):
        super().__init__()

        self.fc = nn.Sequential(
            nn.Linear(n_in + n_cat, n_out),
            nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001),
            nn.Sigmoid()
        )

    def forward(self, x, x_cat):
        x = torch.cat((x, x_cat), dim=1)
        return self.fc(x)

"""
Code below is modified from Xidong Wu (https://github.com/xidongwu), the original algorithm publication 
is available from https://doi.org/10.1109/ICDM54844.2022.00068
"""

class AUPRCSampler(Sampler):

    def __init__(self, labels, batchSize, posNum=1):
        # positive class: minority class
        # negative class: majority class

        self.labels = labels
        self.posNum = posNum
        self.batchSize = batchSize

        self.clsLabelList = np.unique(labels)
        self.dataDict = {}

        for label in self.clsLabelList:
            self.dataDict[str(label)] = []

        for i in range(len(self.labels)):
            self.dataDict[str(self.labels[i])].append(i)

        self.ret = []

    def __iter__(self):
        minority_data_list = self.dataDict[str(1)]
        majority_data_list = self.dataDict[str(0)]

        np.random.shuffle(minority_data_list)
        np.random.shuffle(majority_data_list)

        # In every iteration : sample 1(posNum) positive sample(s), and sample batchSize - 1(posNum) negative samples
        if len(minority_data_list) // self.posNum > len(majority_data_list) // (
                self.batchSize - self.posNum):  # At this case, we go over the all positive samples in every epoch.
            # extend the length of majority_data_list from  len(majority_data_list) to len(minority_data_list)* (batchSize-posNum)
            majority_data_list.extend(np.random.choice(majority_data_list, len(minority_data_list) // self.posNum * (
                        self.batchSize - self.posNum) - len(majority_data_list), replace=True).tolist())

        elif len(minority_data_list) // self.posNum < len(majority_data_list) // (
                self.batchSize - self.posNum):  # At this case, we go over the all negative samples in every epoch.
            # extend the length of minority_data_list from len(minority_data_list) to len(majority_data_list)//(batchSize-posNum) + 1s
            minority_data_list.extend(np.random.choice(minority_data_list, len(majority_data_list) // (
                        self.batchSize - self.posNum) * self.posNum - len(minority_data_list), replace=True).tolist())

        # print("AUPRCSampler333", len(minority_data_list), len(majority_data_list))

        self.ret = []
        for i in range(len(minority_data_list) // self.posNum):
            self.ret.extend(minority_data_list[i * self.posNum:(i + 1) * self.posNum])
            startIndex = i * (self.batchSize - self.posNum)
            endIndex = (i + 1) * (self.batchSize - self.posNum)
            self.ret.extend(majority_data_list[startIndex:endIndex])

        return iter(self.ret)

    def __len__(self):
        return len(self.ret)


class STOAPLOSS(nn.Module):
    def __init__(self, margin, totalPNum, posNum, data_length, device):
        '''
        :param margin: margin for squred hinge loss
        '''
        super(STOAPLOSS, self).__init__()
        self.u_all = torch.tensor([0.0] * data_length).view(-1, 1).to(device)
        self.u_pos = torch.tensor([0.0] * data_length).view(-1, 1).to(device)
        self.threshold = margin
        self.device = device
        self.lmt = min(totalPNum / posNum, 1.5)
        self.alpha = 0.1

    def forward(self, f_ps, f_ns, f_ps_, f_ns_, index_s):
        f_ps = f_ps.view(-1)
        f_ns = f_ns.view(-1)

        f_ps_ = f_ps_.view(-1)
        f_ns_ = f_ns_.view(-1)

        vec_dat = torch.cat((f_ps, f_ns), 0)
        mat_data = vec_dat.repeat(len(f_ps), 1)
        vec_dat_ = torch.cat((f_ps_, f_ns_), 0)
        mat_data_ = vec_dat_.repeat(len(f_ps_), 1)

        f_ps = f_ps.view(-1, 1)
        f_ps_ = f_ps_.view(-1, 1)

        neg_mask = torch.ones_like(mat_data)  # 1
        neg_mask[:, 0:f_ps.size(0)] = 0
        pos_mask = torch.zeros_like(mat_data)  # 0
        pos_mask[:, 0:f_ps.size(0)] = 1

        neg_loss = torch.max(self.threshold - (f_ps - mat_data), torch.zeros_like(mat_data)) ** 2 * neg_mask
        pos_loss = torch.max(self.threshold - (f_ps - mat_data), torch.zeros_like(mat_data)) ** 2 * pos_mask
        neg_loss_ = torch.max(self.threshold - (f_ps_ - mat_data_), torch.zeros_like(mat_data_)) ** 2 * neg_mask
        pos_loss_ = torch.max(self.threshold - (f_ps_ - mat_data_), torch.zeros_like(mat_data_)) ** 2 * pos_mask

        loss = pos_loss + neg_loss
        loss_ = pos_loss_ + neg_loss_

        self.u_pos = (1 - self.alpha) * self.u_pos
        self.u_all = (1 - self.alpha) * self.u_all

        if f_ps.size(0) == 1:
            self.u_pos[index_s] += (pos_loss.mean() - (1 - self.alpha) * pos_loss_.mean()) * self.lmt
            self.u_all[index_s] += (loss.mean() - (1 - self.alpha) * loss_.mean()) * self.lmt
        else:
            self.u_pos[index_s] += (pos_loss.mean(1, keepdim=True) - (1 - self.alpha) * pos_loss_.mean(1,
                                                                                                       keepdim=True).detach_()) * self.lmt
            self.u_all[index_s] += (loss.mean(1, keepdim=True) - (1 - self.alpha) * loss_.mean(1,
                                                                                               keepdim=True).detach_()) * self.lmt

        p = (self.u_pos[index_s] - (self.u_all[index_s]) * pos_mask) / (self.u_all[index_s] ** 2)

        p.detach_()  # detach
        loss = torch.mean(p.detach_() * loss)

        return loss


