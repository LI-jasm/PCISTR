import torch
import pickle
from scipy.sparse import csc_matrix, eye
import numpy as np
from src.Utils import get_model
from src.args import get_citation_args
from src.emb_train import emb_train
from src.Complex_identification_DPPN import get_nets, get_completes

args = get_citation_args()

# krogan
args.dataset = 'krogan14k'
# path ="./data/"+args.dataset+"/" + args.dataset + ".pkl"
path ="./data/"+args.dataset+"/" + args.dataset + "_attr_vector.txt"
data =open(path, 'rb')
# feature_s = csc_matrix(pickle.load(data))

file=open(path)
vector=[]
for i in file:
    if not i.strip(): continue
    node_name=i.split(' ',1)[0]
    node_vector=i.split(' ',1)[1].rstrip('\n')
    node_vector = node_vector.split()
    node_vector = list(map(float, node_vector))
    vector.append(node_vector)

feature_s=torch.tensor(vector)

file.close()



train_data = np.array(get_nets(args.dataset))[np.newaxis, :]
A = train_data
num_prot = train_data[0][0].shape[0]
num_view = train_data.shape[1]
# feature = csc_matrix(eye(num_prot))

encode = np.zeros((num_prot, num_view))
# 遍历每个稀疏矩阵
for i in range(num_view):
    # 假设当前稀疏矩阵为 current_sparse_matrix
    current_sparse_matrix = train_data[0, i]
    # 计算当前稀疏矩阵中每个节点的度，这里使用 axis=1 对行进行求和
    node_degrees = np.array(current_sparse_matrix.sum(axis=1)).flatten()
    # 将结果存储到相应的列
    encode[:, i] = node_degrees
encode = torch.tensor(encode)


def normarlize(H):
    DV = torch.sum(H, dim=1)
    DV += 1e-12
    DV2 = torch.diag(torch.pow(DV, -1))
    G = torch.mm(DV2, H)

    return G

# new_adj=normarlize(new_adj)



# adj, features, labels, idx_train, idx_val, idx_test = load_our_data(args.dataset, args.cuda)

model = get_model(args.model, feature_s.shape[1], A, args.hidden, args.out, args.dropout, False)


emb_train(model, args.dataset, feature_s, A, encode, device=torch.device('cpu'))

