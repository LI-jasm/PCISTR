import numpy as np
import torch
from scipy.sparse import coo_matrix


def coototensor(A):
    """
    Convert a coo_matrix to a torch sparse tensor
    """

    values = A.data
    indices = np.vstack((A.row, A.col))
    i = torch.LongTensor(indices)
    v = torch.FloatTensor(values)
    shape = A.shape

    return torch.sparse.FloatTensor(i, v, torch.Size(shape))

def adj_matrix_weight_merge(A, adj_weight):
    """
    Multiplex Relation Aggregation
    """

    # 您的原始代码
    A_temp = []  # 初始化一个空列表
    for i in range(A.shape[1]):
        n = coototensor(A[0][i].tocoo())
        A_temp.append(n)
    A_t = torch.stack(A_temp, dim=2).to_dense()
    temp = torch.matmul(A_t, adj_weight)
    temp = torch.squeeze(temp, 2)
    return temp + temp.transpose(0, 1)
''' # Alibaba_small
    a = coototensor(A[0][0].tocoo())
    b = coototensor(A[0][1].tocoo())
    c = coototensor(A[0][2].tocoo())
    d = coototensor(A[0][3].tocoo())
    e = coototensor(A[0][4].tocoo())
    f = coototensor(A[0][5].tocoo())
    g = coototensor(A[0][6].tocoo())

    A_t = torch.stack([a, b, c, d, e, f, g], dim=2).to_dense()
'''




def normarlize(H):
    # DV = np.sum(H, axis=1)
    # DV += 1e-12
    # DV2 = np.mat(np.diag(np.power(DV, -1)))
    # G = DV2 * H

    DV = torch.sum(H, dim=1)
    DV += 1e-12

    DE=torch.sum(H,dim=0)
    DE += 1e-12
    DV2 = torch.diag(torch.pow(DV, -1))
    DE2 = torch.diag(torch.pow(DE, -1/2))
    G = torch.mm(DV2,H)
    G = torch.mm(G,DE2)

    return G

def construct_adj(encode, struct_weight):
    weight=torch.diag(struct_weight)
    adjust_encode=torch.mm(encode.to(torch.float32), weight)
    # print(adjust_encode)
    struct_adj=torch.mm(adjust_encode,adjust_encode.t())
    normal_struct_adj=torch.nn.functional.softmax(struct_adj, dim=1)
    return normal_struct_adj

