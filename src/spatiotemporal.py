##############Protein Complex identification algorithm based on DPPN########################
import sys
from math import sqrt

from optparse import OptionParser
from scipy import *
import torch as th
import scipy as sp
from scipy import sparse
from scipy.linalg import *

from numpy import *
import numpy as np

from numpy.linalg import *
import string
import os

from scipy.sparse import coo_matrix, save_npz, vstack


def f_key(a):
    return(a[-1])


def getOptions():
    optparser = OptionParser(usage="%prog [options]\n-h for help")
    
    optparser.add_option("-e", "--Gene_expression_file", dest="EXPRESSION_file", help="Input file of gene expression data",default="./data/series_matrix.txt")
    
    optparser.add_option("-o", "--output", dest="output", help="Output file for writing the identified complexes",default="result.txt",)
   
    optparser.add_option("-c", "--core_thes_parameter", dest="Core_thes_parameter", help="expand_thes_parameter: range from 0 to 1", default="0.09")

    optparser.add_option("-a", "--attach_thes_parameter", dest="Attach_thes_parameter", help="expand_thes_parameter: range from 0 to 1", default="0.05")

    
    
    (options, args) = optparser.parse_args()

    if not options.EXPRESSION_file:
        optparser.error("No input Gene expression data file")
    if not options.output:
        optparser.error("No output file defined")
    return options, args



# if __name__ == "__main__":

def get_nets(dataset):
    protein_out_file = "./data/" + dataset + "/protein.temp"

    options, args = getOptions()
    Dic_map={}
    index=0
    All_node=set([])
    All_node_index=set([])
    go_set_list=[]
    expression_list=[]
    subcellular_list=[]
    All_node_go=set([])
    neighbor_list=[]
    protein_time_list=[]
    time_protein_list = []
    protein_spacial_list=[]
    Seed_edge_list=[]
    Cliques_list=[]
    Core_list=[]
    Complex_list=[]
    Core_protein_set=set([])

    Core_edge_set=set([])
    Protein_SD_weight={}
    
    Node1=[]
    Node2=[]
    
    Node_PPI_1=[]
    Node_PPI_2=[]

    neighbor_PPI_list=[]

    locations = ["Extracellular_space", "Nucleus", "Mitochondrion", "Endosome",
                 "Vacoule", "Peroxisome", "Endoplasmic_Reticulum", "Golgi_apparatus",
                 "Cytosol", "Cytoskeleton", "Plasma_membrane"]
    Spacial_num=11
    Time_num=12
    Times_for_SD=3
    
    SD3_weight=0.99  # 三个活性概率
    SD2_weight=0.95
    SD1_weight=0.68

    
    Seed_weight_thresh=float(options.Core_thes_parameter)  # 两个阈值
    Attach_thresh=float(options.Attach_thes_parameter)
    
   
    ###########input high-throughput PPI data############  # 输入高通量PPI数据
    PPI_file = "./data/" + dataset + "/" + dataset + ".txt"
    f = open(PPI_file,"r")
    f_protein_out=open(protein_out_file,"w")
    

    
    for line in f:
        line = line.strip().split()
        if len(line)==2:
            
            
            if line[0] not in Dic_map:
                
                Dic_map[line[0]]=index  # 建立蛋白质节点索引
                f_protein_out.write(line[0]+"\n")  # 将蛋白质节点名称放到protein.temp文件中
                
                index+=1
                go_set_list.append(set([]))
                expression_list.append([])
                subcellular_list.append([])
                neighbor_list.append(set([]))
                neighbor_PPI_list.append(set([]))
                protein_time_list.append(set([]))
                protein_spacial_list.append(set([]))
                
                
            if line[1] not in Dic_map:
                
                Dic_map[line[1]]=index
                f_protein_out.write(line[1]+"\n")
                
                index+=1
                go_set_list.append(set([]))
                expression_list.append([])
                subcellular_list.append([])
                neighbor_list.append(set([]))
                neighbor_PPI_list.append(set([]))
                protein_time_list.append(set([]))
                protein_spacial_list.append(set([]))

           
                
    Protein_num=index
    f.close()
    f_protein_out.close()

    ###########input Gene expression data##################  # 输入基因表达数据
   
    f=open(options.EXPRESSION_file,"r")
    
    
    for line in f:
        line = line.strip().split()
        if len(line)==38:
            if line [1] in Dic_map:
                
                for i in range(2,Time_num+2):
                    expression_list[Dic_map[line[1]]].append((float(line[i])+float(line[i+12])+float(line[i+24]))/3)
                    # 计算蛋白质节点在12个时间节点上的活性表达值（三个周期的平均值）
                
           
    f.close()


    ###############compute active time attribute for proteins###################  # 计算蛋白质活性时间属性（包括活跃时间点及对应的活性概率）

    for t in range(0,Time_num):
        time_protein_list.append(set([]))  # 建立12个时间集合

    for instance in Dic_map:
        
        if len(expression_list[Dic_map[instance]])>=12:  # 如果该蛋白质有基因表达数据
            Temp_mean_value= 0.
            Temp_sd_value=0.
            Temp_thresh_3SD=0.
            Temp_thresh_2SD=0.
            Temp_thresh_1SD=0.
            
            for j in range(0,Time_num):
                Temp_mean_value+=expression_list[Dic_map[instance]][j]
            Temp_mean_value/=Time_num  # 计算蛋白质节点在12个时间节点上的平均表达值

            for j in range(0,Time_num):
                Temp_sd_value+=(expression_list[Dic_map[instance]][j]-Temp_mean_value)**2
            Temp_sd_value/=(Time_num-1)  # 计算蛋白质节点在12个时间节点上的表达值方差

            k=1

            Temp_thresh_3SD=Temp_mean_value+3*(Temp_sd_value**0.5)*(Temp_sd_value/(1+Temp_sd_value))
            Temp_thresh_2SD=Temp_mean_value+2*(Temp_sd_value**0.5)*(Temp_sd_value/(1+Temp_sd_value))
            Temp_thresh_1SD=Temp_mean_value+1*(Temp_sd_value**0.5)*(Temp_sd_value/(1+Temp_sd_value))
            Temp_thresh_kSD=Temp_mean_value+k*(Temp_sd_value**0.5)*(Temp_sd_value/(1+Temp_sd_value))
            '''for j in range(0,Time_num):
                if expression_list[Dic_map[instance]][j]>=Temp_thresh_3SD:
                    protein_time_list[Dic_map[instance]].add(j)
                    Protein_SD_weight[instance+"&"+str(j)]=SD3_weight
                    
                elif Temp_thresh_3SD>expression_list[Dic_map[instance]][j]>=Temp_thresh_2SD:
                    protein_time_list[Dic_map[instance]].add(j)
                    Protein_SD_weight[instance+"&"+str(j)]=SD2_weight
                    
                elif Temp_thresh_2SD>expression_list[Dic_map[instance]][j]>=Temp_thresh_1SD:
                    protein_time_list[Dic_map[instance]].add(j)  # 记录蛋白质节点的活跃时间节点
                    Protein_SD_weight[instance+"&"+str(j)]=SD1_weight  # 记录蛋白质活跃时间节点对应的活性概率'''

            for j in range(0, Time_num):
                if expression_list[Dic_map[instance]][j] >= Temp_thresh_kSD:
                    protein_time_list[Dic_map[instance]].add(j)  # 记录蛋白质节点的活跃时间节点
                    time_protein_list[j].add(Dic_map[instance])  # 记录各时间活跃的蛋白质
#------------------------------------------------亚细胞位置信息----------------------------------------------
    f_c = open("./data/yeast_compartment_benchmark.txt", "r",encoding='utf-16')

    location_dict = {}
    for i, location in enumerate(locations):
        location_dict[location] = i

    for line in f_c:
        line = line.strip().split()
        if line[0] in Dic_map:
            subcellular_list[Dic_map[line[0]]].append(line[1])  # 记录蛋白质节点的活跃空间位置(位置文字表示)
            protein_spacial_list[Dic_map[line[0]]].add(location_dict[line[1]])  # 记录蛋白质节点的活跃空间位置（位置序号表示）


    #####################时空子图构建########################################
    f = open(PPI_file, "r")
    matrix_list = []
    rows = [[] for _ in range(132)]
    cols = [[] for _ in range(132)]

    for line in f:
        line = line.strip().split()
        if len(line) == 2:
            if len(protein_time_list[Dic_map[line[0]]]&protein_time_list[Dic_map[line[1]]])>0 and \
                    len(protein_spacial_list[Dic_map[line[0]]] & protein_spacial_list[Dic_map[line[1]]]) > 0:
                ctime = protein_time_list[Dic_map[line[0]]] & protein_time_list[Dic_map[line[1]]]
                cspace = protein_spacial_list[Dic_map[line[0]]] & protein_spacial_list[Dic_map[line[1]]]
                ts = [(x+1) * (y+1) for x in ctime for y in cspace]
                for i in ts:
                    rows[i-1].append(Dic_map[line[0]])
                    cols[i-1].append(Dic_map[line[1]])


    for ts in range(0, Time_num * Spacial_num):
        matrix = coo_matrix((np.ones(len(rows[ts]), dtype=float64), (rows[ts], cols[ts])), shape=(len(Dic_map), len(Dic_map)))
        # 将对象添加到列表中
        matrix_list.append(matrix)

    filtered_matrix_list = [matrix for matrix in matrix_list if matrix.nnz > 0]

    f.close()


    ###########input PPI data############  # 输入PPI数据（有共活跃时间的PPI构建的新PPI网络）
    f = open(PPI_file, "r")
    f_out = open("./data/" + dataset + "/new_PPI_file", "w")

    for line in f:
        line = line.strip().split()
        if len(line) == 2:
            if line[0] != line[1]:
                Node_PPI_1.append(line[0])

                Node_PPI_2.append(line[1])

                neighbor_PPI_list[Dic_map[line[0]]].add(line[1])  # 记录原始PPI网络中各节点的邻居节点
                neighbor_PPI_list[Dic_map[line[1]]].add(line[0])

                if len(protein_time_list[Dic_map[line[0]]] & protein_time_list[Dic_map[line[1]]]) > 0 \
                        or len(protein_time_list[Dic_map[line[0]]]) == 0 or len(
                    protein_time_list[Dic_map[line[1]]]) == 0:
                    # 如果在原始PPI中的邻接节点有共同活跃的时间点，将它们加入彼此的邻居节点列表（在新PPI网络中保留）
                    if len(protein_spacial_list[Dic_map[line[0]]]) == 0 or len(
                            protein_spacial_list[Dic_map[line[1]]]) == 0:
                        Node1.append(line[0])

                        Node2.append(line[1])

                        neighbor_list[Dic_map[line[0]]].add(line[1])
                        neighbor_list[Dic_map[line[1]]].add(line[0])

                        f_out.write(line[0] + " " + line[1] + "\n")  # 将过滤后PPI（邻接节点有共活跃时间）记录在新的PPI文件中（new_ppi.txt）

                    elif len(protein_spacial_list[Dic_map[line[0]]] & protein_spacial_list[Dic_map[line[1]]]) > 0:
                        Node1.append(line[0])

                        Node2.append(line[1])

                        neighbor_list[Dic_map[line[0]]].add(line[1])
                        neighbor_list[Dic_map[line[1]]].add(line[0])

                        f_out.write(line[0] + " " + line[1] + "\n")  # 将过滤后PPI（邻接节点有共活跃时间）记录在新的PPI文件中（new_ppi.txt）

    f_out.close()
    f.close()
    return filtered_matrix_list


def get_completes(dataset):

    protein_out_file = "./data/" + dataset + "/protein.temp"

    options, args = getOptions()
    Dic_map = {}
    index = 0
    All_node = set([])
    All_node_index = set([])
    go_set_list = []
    expression_list = []
    subcellular_list = []
    All_node_go = set([])
    neighbor_list = []
    protein_time_list = []
    time_protein_list = []
    protein_spacial_list = []

    tNode1 = []
    tNode2 = []

    sNode1 = []
    sNode2 = []

    Node_PPI_1 = []
    Node_PPI_2 = []

    neighbor_PPI_list = []

    locations = ["Extracellular_space", "Nucleus", "Mitochondrion", "Endosome",
                 "Vacoule", "Peroxisome", "Endoplasmic_Reticulum", "Golgi_apparatus",
                 "Cytosol", "Cytoskeleton", "Plasma_membrane"]
    Spacial_num = 11
    Time_num = 12
    Times_for_SD = 3

    SD3_weight = 0.99  # 三个活性概率
    SD2_weight = 0.95
    SD1_weight = 0.68

    Seed_weight_thresh = float(options.Core_thes_parameter)  # 两个阈值
    Attach_thresh = float(options.Attach_thes_parameter)

    ###########input high-throughput PPI data############  # 输入高通量PPI数据
    PPI_file = "./data/" + dataset + "/" + dataset + ".txt"
    f = open(PPI_file, "r")
    f_protein_out = open(protein_out_file, "w")

    for line in f:
        line = line.strip().split()
        if len(line) == 2:

            if line[0] not in Dic_map:
                Dic_map[line[0]] = index  # 建立蛋白质节点索引
                f_protein_out.write(line[0] + "\n")  # 将蛋白质节点名称放到protein.temp文件中

                index += 1
                go_set_list.append(set([]))
                expression_list.append([])
                subcellular_list.append([])
                neighbor_list.append(set([]))
                neighbor_PPI_list.append(set([]))
                protein_time_list.append(set([]))
                protein_spacial_list.append(set([]))

            if line[1] not in Dic_map:
                Dic_map[line[1]] = index
                f_protein_out.write(line[1] + "\n")

                index += 1
                go_set_list.append(set([]))
                expression_list.append([])
                subcellular_list.append([])
                neighbor_list.append(set([]))
                neighbor_PPI_list.append(set([]))
                protein_time_list.append(set([]))
                protein_spacial_list.append(set([]))

    Protein_num = index
    f.close()
    f_protein_out.close()

    ###########input Gene expression data##################  # 输入基因表达数据

    f = open(options.EXPRESSION_file, "r")

    for line in f:
        line = line.strip().split()
        if len(line) == 38:
            if line[1] in Dic_map:

                for i in range(2, Time_num + 2):
                    expression_list[Dic_map[line[1]]].append(
                        (float(line[i]) + float(line[i + 12]) + float(line[i + 24])) / 3)
                    # 计算蛋白质节点在12个时间节点上的活性表达值（三个周期的平均值）

    f.close()

    ###############compute active time attribute for proteins###################  # 计算蛋白质活性时间属性（包括活跃时间点及对应的活性概率）

    for t in range(0, Time_num):
        time_protein_list.append(set([]))  # 建立12个时间集合

    for instance in Dic_map:

        if len(expression_list[Dic_map[instance]]) >= 12:  # 如果该蛋白质有基因表达数据
            Temp_mean_value = 0.
            Temp_sd_value = 0.
            Temp_thresh_3SD = 0.
            Temp_thresh_2SD = 0.
            Temp_thresh_1SD = 0.

            for j in range(0, Time_num):
                Temp_mean_value += expression_list[Dic_map[instance]][j]
            Temp_mean_value /= Time_num  # 计算蛋白质节点在12个时间节点上的平均表达值

            for j in range(0, Time_num):
                Temp_sd_value += (expression_list[Dic_map[instance]][j] - Temp_mean_value) ** 2
            Temp_sd_value /= (Time_num - 1)  # 计算蛋白质节点在12个时间节点上的表达值方差

            k = 1

            Temp_thresh_3SD = Temp_mean_value + 3 * (Temp_sd_value ** 0.5) * (Temp_sd_value / (1 + Temp_sd_value))
            Temp_thresh_2SD = Temp_mean_value + 2 * (Temp_sd_value ** 0.5) * (Temp_sd_value / (1 + Temp_sd_value))
            Temp_thresh_1SD = Temp_mean_value + 1 * (Temp_sd_value ** 0.5) * (Temp_sd_value / (1 + Temp_sd_value))
            Temp_thresh_kSD = Temp_mean_value + k * (Temp_sd_value ** 0.5) * (Temp_sd_value / (1 + Temp_sd_value))
            '''for j in range(0,Time_num):
                if expression_list[Dic_map[instance]][j]>=Temp_thresh_3SD:
                    protein_time_list[Dic_map[instance]].add(j)
                    Protein_SD_weight[instance+"&"+str(j)]=SD3_weight

                elif Temp_thresh_3SD>expression_list[Dic_map[instance]][j]>=Temp_thresh_2SD:
                    protein_time_list[Dic_map[instance]].add(j)
                    Protein_SD_weight[instance+"&"+str(j)]=SD2_weight

                elif Temp_thresh_2SD>expression_list[Dic_map[instance]][j]>=Temp_thresh_1SD:
                    protein_time_list[Dic_map[instance]].add(j)  # 记录蛋白质节点的活跃时间节点
                    Protein_SD_weight[instance+"&"+str(j)]=SD1_weight  # 记录蛋白质活跃时间节点对应的活性概率'''

            for j in range(0, Time_num):
                if expression_list[Dic_map[instance]][j] >= Temp_thresh_kSD:
                    protein_time_list[Dic_map[instance]].add(j)  # 记录蛋白质节点的活跃时间节点
                    time_protein_list[j].add(Dic_map[instance])  # 记录各时间活跃的蛋白质
    # ------------------------------------------------亚细胞位置信息----------------------------------------------
    f_c = open("./data/yeast_compartment_benchmark.txt", "r", encoding='utf-16')

    location_dict = {}
    for i, location in enumerate(locations):
        location_dict[location] = i

    for line in f_c:
        line = line.strip().split()
        if line[0] in Dic_map:
            subcellular_list[Dic_map[line[0]]].append(line[1])  # 记录蛋白质节点的活跃空间位置(位置文字表示)
            protein_spacial_list[Dic_map[line[0]]].add(location_dict[line[1]])  # 记录蛋白质节点的活跃空间位置（位置序号表示）

    # ---------------------------PRP-----------------------------

    proteinRprotein = sparse.load_npz('./data/' + dataset + '/protein-R-protein.npz')
    proteinRprotein = sparse_mx_to_torch_sparse_tensor(proteinRprotein).coalesce().indices()
    rNode1 = proteinRprotein[0]
    rNode2 = proteinRprotein[1]

    ###########input PPI data############  # 输入PPI数据（有共活跃时间或空间的PPI构建的新PPI网络）
    f = open(PPI_file, "r")
    matrix_list = []

    for line in f:
        line = line.strip().split()
        if len(line) == 2:
            if len(protein_time_list[Dic_map[line[0]]] & protein_time_list[Dic_map[line[1]]]) > 0:
                # or len(protein_time_list[Dic_map[line[0]]]) == 0 or len(protein_time_list[Dic_map[line[1]]]) == 0:
                tNode1.append(Dic_map[line[0]])
                tNode2.append(Dic_map[line[1]])
            elif len(protein_spacial_list[Dic_map[line[0]]]) == 0 or len(protein_spacial_list[Dic_map[line[1]]]) == 0:
                # or len(protein_spacial_list[Dic_map[line[0]]] & protein_spacial_list[Dic_map[line[1]]]) > 0:
                sNode1.append(Dic_map[line[0]])
                sNode2.append(Dic_map[line[1]])
    tmatrix = coo_matrix((np.ones(len(tNode1), dtype=float64), (tNode1, tNode2)), shape=(len(Dic_map), len(Dic_map)))
    matrix_list.append(tmatrix)
    smatrix = coo_matrix((np.ones(len(sNode1), dtype=float64), (sNode1, sNode2)), shape=(len(Dic_map), len(Dic_map)))
    matrix_list.append(smatrix)
    rmatrix = coo_matrix((np.ones(len(rNode1), dtype=float64), (rNode1, rNode2)), shape=(len(Dic_map), len(Dic_map)))
    matrix_list.append(rmatrix)

    filtered_matrix_list = [matrix for matrix in matrix_list if matrix.nnz > 0]

    f.close()
    return filtered_matrix_list

def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = th.from_numpy(
        np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = th.from_numpy(sparse_mx.data)
    shape = th.Size(sparse_mx.shape)
    return th.sparse.FloatTensor(indices, values, shape)