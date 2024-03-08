# import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
#calculate the overlap number of the complex i get with the known conplex
#file=open("dip_complex_degree.txt")
dataset = "biogrid"
p = '../data' + '/' + dataset + '/'
str='../data/result/final_' + dataset + '_attr_output'
# str="result/result.txt"
file=open(str)
file1=open("../data/Form_CYC20083.txt")
# file1=open("../data/golden_standard.txt")

#g=nx.Graph()
predicted_num=len(file.readlines())
reference_num=len(file1.readlines())
file.close()
file1.close()

file=open(str)
file1=open("../data/Form_CYC20083.txt")
# file1=open("../data/golden_standard.txt")
reference_complex=[]
for j in file1:
    j = j.rstrip()
    j = j.rstrip('\n')
    complex_list = j.split('	')
    reference_complex.append(complex_list)

predicted_complex=[]
for i in file:
    i = i.rstrip()
    i = i.rstrip('\n')
    node_list = i.split(' ')
    node_list.remove(node_list[0])

    predicted_complex.append(node_list)
#precision
number=0
c_number=0
row=1
for i in predicted_complex:
    overlapscore=0.0
    for j in reference_complex:
        set1 = set(i)
        set2 = set(j)
        overlap = set1 & set2
        score = float((pow(len(overlap), 2))) / float((len(set1) * len(set2)))
        if (score > overlapscore):
            overlapscore = score
    if (overlapscore > 0.20):
            number = number + 1
            print(row,)
            print(" ",)
    if len(i)==5:
        print(i)
    row=row+1
#recall
for i in reference_complex:
    overlapscore=0.0
    for j in predicted_complex:
        set1 = set(i)
        set2 = set(j)
        overlap = set1 & set2
        score = float((pow(len(overlap), 2))) / float((len(set1) * len(set2)))
        if(score > overlapscore):
            overlapscore = score
    if(overlapscore > 0.20):
        c_number=c_number+1
#sn
T_sum1=0.0
N_sum=0.0
for i in reference_complex:
    max=0.0
    for j in predicted_complex:
        set1=set(i)
        set2=set(j)
        overlap=set1&set2
        if len(overlap)>max:
            max=len(overlap)
    T_sum1=T_sum1+max
    N_sum=N_sum+len(set1)
#ppv
T_sum2=0.0
T_sum=0.0
for i in predicted_complex:
    max=0.0
    for j in reference_complex:
        set1=set(i)
        set2=set(j)
        overlap=set1&set2
        T_sum=T_sum+len(overlap)
        if len(overlap)>max:
            max=len(overlap)
    T_sum2=T_sum2+max

print("\n")
print(number,predicted_num)# matched predicted complex number
#print c_number,reference_num# matched reference complex number
precision=float(number/float(predicted_num))
recall=float(c_number/float(reference_num))
F1=float((2*precision*recall)/(precision+recall))
Sn=float(T_sum1)/float(N_sum)
PPV=float(T_sum2)/float(T_sum)
Acc=pow(float(Sn*PPV),0.5)
print("precision",)
print(precision)
print("recall",)
print(recall)
print("F1",)
print(F1)
print("Sn",)
print(Sn)
print("PPV",)
print(PPV)
print("Acc",)
print(Acc)
file.close()
file1.close()

#nx.draw(g)
#plt.savefig("wuxiangtu.png")
#plt.show()