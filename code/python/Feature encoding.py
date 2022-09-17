# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 18:43:40 2021

@author: wenjin
"""

import pandas as pd
import numpy as np
from collections import Counter
import math
import re
import itertools
import os
from numpy import *
from pandas.core.frame import DataFrame 
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import KFold





def compositiontransition(aminoseq):
    xname=pd.Series(list(aminoseq))
    aminocomposition=xname.value_counts(normalize=True)
    name=list(aminocomposition.index)
    tdata=aminocomposition.values
    oname=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    data=np.zeros(20)
    for i in range(len(name)):
        j=oname.index(name[i])
        data[j]=tdata[i]
    return data

def vandertransition(aminoseq): 
    #氨基酸序列被字母P,N,H 代替
    x1=aminoseq.replace("G","P")
    x2=x1.replace("A","P")
    x3=x2.replace("S","P")
    x4=x3.replace("C","P")
    x5=x4.replace("T","P")
    x6=x5.replace("P","P")
    x7=x6.replace("D","P")
    ###################   
    x8=x7.replace("N","N")
    x9=x8.replace("V","N")
    x10=x9.replace("E","N")
    x11=x10.replace("Q","N")
    x12=x11.replace("I","N")
    x13=x12.replace("L","N")
    
    x14=x13.replace("M","H")
    x15=x14.replace("H","H")
    x16=x15.replace("K","H")
    x17=x16.replace("F","H")
    x18=x17.replace("R","H")
    x19=x18.replace("Y","H")
    x20=x19.replace("W","H")
    ####计算distribution
    p1=x20.count("P")
    li=np.array([i.start() for i in re.finditer('P',x20)])+1
    if p1!=0:
        ap=math.floor(0.25*p1)-1
        if ap==0:
            ap=0
        bp=math.floor(0.50*p1)-1
        if bp==0:
            bp=0
        cp=math.floor(0.75*p1)-1
        if cp==0:
            cp=0
        c=[0,ap,bp,cp,-1]
        dp1=li[c]/len(x20)
    else:
        dp1=np.array([0,0,0,0,0])
#
#提取N的位置distribution
    N1=x20.count("N")
    Nli=np.array([i.start() for i in re.finditer('N',x20)])+1
    if N1!=0:
        aN=math.floor(0.25*N1)-1
        if aN==0:
            aN=0
        bN=math.floor(0.50*N1)-1
        if bN==0:
            bN=0
        cN=math.floor(0.75*N1)-1
        if cN==0:
            cN=0
        c=[0,aN,bN,cN,-1]
        dN1=Nli[c]/len(x20)
    else:
        dN1=np.array([0,0,0,0,0])    
        
#提取H的位置distribution
    H1=x20.count("H")
    Hli=np.array([i.start() for i in re.finditer('H',x20)])+1
    if H1!=0:
        aH=math.floor(0.25*H1)-1
        if aH==0:
            aH=0
        bH=math.floor(0.50*H1)-1
        if bH==0:
            bH=0
        cH=math.floor(0.75*H1)-1
        if cH==0:
            cH=0
        c=[0,aH,bH,cH,-1]
        dH1=Hli[c]/len(x20)
    else:
        dH1=np.array([0,0,0,0,0]) 
    alld=np.hstack((dp1,dN1,dH1))
    ##计算组成composition
    cp1=p1/len(x20)
    cN1=N1/len(x20)
    cH1=H1/len(x20)
    allc=np.array([cp1,cN1,cH1])
    
    #计算转换数据transition
    x22=x20.replace("P","2")
    x23=x22.replace("N","0")
    x24=x23.replace("H","1")
    x25=np.array([int(x) for x in list(x24)])
    x26=np.delete(x25,0)
    x27=np.delete(x25,-1)
    tx=np.where(x26!=x27)[0]
    tx1=x25[tx]+x25[tx+1]
    tpn=Counter(tx1)[2]/len(x20)
    tph=Counter(tx1)[3]/len(x20)
    tnh=Counter(tx1)[1]/len(x20)
    allt=[tpn,tph,tnh]
    data=np.hstack([allc,allt,alld])
    return data

def hydoptransition(aminoseq):
    #氨基酸序列被字母P,N,H 代替
    x1=aminoseq.replace("G","Z")
    x2=x1.replace("A","Z")
    x3=x2.replace("S","Z")
    x4=x3.replace("P","Z")
    x5=x4.replace("H","Z")
    x6=x5.replace("Y","Z")
    x7=x6.replace("T","Z")
    ###################   
    x8=x7.replace("R","P")
    x9=x8.replace("K","P")
    x10=x9.replace("E","P")
    x11=x10.replace("Q","P")
    x12=x11.replace("N","P")
    x13=x12.replace("D","P")
    
    x14=x13.replace("C","H")
    x15=x14.replace("V","H")
    x16=x15.replace("L","H")
    x17=x16.replace("I","H")
    x18=x17.replace("M","H")
    x19=x18.replace("F","H")
    x20=x19.replace("W","H")
    x21=x20.replace("Z","N")
    ####计算distribution
    p1=x21.count("P")#提取P的位置
    li=np.array([i.start() for i in re.finditer('P',x21)])+1
    if p1!=0:
        ap=math.floor(0.25*p1)-1
        if ap==0:
            ap=0
        bp=math.floor(0.50*p1)-1
        if bp==0:
            bp=0
        cp=math.floor(0.75*p1)-1
        if cp==0:
            cp=0
        c=[0,ap,bp,cp,-1]
        dp1=li[c]/len(x21)
    else:
        dp1=np.array([0,0,0,0,0])
#
#提取N的位置distribution
    N1=x21.count("N")
    Nli=np.array([i.start() for i in re.finditer('N',x21)])+1
    if N1!=0:
        aN=math.floor(0.25*N1)-1
        if aN==0:
            aN=0
        bN=math.floor(0.50*N1)-1
        if bN==0:
            bN=0
        cN=math.floor(0.75*N1)-1
        if cN==0:
            cN=0
        c=[0,aN,bN,cN,-1]
        dN1=Nli[c]/len(x21)
    else:
        dN1=np.array([0,0,0,0,0])    
        
#提取H的位置distribution
    H1=x21.count("H")
    Hli=np.array([i.start() for i in re.finditer('H',x21)])+1
    if H1!=0:
        aH=math.floor(0.25*H1)-1
        if aH==0:
            aH=0
        bH=math.floor(0.50*H1)-1
        if bH==0:
            bH=0
        cH=math.floor(0.75*H1)-1
        if cH==0:
            cH=0
        c=[0,aH,bH,cH,-1]
        dH1=Hli[c]/len(x21)
    else:
        dH1=np.array([0,0,0,0,0]) 
    alld=np.hstack((dp1,dN1,dH1))
    ##计算组成composition
    cp1=p1/len(x21)
    cN1=N1/len(x21)
    cH1=H1/len(x21)
    allc=np.array([cp1,cN1,cH1])
    
    #计算转换数据transition
    x22=x21.replace("P","2")
    x23=x22.replace("N","0")
    x24=x23.replace("H","1")
    x25=np.array([int(x) for x in list(x24)])
    x26=np.delete(x25,0)
    x27=np.delete(x25,-1)
    tx=np.where(x26!=x27)[0]
    tx1=x25[tx]+x25[tx+1]
    tpn=Counter(tx1)[2]/len(x21)
    tph=Counter(tx1)[3]/len(x21)
    tnh=Counter(tx1)[1]/len(x21)
    allt=[tpn,tph,tnh]
    data=np.hstack([allc,allt,alld])
    return data

def polaritytransiton(aminoseq):
    #氨基酸序列被字母P,N,H 代替
    x1=aminoseq.replace("L","Z")
    x2=x1.replace("I","Z")
    x3=x2.replace("F","Z")
    x4=x3.replace("W","Z")
    x5=x4.replace("C","Z")
    x6=x5.replace("M","Z")
    x7=x6.replace("Y","Z") 
    x8=x7.replace("V","Z")

    ###################  
    x9=x8.replace("P","X")
    x10=x9.replace("A","X")
    x11=x10.replace("T","X")
    x12=x11.replace("G","X")
    x13=x12.replace("S","X")
    
    x14=x13.replace("H","H")
    x15=x14.replace("Q","H")
    x16=x15.replace("R","H")
    x17=x16.replace("K","H")
    x18=x17.replace("N","H")
    x19=x18.replace("E","H")
    x20=x19.replace("D","H")
    x201=x20.replace("Z","P")
    x21=x201.replace("X","N")

    ####计算distribution
    p1=x21.count("P")#提取P的位置
    li=np.array([i.start() for i in re.finditer('P',x21)])+1
    if p1!=0:
        ap=math.floor(0.25*p1)-1
        if ap==0:
            ap=0
        bp=math.floor(0.50*p1)-1
        if bp==0:
            bp=0
        cp=math.floor(0.75*p1)-1
        if cp==0:
            cp=0
        c=[0,ap,bp,cp,-1]
        dp1=li[c]/len(x21)
    else:
        dp1=np.array([0,0,0,0,0])
#
#提取N的位置distribution
    N1=x21.count("N")
    Nli=np.array([i.start() for i in re.finditer('N',x21)])+1
    if N1!=0:
        aN=math.floor(0.25*N1)-1
        if aN==0:
            aN=0
        bN=math.floor(0.50*N1)-1
        if bN==0:
            bN=0
        cN=math.floor(0.75*N1)-1
        if cN==0:
            cN=0
        c=[0,aN,bN,cN,-1]
        dN1=Nli[c]/len(x21)
    else:
        dN1=np.array([0,0,0,0,0])    
        
#提取H的位置distribution
    H1=x21.count("H")
    Hli=np.array([i.start() for i in re.finditer('H',x21)])+1
    if H1!=0:
        aH=math.floor(0.25*H1)-1
        if aH==0:
            aH=0
        bH=math.floor(0.50*H1)-1
        if bH==0:
            bH=0
        cH=math.floor(0.75*H1)-1
        if cH==0:
            cH=0
        c=[0,aH,bH,cH,-1]
        dH1=Hli[c]/len(x21)
    else:
        dH1=np.array([0,0,0,0,0]) 
    alld=np.hstack((dp1,dN1,dH1))
    ##计算组成composition
    cp1=p1/len(x21)
    cN1=N1/len(x21)
    cH1=H1/len(x21)
    allc=np.array([cp1,cN1,cH1])
    
    #计算转换数据transition
    x22=x21.replace("P","2")
    x23=x22.replace("N","0")
    x24=x23.replace("H","1")
    x25=np.array([int(x) for x in list(x24)])
    x26=np.delete(x25,0)
    x27=np.delete(x25,-1)
    tx=np.where(x26!=x27)[0]
    tx1=x25[tx]+x25[tx+1]
    tpn=Counter(tx1)[2]/len(x21)
    tph=Counter(tx1)[3]/len(x21)
    tnh=Counter(tx1)[1]/len(x21)
    allt=[tpn,tph,tnh]
    data=np.hstack([allc,allt,alld])
    return data

def polarizability(aminoseq):
    #氨基酸序列被字母P,N,H 代替
    x1=aminoseq.replace("G","Z")
    x2=x1.replace("A","Z")
    x3=x2.replace("S","Z")
    x4=x3.replace("T","Z")
    x5=x4.replace("D","Z")

    x6=x5.replace("C","N")
    x7=x6.replace("P","N") 
    x8=x7.replace("N","N")
    x9=x8.replace("V","N")
    x10=x9.replace("E","N")
    x11=x10.replace("Q","N")
    x12=x11.replace("I","N")
    x13=x12.replace("L","N")
    
    x14=x13.replace("K","H")
    x15=x14.replace("M","H")
    x16=x15.replace("H","H")
    x17=x16.replace("F","H")
    x18=x17.replace("R","H")
    x19=x18.replace("Y","H")
    x20=x19.replace("W","H")
    x21=x20.replace("Z","P")


    ####计算distribution
    p1=x21.count("P")#提取P的位置
    li=np.array([i.start() for i in re.finditer('P',x21)])+1
    if p1!=0:
        ap=math.floor(0.25*p1)-1
        if ap==0:
            ap=0
        bp=math.floor(0.50*p1)-1
        if bp==0:
            bp=0
        cp=math.floor(0.75*p1)-1
        if cp==0:
            cp=0
        c=[0,ap,bp,cp,-1]
        dp1=li[c]/len(x21)
    else:
        dp1=np.array([0,0,0,0,0])
#
#提取N的位置distribution
    N1=x21.count("N")
    Nli=np.array([i.start() for i in re.finditer('N',x21)])+1
    if N1!=0:
        aN=math.floor(0.25*N1)-1
        if aN==0:
            aN=0
        bN=math.floor(0.50*N1)-1
        if bN==0:
            bN=0
        cN=math.floor(0.75*N1)-1
        if cN==0:
            cN=0
        c=[0,aN,bN,cN,-1]
        dN1=Nli[c]/len(x21)
    else:
        dN1=np.array([0,0,0,0,0])    
        
#提取H的位置distribution
    H1=x21.count("H")
    Hli=np.array([i.start() for i in re.finditer('H',x21)])+1
    if H1!=0:
        aH=math.floor(0.25*H1)-1
        if aH==0:
            aH=0
        bH=math.floor(0.50*H1)-1
        if bH==0:
            bH=0
        cH=math.floor(0.75*H1)-1
        if cH==0:
            cH=0
        c=[0,aH,bH,cH,-1]
        dH1=Hli[c]/len(x21)
    else:
        dH1=np.array([0,0,0,0,0]) 
    alld=np.hstack((dp1,dN1,dH1))
    ##计算组成composition
    cp1=p1/len(x21)
    cN1=N1/len(x21)
    cH1=H1/len(x21)
    allc=np.array([cp1,cN1,cH1])
    
    #计算转换数据transition
    x22=x21.replace("P","2")
    x23=x22.replace("N","0")
    x24=x23.replace("H","1")
    x25=np.array([int(x) for x in list(x24)])
    x26=np.delete(x25,0)
    x27=np.delete(x25,-1)
    tx=np.where(x26!=x27)[0]
    tx1=x25[tx]+x25[tx+1]
    tpn=Counter(tx1)[2]/len(x21)
    tph=Counter(tx1)[3]/len(x21)
    tnh=Counter(tx1)[1]/len(x21)
    allt=[tpn,tph,tnh]
    data=np.hstack([allc,allt,alld])
    return data

def find(sequence,dul):
    k=sequence
    single=np.zeros(420)
    one=[a+b for a, b in zip(k[0:len(k)-1],k[1:len(k)])]
    k1=pd.Series(one)
    k2=pd.Series(list(k))
    kk=k1.value_counts(normalize=True)
    kk2=k2.value_counts(normalize=True)
    site=[dul.index(i) for i in kk.index]
    site2=[dul.index(i) for i in kk2.index]
    np.put(single,site,kk.values)
    np.put(single,site2,kk2.values)
    return single

############获得特征矩阵########
def trainmatrix(data,dul):
    aa=data['AA']
    xname=data['Entryname'].values 
    matrix84=np.zeros((len(aa),84))
    j=0
    for i in aa:
        np.put(matrix84[j],np.arange(0,21),hydoptransition(i))
        np.put(matrix84[j],np.arange(21,42),vandertransition(i))
        np.put(matrix84[j],np.arange(42,63),polaritytransiton(i))
        np.put(matrix84[j],np.arange(63,84),polarizability(i))
        j=j+1
    matrix420=np.zeros((len(aa),len(dul)))
    b=0
    for i in aa:
        np.put(matrix420[b],np.arange(0,len(dul)),find(i,dul))
        b=b+1 
    matrix=pd.DataFrame(np.c_[matrix420,matrix84],columns=nameall)
    return matrix

############获得特征矩阵########
def trainmatrix420(data,dul):
    aa=data['AA']
    xname=data['Entryname'].values   
    matrix420=np.zeros((len(aa),len(dul)))
    b=0
    for i in aa:
        np.put(matrix420[b],np.arange(0,len(dul)),find(i,dul))
        b=b+1  
    matrix=pd.DataFrame(matrix420,columns=dul,index=xname)
    return matrix

def number():
    path='F:/zhulin/1datasz1/RNA/predict/predictRBPspot/outputs'
    files=os.listdir(path)
    s=[]
    for file in files:
        a=file[8:]
        s.append(int(a[:-7]))
    return s

#####fasta 数据
def read_data():
    numb=number()
    protein_name=[]
    protein_seq=[]
    path='F:/zhulin/1datasz1/RNA/predict/predictRBP'
    for num in numb:
        with open(path+'/'+'predict_{}.fasta'.format(num),'r') as f:
            lines=f.readlines()
            head=lines[0].strip()
            protein_name.append(head[1:])
            seq=lines[1].strip()
            protein_seq.append(seq[0:])
    dic={"Entryname":protein_name,
         "AA":protein_seq}
    data=pd.DataFrame(dic)
    return data

###统计amino acid length
data=read_data()
aa=[len(list(data['AA'][i])) for i in range(0,421)]
seq_len=pd.DataFrame(aa)
outpath="F:/zhulin/1datasz1/RNA/predict/data"
seq_len.to_csv(os.path.join(outpath,"predict_sequence_length.csv"))

###获取504的feature
name182=["x"+str(i) for i in list(range(1,85))]
acid="ACDEFGHIKLMNPQRSTVWY"
dacid=[]
for i in itertools.product(acid,repeat=2):##产生400个变量特征的名字
    dacid.append("".join(i))
dul=dacid+list(acid)
nameall=dul+name182
data=read_data()
matrix=trainmatrix(data,dul)


##structure feature matrix
struc_data=pd.read_csv('F:/zhulin/1datasz1/RNA/predict/data/predict_structure_feature.csv')
structure_x=struc_data.iloc[:,1:]
##把原来的feature与structure feature 合并为700dim
all_x=pd.concat([matrix,structure_x],axis=1)

outpath="F:/zhulin/1datasz1/RNA/predict/data"
all_x.to_csv(os.path.join(outpath,"all_x_matrix_predict.csv"))



###
#把csv文件里的序列转化为fasta文件
data_pre=pd.read_csv("F:/zhulin/1datasz1/RNA/predictRBP.csv")
def fast(data):
    name=np.array(data['Entryname'])
    seq=np.array(data['AA'])
    i=1
    for d in range(0,len(seq)):
        fw=open('F:/zhulin/1datasz1/RNA/predictRBP/predict_{}.fasta'.format(i),'w')
        fw.write('>'+name[d]+'\n')
        fw.write(seq[d]+'\n')
        i+=1
        fw.close()
fast(data_pre)
