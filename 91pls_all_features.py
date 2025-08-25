# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 08:50:42 2021

@author: wenjin
"""

import pandas as pd
import numpy as np
from collections import Counter
import math
import re
import itertools
from numpy import *
from pandas.core.frame import DataFrame 
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import KFold
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import feature_selection
from mrmr import mrmr_classif


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

def performance(y_predict,y_test):
   #定义RMSE,NRMSE,R2,Pearson,Spearman
   y_predict=y_predict.flatten()
   SS_R=sum((y_test-y_predict)**2)
   SS_T=sum((y_test-np.mean(y_test))**2)
   R2=1-(float(SS_R))/SS_T
   R2=r2_score(y_test,y_predict)
   rmse=np.sqrt(mean_squared_error(y_test, y_predict))
   nrmse=rmse/np.std(y_test)
   pear=pearsonr(y_test, y_predict)[0]
   spear=spearmanr(y_test,y_predict)[0]
   data=[R2,rmse,nrmse,pear,spear]
   return data

############获得特征矩阵########
def trainmatrix(data,dul):
    aa=data['AA']
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

#随机划分的PLSregression########  
def plsmodel(x,y,component):
    ######pls##########
    ##划分数据集为train和test
#    x_train,x_test,y_train,y_test=train_test_split(x,y,test_size = 0.2)
    ##training set 建立模型
    pls_model_set=PLSRegression(n_components=component,scale=True)###scale=True默认标准化数据
    pls_model=pls_model_set.fit(x_train,y_train)
    ####
    y_predict=pls_model.predict(x_train)
    #plt.scatter(y_predict,y_train)
    train_outdata=performance(y_predict,y_train)
    ###testing数据集
    y_test_predict=pls_model.predict(x_test)
    test_outdata=performance(y_test_predict,y_test)
    ###合并training和testing结果
    outdata=train_outdata+test_outdata
    return outdata 
        
####随即划分重复repeat次的均值输出
def output(x,y,repeat,component):
    ###重复训练次
    outdata=pd.DataFrame(columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))
    for i in range(0,repeat):
        outdata.loc[i]=plsmodel(x,y,component)
    out=np.array(outdata.mean(axis=0))
#    outpath="F:/zhulin/1datasz1/RNA/outdata/"
#    outdata.to_csv(os.path.join(outpath,"psl_component{}.csv".format(component)))
    return out



####读取数据###all feature
data1=pd.read_csv("F:/zhulin/1datasz1/RNA/plsvalidedataset.csv")
data91=data1.drop([4,44,45,56]).reset_index(drop=True) 
data85=data91.drop([32,33,35,36,38,41]).reset_index(drop=True)
Entryname=list(data85['Entryname'])

length=[len(data85['AA'][i]) for i in range(0,85)]
seq_len=pd.DataFrame(length)
outpath="F:/zhulin/1datasz1/RNA/putout/"
seq_len.to_csv(os.path.join(outpath,"all_x_seq_length.csv"))

##structure feature 
structure_feature=pd.read_csv('F:/zhulin/1datasz1/RNA/outdata/structure_feature.csv')
structure_x=structure_feature.iloc[:,1:-1]###

###设置feature name
##
name182=["x"+str(i) for i in list(range(1,85))]
acid="ACDEFGHIKLMNPQRSTVWY"
dacid=[]
for i in itertools.product(acid,repeat=2):##产生400个变量特征的名字
    dacid.append("".join(i))
dul=dacid+list(acid)
nameall=dul+name182

##把原来的feature与structure feature 合并为700dim
matrix=trainmatrix(data85,dul)
all_x=pd.concat([matrix,structure_x],axis=1)

outpath="F:/zhulin/1datasz1/RNA/putout/"
all_x.to_csv(os.path.join(outpath,"all_x_matrix.csv"))


y=np.log10(data85['fold']).values
x_train,x_test,y_train,y_test=train_test_split(all_x,y,test_size = 0.2)

####参照论文，使用论文中的123个feature做pls model
#######读取123个特征
feature123=pd.read_csv("F:/zhulin/1datasz1/RNA/features.csv")
feature123=list(feature123['Feature name'])
x_123=x_train

#使用700个特征做spls 模型
x=all_x
y=np.log10(data85['fold']).values
repeat=100
outdata700=[]
for i in range(1, 11):
    outdata700.append(output(x,y,repeat,i))
outdata700=pd.DataFrame(outdata700,columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))
outpath="F:/zhulin/1datasz1/RNA/putout/"
outdata700.to_csv(os.path.join(outpath,"feature700_component1_10.csv"))




###feature selection：
##1: 先删除全部为0的特征
mask=(all_x==0).all(0)
col_index=np.where(mask)[0]
all_x.drop(all_x.columns[col_index],axis=1,inplace=True)

##2：移除掉方差小于0.001的feature,剩余219个特征
var_feature= all_x.var().sort_values()
aa=np.where(var_feature.values>=0.001)
feature219=all_x[var_feature.index[aa]]

#
#先使用219feature做pls模型
###feature的pls 模型：
repeat=100
outdata=[]
x=feature219
for i in range(1,10):
    out=output(x,y,repeat,i)
    outdata.append(out)
outdata=pd.DataFrame(outdata,columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))
outpath="F:/zhulin/1datasz1/RNA/putout/"
outdata.to_csv(os.path.join(outpath,"feature219_component1_10.csv"))

###使用backward 选择特征：
def backward(feature219,y,repeat,component):
    repeat=100
    component=2
    featurename=list(feature219.columns)
    result=[]    
    com=featurename
    re_feature=[]
    x=feature219[featurename]
    result.append(output(x,y,repeat,component)[-1])
    for i in range(1,200):
        re_one=0
        ll=len(com)
        for d in range(0,ll):
            a=np.r_[0:i,i+1:ll]
            x=feature219.iloc[:,a]
            one=output(x,y,repeat,component)[-1]
            if one>re_one:
                re_one=one
                re_feature1=d
        result.append(re_one)
        com.pop(re_feature1)
        re_feature.append(re_feature1)
    return result, re_feature



####进一步做mrmr特征的筛选
outpath="F:/zhulin/1datasz1/RNA/putout/"
feature219.to_csv(os.path.join(outpath,"feature219.csv"))

##mrmr特征筛选,component=2

all_mrmr=pd.read_csv('F:/zhulin/1datasz1/RNA/outdata/all_feature_mrmr.csv')###读取mrmr数据
FeatureName=list(all_mrmr['FeatureName'])##提取mrmr order feature name
FeatureName[641]='NA'
repeat=100
outdata=[]
y=structure_feature['fold']
for i in range(10,700,10):
    name_mrmr=FeatureName[0:i]
    x=all_x[name_mrmr]
    outdata.append(output(x,y,repeat,2))
out=pd.DataFrame(outdata,columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea')) 
outpath="F:/zhulin/1datasz1/RNA/outdata/"
out.to_csv(os.path.join(outpath,"all_feature_mrmr10_700.csv"))




 