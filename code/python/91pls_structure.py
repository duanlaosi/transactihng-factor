# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 19:52:26 2021

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


####5折交叉验证的y_predictor ,y_test
def Y_Predict_cross(data,dul,featuremrmr,component):
    dataset=trainmatrix(data,dul,featuremrmr)##获得特征矩阵
    outdata=pd.DataFrame()
    kf=KFold(n_splits=5)
    for train, test in kf.split(dataset):
        df=return_predict(dataset.iloc[train],dataset.iloc[test],component)
        df['site']=np.concatenate((train, test), axis=0)
        outdata=pd.concat([df,outdata],axis=1,ignore_index=True)    
    outpath="F:/zhulin/1datasz1/RNA/outdata/"
    outdata.to_csv(os.path.join(outpath,"pls5fold_predict.csv")) 
    return outdata


def performance(y_predict,y_test):
   #定义RMSE,NRMSE,R2,Pearson,Spearman
   y_predict=y_predict.flatten()
   R2=r2_score(y_test,y_predict)
   rmse=np.sqrt(mean_squared_error(y_test, y_predict))
   nrmse=rmse/np.std(y_test)
   pear=pearsonr(y_test, y_predict)[0]
   spear=spearmanr(y_test,y_predict)[0]
   data=[R2,rmse,nrmse,pear,spear]
   return data

#随机划分的PLSregression########  
def plsmodel(x,y,component):
    ######pls##########
    ##划分数据集为train和test
    x_train,x_test,y_train,y_test=train_test_split(x,y,test_size = 0.2)
    ##training set 建立模型
    pls_model_set=PLSRegression(n_components=component)###scale=True默认标准化数据
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
def output(structure_x,fold85,repeat,component):
    ###重复训练1000次
    outdata=pd.DataFrame(columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))
    for i in range(0,repeat):
        outdata.loc[i]=plsmodel(structure_x,fold85,component)
    out=np.array(outdata.mean(axis=0))
#    outpath="F:/zhulin/1datasz1/RNA/outdata/"
#    outdata.to_csv(os.path.join(outpath,"psl_component{}.csv".format(component)))
    return out

structure_feature=pd.read_csv('F:/zhulin/1datasz1/RNA/outdata/structure_feature.csv')
structure_x=structure_feature.iloc[:,1:-1]###
structure_y=structure_feature['fold']
#structure_x.index=np.array(structure_feature['name'])###X--structure feature

####196个特征，全部structure feature 未做feature selection
x=structure_x
y=np.log10(structure_y).values
repeat=100
feature196_component10=[]
component=list(range(1,11))
for i in component:
    outdata=output(x,y,repeat,i)
    feature196_component10.append(outdata)
print(feature196_component10)    
out=pd.DataFrame(feature196_component10,columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))    
outpath="F:/zhulin/1datasz1/RNA/outdata/"
out.to_csv(os.path.join(outpath,"feature196_component1_10.csv"))

####
##mrmr特征筛选,component=2
structure_mrmr=pd.read_csv('F:/zhulin/1datasz1/RNA/outdata/structure_196_mrm.csv')###读取mrmr数据
FeatureName=list(structure_mrmr['FeatureName'])##提取mrmr order feature name
repeat=100
outdata=[]
for i in range(10,196,10):
    name_mrmr=FeatureName[0:i]
    x=structure_x[name_mrmr]
    outdata.append(output(x,y,repeat,2))
print(outdata)
out=pd.DataFrame(outdata,columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea')) 
outpath="F:/zhulin/1datasz1/RNA/outdata/"
out.to_csv(os.path.join(outpath,"feature196_mrmr10_190.csv"))


###
#intronic data的structure feature
intronic=pd.read_csv('F:/zhulin/1datasz1/RNA/Intronic_Splicing_Model_(GUR).csv')
intronic=intronic.drop([67])
intronic_name=list(intronic['Name-New'])  
name91=list(structure_feature['name'])
intronic_fold=np.log10(intronic['Fold.1']).values
site=[]
for i in intronic_name:
    if i in name91:
        site.append(name91.index(i))
intronic_feature=structure_feature.iloc[site,:].reset_index()
site1=[intronic_name.index(i) for i in list(intronic_feature['name'])]
intro_fold=intronic_fold[site1]
intronic_feature['fold']=list(intro_fold)
intronic_feature.drop(['index'],axis=1,inplace=True)
###把可以做mrmr的intronic的structure feature 保存
outpath="F:/zhulin/1datasz1/RNA/outdata/"
intronic_feature.to_csv(os.path.join(outpath,"intronic_structure_feature.csv"))



def return_predict(data_train,data_test,component):####返回预测的值y_predict,y_test
     ######pls##########
    ##划分数据集为train和test
    data_train,data_test=train_test_split(dataset,test_size = 0.2)
    ##training set 建立模型
    x_trains=np.array(data_train.loc[:,featuremrmr])
    y_trains=np.array(data_train.loc[:,'class'])
    ####数据标准化处理

    pls_model_set=PLSRegression(n_components=component)
    pls_model=pls_model_set.fit(x_train,y_train)
    ####
    y_predict=pls_model.predict(x_train)
    yy=list(y_predict.reshape(len(y_trains),))
    datatrain=pd.DataFrame({"predict":yy,"measured":list(y_trains),"class":np.repeat('train',len(yy))})
    #plt.scatter(y_predict,y_train)
    ###testing数据集
    x_tests=np.array(data_test.loc[:,featuremrmr])
    y_tests=np.array(data_test.loc[:,'class'])
    ###数据标准化处理
    x_test=preprocessing.scale(x_tests)
    y_test_predict=pls_model.predict(x_test)
    yy_test=list(y_test_predict.reshape(len(y_tests),))
    datatest=pd.DataFrame({'predict':yy_test,'measured':list(y_tests),'class':np.repeat('test',len(yy_test))})
    datafold=pd.concat([datatrain,datatest],axis=0,ignore_index=True)
    sns.lmplot(data = datafold,x='predict',y='measured',fit_reg=False,hue='class')
    return datafold 

#交叉验证的PLSregression########
def plsmodel(data_train,data_test,component):
    ######pls##########
    ##划分数据集为train和test
    ##training set 建立模型
    x_trains=np.array(data_train.loc[:,featuremrmr])
    y_trains=np.array(data_train.loc[:,'class'])
    ####数据标准化处理
    x_train=preprocessing.scale(x_trains)
    y_train=preprocessing.scale(y_trains)
    pls_model_set=PLSRegression(n_components=component)
    pls_model=pls_model_set.fit(x_train,y_train)
    ####
    y_predict=pls_model.predict(x_train)
    #plt.scatter(y_predict,y_train)
    train_outdata=performance(y_predict,y_train)
    ###testing数据集
    x_tests=np.array(data_test.loc[:,featuremrmr])
    y_tests=np.array(data_test.loc[:,'class'])
    ###数据标准化处理
    x_test=preprocessing.scale(x_tests)
    y_test=preprocessing.scale(y_tests)
    y_test_predict=pls_model.predict(x_test)
    test_outdata=performance(y_test_predict,y_test)
    ###合并training和testing结果
    outdata=train_outdata+test_outdata
    return outdata  

def allmean(outdata):
    outdata.loc['mean']=np.array(outdata.mean(axis=0))
    outdata.loc['var']=np.array(outdata.var(axis=0))
    outdata.loc['std']=np.array(outdata.std(axis=0))
    return outdata

##5折交叉验证的performance结果输出
def coss(data,dul,featuremrmr,component):
    dataset=trainmatrix(data,dul,featuremrmr)##获得特征矩阵
    outdata=[]
    kf=KFold(n_splits=5)
    for train, test in kf.split(dataset):
        outdata.append(plsmodel(dataset.iloc[train],dataset.iloc[test],component))
    df_koutdata=pd.DataFrame(outdata)
    df_koutdata.columns=['train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea']
    k_outdata=allmean(df_koutdata)  
    return k_outdata
k_outdata=coss(data,dul,featuremrmr,component)
outpath="F:/zhulin/1datasz1/RNA/outdata/"
k_outdata.to_csv(os.path.join(outpath,"pls-10-fold.csv")) 
    
##intronic       
####196个特征，全部structure feature 未做feature selection
x=intronic_feature.iloc[:,1:-1]
y=intronic_feature.iloc[:,0]
repeat=100
intronic196_component10=[]
component=list(range(1,11))
for i in component:
    outdata=output(x,y,repeat,i)
    intronic196_component10.append(outdata)
print(intronic196_component10)    
out=pd.DataFrame(intronic196_component10,columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))    
outpath="F:/zhulin/1datasz1/RNA/outdata/"
out.to_csv(os.path.join(outpath,"intronic196_component1_10.csv"))

##
##mrmr特征筛选,component=2
intronic_structure_mrmr=pd.read_csv('F:/zhulin/1datasz1/RNA/outdata/intronic_structure_mrmr.csv')###读取mrmr数据
intronic_FeatureName=list(intronic_structure_mrmr['FeatureName'])##提取mrmr order feature name
repeat=100
outdata=[]
for i in range(10,196,10):
    name_mrmr=FeatureName[0:i]
    x=intronic_feature[name_mrmr]
    outdata.append(output(x,y,repeat,10))
print(outdata)
out=pd.DataFrame(outdata,columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea')) 
outpath="F:/zhulin/1datasz1/RNA/outdata/"
out.to_csv(os.path.join(outpath,"intronic_structure_mrmr_10_190.csv"))

seq=[]
with open("C:/Users/wenjin/Desktop/91pls.fasta",'r') as f:
    for line in f.readlines():
        if '>' in line:
            pass
        else:
            new_line=line.replace('\n',"")
            seq.append(list(new_line))
#            seq.append([line.replace('\n',"")])
seq1=[list(i) for i in seq ]
for i in seq:
    d=list(i)


   