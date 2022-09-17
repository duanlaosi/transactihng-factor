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

def plsmodel(x,y,component):
    ######pls##########
    ##划分数据集为train和test
    x=preprocessing.scale(x)
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

############获得特征矩阵########
def trainmatrix(data,dul):
    aa=data['AA']
    xname=data['Entryname'].values   
    matrix420=np.zeros((len(aa),len(dul)))
    b=0
    for i in aa:
        np.put(matrix420[b],np.arange(0,len(dul)),find(i,dul))
        b=b+1  
    matrix=pd.DataFrame(matrix420,columns=dul,index=xname)
    return matrix

####随即划分重复1000次的输出
def output(x,y,repeat,component):
    ###重复训练1000次
    outdata=pd.DataFrame(columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))
    for i in range(0,repeat):
        outdata.loc[i]=plsmodel(x,y,component)
    out=np.array(outdata.mean(axis=0))
#    outpath="F:/zhulin/1datasz1/RNA/outdata/"
#    outdata.to_csv(os.path.join(outpath,"psl_component{}.csv".format(component)))
    return out

def backward(data,dul,featuremrmr):
    dataset=trainmatrix(data,dul,featuremrmr)
    result=[]    
    com=featuremrmr
    re_com=[]
    result[0]=plsmodel(dataset[list(com)+['class']],component)[0]
    for i in range(1,182):
        re_one=0
        for d in com:
            a=com
            a.remove(d)
            one=plsmodel(dataset[list(a)+['class']],component)[0]
            if one>re_one:
                re_one=one
                re_com1=d
        result.append(re_one)
        com.remove(re_com1)
        re_com.append(re_com1)
    return result, re_com


###读取数据
data1=pd.read_csv("F:/zhulin/1datasz1/RNA/plsvalidedataset.csv")
data91=data1.drop([4,44,45,56]).reset_index(drop=True)
data85=data91.drop([32,33,35,36,38,41]).reset_index(drop=True)   
 #######读取132个特征
feature123=pd.read_csv("F:/zhulin/1datasz1/RNA/features.csv")
feature123=list(feature123['Feature name'])
###设置feature name
acid="ACDEFGHIKLMNPQRSTVWY"
dacid=[]
for i in itertools.product(acid,repeat=2):##产生400个变量特征的名字
    dacid.append("".join(i))
dul=dacid+list(acid)

###
#
x=trainmatrix(data85,dul)[feature123]
y=np.log10(data85['fold']).values
repeat=100
out123=[]
for i in range(1,10):
    out123.append(output(x,y,repeat,i))
outdata=pd.DataFrame(out123,columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))
outpath="F:/zhulin/1datasz1/RNA/putout/"
outdata.to_csv(os.path.join(outpath,"out123_1to10.csv"))   
    


   



        

    
    
    
    
    
    
   
        