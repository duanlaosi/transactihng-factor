# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 21:02:22 2021

@author: wenjin
"""

def return_predict(dataset,component):####返回预测的值y_predict,y_test
     ######pls##########
    ##划分数据集为train和test
    data_train,data_test=train_test_split(dataset,test_size = 0.2)
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
    yy=list(y_predict.reshape(len(y_trains),))
    datatrain=pd.DataFrame({"train_predict":yy,"train_measured":list(y_trains)})
    #plt.scatter(y_predict,y_train)
    ###testing数据集
    x_tests=np.array(data_test.loc[:,featuremrmr])
    y_tests=np.array(data_test.loc[:,'class'])
    ###数据标准化处理
    x_test=preprocessing.scale(x_tests)
    y_test_predict=pls_model.predict(x_test)
    yy_test=list(y_test_predict.reshape(len(y_tests),))
    datatest=pd.DataFrame({'test_predict':yy_test,'test_measured':list(y_tests)})
    datafold=datatrain.append(datatest,ignore_index=True)
    return datafold 
#获得91个数据504个特征的矩阵数据
def output(data,dul,featuremrmr,repeat,component):
    dataset=trainmatrix(data,dul,featuremrmr)##获得特征矩阵
    ###重复训练1000次
    outdata=pd.DataFrame(columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))
    for i in range(0,repeat):
        outdata.loc[i]=plsmodel(dataset,component)
    out=np.array(outdata.mean(axis=0))
#    outpath="F:/zhulin/1datasz1/RNA/outdata/"
#    outdata.to_csv(os.path.join(outpath,"psl_component{}.csv".format(component)))
    return out

########
####循环调用代码 每次重复100次，然后重复10次
repeat=100
all=pd.DataFrame(columns=('train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea'))    
for i in range(0,10):
    all.loc[i]=(output(data,dul,featuremrmr,repeat,component))
outpath="F:/zhulin/1datasz1/RNA/outdata/"
all.to_csv(os.path.join(outpath,"pls_10_100.csv")) 
##component 2,repeat 100 to 5000, 
mean=[]
for i in list(range(100,5000,300)):
    mean.append(output(data,dul,featuremrmr,i,2))
df_mean=pd.DataFrame(mean)
df_mean.columns=['train_R^2','train_RMSE','train_NRMSE','train_Pear','train_Spea','test_R^2','test_RMSE','test_NRMSE','test_Pear','test_Spea']
plt.plot(df_mean['test_Pear'])
outpath="F:/zhulin/1datasz1/RNA/outdata/"
df_mean.to_csv(os.path.join(outpath,"psl_comp2_100_to_5000.csv"))


###读取数据
data1=pd.read_csv("F:/zhulin/1datasz1/RNA/plsvalidedataset.csv")
data=data1.drop([4,44,45,56])

#把csv文件里的序列转化为fasta文件
def fast(data):
    name=np.array(data['Entryname'])
    seq=np.array(data['AA'])
    i=1
    for d in range(0,91):
        fw=open('F:/zhulin/1datasz1/RNA/91pls/91pls_{}.fasta'.format(i),'w')
        fw.write('>'+name[d]+'\n')
        fw.write(seq[d]+'\n')
        i+=1
        fw.close()