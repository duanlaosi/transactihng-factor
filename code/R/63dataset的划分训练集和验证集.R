library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)
library(caTools)
library(scorecard)
##读取数据
data=read.table("F:/zhulin/1datasz1/RNA/trainingdatasets.csv",header=TRUE,sep=",")
s=c(5,45,46,57)
data=data[-s,]

re=10
n=0.8
pre.value=rep(0,re)
for(i in 1:re)
{pre.value[i]=onetest63(n,data)
}

onetest63=function(n,data)
{
nr=nrow(data)
trainindex=sample(1:nr,n*nr)
train=data[trainindex,]
test=data[-trainindex,]
train.aa=train$AA  
train.xname=train$Name.1                    
train.class2=train$class2
##############################
test.aa=test$AA  
test.xname=test$Name.1                    
test.class2=test$class2
acid=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
dacid=c("AA","CC","DD","EE","FF","GG","HH","II","KK","LL","MM","NN","PP","QQ","RR","SS","TT","VV","WW","YY")
amino.acid=combn(acid,2)
acid.1=paste(c("^"),acid,c("$"),sep="")
dacid.1=paste(c("^"),dacid,c("$"),sep="")
dual.s=c(paste(c("^"),amino.acid[2,],amino.acid[1,],c("$"),sep=""),paste(c("^"),amino.acid[1,],amino.acid[2,],c("$"),sep=""),dacid.1,acid.1)
dual.sing=c(paste(amino.acid[2,],amino.acid[1,],sep=""),paste(amino.acid[1,],amino.acid[2,],sep=""),dacid,acid)
##sliding windows 计算频率
find=function(sequence){
k=sequence
str.one=str_split(k,"")[[1]]
kk=c(paste(str.one[-length(str.one)],str.one[-1],sep=""),str.one)
single=rep(0,420)
for(i in 1:420){
single[i]=length(grep(kk,pattern=dual.s[i]))
}
single=single/sum(single)
return(single)
}
train.matrix=matrix(0,length(train.aa),420)
for(i in 1:length(train.aa)){
train.matrix[i,]=find(train.aa[i])
}
rownames(train.matrix)=train.xname
colnames(train.matrix)=dual.sing
#############################################
##test.aa   sliding windows 计算频率
find=function(sequence){
k=sequence
str.one=str_split(k,"")[[1]]
kk=c(paste(str.one[-length(str.one)],str.one[-1],sep=""),str.one)
single=rep(0,420)
for(i in 1:420){
single[i]=length(grep(kk,pattern=dual.s[i]))
}
single=single/sum(single)
return(single)
}
test.matrix=matrix(0,length(test.aa),420)
for(i in 1:length(test.aa)){
train.matrix[i,]=find(test.aa[i])
}
rownames(test.matrix)=test.xname
colnames(test.matrix)=dual.sing
###################################
##直接使用123个特征
featuredata=read.table("F:/zhulin/1datasz1/RNA/features.csv",header=TRUE,sep=",")
featurename=featuredata[,1][-1]
site=match(featurename,dual.sing)
feature.name=dual.sing[site]
train.matrix123=train.matrix[,site]
test.matrix123=test.matrix[,site]

##PLD-DA分析
spls.train.matrix=splsda(train.matrix123,train.class2,ncomp=3,scale=TRUE)
#plotIndiv(spls.train.matrix,group=train.class2,col=c("blue","grey","red"),lty=2,ind.names=TRUE,ellipse=TRUE,legend=TRUE,X.label='Component 1',
#Y.label='Component 2',title='train-63-123',plot.title = element_text(hjust = 0))
#perf.train.plsda<- perf(spls.train.matrix, validation = "Mfold", folds =3,progressBar = TRUE, auc = TRUE, nrepeat = 10)
#plot(perf.train.plsda)
test.predict=predict(spls.train.matrix,test.matrix123,dist="max.dist")
prediction=test.predict$class$max.dist[,2]
confusion.mat=get.confusion_matrix(truth=test.class2,predict=prediction)
return(get.BER(confusion.mat))
}


